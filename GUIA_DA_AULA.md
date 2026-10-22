# Guia da aula — Pipeline para análise de variantes genéticas

**Bioinformática aplicada ao estudo do câncer · IEC – PUC Minas · Prof. Dr. Fábio Ribeiro Queiroz**

Este guia segue as mesmas etapas do roteiro da aula, já com os caminhos dos arquivos deste repositório. Copie e cole os comandos no terminal, **na ordem**, um bloco de cada vez.

> **Sobre os dados.** São 12 pacientes **simulados** (não são dados reais), sequenciados com um painel de captura de 13 genes: `BRCA1`, `BRCA2`, `PALB2`, `TP53`, `MLH1`, `MSH2`, `APC`, `PTEN`, `CDH1`, `KRAS`, `BRAF`, `PIK3CA` e `EGFR`. As leituras (Illumina 2×150 pb) foram geradas a partir das sequências reais do genoma humano **GRCh38** e carregam variantes reais, descritas no ClinVar — algumas patogênicas, outras benignas. Os dados clínicos estão em `dados/dados_clinicos.tsv`.
>
> **Sobre o genoma.** Em vez do genoma humano inteiro (~3 Gb), usamos `dados/referencia/genoma_aula.fa`: só as regiões dos 13 genes (~1,2 Mb), cada uma como um "contig" com o nome do gene. Isso faz a indexação levar segundos em vez de horas — mas exige uma etapa extra (7b) para devolver as variantes às coordenadas reais do GRCh38 antes do VEP.

---

## Antes de começar

> **Como copiar os comandos:** copie só o texto de dentro dos quadros cinza. Nunca copie as linhas ```` ```bash ```` e ```` ``` ```` que aparecem quando o arquivo é aberto como texto. Se colar essas linhas por engano, o terminal fica "preso": `ls` e os outros comandos param de mostrar resultados. Nesse caso, feche o terminal (ícone de lixeira 🗑) e abra outro (menu ☰ → *Terminal* → *New Terminal*).
>
> Para ver o guia formatado: botão direito em `GUIA_DA_AULA.md` → **Open Preview**.

Confira se o ambiente está pronto:

```bash
bash scripts/verificar_ambiente.sh
```

Escolha a amostra com que você vai trabalhar (o professor vai indicar) e crie as pastas de resultados:

```bash
AMOSTRA=P01
mkdir -p resultados/01_fastqc_bruto resultados/02_trimagem resultados/03_fastqc_trimado \
         resultados/04_alinhamento resultados/07_variantes resultados/08_para_o_VEP resultados/09_maf
```

> ⚠️ Se você fechar o terminal, rode de novo `AMOSTRA=P01` (com a sua amostra) antes de continuar.

---

## 1ª Etapa — Análise de qualidade dos dados brutos (FastQC)

```bash
fastqc --outdir resultados/01_fastqc_bruto dados/fastq/${AMOSTRA}_R1.fastq.gz dados/fastq/${AMOSTRA}_R2.fastq.gz
```

Abra os relatórios `resultados/01_fastqc_bruto/*_fastqc.html` (no Codespaces: clique com o botão direito no arquivo → **Show Preview**; ou → **Download** e abra no seu navegador).

**O que observar:** compare R1 e R2 em *Per base sequence quality* (a qualidade cai no fim das leituras, sobretudo no R2) e em *Adapter Content* (adaptador Illumina aparecendo no fim das leituras — sinal de fragmentos de DNA mais curtos que 150 pb, comum em amostras tumorais FFPE).

---

## 2ª Etapa — Remoção de adaptadores e bases de baixa qualidade (Trimmomatic)

No ambiente da aula, `java -jar trimmomatic.jar` é simplesmente `trimmomatic`:

```bash
trimmomatic PE -phred33 \
  dados/fastq/${AMOSTRA}_R1.fastq.gz dados/fastq/${AMOSTRA}_R2.fastq.gz \
  resultados/02_trimagem/${AMOSTRA}_R1_paired.fastq.gz resultados/02_trimagem/${AMOSTRA}_R1_unpaired.fastq.gz \
  resultados/02_trimagem/${AMOSTRA}_R2_paired.fastq.gz resultados/02_trimagem/${AMOSTRA}_R2_unpaired.fastq.gz \
  ILLUMINACLIP:dados/adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

Anote a linha final (`Input Read Pairs ... Both Surviving ...`): ela é a base da tabela de aproveitamento.

> **Curiosidade:** por que muitos pares aparecem como *Forward Only Surviving*? Quando o fragmento é menor que a leitura, R1 e R2 leem exatamente o mesmo trecho de DNA (em fitas opostas). No modo *palíndromo* do `ILLUMINACLIP`, o Trimmomatic remove o adaptador e descarta o R2, que seria redundante.

---

## 3ª Etapa — Qualidade após a trimagem

```bash
fastqc --outdir resultados/03_fastqc_trimado resultados/02_trimagem/${AMOSTRA}_R1_paired.fastq.gz resultados/02_trimagem/${AMOSTRA}_R2_paired.fastq.gz
```

**Tabela de dados brutos × dados trimados** (número de leituras = linhas do FASTQ ÷ 4):

```bash
for f in dados/fastq/${AMOSTRA}_R1.fastq.gz resultados/02_trimagem/${AMOSTRA}_R1_paired.fastq.gz; do
  echo -e "$f\t$(( $(zcat $f | wc -l) / 4 )) leituras"
done
```

---

## 4ª Etapa — Alinhamento ao genoma de referência (BWA)

**Comando 1** — indexar o genoma (só uma vez para toda a turma; leva ~1 segundo):

```bash
bwa index dados/referencia/genoma_aula.fa
samtools faidx dados/referencia/genoma_aula.fa
```

**Comando 2** — alinhar. Acrescentamos `-R` (*read group*), que grava o nome da amostra no alinhamento. Sem ele, o VCF sairia com o nome do arquivo no lugar do nome da amostra — e não conseguiríamos juntar os 12 pacientes depois.

```bash
bwa mem -R "@RG\tID:${AMOSTRA}\tSM:${AMOSTRA}\tPL:ILLUMINA" dados/referencia/genoma_aula.fa \
  resultados/02_trimagem/${AMOSTRA}_R1_paired.fastq.gz resultados/02_trimagem/${AMOSTRA}_R2_paired.fastq.gz \
  > resultados/04_alinhamento/${AMOSTRA}.sam
```

Dê uma olhada no SAM (cabeçalho e primeiras leituras):

```bash
head -n 20 resultados/04_alinhamento/${AMOSTRA}.sam | cut -c 1-200
```

---

## Etapa 5 — Converter SAM em BAM

```bash
samtools view -S -b resultados/04_alinhamento/${AMOSTRA}.sam > resultados/04_alinhamento/${AMOSTRA}.bam
ls -lh resultados/04_alinhamento/
```

Compare o tamanho do `.sam` e do `.bam`.

---

## Etapa 6 — Ordenar, indexar e marcar duplicatas

```bash
samtools sort resultados/04_alinhamento/${AMOSTRA}.bam -o resultados/04_alinhamento/${AMOSTRA}_sorted.bam
samtools index resultados/04_alinhamento/${AMOSTRA}_sorted.bam
```

Marcar duplicatas de PCR (recomendado para dados de câncer):

```bash
picard MarkDuplicates \
  I=resultados/04_alinhamento/${AMOSTRA}_sorted.bam \
  O=resultados/04_alinhamento/${AMOSTRA}_dedup.bam \
  M=resultados/04_alinhamento/${AMOSTRA}_dup_metrics.txt
samtools index resultados/04_alinhamento/${AMOSTRA}_dedup.bam
```

Qual a porcentagem de duplicatas? (coluna `PERCENT_DUPLICATION`)

```bash
grep -A1 "^LIBRARY" resultados/04_alinhamento/${AMOSTRA}_dup_metrics.txt | cut -f 3,7,9
samtools flagstat resultados/04_alinhamento/${AMOSTRA}_dedup.bam
```

---

## Etapa 7 — Chamada e filtragem de variantes (BCFtools)

```bash
bcftools mpileup -Ou -f dados/referencia/genoma_aula.fa resultados/04_alinhamento/${AMOSTRA}_dedup.bam \
  | bcftools call -mv -Ov -o resultados/07_variantes/${AMOSTRA}_variants.vcf

bcftools filter -e 'QUAL<30 || DP<10' -Ov -o resultados/07_variantes/${AMOSTRA}_filtered_variants.vcf \
  resultados/07_variantes/${AMOSTRA}_variants.vcf
```

Quantas variantes antes e depois do filtro? E onde elas estão?

```bash
grep -vc "^#" resultados/07_variantes/${AMOSTRA}_variants.vcf
grep -vc "^#" resultados/07_variantes/${AMOSTRA}_filtered_variants.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT\t%AD]\n' resultados/07_variantes/${AMOSTRA}_filtered_variants.vcf
```

Repare na primeira coluna: `TP53`, `BRCA1`... e posições pequenas. São coordenadas **dentro do recorte** de cada gene, não do cromossomo.

---

## Etapa 7b — Restaurar as coordenadas reais do GRCh38

*(Etapa extra, necessária apenas porque usamos o genoma reduzido.)*

```bash
bash scripts/restaurar_coordenadas.sh resultados/07_variantes/${AMOSTRA}_filtered_variants.vcf \
  resultados/07_variantes/${AMOSTRA}_GRCh38.vcf

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD]\n' resultados/07_variantes/${AMOSTRA}_GRCh38.vcf
```

Agora `TP53  8480` virou `chr17  7675088`. O script soma o início de cada recorte (tabela `dados/referencia/regioes.tsv`) e também normaliza os indels (`bcftools norm`), como se faz antes de qualquer anotação.

---

## Processar todos os pacientes de uma vez

Para analisar a coorte no maftools precisamos das 12 amostras. O script abaixo repete as Etapas 1 a 7b para cada paciente (mostrando cada comando) e junta tudo em um único VCF (~4 minutos):

```bash
bash scripts/rodar_pipeline.sh
```

Ao final: `resultados/resumo_qc.tsv` (tabela de aproveitamento, duplicatas e profundidade de todas as amostras) e **`resultados/08_para_o_VEP/coorte_GRCh38.vcf`**, o arquivo que vai para o VEP.

```bash
cat resultados/resumo_qc.tsv
```

> Se você travou em alguma etapa individual, este script também serve de "resgate": ele gera tudo de que você precisa para continuar.

---

## Etapa 8 — Anotação das variantes com o VEP (online)

1. **Baixe** `resultados/08_para_o_VEP/coorte_GRCh38.vcf` (botão direito → **Download**).
2. Acesse **https://www.ensembl.org/Tools/VEP** → *Species*: **Human (Homo sapiens) – GRCh38**.
3. Em *Upload file*, envie o VCF.
4. Opções (as que não forem citadas podem ficar no padrão):
   - **Identifiers**: marque **HGVS**.
   - **Variants and frequency data**: mantenha *Find co-located known variants* e marque **gnomAD (exomes) allele frequencies**.
   - **Predictions**: SIFT e PolyPhen (padrão).
   - **Filtering options → Restrict results**: **Show one selected consequence** (equivale ao `--pick` do roteiro).
   - ⚠️ **Não** marque *Exclude common variants* — queremos ver os polimorfismos.
5. **Run**. Quando terminar, explore os resultados na própria página (filtros, gráficos, links para ClinVar/gnomAD).
6. Baixe o resultado em **VCF**, renomeie para `coorte_anotada_VEP.vcf` e **arraste o arquivo** para a pasta `resultados/08_para_o_VEP/` no explorador de arquivos do Codespaces.

Equivalente em linha de comando (precisa do cache do VEP, ~20 GB — por isso usamos a versão online):
`vep -i coorte_GRCh38.vcf -o coorte_anotada_VEP.vcf --vcf --pick --cache --assembly GRCh38 --hgvs --af_gnomade`

> **Plano B** — se o site do VEP estiver fora do ar ou lento, use o arquivo já anotado (gerado com a API do próprio Ensembl VEP a partir deste mesmo VCF):
> ```bash
> cp dados/plano_b_VEP/coorte_anotada_VEP.vcf resultados/08_para_o_VEP/
> ```

---

## Etapa 9 — Converter o VCF anotado em MAF

**Comando 1 — extrair os campos do CSQ com o `split-vep`.** Em arquivos com várias amostras, os campos por amostra (`%SAMPLE`, `%GT`) precisam ficar entre colchetes `[ ]`; `-i 'GT="alt"'` mantém só as amostras que carregam a variante:

```bash
bcftools +split-vep -s worst -d -i 'GT="alt"' \
  -f '[%CHROM\t%POS\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%SAMPLE\n]' \
  resultados/08_para_o_VEP/coorte_anotada_VEP.vcf > resultados/09_maf/variantes_extraidas.tsv

head resultados/09_maf/variantes_extraidas.tsv
```

**Comando 2 — converter para MAF.** O script faz a extração acima com mais campos (HGVS, ClinVar, frequência no gnomAD, profundidade) e usa `awk` para traduzir as consequências do VEP para as categorias do maftools (`Variant_Classification`), classificar o tipo de variante (`Variant_Type`) e criar a coluna `HGVSp_Short` (p.Arg175His → p.R175H):

```bash
bash scripts/vcf_para_maf.sh resultados/08_para_o_VEP/coorte_anotada_VEP.vcf resultados/09_maf/coorte.maf
```

---

## Etapa 10 — Análise do perfil das variantes com maftools (R)

Rode tudo de uma vez (as figuras vão para `resultados/10_maftools/`):

```bash
Rscript scripts/analise_maftools.R resultados/09_maf/coorte.maf
```

Ou, para ir passo a passo, abra o R digitando `R` no terminal e rode os blocos do arquivo `scripts/analise_maftools.R` (saia do R com `q()`).

**Para discutir:** compare `02_oncoplot_todas_variantes.png` com `04_oncoplot_variantes_raras.png`. Por que *BRCA2* aparece "mutado" em 100% dos pacientes no primeiro gráfico? E por que o paciente P05 (endométrio, história familiar de câncer colorretal) não mostra alteração em *MSH2* no oncoplot?

---

## Extra — Ver os alinhamentos no IGV-Web

1. Baixe `dados/referencia/genoma_aula.fa` e `genoma_aula.fa.fai`, e `resultados/04_alinhamento/P01_dedup.bam` e `.bai`.
2. Acesse **https://igv.org/app** → *Genome* → *Local File* → selecione o `.fa` e o `.fai` juntos.
3. *Tracks* → *Local File* → selecione o `.bam` e o `.bai` juntos.
4. Digite `TP53:8400-8560` na busca para ver a variante R175H do P01.

---

## Problemas comuns

| Sintoma | Solução |
|---|---|
| `ls` e os comandos não mostram nada (só as mensagens de erro aparecem) | As linhas ```` ```bash ```` foram coladas junto com os comandos e o terminal ficou preso. Feche o terminal (🗑), abra outro e rode `AMOSTRA=P01` de novo |
| `command not found: fastqc` (ou bwa, samtools...) | O ambiente não está ativo: `micromamba activate aula-variantes` |
| `No such file or directory` com `_R1.fastq.gz` | A variável sumiu (terminal novo). Rode `AMOSTRA=P01` de novo |
| `[E::bwa_idx_load_from_disk] fail to locate the index` | Faltou o `bwa index` da Etapa 4 |
| `bcftools merge`: *Duplicate sample names* | O BAM foi gerado sem `-R` no `bwa mem` |
| VEP: *chromosome not found* / variantes nos lugares errados | Faltou a Etapa 7b (o VCF ainda tem `TP53`, `BRCA1`... na 1ª coluna) |
| `split-vep`: *no such tag defined in the VCF header: INFO/SAMPLE* | `%SAMPLE` precisa estar entre colchetes `[ ]` (veja a Etapa 9) |
| `split-vep`: *The CSQ tag is not present* | Esse VCF não foi anotado pelo VEP — use o arquivo baixado do VEP |
