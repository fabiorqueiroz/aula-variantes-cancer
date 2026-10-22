# Guia da aula — Pipeline para análise de variantes genéticas

**Bioinformática aplicada ao estudo do câncer · IEC – PUC Minas · Prof. Dr. Fábio Ribeiro Queiroz**

Este guia segue as mesmas etapas do roteiro da aula, já com os caminhos dos arquivos deste repositório. Copie e cole os comandos no terminal, **na ordem**, um bloco de cada vez. Todos os alunos processam os **12 pacientes**, para que o maftools tenha uma coorte inteira para analisar.

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

Crie as pastas de resultados:

```bash
mkdir -p resultados/01_fastqc_bruto resultados/02_trimagem resultados/03_fastqc_trimado \
         resultados/04_alinhamento resultados/07_variantes resultados/08_para_o_VEP resultados/09_maf
```

**Todos processam os 12 pacientes.** A lista das amostras está no arquivo `dados/amostras.txt`:

```bash
cat dados/amostras.txt
```

> **Como funciona o `for`.** Para não digitar o mesmo comando 12 vezes, usamos um laço:
> `for A in $(cat dados/amostras.txt); do ... ; done`
> O bash lê a lista e repete o que está entre `do` e `done` para cada amostra. Em cada volta, `$A` vale `P01`, depois `P02`, e assim por diante. Por isso o nome dos arquivos aparece como `${A}_R1.fastq.gz`. Cole o bloco inteiro e espere o `$` voltar.

---

## 1ª Etapa — Análise de qualidade dos dados brutos (FastQC)

O FastQC aceita vários arquivos de uma vez (`*` quer dizer "todos"); `-t 2` analisa dois por vez:

```bash
fastqc -t 2 --outdir resultados/01_fastqc_bruto dados/fastq/*.fastq.gz
```

Abra os relatórios `resultados/01_fastqc_bruto/*_fastqc.html` (no Codespaces: botão direito no arquivo → **Show Preview**; ou → **Download** e abra no seu navegador).

**O que observar:** compare R1 e R2 em *Per base sequence quality* (a qualidade cai no fim das leituras, sobretudo no R2) e em *Adapter Content* (adaptador Illumina aparecendo no fim das leituras — sinal de fragmentos de DNA mais curtos que 150 pb, comum em amostras tumorais FFPE).

---

## 2ª Etapa — Remoção de adaptadores e bases de baixa qualidade (Trimmomatic)

No ambiente da aula, `java -jar trimmomatic.jar` é simplesmente `trimmomatic`. O relatório de cada amostra vai para um arquivo `_trimmomatic.log` (~1 minuto no total):

```bash
for A in $(cat dados/amostras.txt); do
  echo "== Trimando $A"
  trimmomatic PE -phred33 \
    dados/fastq/${A}_R1.fastq.gz dados/fastq/${A}_R2.fastq.gz \
    resultados/02_trimagem/${A}_R1_paired.fastq.gz resultados/02_trimagem/${A}_R1_unpaired.fastq.gz \
    resultados/02_trimagem/${A}_R2_paired.fastq.gz resultados/02_trimagem/${A}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:dados/adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
    2> resultados/02_trimagem/${A}_trimmomatic.log
  grep "Input Read Pairs" resultados/02_trimagem/${A}_trimmomatic.log
done
```

> **Curiosidade:** por que muitos pares aparecem como *Forward Only Surviving*? Quando o fragmento é menor que a leitura, R1 e R2 leem exatamente o mesmo trecho de DNA (em fitas opostas). No modo *palíndromo* do `ILLUMINACLIP`, o Trimmomatic remove o adaptador e descarta o R2, que seria redundante.

---

## 3ª Etapa — Qualidade após a trimagem

```bash
fastqc -t 2 --outdir resultados/03_fastqc_trimado resultados/02_trimagem/*_paired.fastq.gz
```

**Tabela de dados brutos × dados trimados**, montada a partir dos relatórios do Trimmomatic:

```bash
echo -e "amostra\tpares_brutos\tpares_mantidos\tpct_aproveitamento" > resultados/02_trimagem/tabela_aproveitamento.tsv
for A in $(cat dados/amostras.txt); do
  grep "Input Read Pairs" resultados/02_trimagem/${A}_trimmomatic.log \
    | awk -v a=$A '{printf "%s\t%s\t%s\t%.1f\n", a, $4, $7, 100*$7/$4}'
done >> resultados/02_trimagem/tabela_aproveitamento.tsv

cat resultados/02_trimagem/tabela_aproveitamento.tsv
```

---

## 4ª Etapa — Alinhamento ao genoma de referência (BWA)

**Comando 1** — indexar o genoma (uma vez só; leva ~1 segundo):

```bash
bwa index dados/referencia/genoma_aula.fa
samtools faidx dados/referencia/genoma_aula.fa
```

**Comando 2** — alinhar as 12 amostras (~1 minuto). Acrescentamos `-R` (*read group*), que grava o nome da amostra no alinhamento. Sem ele, o VCF sairia com o nome do arquivo no lugar do nome da amostra — e não conseguiríamos juntar os 12 pacientes depois.

```bash
for A in $(cat dados/amostras.txt); do
  echo "== Alinhando $A"
  bwa mem -t 2 -R "@RG\tID:${A}\tSM:${A}\tPL:ILLUMINA" dados/referencia/genoma_aula.fa \
    resultados/02_trimagem/${A}_R1_paired.fastq.gz resultados/02_trimagem/${A}_R2_paired.fastq.gz \
    > resultados/04_alinhamento/${A}.sam 2> resultados/04_alinhamento/${A}_bwa.log
done
```

Dê uma olhada no SAM de um paciente (cabeçalho e primeiras leituras):

```bash
head -n 20 resultados/04_alinhamento/P01.sam | cut -c 1-200
```

---

## Etapa 5 — Converter SAM em BAM

```bash
for A in $(cat dados/amostras.txt); do
  echo "== Convertendo $A"
  samtools view -S -b resultados/04_alinhamento/${A}.sam > resultados/04_alinhamento/${A}.bam
done

ls -lh resultados/04_alinhamento/P01.sam resultados/04_alinhamento/P01.bam
```

Compare o tamanho do `.sam` e do `.bam`.

---

## Etapa 6 — Ordenar, indexar e marcar duplicatas

```bash
for A in $(cat dados/amostras.txt); do
  echo "== Ordenando e indexando $A"
  samtools sort resultados/04_alinhamento/${A}.bam -o resultados/04_alinhamento/${A}_sorted.bam
  samtools index resultados/04_alinhamento/${A}_sorted.bam
done
```

Marcar duplicatas de PCR (recomendado para dados de câncer; ~1 minuto). As mensagens do Picard vão para um arquivo `_picard.log`:

```bash
for A in $(cat dados/amostras.txt); do
  echo "== Marcando duplicatas em $A"
  picard MarkDuplicates \
    I=resultados/04_alinhamento/${A}_sorted.bam \
    O=resultados/04_alinhamento/${A}_dedup.bam \
    M=resultados/04_alinhamento/${A}_dup_metrics.txt \
    2> resultados/04_alinhamento/${A}_picard.log
  samtools index resultados/04_alinhamento/${A}_dedup.bam
done
```

Qual a porcentagem de duplicatas de cada paciente? (coluna `PERCENT_DUPLICATION` das métricas do Picard)

```bash
for A in $(cat dados/amostras.txt); do
  echo -e "$A\t$(grep -A1 '^LIBRARY' resultados/04_alinhamento/${A}_dup_metrics.txt | tail -1 | cut -f9)"
done

samtools flagstat resultados/04_alinhamento/P01_dedup.bam
```

---

## Etapa 7 — Chamada e filtragem de variantes (BCFtools)

```bash
for A in $(cat dados/amostras.txt); do
  echo "== Chamando variantes de $A"
  bcftools mpileup -Ou -f dados/referencia/genoma_aula.fa resultados/04_alinhamento/${A}_dedup.bam 2> /dev/null \
    | bcftools call -mv -Ov -o resultados/07_variantes/${A}_variants.vcf 2> /dev/null
  bcftools filter -e 'QUAL<30 || DP<10' -Ov -o resultados/07_variantes/${A}_filtered_variants.vcf \
    resultados/07_variantes/${A}_variants.vcf
done
```

Quantas variantes cada paciente tem antes e depois do filtro?

```bash
echo -e "amostra\tbrutas\tfiltradas"
for A in $(cat dados/amostras.txt); do
  echo -e "$A\t$(grep -vc '^#' resultados/07_variantes/${A}_variants.vcf)\t$(grep -vc '^#' resultados/07_variantes/${A}_filtered_variants.vcf)"
done
```

E onde elas estão? Veja as variantes de um paciente:

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT\t%AD]\n' resultados/07_variantes/P01_filtered_variants.vcf
```

Repare na primeira coluna: `TP53`, `BRCA1`... e posições pequenas. São coordenadas **dentro do recorte** de cada gene, não do cromossomo.

---

## Etapa 7b — Restaurar as coordenadas reais do GRCh38

*(Etapa extra, necessária apenas porque usamos o genoma reduzido.)*

```bash
for A in $(cat dados/amostras.txt); do
  bash scripts/restaurar_coordenadas.sh resultados/07_variantes/${A}_filtered_variants.vcf \
    resultados/07_variantes/${A}_GRCh38.vcf
done

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD]\n' resultados/07_variantes/P01_GRCh38.vcf
```

Agora `TP53  8480` virou `chr17  7675088`. O script soma o início de cada recorte (tabela `dados/referencia/regioes.tsv`) e também normaliza os indels (`bcftools norm`), como se faz antes de qualquer anotação.

---

## Etapa 7c — Juntar os 12 pacientes em um único VCF

O VEP online recebe um arquivo só. O `bcftools merge` junta os VCFs lado a lado: uma linha por variante, uma coluna por paciente. Antes, cada VCF precisa ser compactado (`bgzip`) e indexado:

```bash
for A in $(cat dados/amostras.txt); do
  bgzip -f -k resultados/07_variantes/${A}_GRCh38.vcf
  bcftools index -f resultados/07_variantes/${A}_GRCh38.vcf.gz
done

bcftools merge -Ov -o resultados/08_para_o_VEP/coorte_GRCh38.vcf resultados/07_variantes/*_GRCh38.vcf.gz

grep -vc "^#" resultados/08_para_o_VEP/coorte_GRCh38.vcf
```

**`resultados/08_para_o_VEP/coorte_GRCh38.vcf`** é o arquivo que vai para o VEP.

> **Se você se perdeu em alguma etapa:** o script abaixo refaz as Etapas 1 a 7c para os 12 pacientes (~4 minutos), mostrando cada comando, e deixa tudo pronto para a Etapa 8. Ele também gera `resultados/resumo_qc.tsv`, com aproveitamento, duplicatas e profundidade de todas as amostras.
> ```bash
> bash scripts/rodar_pipeline.sh
> ```

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
| `ls` e os comandos não mostram nada (só as mensagens de erro aparecem) | As linhas ```` ```bash ```` foram coladas junto com os comandos e o terminal ficou preso. Feche o terminal (🗑) e abra outro |
| `command not found: fastqc` (ou bwa, samtools...) | O ambiente não está ativo: `micromamba activate aula-variantes` |
| `No such file or directory` dentro de um `for` | Uma etapa anterior não terminou para todas as amostras. Volte à etapa anterior ou rode `bash scripts/rodar_pipeline.sh` |
| `cat: dados/amostras.txt: No such file or directory` | O terminal não está na pasta do projeto. Rode `cd /workspaces/aula-variantes-cancer` |
| `[E::bwa_idx_load_from_disk] fail to locate the index` | Faltou o `bwa index` da Etapa 4 |
| `bcftools merge`: *Duplicate sample names* | O BAM foi gerado sem `-R` no `bwa mem` |
| VEP: *chromosome not found* / variantes nos lugares errados | Faltou a Etapa 7b (o VCF ainda tem `TP53`, `BRCA1`... na 1ª coluna) |
| `split-vep`: *no such tag defined in the VCF header: INFO/SAMPLE* | `%SAMPLE` precisa estar entre colchetes `[ ]` (veja a Etapa 9) |
| `split-vep`: *The CSQ tag is not present* | Esse VCF não foi anotado pelo VEP — use o arquivo baixado do VEP |
