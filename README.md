# Pipeline para análise de variantes genéticas — aula prática

**Bioinformática aplicada ao estudo do câncer · Pós-graduação em Bioinformática aplicada à saúde · IEC – PUC Minas**
**Prof. Dr. Fábio Ribeiro Queiroz**

Material da aula prática de identificação e anotação de variantes em genes de câncer: do FASTQ ao oncoplot, com FastQC, Trimmomatic, BWA, SAMtools, Picard, BCFtools, VEP e maftools.

👉 **O passo a passo da aula está em [GUIA_DA_AULA.md](GUIA_DA_AULA.md).**

---

## Como abrir (escolha uma opção)

### Opção 1 — No navegador, sem instalar nada (recomendado)

[![Abrir no GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/fabiorqueiroz/aula-variantes-cancer?quickstart=1)

1. Faça login no GitHub (crie uma conta gratuita em https://github.com/signup, se não tiver).
2. Clique no botão acima → **Create codespace**.
3. Aguarde a instalação automática das ferramentas (3 a 8 minutos na primeira vez — acompanhe pelo terminal). Ao final aparece **"Ambiente pronto!"**.
4. Abra um terminal novo (menu ☰ → *Terminal* → *New Terminal*) e confira:
   ```bash
   bash scripts/verificar_ambiente.sh
   ```

O Codespaces é gratuito até 60 horas/mês em contas pessoais. Ao terminar a aula, pare ou apague o codespace em https://github.com/codespaces.

### Opção 2 — No seu computador (Linux, macOS ou Windows com WSL)

```bash
git clone https://github.com/fabiorqueiroz/aula-variantes-cancer.git
cd aula-variantes-cancer
bash instalar.sh          # instala tudo em um ambiente isolado (não mexe no seu conda)
source ~/.bashrc          # ou abra um terminal novo
bash scripts/verificar_ambiente.sh
```

No Windows, use o terminal do **WSL (Ubuntu)** — instruções em https://learn.microsoft.com/windows/wsl/install.

### Atalho — rodar o pipeline inteiro automaticamente

```bash
bash scripts/rodar_pipeline.sh
```

Executa as Etapas 1 a 7 para os 12 pacientes (mostrando cada comando na tela) e gera o VCF pronto para o VEP em `resultados/08_para_o_VEP/coorte_GRCh38.vcf`.

---

## O que tem aqui

```
├── GUIA_DA_AULA.md                 ← roteiro passo a passo
├── dados/
│   ├── fastq/                      12 pacientes simulados (P01…P12), leituras pareadas 2×150 pb
│   ├── referencia/
│   │   ├── genoma_aula.fa          genoma reduzido: regiões de 13 genes do GRCh38 (~1,2 Mb)
│   │   ├── regioes.tsv             onde cada recorte fica no GRCh38 (usado na Etapa 7b)
│   │   └── alvos_painel.bed        regiões capturadas pelo painel (éxons codificantes ± 50 pb)
│   ├── adapters.fasta              adaptadores Illumina TruSeq (Trimmomatic)
│   ├── dados_clinicos.tsv          dados clínicos simulados (para o maftools)
│   └── plano_b_VEP/                VCF já anotado pelo VEP, caso o site esteja fora do ar
├── scripts/
│   ├── rodar_pipeline.sh           pipeline completo (Etapas 1–7b), todas as amostras
│   ├── restaurar_coordenadas.sh    Etapa 7b: genoma reduzido → coordenadas GRCh38
│   ├── vcf_para_maf.sh             Etapa 9: VCF anotado pelo VEP → MAF
│   ├── analise_maftools.R          Etapa 10: oncoplot, Ti/Tv, lollipop...
│   └── verificar_ambiente.sh       confere se as ferramentas estão instaladas
├── environment.yml                 versões das ferramentas
└── instalar.sh                     instalação local (micromamba)
```

## Sobre os dados

- **Pacientes simulados.** Nenhum dado é de paciente real. As leituras foram geradas por simulação a partir da sequência de referência GRCh38, imitando um painel de captura sequenciado em Illumina (qualidade binada estilo NovaSeq, adaptadores em fragmentos curtos, caudas de baixa qualidade e ~18% de duplicatas de PCR).
- **Variantes reais.** As variantes inseridas existem em bancos públicos (dbSNP, ClinVar, COSMIC) e foram conferidas no Ensembl VEP: variantes patogênicas germinativas em genes de alta penetrância, mutações somáticas *driver* e polimorfismos comuns. Como as coordenadas são as do GRCh38, o VCF da aula pode ser anotado no VEP online normalmente.
- **Genoma reduzido.** Cada gene é um contig (`>TP53 chr17:7666609-7679594 GRCh38`). A Etapa 7b converte as posições de volta para o cromossomo real.
- **Dados clínicos.** Idade, sexo, diagnóstico, história familiar e sobrevida são fictícios e coerentes com cada caso, apenas para fins didáticos.

Ferramentas: FastQC 0.13 · Trimmomatic 0.40 · BWA 0.7.19 · SAMtools/BCFtools/HTSlib 1.24 · Picard 3.5 · R 4.4 + maftools 2.22.
