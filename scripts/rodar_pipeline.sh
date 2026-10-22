#!/usr/bin/env bash
# =============================================================================
#  Pipeline de análise de variantes — Bioinformática aplicada ao estudo do câncer
#  IEC – PUC Minas  ·  Prof. Dr. Fábio Ribeiro Queiroz
#
#  Executa automaticamente as Etapas 1 a 7 do roteiro da aula (FastQC →
#  Trimmomatic → FastQC → BWA → SAMtools → Picard → BCFtools), restaura as
#  coordenadas GRCh38 e junta todas as amostras em um único VCF pronto para
#  ser enviado ao VEP online (Etapa 8).
#
#  Cada comando é exibido na tela ANTES de ser executado — é exatamente o que
#  você digitaria no terminal seguindo o roteiro.
#
#  USO:
#     bash scripts/rodar_pipeline.sh              # todas as amostras
#     bash scripts/rodar_pipeline.sh P01 P04      # só as amostras indicadas
# =============================================================================
set -euo pipefail

DIR_REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$DIR_REPO"

REF=dados/referencia/genoma_aula.fa
ALVOS=dados/referencia/alvos_painel.bed
ADAPT=dados/adapters.fasta
OUT=resultados
THREADS=$(( $(nproc 2>/dev/null || echo 2) > 4 ? 4 : $(nproc 2>/dev/null || echo 2) ))

AZUL='\033[1;34m'; VERDE='\033[0;32m'; AMARELO='\033[1;33m'; VERMELHO='\033[1;31m'; SEM='\033[0m'
etapa()  { echo -e "\n${AZUL}━━━ $* ━━━${SEM}"; }
LOG=/dev/null
run() {
    echo -e "${VERDE}\$ $1${SEM}"
    if ! bash -c "$1" 2>>"$LOG"; then
        echo -e "${VERMELHO}ERRO ao executar o comando acima. Últimas linhas do log ($LOG):${SEM}"
        tail -n 15 "$LOG"
        exit 1
    fi
}

# ---------------------------------------------------------------- verificações
for prog in fastqc trimmomatic bwa samtools picard bcftools bgzip; do
    if ! command -v "$prog" >/dev/null 2>&1; then
        echo -e "${VERMELHO}ERRO: '$prog' não encontrado.${SEM}"
        echo "O ambiente da aula não está ativo. Rode:  micromamba activate aula-variantes"
        echo "(ou veja a seção 'Instalação' do README.md)"
        exit 1
    fi
done

if [ $# -gt 0 ]; then
    AMOSTRAS=("$@")
else
    AMOSTRAS=()
    for f in dados/fastq/*_R1.fastq.gz; do AMOSTRAS+=("$(basename "$f" _R1.fastq.gz)"); done
fi
for a in "${AMOSTRAS[@]}"; do
    [ -s "dados/fastq/${a}_R1.fastq.gz" ] || { echo "ERRO: amostra '$a' não existe em dados/fastq/"; exit 1; }
done

mkdir -p $OUT/logs $OUT/01_fastqc_bruto $OUT/02_trimagem $OUT/03_fastqc_trimado \
         $OUT/04_alinhamento $OUT/07_variantes $OUT/08_para_o_VEP

echo -e "${AMARELO}Amostras: ${AMOSTRAS[*]}   ·   threads: $THREADS${SEM}"

# ---------------------------------------------------- indexação da referência
if [ ! -s "$REF.bwt" ]; then
    etapa "Etapa 4 (preparo) — Indexar o genoma de referência da aula (só uma vez)"
    LOG=$OUT/logs/indexacao.log
    run "bwa index $REF"
    run "samtools faidx $REF"
fi

RESUMO=$OUT/resumo_qc.tsv
[ -s "$RESUMO" ] || echo -e "amostra\tpares_brutos\tpares_apos_trimagem\tpct_aproveitamento\tpct_duplicatas\tprofundidade_media_alvos\tvariantes_brutas\tvariantes_filtradas" > "$RESUMO"

for A in "${AMOSTRAS[@]}"; do
    LOG=$OUT/logs/$A.log; : > "$LOG"
    R1=dados/fastq/${A}_R1.fastq.gz
    R2=dados/fastq/${A}_R2.fastq.gz
    T=$OUT/02_trimagem/$A
    B=$OUT/04_alinhamento/$A
    V=$OUT/07_variantes/$A

    echo -e "\n${AMARELO}████████████████████  AMOSTRA $A  ████████████████████${SEM}"

    etapa "Etapa 1 — Qualidade dos dados brutos (FastQC)"
    run "fastqc -q -t 2 --outdir $OUT/01_fastqc_bruto $R1 $R2"

    etapa "Etapa 2 — Remoção de adaptadores e bases de baixa qualidade (Trimmomatic)"
    run "trimmomatic PE -threads $THREADS -phred33 $R1 $R2 ${T}_R1_paired.fastq.gz ${T}_R1_unpaired.fastq.gz ${T}_R2_paired.fastq.gz ${T}_R2_unpaired.fastq.gz ILLUMINACLIP:$ADAPT:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
    grep "Input Read Pairs" "$LOG" | tail -1

    etapa "Etapa 3 — Qualidade após a trimagem (FastQC)"
    run "fastqc -q -t 2 --outdir $OUT/03_fastqc_trimado ${T}_R1_paired.fastq.gz ${T}_R2_paired.fastq.gz"

    etapa "Etapa 4 — Alinhamento ao genoma de referência (BWA-MEM)"
    run "bwa mem -t $THREADS -R '@RG\tID:$A\tSM:$A\tPL:ILLUMINA' $REF ${T}_R1_paired.fastq.gz ${T}_R2_paired.fastq.gz > $B.sam"

    etapa "Etapa 5 — Conversão SAM → BAM (SAMtools)"
    run "samtools view -S -b $B.sam > $B.bam"

    etapa "Etapa 6 — Ordenar, indexar e marcar duplicatas (SAMtools + Picard)"
    run "samtools sort $B.bam -o ${B}_sorted.bam"
    run "samtools index ${B}_sorted.bam"
    run "picard MarkDuplicates I=${B}_sorted.bam O=${B}_dedup.bam M=${B}_dup_metrics.txt"
    run "samtools index ${B}_dedup.bam"

    etapa "Etapa 7 — Chamada e filtragem de variantes (BCFtools)"
    run "bcftools mpileup -Ou -f $REF ${B}_dedup.bam | bcftools call -mv -Ov -o ${V}_variants.vcf"
    run "bcftools filter -e 'QUAL<30 || DP<10' -Ov -o ${V}_filtered_variants.vcf ${V}_variants.vcf"

    etapa "Etapa 7b — Restaurar coordenadas GRCh38 (necessário por causa do genoma reduzido)"
    run "bash scripts/restaurar_coordenadas.sh ${V}_filtered_variants.vcf ${V}_GRCh38.vcf"

    # ---- resumo de QC da amostra
    brutos=$(grep "Input Read Pairs" "$LOG" | tail -1 | awk '{print $4}')
    ambos=$(grep "Input Read Pairs" "$LOG" | tail -1 | awk '{print $7}')
    pct=$(awk -v a="$ambos" -v b="$brutos" 'BEGIN{printf "%.1f", 100*a/b}')
    dup=$(awk -F'\t' '/^LIBRARY/{getline; printf "%.1f", 100*$9}' ${B}_dup_metrics.txt)
    prof=$(samtools depth -a -b $ALVOS ${B}_dedup.bam 2>/dev/null | awk '{s+=$3;n++} END{printf "%.0f", s/n}')
    nbr=$(grep -vc '^#' ${V}_variants.vcf || true)
    nfi=$(grep -vc '^#' ${V}_filtered_variants.vcf || true)
    awk -v a="$A" 'BEGIN{FS=OFS="\t"} $1!=a' "$RESUMO" > "$RESUMO.tmp" && mv "$RESUMO.tmp" "$RESUMO"
    echo -e "$A\t$brutos\t$ambos\t$pct\t$dup\t$prof\t$nbr\t$nfi" >> "$RESUMO"
    echo -e "${AMARELO}$A: $ambos/$brutos pares mantidos ($pct%) · duplicatas $dup% · profundidade média ${prof}x · $nfi variantes após o filtro${SEM}"
done

# ------------------------------------------------- juntar amostras para o VEP
etapa "Etapa 8 (preparo) — Juntar todas as amostras em um único VCF para o VEP online"
LOG=$OUT/logs/juntar.log; : > "$LOG"
LISTA=()
for f in $OUT/07_variantes/*_GRCh38.vcf; do
    run "bgzip -f -k $f && bcftools index -f $f.gz"
    LISTA+=("$f.gz")
done
if [ ${#LISTA[@]} -gt 1 ]; then
    run "bcftools merge -Ov -o $OUT/08_para_o_VEP/coorte_GRCh38.vcf ${LISTA[*]}"
else
    run "bcftools view -Ov -o $OUT/08_para_o_VEP/coorte_GRCh38.vcf ${LISTA[0]}"
fi

echo -e "\n${AZUL}━━━ Concluído ━━━${SEM}"
echo "Resumo de QC por amostra:            $RESUMO"
echo "VCF de cada amostra (GRCh38):        $OUT/07_variantes/<amostra>_GRCh38.vcf"
echo -e "${AMARELO}Arquivo para enviar ao VEP online:   $OUT/08_para_o_VEP/coorte_GRCh38.vcf${SEM}"
echo "   → https://www.ensembl.org/Tools/VEP  (veja o passo a passo no GUIA_DA_AULA.md, Etapa 8)"
