#!/usr/bin/env bash
# =============================================================================
# Restaura as coordenadas reais do GRCh38 em um VCF chamado contra o genoma
# reduzido da aula (dados/referencia/genoma_aula.fa).
#
# POR QUE ISSO É NECESSÁRIO?
#   O genoma da aula contém apenas trechos do GRCh38 (um "contig" por gene).
#   Por isso, no VCF as variantes aparecem como, por exemplo,  TP53  8480 —
#   ou seja, posição 8480 DENTRO do recorte do TP53. O VEP (e qualquer banco
#   de dados) só entende coordenadas cromossômicas reais:  chr17  7675088.
#   Este script soma o deslocamento de cada recorte (tabela regioes.tsv) e
#   troca o nome do contig pelo cromossomo verdadeiro.
#
#   Antes disso, as variantes são normalizadas (bcftools norm): indels são
#   alinhados à esquerda e alelos multialélicos são separados — boa prática
#   antes de qualquer anotação.
#
# USO:
#   bash scripts/restaurar_coordenadas.sh  entrada.vcf  saida_GRCh38.vcf
# =============================================================================
set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Uso: bash scripts/restaurar_coordenadas.sh  entrada.vcf  saida_GRCh38.vcf" >&2
    exit 1
fi

ENTRADA="$1"
SAIDA="$2"
DIR_REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF="$DIR_REPO/dados/referencia/genoma_aula.fa"
REGIOES="$DIR_REPO/dados/referencia/regioes.tsv"

[ -s "$ENTRADA" ] || { echo "ERRO: arquivo '$ENTRADA' não encontrado ou vazio." >&2; exit 1; }
[ -s "$REF.fai" ] || samtools faidx "$REF"

bcftools norm -f "$REF" -m -any -Ov "$ENTRADA" 2>/dev/null \
| awk -v OFS='\t' -v regioes="$REGIOES" '
    BEGIN {
        while ((getline linha < regioes) > 0) {
            if (linha ~ /^#/) continue
            split(linha, c, "\t")
            crom[c[1]] = c[2]; desloc[c[1]] = c[3] - 1; tam[c[2]] = c[5]
            if (!(c[2] in visto)) { visto[c[2]] = 1; ordem[++n] = c[2] }
        }
    }
    /^##contig=/   { next }                       # contigs do recorte saem...
    /^##reference=/ { next }
    /^#CHROM/ {                                    # ...e entram os cromossomos reais
        print "##reference=GRCh38"
        for (i = 1; i <= n; i++)
            print "##contig=<ID=" ordem[i] ",length=" tam[ordem[i]] ",assembly=GRCh38>"
        print "##restaurar_coordenadas=coordenadas convertidas do genoma_aula.fa para GRCh38"
        print; next
    }
    /^#/ { print; next }
    {
        if (!($1 in crom)) { print "ERRO: contig desconhecido: " $1 > "/dev/stderr"; exit 1 }
        $2 = $2 + desloc[$1]
        $1 = crom[$1]
        print
    }' \
| bcftools sort -Ov -o "$SAIDA" 2>/dev/null

echo "Coordenadas restauradas: $SAIDA ($(grep -vc '^#' "$SAIDA") variantes)"
