#!/usr/bin/env bash
# =============================================================================
# Etapa 9 — Converter o VCF anotado pelo VEP em MAF (para o maftools)
#
#  1) bcftools +split-vep extrai os campos do CSQ (anotação do VEP) em tabela,
#     uma linha por amostra portadora da variante (genótipo com alelo alternativo);
#  2) awk traduz os termos de consequência do VEP (Sequence Ontology) para as
#     categorias do maftools (Variant_Classification), classifica o tipo de
#     variante (Variant_Type), ajusta indels ao padrão MAF e gera HGVSp_Short
#     (ex.: p.Arg175His → p.R175H), usado no lollipopPlot.
#
# USO:
#   bash scripts/vcf_para_maf.sh  anotado_pelo_VEP.vcf  saida.maf
#
# No VEP online, marque: "Identifiers and frequency data" (HGVS, frequências do
# gnomAD, significância clínica) e, em "Filtering options", "Restrict results"
# → "Show one selected consequence per variant" (equivalente ao --pick).
# Baixe o resultado no formato VCF.
# =============================================================================
set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Uso: bash scripts/vcf_para_maf.sh  anotado_pelo_VEP.vcf  saida.maf" >&2
    exit 1
fi
ENTRADA="$1"; SAIDA="$2"
[ -s "$ENTRADA" ] || { echo "ERRO: '$ENTRADA' não encontrado ou vazio." >&2; exit 1; }

if ! grep -q 'ID=CSQ' <( (gzip -cd "$ENTRADA" 2>/dev/null || cat "$ENTRADA") | head -500 ); then
    echo "ERRO: '$ENTRADA' não tem a anotação CSQ do VEP. Você enviou o arquivo de saída do VEP em formato VCF?" >&2
    exit 1
fi

# Campos disponíveis no CSQ (variam conforme as opções marcadas no VEP)
CAMPOS=$(bcftools +split-vep -l "$ENTRADA" | cut -f2)
tem() { grep -qx "$1" <<< "$CAMPOS"; }
campo() { if tem "$1"; then echo "%$1"; else echo "."; fi; }

for obrig in SYMBOL Consequence; do
    tem "$obrig" || { echo "ERRO: o campo '$obrig' não está no CSQ do VEP." >&2; exit 1; }
done
FREQ="."
for f in gnomADe_AF gnomAD_AF gnomADg_AF MAX_AF AF; do
    if tem "$f"; then FREQ="%$f"; NOME_FREQ="$f"; break; fi
done
[ "$FREQ" = "." ] && echo "AVISO: nenhuma frequência populacional (gnomAD) no CSQ — a coluna gnomAD_AF ficará vazia." >&2

FORMATO="[%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t$(campo IMPACT)\t$(campo HGVSc)\t$(campo HGVSp)\t$(campo Existing_variation)\t$(campo CLIN_SIG)\t$FREQ\t%SAMPLE\t%GT\t%AD\n]"

bcftools +split-vep -s worst -d -i 'GT="alt"' -f "$FORMATO" "$ENTRADA" 2>/dev/null \
| awk 'BEGIN { FS = OFS = "\t"
    # ordem de gravidade (a primeira encontrada define a classificação)
    n = split("transcript_ablation splice_acceptor_variant splice_donor_variant stop_gained frameshift_variant stop_lost start_lost initiator_codon_variant inframe_insertion inframe_deletion missense_variant protein_altering_variant splice_donor_5th_base_variant splice_region_variant splice_donor_region_variant splice_polypyrimidine_tract_variant incomplete_terminal_codon_variant start_retained_variant stop_retained_variant synonymous_variant coding_sequence_variant 5_prime_UTR_variant 3_prime_UTR_variant non_coding_transcript_exon_variant intron_variant upstream_gene_variant downstream_gene_variant intergenic_variant", ordem, " ")
    split("Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val Ter Sec Pyl Xaa", tres, " ")
    split("A R N D C Q E G H I L K M F P S T W Y V * U O X", um, " ")
    print "Hugo_Symbol","Entrez_Gene_Id","Center","NCBI_Build","Chromosome","Start_Position","End_Position","Strand","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","dbSNP_RS","Tumor_Sample_Barcode","HGVSc","HGVSp","HGVSp_Short","Consequence","IMPACT","CLIN_SIG","gnomAD_AF","t_depth","t_ref_count","t_alt_count","VAF"
}
function classifica(cons, tipo, inframe,    i) {
    for (i = 1; i <= n; i++) if (index("&" cons "&", "&" ordem[i] "&")) break
    t = ordem[i]
    if (t ~ /^(transcript_ablation|splice_acceptor_variant|splice_donor_variant)$/) return "Splice_Site"
    if (t == "stop_gained")        return "Nonsense_Mutation"
    if (t == "frameshift_variant") return (tipo == "DEL") ? "Frame_Shift_Del" : "Frame_Shift_Ins"
    if (t == "stop_lost")          return "Nonstop_Mutation"
    if (t ~ /^(start_lost|initiator_codon_variant)$/) return "Translation_Start_Site"
    if (t == "inframe_insertion")  return "In_Frame_Ins"
    if (t == "inframe_deletion")   return "In_Frame_Del"
    if (t ~ /^(missense_variant|coding_sequence_variant)$/) return "Missense_Mutation"
    if (t == "protein_altering_variant") return inframe ? ((tipo == "DEL") ? "In_Frame_Del" : "In_Frame_Ins") : "Missense_Mutation"
    if (t ~ /^splice_(donor_5th_base|region|donor_region|polypyrimidine_tract)_variant$/) return "Splice_Region"
    if (t ~ /^(synonymous_variant|stop_retained_variant|start_retained_variant|incomplete_terminal_codon_variant)$/) return "Silent"
    if (t == "5_prime_UTR_variant") return "5'"'"'UTR"
    if (t == "3_prime_UTR_variant") return "3'"'"'UTR"
    if (t == "non_coding_transcript_exon_variant") return "RNA"
    if (t == "intron_variant")      return "Intron"
    if (t == "upstream_gene_variant")   return "5'"'"'Flank"
    if (t == "downstream_gene_variant") return "3'"'"'Flank"
    return "IGR"
}
function curto(hgvsp,    p, i) {
    if (hgvsp == "." || hgvsp == "") return "."
    p = hgvsp; sub(/^[^:]*:/, "", p); gsub(/%3D/, "=", p)
    for (i = 1; i <= 24; i++) gsub(tres[i], um[i], p)
    return p
}
{
    crom = $1; pos = $2; ref = $3; alt = $4
    if (alt == "*" || alt == ".") next
    # tipo de variante e coordenadas no padrão MAF
    if (length(ref) == length(alt)) {
        tipo = (length(ref) == 1) ? "SNP" : (length(ref) == 2) ? "DNP" : (length(ref) == 3) ? "TNP" : "ONP"
        ini = pos; fim = pos + length(ref) - 1; r = ref; a = alt
    } else if (length(ref) > length(alt)) {          # deleção: remove a base âncora
        tipo = "DEL"; r = substr(ref, length(alt) + 1); a = "-"
        ini = pos + length(alt); fim = ini + length(r) - 1
    } else {                                         # inserção
        tipo = "INS"; r = "-"; a = substr(alt, length(ref) + 1)
        ini = pos + length(ref) - 1; fim = ini + 1
    }
    inframe = ((length(r) - length(a)) % 3 == 0)
    vc = classifica($6, tipo, inframe)
    split($15, ad, ",")
    dp = ad[1] + ad[2]
    vaf = (dp > 0) ? sprintf("%.2f", ad[2] / dp) : "."
    rs = "."
    if ($10 != ".") { m = split($10, ids, "&"); for (i = 1; i <= m; i++) if (ids[i] ~ /^rs/) { rs = ids[i]; break } }
    sym = ($5 == "." || $5 == "") ? "Unknown" : $5
    print sym, 0, ".", "GRCh38", crom, ini, fim, "+", vc, tipo, r, r, a, rs, $13, $8, $9, curto($9), $6, $7, $11, $12, dp, ad[1], ad[2], vaf
}' > "$SAIDA"

NV=$(( $(wc -l < "$SAIDA") - 1 ))
NA=$(awk -F'\t' 'NR>1{print $15}' "$SAIDA" | sort -u | wc -l)
echo "MAF gerado: $SAIDA  ($NV variantes em $NA amostras)"
[ "${NOME_FREQ:-}" != "" ] && echo "Frequência populacional usada na coluna gnomAD_AF: $NOME_FREQ"
exit 0
