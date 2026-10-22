#!/usr/bin/env Rscript
# =============================================================================
# Etapa 10 — Análise do perfil das variantes com maftools (R)
#
# USO (no terminal):
#   Rscript scripts/analise_maftools.R  resultados/09_maf/coorte.maf
#
# Ou abra o R (digite R no terminal) e rode os blocos abaixo um a um.
# As figuras são salvas em resultados/10_maftools/ (abra os .png pelo
# explorador de arquivos à esquerda).
# =============================================================================
suppressPackageStartupMessages({
  library(maftools)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
arquivo_maf <- if (length(args) >= 1) args[1] else "resultados/09_maf/coorte.maf"
dados_clin  <- if (length(args) >= 2) args[2] else "dados/dados_clinicos.tsv"
saida <- "resultados/10_maftools"
dir.create(saida, recursive = TRUE, showWarnings = FALSE)
figura <- function(nome, largura = 1800, altura = 1100) {
  png(file.path(saida, nome), width = largura, height = altura, res = 150)
}

# ---------------------------------------------------------------------------
# 1) Carregar o MAF (dados clínicos são opcionais, mas enriquecem a análise)
# ---------------------------------------------------------------------------
maf <- read.maf(maf = arquivo_maf, clinicalData = dados_clin, verbose = FALSE)

cat("\n========== Resumo geral ==========\n");  print(maf)
cat("\n========== Por amostra ==========\n");   print(getSampleSummary(maf))
cat("\n========== Por gene ==========\n");      print(getGeneSummary(maf))

figura("01_resumo_todas_variantes.png")
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = "median", dashboard = TRUE)
dev.off()

figura("02_oncoplot_todas_variantes.png")
oncoplot(maf = maf, top = 50, clinicalFeatures = "Tipo_tumoral", sortByAnnotation = TRUE)
dev.off()

# ---------------------------------------------------------------------------
# 2) Separar polimorfismos comuns das variantes raras
#    Polimorfismos (frequência populacional >= 1% no gnomAD) não explicam
#    risco de câncer de alta penetrância nem são mutações "driver".
# ---------------------------------------------------------------------------
tabela <- read.delim(arquivo_maf, stringsAsFactors = FALSE, na.strings = c(".", "", "NA"))
tabela$gnomAD_AF <- suppressWarnings(as.numeric(sub("&.*", "", tabela$gnomAD_AF)))

cat("\n========== Variantes comuns (gnomAD_AF >= 0.01) ==========\n")
print(tabela %>% filter(!is.na(gnomAD_AF) & gnomAD_AF >= 0.01) %>%
        count(Hugo_Symbol, HGVSp_Short, gnomAD_AF, name = "n_amostras") %>%
        arrange(Hugo_Symbol))

raras <- tabela %>% filter(is.na(gnomAD_AF) | gnomAD_AF < 0.01)
cat("\n========== Variantes raras (gnomAD_AF < 0.01 ou ausentes do gnomAD) ==========\n")
print(raras %>% select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification,
                       HGVSp_Short, HGVSc, CLIN_SIG, gnomAD_AF, VAF) %>%
        arrange(Tumor_Sample_Barcode))

maf_raras <- read.maf(maf = raras, clinicalData = dados_clin, verbose = FALSE)

figura("03_resumo_variantes_raras.png")
plotmafSummary(maf = maf_raras, rmOutlier = TRUE, addStat = "median", dashboard = TRUE)
dev.off()

figura("04_oncoplot_variantes_raras.png")
oncoplot(maf = maf_raras, top = 50, clinicalFeatures = "Tipo_tumoral", sortByAnnotation = TRUE)
dev.off()

# ---------------------------------------------------------------------------
# 3) Ocultar/tratar "Multi_Hit" no oncoplot
#    Multi_Hit = mais de uma mutação no mesmo gene e na mesma amostra.
#    Opção 2 do roteiro: manter só a mutação mais deletéria por gene/amostra.
# ---------------------------------------------------------------------------
maf_dedup <- maf_raras@data %>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  slice_min(order_by = match(Variant_Classification,
              c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
                "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation",
                "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Silent")),
            n = 1, with_ties = FALSE) %>%
  ungroup()
maf_sem_multihit <- read.maf(maf = maf_dedup, clinicalData = dados_clin, verbose = FALSE)

figura("05_oncoplot_sem_multihit.png")
oncoplot(maf = maf_sem_multihit, top = 50, clinicalFeatures = "Tipo_tumoral", sortByAnnotation = TRUE)
dev.off()

# ---------------------------------------------------------------------------
# 4) Transições/transversões (Ti/Tv) — controle de qualidade
# ---------------------------------------------------------------------------
titv_res <- titv(maf = maf, plot = FALSE, useSyn = TRUE)
figura("06_titv.png")
plotTiTv(res = titv_res)
dev.off()

# ---------------------------------------------------------------------------
# 5) Distribuição das mutações ao longo do TP53 (lollipop)
# ---------------------------------------------------------------------------
figura("07_lollipop_TP53.png", 1800, 900)
lollipopPlot(maf = maf_raras, gene = "TP53", AACol = "HGVSp_Short",
             showMutationRate = TRUE, labelPos = "all")
dev.off()

# ---------------------------------------------------------------------------
# 6) (Opcional) Associação com sobrevida — dados clínicos SIMULADOS
#    Com 12 pacientes, serve apenas para mostrar o comando.
# ---------------------------------------------------------------------------
# mafSurvival(maf = maf_raras, genes = "TP53", time = "tempo_meses",
#             Status = "status", isTCGA = FALSE)

cat("\nFiguras salvas em", saida, ":\n")
print(list.files(saida))
