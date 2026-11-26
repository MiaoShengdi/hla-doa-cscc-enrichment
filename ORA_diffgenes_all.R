#######################################################################
# Over-representation analysis of diff genes
# Collections: GO + KEGG + Hallmark + C7 + C2-Reactome
# ---------------------------------------------------------------
# Purpose: perform hypergeometric over-representation analysis (ORA) on
#          sigDiff.txt, which contains significant differentially regulated
#          genes obtained after network-level gain-of-function perturbation
#          of HLA-DOA. This script identifies enriched pathways but does not
#          determine direction. For direction associated with endogenous
#          HLA-DOA expression, use GSVA_direction_from_ORA.R.
# Gene sets and backgrounds:
#   - GO (BP/CC/MF): enrichGO; all annotated genes in org.Hs.eg.db
#   - KEGG: enrichKEGG; all KEGG-annotated genes
#   - H / C7 / C2-Reactome: enricher; union of the selected MSigDB gene sets
# Outputs (ORA_out/): ORA_GO.txt / ORA_KEGG.txt / ORA_H.txt /
#                     ORA_C7.txt / ORA_C2Reactome.txt, bar plots, and dot plots
# Filtering: nominal P < 0.05, with no restriction on p.adjust
#
# [Revisions, 2026-11-02]
#   - Added tryCatch handling for enrichKEGG because failures of the online
#     KEGG API can otherwise interrupt the script.
#   - Added NA filtering during KEGG ENTREZID-to-SYMBOL conversion to prevent
#     "NA" strings from appearing in geneID.
#   - All other logic, column names, and output filenames remain unchanged;
#     the GSVA script depends on this output format.
# Save this file with UTF-8 encoding
# (RStudio: File -> Save with Encoding -> UTF-8).
#######################################################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(msigdbr)
library(dplyr)
# ===================== 0. Settings =====================
sigDiffFile    <- "sigDiff.txt"         # Input: table of significant diff genes; first column contains gene symbols
outDir         <- "ORA_out"             # Output directory
species        <- "Homo sapiens"        # Species: Homo sapiens / Mus musculus
pvalueFilter   <- 0.05                  # Nominal P-value threshold
p.adjustFilter <- 1                     # Adjusted P-value threshold; 1 means unrestricted, so filtering uses nominal P only
minGSSize      <- 5                     # Minimum gene-set size for enricher; smaller sets are excluded
topN           <- 30                    # Number of pathways displayed in plots
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
colorSel <- if (p.adjustFilter > 0.05) "pvalue" else "p.adjust"
# ===================== 1. Read diff genes =====================
rt <- read.table(sigDiffFile, header = TRUE, sep = "\t", check.names = FALSE)
genes <- unique(as.vector(rt[, 1]))
cat("输入显著基因数:", length(genes), "\n")
# Convert gene symbols to ENTREZID, as required by enrichGO and enrichKEGG
entrez <- mapIds(org.Hs.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL")
geneEntrez <- as.character(entrez[!is.na(entrez)])
cat("转换到 ENTREZID 的基因数:", length(geneEntrez), "\n")
# ===================== 2. GO and KEGG ORA using whole-genome backgrounds =====================
# ---- GO (BP/CC/MF) ----
kkGO <- enrichGO(gene = geneEntrez, OrgDb = org.Hs.eg.db,
                 pvalueCutoff = 1, qvalueCutoff = 1, ont = "all", readable = TRUE)
GO <- as.data.frame(kkGO)
GO <- GO[GO$pvalue < pvalueFilter & GO$p.adjust < p.adjustFilter, ]
write.table(GO, file = file.path(outDir, "ORA_GO.txt"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("GO 富集通路数(p<0.05):", nrow(GO), "\n")
# ---- KEGG: use tryCatch so that KEGG API/network failures do not interrupt the entire script ----
kkKEGG <- tryCatch(
  enrichKEGG(gene = geneEntrez, organism = "hsa",
             pvalueCutoff = 1, qvalueCutoff = 1),
  error = function(e) {
    cat("[警告] enrichKEGG 失败（多为网络/KEGG API 问题）:", e$message, "\n")
    NULL
  }
)
KEGG <- if (!is.null(kkKEGG)) as.data.frame(kkKEGG) else data.frame()
KEGG <- KEGG[KEGG$pvalue < pvalueFilter & KEGG$p.adjust < p.adjustFilter, ]
if (nrow(KEGG) > 0) {
  eg2sym <- mapIds(org.Hs.eg.db, keys = unique(unlist(strsplit(KEGG$geneID, "/"))),
                   column = "SYMBOL", keytype = "ENTREZID")
  KEGG$geneSym <- sapply(strsplit(KEGG$geneID, "/"), function(x) {
    s <- eg2sym[x]
    s <- s[!is.na(s)]                      # Remove ENTREZIDs without a gene-symbol mapping
    paste(s, collapse = "/")
  })
  KEGG$geneID <- KEGG$geneSym
  KEGG <- KEGG[KEGG$geneID != "", , drop = FALSE]
}
write.table(KEGG, file = file.path(outDir, "ORA_KEGG.txt"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("KEGG 富集通路数(p<0.05):", nrow(KEGG), "\n")
# ===================== 3. ORA for H, C7, and C2-Reactome using predefined MSigDB gene sets =====================
# Note: enricher requires gene symbols matching the gene_symbol column in
# msigdbr; do not use the ENTREZIDs generated above.
runMSigORA <- function(nm, t2g) {
  cat("\n########## 正在处理:", nm, "##########\n")
  t2g <- as.data.frame(t2g)
  kk <- tryCatch(
    enricher(gene = genes, TERM2GENE = t2g,
             pvalueCutoff = 1, qvalueCutoff = 1,
             minGSSize = minGSSize),
    error = function(e) { cat("  [警告]", nm, "enricher 失败:", e$message, "\n"); NULL }
  )
  if (is.null(kk) || nrow(as.data.frame(kk)) == 0) {
    cat("  ", nm, "无任何富集结果\n")
    return(NULL)
  }
  res <- as.data.frame(kk)
  res <- res[res$pvalue < pvalueFilter & res$p.adjust < p.adjustFilter, ]
  if (nrow(res) > 0) {
    write.table(res, file = file.path(outDir, paste0("ORA_", nm, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  cat("富集通路数(p<0.05):", nrow(res), "\n")
  list(kk = kk, res = res)
}
# ---- Hallmark: all 50 gene sets ----
h_gs <- msigdbr(species = species, collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)
cat("\nHallmark 基因集数:", length(unique(h_gs$gs_name)), "\n")
outH <- runMSigORA("H", h_gs)
# ---- C7 immune signatures: selected M1/M2, antigen processing/presentation, MHC, phagocytosis, and lysosome sets ----
c7_all <- msigdbr(species = species, collection = "C7") %>%
  dplyr::select(gs_name, gene_symbol)
c7_keep <- unique(c(
  grep("GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE",
       unique(c7_all$gs_name), value = TRUE),
  grep("ANTIGEN_PROCESSING|ANTIGEN_PRESENTATION|CROSS_PRESENTATION|MHC|PHAGOCYTOSIS|LYSOSOME",
       unique(c7_all$gs_name), value = TRUE)
))
cat("C7 精选免疫签名数:", length(c7_keep), "\n")
print(c7_keep)
outC7 <- runMSigORA("C7", c7_all[c7_all$gs_name %in% c7_keep, ])
# ---- C2-Reactome: selected lysosome, phagosome, and antigen processing/presentation pathways ----
c2r_all <- msigdbr(species = species, collection = "C2", subcollection = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol)
c2r_keep <- grep("LYSOSOME|PHAGOSOME|ANTIGEN_PROCESSING|MHC_CLASS",
                 unique(c2r_all$gs_name), value = TRUE)
cat("C2-Reactome 精选通路数:", length(c2r_keep), "\n")
print(c2r_keep)
outC2R <- runMSigORA("C2Reactome", c2r_all[c2r_all$gs_name %in% c2r_keep, ])
# ===================== 4. Generate bar plots and dot plots =====================
plotORA <- function(kk, nm) {
  if (is.null(kk) || nrow(as.data.frame(kk)) == 0) {
    cat("  ", nm, "无显著通路，跳过绘图\n")
    return(invisible(NULL))
  }
  nshow <- min(topN, nrow(as.data.frame(kk)))
  pdf(file.path(outDir, paste0("ORA_", nm, ".barplot.pdf")), width = 10, height = 7)
  print(barplot(kk, drop = TRUE, showCategory = nshow, label_format = 130, color = colorSel))
  dev.off()
  pdf(file.path(outDir, paste0("ORA_", nm, ".bubble.pdf")), width = 10, height = 7)
  print(dotplot(kk, showCategory = nshow, orderBy = "GeneRatio",
                label_format = 130, color = colorSel))
  dev.off()
  cat("  图表已保存:", nm, "\n")
}
plotORA(kkGO, "GO")
plotORA(if (!is.null(kkKEGG)) kkKEGG else NULL, "KEGG")
plotORA(if (!is.null(outH)) outH$kk else NULL, "H")
plotORA(if (!is.null(outC7)) outC7$kk else NULL, "C7")
plotORA(if (!is.null(outC2R)) outC2R$kk else NULL, "C2Reactome")
cat("\n全部完成，输出目录:", outDir, "\n")
#######################################################################
# Output description
# - ORA_GO.txt / ORA_KEGG.txt / ORA_H.txt / ORA_C7.txt /
#   ORA_C2Reactome.txt: enrichment tables for each collection (nominal P < 0.05)
# - ORA_*.barplot.pdf / ORA_*.bubble.pdf: bar plots and dot plots
# - This script identifies enriched pathways without assigning direction.
#   To evaluate pathway scores in relation to endogenous HLA-DOA-high versus
#   HLA-DOA-low expression, use GSVA_direction_from_ORA.R.
# - If KEGG fails, an empty table is written or plotting is skipped without
#   affecting the GO and MSigDB analyses.
#######################################################################
