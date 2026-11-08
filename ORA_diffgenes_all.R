#######################################################################
# ORA enrichment analysis (differentially expressed genes)
#   -- GO + KEGG + Hallmark (H) + C7 + C2-Reactome
# ---------------------------------------------------------------
# Purpose: over-representation analysis (ORA) of sigDiff.txt, i.e. the
#          significant genes obtained after virtual overexpression of the
#          target gene. It answers "which pathways do the diff genes
#          concentrate in" (presence only, no direction; pathway
#          direction is handled by GSVA_direction_from_ORA.R).
# Gene sets and backgrounds:
#   - GO (BP/CC/MF)        : enrichGO, background = all annotated genes
#                            in org.Hs.eg.db
#   - KEGG                  : enrichKEGG, background = all KEGG-annotated genes
#   - H / C7 / C2-Reactome  : enricher, background = union of the selected
#                            MSigDB gene sets
# Output (ORA_out/): ORA_GO.txt / ORA_KEGG.txt / ORA_H.txt /
#                    ORA_C7.txt / ORA_C2Reactome.txt,
#                    plus bar plots and dot plots
# Filtering: nominal p < 0.05 (p.adjust is left unrestricted, consistent
#            with the analysis reported in the manuscript)
# Note: keep this file in UTF-8 encoding
#######################################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(msigdbr)
library(dplyr)

# ===================== 0. Settings =====================
sigDiffFile    <- "sigDiff.txt"         # input: table of significant genes (1st column = gene symbol)
outDir         <- "ORA_out"             # output directory
species        <- "Homo sapiens"        # species: Homo sapiens / Mus musculus
pvalueFilter   <- 0.05                  # nominal P cutoff
p.adjustFilter <- 1                     # adjusted P cutoff (1 = unrestricted, filter on nominal P only)
minGSSize      <- 5                     # minimum gene-set size for enricher (smaller sets are dropped)
topN           <- 30                    # number of pathways shown in the plots

dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
colorSel <- if (p.adjustFilter > 0.05) "pvalue" else "p.adjust"

# ===================== 1. Read the diff genes =====================
rt <- read.table(sigDiffFile, header = TRUE, sep = "\t", check.names = FALSE)
genes <- unique(as.vector(rt[, 1]))
cat("输入显著基因数:", length(genes), "\n")

# symbol -> ENTREZID (required by enrichGO / enrichKEGG)
entrez <- mapIds(org.Hs.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL")
geneEntrez <- as.character(entrez[!is.na(entrez)])
cat("转换到 ENTREZID 的基因数:", length(geneEntrez), "\n")

# ===================== 2. GO / KEGG ORA (whole-genome background) =====================
# ---- GO (BP/CC/MF) ----
kkGO <- enrichGO(gene = geneEntrez, OrgDb = org.Hs.eg.db,
                 pvalueCutoff = 1, qvalueCutoff = 1, ont = "all", readable = TRUE)
GO <- as.data.frame(kkGO)
GO <- GO[GO$pvalue < pvalueFilter & GO$p.adjust < p.adjustFilter, ]
write.table(GO, file = file.path(outDir, "ORA_GO.txt"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("GO 富集通路数(p<0.05):", nrow(GO), "\n")

# ---- KEGG ----
kkKEGG <- enrichKEGG(gene = geneEntrez, organism = "hsa",
                     pvalueCutoff = 1, qvalueCutoff = 1)
KEGG <- as.data.frame(kkKEGG)
KEGG <- KEGG[KEGG$pvalue < pvalueFilter & KEGG$p.adjust < p.adjustFilter, ]
if (nrow(KEGG) > 0) {
  eg2sym <- mapIds(org.Hs.eg.db, keys = unique(unlist(strsplit(KEGG$geneID, "/"))),
                   column = "SYMBOL", keytype = "ENTREZID")
  KEGG$geneSym <- sapply(strsplit(KEGG$geneID, "/"), function(x)
    paste(eg2sym[x], collapse = "/"))
  KEGG$geneID <- KEGG$geneSym
}
write.table(KEGG, file = file.path(outDir, "ORA_KEGG.txt"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("KEGG 富集通路数(p<0.05):", nrow(KEGG), "\n")

# ===================== 3. H / C7 / C2-Reactome ORA (pre-defined MSigDB sets) =====================
# Note: enricher takes gene symbols (the gene_symbol column of msigdbr),
#       not the ENTREZID used above.
runMSigORA <- function(nm, t2g) {
  cat("\n########## 正在处理:", nm, "##########\n")
  t2g <- as.data.frame(t2g)
  kk <- enricher(gene = genes, TERM2GENE = t2g,
                 pvalueCutoff = 1, qvalueCutoff = 1,
                 minGSSize = minGSSize)
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

# ---- H (Hallmark, all 50 sets) ----
h_gs <- msigdbr(species = species, collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)
cat("\nHallmark 基因集数:", length(unique(h_gs$gs_name)), "\n")
outH <- runMSigORA("H", h_gs)

# ---- C7 immune signatures (selected: M1/M2 plus antigen processing /
#      presentation, MHC, phagocytosis, lysosome) ----
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

# ---- C2-Reactome (selected: lysosome / phagosome / antigen processing
#      and presentation) ----
c2r_all <- msigdbr(species = species, collection = "C2", subcollection = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol)
c2r_keep <- grep("LYSOSOME|PHAGOSOME|ANTIGEN_PROCESSING|MHC_CLASS",
                 unique(c2r_all$gs_name), value = TRUE)
cat("C2-Reactome 精选通路数:", length(c2r_keep), "\n")
print(c2r_keep)
outC2R <- runMSigORA("C2Reactome", c2r_all[c2r_all$gs_name %in% c2r_keep, ])

# ===================== 4. Plots (bar plot + dot plot) =====================
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
plotORA(kkKEGG, "KEGG")
plotORA(if (!is.null(outH)) outH$kk else NULL, "H")
plotORA(if (!is.null(outC7)) outC7$kk else NULL, "C7")
plotORA(if (!is.null(outC2R)) outC2R$kk else NULL, "C2Reactome")

cat("\n全部完成，输出目录:", outDir, "\n")

#######################################################################
# Output description
# - ORA_GO.txt / ORA_KEGG.txt / ORA_H.txt / ORA_C7.txt /
#   ORA_C2Reactome.txt : enrichment tables for each collection (nominal p < 0.05)
# - ORA_*.barplot.pdf / ORA_*.bubble.pdf : bar plots / dot plots
# - This script only reports "where the diff genes are enriched" (no
#   direction). Pathway direction under the overexpressed state is
#   obtained with GSVA_direction_from_ORA.R.
#######################################################################
