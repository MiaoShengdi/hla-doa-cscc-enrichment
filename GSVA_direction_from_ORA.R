#######################################################################
# GSVA direction analysis (on the ORA-enriched pathways)
#   -- assigns a sign (up / down) to the pathways enriched by the diff genes
# ---------------------------------------------------------------
# Purpose: read the ORA tables (ORA_GO.txt / ORA_KEGG.txt / ORA_H.txt /
#          ORA_C7.txt / ORA_C2Reactome.txt, produced by ORA_diffgenes_all.R),
#          build a gene module from the diff-gene members of every
#          significant pathway, then score each module by GSVA and compare
#          groups with limma, so that each pathway receives an up / down
#          direction in the "target gene high" (virtual overexpression)
#          group.
# Direction interpretation:
#   - logFC > 0 : the pathway is more active in the HLA-DOA-high group
#                 (positive / up-regulated)
#   - logFC < 0 : the pathway is more active in the HLA-DOA-low group
#                 (negative / down-regulated)
# Fallback: when an ORA table is empty or missing (C7 usually has no
#          significant pathway), the script switches to the pre-defined
#          MSigDB sets (all H / selected C7 / selected C2-Reactome) for
#          GSVA. The switch is printed to the console, so every collection
#          still yields a direction table.
# Note: keep this file in UTF-8 encoding
#######################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(GSVA)
  library(GSEABase)
  library(limma)
  library(edgeR)
  library(msigdbr)
  library(pheatmap)
  library(ggplot2)
  library(dplyr)
})

set.seed(123)

# ===================== 0. Settings =====================
targetGene        <- "HLA-DOA"      # target gene (used for expression-based grouping, i.e. the virtual overexpression state)
cellName          <- "macrophages"  # cell type analysed
keepTypes         <- c("Control", "Treat") # sample groups kept for analysis
dataFile          <- "seurat.Rdata"  # input 1: Seurat object (must be named `pbmc`)
oraDir            <- "ORA_out"       # input 2: directory holding the ORA tables
sigDiffFile       <- "sigDiff.txt"   # input 3 (optional): diff-gene table (whole-module direction)
outDir            <- "GSVA_out"      # output directory
species           <- "Homo sapiens"
groupMethod       <- "median"        # "median" = median split; "expressed" = expressed (>0) vs not expressed
minCellsPerSample <- 10
pairedBySample    <- TRUE
minModuleSize     <- 5               # minimum gene-module size (smaller modules are dropped)
topN              <- 30
useFallback       <- TRUE            # TRUE = fall back to pre-defined MSigDB sets when ORA is empty

dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# ===================== Helper functions (plotting) =====================
cleanName <- function(x) {
  x <- gsub("^HALLMARK_|^REACTOME_|^KEGG_|^GOBP_|^GOCC_|^GOMF_", "", x)
  x <- gsub("_", " ", x)
  x
}
strwrap2 <- function(x, width = 45) {
  vapply(as.character(x), function(s)
    paste(strwrap(s, width = width), collapse = "\n"), character(1))
}
sigStar <- function(p, padj) ifelse(padj < 0.05, "**", ifelse(p < 0.05, "*", ""))

# ===================== 1. Read the inputs =====================
load(dataFile)
stopifnot(exists("pbmc"))

if ("Type" %in% colnames(pbmc@meta.data) && !is.null(keepTypes)) {
  pbmc <- subset(pbmc, subset = Type %in% keepTypes)
  cat("保留样本组:", paste(keepTypes, collapse = ", "), "\n")
}
pbmc <- subset(pbmc, subset = cellType == cellName)
cat("剩余细胞数:", ncol(pbmc), "\n")
stopifnot(targetGene %in% rownames(pbmc))

geneExpr <- GetAssayData(pbmc, assay = "RNA", layer = "data")[targetGene, ]
if (groupMethod == "median") {
  pbmc$gene_group <- ifelse(geneExpr > median(geneExpr), "HLA-DOA-high", "HLA-DOA-low")
} else {
  pbmc$gene_group <- ifelse(geneExpr > 0, "HLA-DOA-high", "HLA-DOA-low")
}
print(table(pbmc$gene_group))

if (!"orig.ident" %in% colnames(pbmc@meta.data)) pbmc$orig.ident <- "sample"
pbmc$grp <- if ("Type" %in% colnames(pbmc@meta.data))
  paste(pbmc$Type, pbmc$orig.ident, pbmc$gene_group, sep = ".") else
  paste(pbmc$orig.ident, pbmc$gene_group, sep = ".")

ncells <- table(pbmc$grp)
cat("各伪样本细胞数:\n"); print(ncells)

pseudo <- AggregateExpression(pbmc, group.by = "grp", assays = "RNA",
                              return.seurat = FALSE)$RNA
keepCols  <- names(ncells)[ncells >= minCellsPerSample]
keepCols  <- intersect(keepCols, colnames(pseudo))
if (length(keepCols) < 2) stop("伪样本不足 2 个：请调低 minCellsPerSample 或检查分组")
pseudo    <- pseudo[, keepCols, drop = FALSE]
keepRows  <- rowSums(pseudo > 0) >= 2
pseudo    <- pseudo[keepRows, , drop = FALSE]
pseudo    <- cpm(pseudo, log = TRUE)          # log2-CPM over all genes
cat("伪样本:", colnames(pseudo), "\n")

pseudoType <- unique(gsub("\\..*$", "", colnames(pseudo)))
typeLabel  <- if (length(pseudoType) == 1) {
  if (pseudoType == "Control") "Control group (normal)" else "Treat group (disease)"
} else paste(paste(pseudoType, collapse = " + "), "groups")
footNote <- paste0(typeLabel, "  |  Cell type: ", cellName)

# ===================== 2. Read the ORA tables and build the modules =====================
readOraTable <- function(f) {
  p <- file.path(oraDir, f)
  if (!file.exists(p)) { cat("  [提示] 未找到", f, "，视为空\n"); return(NULL) }
  df <- read.delim(p, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(NULL)
  df
}

# Build one gene module per significant pathway from its diff-gene members
# (geneID column, "/"-separated symbols)
oraModules <- function(ora) {
  gsList <- lapply(seq_len(nrow(ora)), function(i) {
    genes <- unique(strsplit(ora$geneID[i], "/")[[1]])
    genes <- genes[genes %in% rownames(pseudo)]
    if (length(genes) < minModuleSize) return(NULL)
    GeneSet(genes, setName = ora$ID[i])
  })
  gsList <- Filter(Negate(is.null), gsList)
  if (length(gsList) == 0) return(NULL)
  GeneSetCollection(gsList)
}

# Pre-defined MSigDB sets (used as fallback)
df2gsc <- function(df) {
  gs <- split(df$gene_symbol, df$gs_name)
  gs <- lapply(gs, unique)
  gs <- gs[sapply(gs, function(g) sum(g %in% rownames(pseudo)) >= minModuleSize)]
  if (length(gs) == 0) return(NULL)
  GeneSetCollection(lapply(names(gs), function(nm) GeneSet(gs[[nm]], setName = nm)))
}
buildFallback <- function(nm) {
  if (nm == "H") {
    h_df <- msigdbr(species = species, collection = "H") %>%
      dplyr::select(gs_name, gene_symbol)
    df2gsc(h_df)
  } else if (nm == "C7") {
    c7_all <- msigdbr(species = species, collection = "C7") %>%
      dplyr::select(gs_name, gene_symbol)
    c7_keep <- unique(c(
      grep("GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE",
           unique(c7_all$gs_name), value = TRUE),
      grep("ANTIGEN_PROCESSING|ANTIGEN_PRESENTATION|CROSS_PRESENTATION|MHC|PHAGOCYTOSIS|LYSOSOME",
           unique(c7_all$gs_name), value = TRUE)
    ))
    df2gsc(c7_all[c7_all$gs_name %in% c7_keep, ])
  } else {  # C2Reactome
    c2r_all <- msigdbr(species = species, collection = "C2", subcollection = "CP:REACTOME") %>%
      dplyr::select(gs_name, gene_symbol)
    c2r_keep <- grep("LYSOSOME|PHAGOSOME|ANTIGEN_PROCESSING|MHC_CLASS",
                     unique(c2r_all$gs_name), value = TRUE)
    df2gsc(c2r_all[c2r_all$gs_name %in% c2r_keep, ])
  }
}

collections <- list(
  GO         = list(file = "ORA_GO.txt",         fallback = FALSE),
  KEGG       = list(file = "ORA_KEGG.txt",       fallback = FALSE),
  H          = list(file = "ORA_H.txt",          fallback = TRUE),
  C7         = list(file = "ORA_C7.txt",         fallback = TRUE),
  C2Reactome = list(file = "ORA_C2Reactome.txt", fallback = TRUE)
)

gscList <- list()
displayList <- list()
for (nm in names(collections)) {
  ora <- readOraTable(collections[[nm]]$file)
  gsc <- NULL; disp <- NULL
  if (!is.null(ora)) {
    gsc  <- oraModules(ora)
    disp <- setNames(ora$Description, ora$ID)   # pathway names shown in the plots
    cat(nm, ": ORA 显著通路", nrow(ora), "条 → 可用模块",
        ifelse(is.null(gsc), 0, length(gsc)), "个\n")
  }
  if (is.null(gsc) && collections[[nm]]$fallback && useFallback) {
    gsc <- buildFallback(nm)
    cat("  !! ", nm, " ORA 无可用通路 → 回退 MSigDB 预设基因集（",
        ifelse(is.null(gsc), 0, length(gsc)), "个），方向结果来自预设基因集\n")
  }
  gscList[[nm]]     <- gsc
  displayList[[nm]] <- disp
}

# Whole diff-gene module (optional): all diff genes as a single set; direction only, no plot
gscAll <- NULL
if (file.exists(sigDiffFile)) {
  sig <- read.delim(sigDiffFile, stringsAsFactors = FALSE)
  sg  <- unique(as.vector(sig[, 1]))
  sg  <- intersect(sg, rownames(pseudo))
  if (length(sg) >= minModuleSize)
    gscAll <- GeneSet(sg, setName = paste0(targetGene, "_diff_all"))
  cat("整体 diff genes 模块基因数:", length(sg), "\n")
}

# ===================== 3. Per collection: GSVA + limma direction + plots =====================
runSet <- function(nm, gsc, display = NULL) {
  cat("\n########## 正在处理:", nm, "##########\n")
  if (is.null(gsc)) { message("  ", nm, " 无可用基因集，跳过"); return(NULL) }
  if (!is.null(gscAll)) gsc <- GeneSetCollection(c(gsc, gscAll))

  # ---- GSVA scoring (the expression matrix always keeps all genes) ----
  gp <- gsvaParam(exprData = pseudo, geneSets = gsc,
                  kcdf = "Gaussian", maxDiff = TRUE, tau = 1)
  es <- gsva(gp)
  write.csv(es, file.path(outDir, paste0("GSVA_scores_", nm, ".csv")))

  # ---- limma: HLA-DOA-high vs HLA-DOA-low ----
  grp  <- factor(gsub(".*\\.", "", colnames(es)), levels = c("HLA-DOA-low", "HLA-DOA-high"))
  type <- gsub("\\..*$", "", colnames(es))
  if (length(unique(type)) > 1) {
    type <- factor(type)
    design <- model.matrix(~ type + grp)
    coefName <- "grpHLA-DOA-high"
    cat("  设计矩阵: ~ Type + HLA-DOA_group（校正", paste(levels(type), collapse = "/"),
        "后比较 HLA-DOA 高低）\n")
  } else {
    design <- model.matrix(~ grp)
    coefName <- "grpHLA-DOA-high"
  }
  if (pairedBySample) {
    sampleID <- gsub("\\.HLA-DOA-.*$", "", colnames(es))
    corfit   <- duplicateCorrelation(es, design, block = sampleID)
    fit      <- lmFit(es, design, block = sampleID, correlation = corfit$consensus)
  } else fit <- lmFit(es, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = coefName, number = Inf, sort.by = "none")
  res <- res[order(res$logFC, decreasing = TRUE), ]
  res$direction <- ifelse(res$logFC > 0, "HLA-DOA_high(正向)", "HLA-DOA_low(负向)")
  if (!is.null(display)) res$Description <- display[rownames(res)]
  write.csv(res, file.path(outDir, paste0("GSVA_direction_", nm, ".csv")))
  cat("基因集数:", nrow(res), "| 正向:", sum(res$logFC > 0), "| 负向:", sum(res$logFC < 0), "\n")

  # Direction of the whole diff-gene module (printed only, not plotted)
  sigRow <- grep(paste0("^", targetGene, "_diff_all$"), rownames(res))
  if (length(sigRow))
    cat("  [整体模块] diff genes 全部 logFC =", round(res$logFC[sigRow], 3),
        " P =", round(res$P.Value[sigRow], 4), " ->", res$direction[sigRow], "\n")

  # ---- Direction bar plot (topN, cleaned names + significance stars) ----
  resMain <- res[!grepl(paste0("^", targetGene, "_diff_all$"), rownames(res)), , drop = FALSE]
  plotDF <- head(resMain[order(abs(resMain$logFC), decreasing = TRUE), ], topN)
  plotDF$label <- if (!is.null(display)) plotDF$Description else cleanName(rownames(plotDF))
  plotDF$star  <- sigStar(plotDF$P.Value, plotDF$adj.P.Val)
  plotDF$label <- strwrap2(plotDF$label, 48)
  plotDF$label <- factor(plotDF$label, levels = plotDF$label[order(plotDF$logFC)])
  p <- ggplot(plotDF, aes(x = label, y = logFC, fill = logFC > 0)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = star, hjust = ifelse(logFC >= 0, -0.35, 1.35)), size = 4.5) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
    labs(title = paste0(nm, ": GSVA direction (HLA-DOA_high vs HLA-DOA_low)"),
         subtitle = "logFC>0（红）= HLA-DOA 高表达组活性更高；* P<0.05；** adj.P<0.05",
         x = "", y = "GSVA logFC", caption = footNote) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          plot.caption = element_text(size = 9, colour = "grey30"))
  h1 <- max(6, 0.32 * nrow(plotDF) + 1.5)
  ggsave(file.path(outDir, paste0("GSVA_direction_barplot_", nm, ".pdf")), p,
         width = 9.5, height = h1)
  ggsave(file.path(outDir, paste0("GSVA_direction_barplot_", nm, ".png")), p,
         width = 9.5, height = h1, dpi = 300)

  # ---- Heatmap (topN, column annotation: Type + HLA-DOA group) ----
  topIds <- rownames(resMain)[order(abs(resMain$logFC), decreasing = TRUE)][1:min(topN, nrow(resMain))]
  esPlot <- es[topIds, , drop = FALSE]
  colType <- gsub("\\..*$", "", colnames(esPlot))
  colnames(esPlot) <- gsub("^(Control|Treat)\\.", "", gsub("\\.HLA-DOA-", "-", colnames(esPlot)))
  rowLab <- if (!is.null(display)) display[topIds] else cleanName(topIds)
  rownames(esPlot) <- strwrap2(rowLab, 34)
  annoCol <- data.frame(Type = colType,
                        `HLA-DOA` = ifelse(grepl("-high$", colnames(esPlot)), "high", "low"),
                        row.names = colnames(esPlot), check.names = FALSE)
  annColors <- list(Type = c(Control = "#3C5488", Treat = "#E64B35"),
                    `HLA-DOA` = c(high = "#F39B7F", low = "#4DBBD5"))
  h2 <- max(6.5, 0.32 * nrow(esPlot) + 3)
  hmArgs <- list(esPlot, scale = "row", cluster_cols = FALSE,
                 clustering_method = "ward.D2",
                 annotation_col = annoCol,
                 annotation_colors = annColors, angle_col = 45,
                 fontsize_row = 8, fontsize_col = 10,
                 main = paste0(nm, ": GSVA scores (row-scaled)"))
  pdf(file.path(outDir, paste0("GSVA_heatmap_", nm, ".pdf")), width = 9, height = h2)
  do.call(pheatmap, hmArgs)
  grid::grid.text(footNote, x = 0.98, y = 0.03, just = "right",
                  gp = grid::gpar(fontsize = 9, col = "grey30"))
  dev.off()
  png(file.path(outDir, paste0("GSVA_heatmap_", nm, ".png")),
      width = 9, height = h2, units = "in", res = 300)
  do.call(pheatmap, hmArgs)
  grid::grid.text(footNote, x = 0.98, y = 0.03, just = "right",
                  gp = grid::gpar(fontsize = 9, col = "grey30"))
  dev.off()
  cat("  图表已保存:", nm, "\n")
  list(res = res, es = es)
}

for (nm in names(gscList)) {
  out <- runSet(nm, gscList[[nm]], displayList[[nm]])
  assign(paste0("out", nm), out)
}

cat("\n全部完成，输出目录:", outDir, "\n")

#######################################################################
# Output description
# - GSVA_direction_*.csv : per-collection direction tables
#                          (logFC > 0 = positive / more active in the high group)
# - GSVA_direction_barplot_*.pdf/png : direction bar plots (with significance stars)
# - GSVA_heatmap_*.pdf/png : GSVA scores of each pathway across pseudobulk samples
# - GSVA_scores_*.csv : raw GSVA score matrices
# - This script answers whether the ORA-enriched pathways are up- or
#   down-regulated under the overexpressed state; "where the diff genes
#   are enriched" is answered by ORA_diffgenes_all.R. The two scripts are
#   designed to be run in this order.
#######################################################################
