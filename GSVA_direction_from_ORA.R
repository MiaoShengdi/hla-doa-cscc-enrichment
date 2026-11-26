#######################################################################
# GSVA 方向性分析（基于 ORA 富集通路）—— 给"diff genes 富集到的通路"加正负方向
# ---------------------------------------------------------------
# 作用：读取 ORA 富集表（ORA_GO.txt / ORA_KEGG.txt / ORA_H.txt /
#       ORA_C7.txt / ORA_C2Reactome.txt，来自 ORA_diffgenes_all.R 或
#       你已有的 GO/KEGG、H/C7/C2 富集脚本输出），把每条显著通路的
#       diff 基因成员构造成基因模块，再做 GSVA 打分 + limma 组间比较，
#       得到每条通路在内源性 HLA-DOA 高/低表达组之间的评分差异。
# 方向解读：
#   - logFC>0 = HLA-DOA 高表达组该通路的 GSVA 评分更高
#   - logFC<0 = HLA-DOA 低表达组该通路的 GSVA 评分更高
# 回退机制：某集合 ORA 表为空或缺失时（如 C7 常无显著通路），
#       自动改用 MSigDB 预设基因集（H 全部 / C7 精选 / C2-Reactome 精选）
#       做 GSVA，控制台会打印注明，保证每个集合都有方向结果。
#
# 【2026-11-02 偏倚修正（重要，务必阅读）】
# 排查发现原脚本存在分组系统性偏倚，导致方向结果几乎全正向（GO 434/436）：
#   ① 62.1% 细胞 HLA-DOA 零表达 → 全细胞中位数 = 0，median 分组实际是
#      "表达 vs 不表达"，而非"高 vs 低表达"；
#   ② "high 组"（表达 HLA-DOA 的细胞）nCount/nFeature 为 low 组的 2.6/1.8 倍，
#      属于更活化/更高测序深度的细胞，GSVA 分数因此整体上移。
#   随机基因集阴性对照（100 个，正向显著 3/100）未提示管线本身存在系统性方向偏倚。
# 本次修正：
#   A. 分组改为"仅在 HLA-DOA 表达>0 的细胞内按中位数切 high/low"；
#   B. limma 设计矩阵加入 sizeFactor（伪样本 log 总表达量）协变量，
#      校正测序深度/活化差异；
#   C. 修正 coefName 提取（model.matrix 会把连字符转为点，原固定字符串会报错）。
# 修正后方向结果方可用于论文；论文方法学应写明上述分组与校正策略。
# 注意：请用 UTF-8 编码保存（RStudio: File -> Save with Encoding -> UTF-8）
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
# ===================== 0. 参数设置（按你的数据修改） =====================
targetGene        <- "HLA-DOA"      # 目标基因（用于内源性表达分组，不等同于虚拟网络扰动）
cellName          <- "macrophages" # 细胞类型
keepTypes         <- c("Treat")     # 论文分析仅纳入 CSCC（Treat）样本组
dataFile          <- "seurat.Rdata"  # 输入1：Seurat 对象
oraDir            <- "ORA_out"       # 输入2：ORA 富集表所在目录
sigDiffFile       <- "sigDiff.txt"   # 输入3（可选）：差异基因表（整体模块）
outDir            <- "GSVA_out"      # 输出目录
species           <- "Homo sapiens"
groupMethod       <- "median"        # "median"=中位数分组；"expressed"=按是否表达(>0)分组
minCellsPerSample <- 10
pairedBySample    <- TRUE
minModuleSize     <- 5               # 基因模块最少基因数（不足剔除）
topN              <- 30
useFallback       <- TRUE            # TRUE=ORA 为空时回退 MSigDB 预设基因集
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
# ===================== 辅助函数（作图用） =====================
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
# ===================== 1. 读取输入 =====================
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
# ----【修正 A】分组：仅在表达>0 的细胞内按中位数切组，避免零膨胀退化成"表达 vs 不表达" ----
if (groupMethod == "median") {
  posCells <- geneExpr > 0
  cat("HLA-DOA 表达>0 细胞比例:", round(mean(posCells) * 100, 1), "%\n")
  pbmc <- pbmc[, posCells]                                  # 剔除零表达细胞
  geneExpr <- GetAssayData(pbmc, assay = "RNA", layer = "data")[targetGene, ]
  pbmc$gene_group <- ifelse(geneExpr > median(geneExpr), "HLA-DOA-high", "HLA-DOA-low")
} else {
  pbmc$gene_group <- ifelse(geneExpr > 0, "HLA-DOA-high", "HLA-DOA-low")
}
print(table(pbmc$gene_group))
# ----【修正 A 复核】分组修正后 high/low 细胞质量是否仍失衡（期望差异大幅缩小） ----
cat("\n【分组修正后质量复核】\n")
for (m in c("nCount_RNA", "nFeature_RNA")) {
  hi <- pbmc[[m]][pbmc$gene_group == "HLA-DOA-high", 1]
  lo <- pbmc[[m]][pbmc$gene_group == "HLA-DOA-low",  1]
  w  <- suppressWarnings(wilcox.test(hi, lo))
  cat(m, "| high中位:", round(median(hi), 1),
      " low中位:", round(median(lo), 1),
      " | wilcoxon P:", format.pval(w$p.value, digits = 3), "\n")
}
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
pseudo    <- cpm(pseudo, log = TRUE)          # 全基因 log2-CPM
cat("伪样本:", colnames(pseudo), "\n")
# ----【修正 B】sizeFactor：每个伪样本的 log 总表达量（测序深度/活化协变量） ----
aggCounts  <- AggregateExpression(pbmc, group.by = "grp", assays = "RNA",
                                  return.seurat = FALSE)$RNA
sizeFactor <- log(colSums(aggCounts[, keepCols, drop = FALSE]))
names(sizeFactor) <- keepCols
cat("sizeFactor（log 总表达量）:\n"); print(round(sizeFactor, 2))
pseudoType <- unique(gsub("\\..*$", "", colnames(pseudo)))
typeLabel  <- if (length(pseudoType) == 1) {
  if (pseudoType == "Control") "Control group (normal)" else "Treat group (disease)"
} else paste(paste(pseudoType, collapse = " + "), "groups")
footNote <- paste0(typeLabel, "  |  Cell type: ", cellName)
# ===================== 2. 读取 ORA 表，构建模块基因集 =====================
readOraTable <- function(f) {
  p <- file.path(oraDir, f)
  if (!file.exists(p)) { cat("  [提示] 未找到", f, "，视为空\n"); return(NULL) }
  df <- read.delim(p, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(NULL)
  df
}
# 把 ORA 显著通路的 diff 基因成员（geneID，"/"分隔，symbol）构造成基因模块
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
# MSigDB 预设基因集（回退用）
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
  GO         = list(file = "ORA_GO.txt",        fallback = FALSE),
  KEGG       = list(file = "ORA_KEGG.txt",      fallback = FALSE),
  H          = list(file = "ORA_H.txt",         fallback = TRUE),
  C7         = list(file = "ORA_C7.txt",        fallback = TRUE),
  C2Reactome = list(file = "ORA_C2Reactome.txt", fallback = TRUE)
)
gscList <- list()
displayList <- list()
for (nm in names(collections)) {
  ora <- readOraTable(collections[[nm]]$file)
  gsc <- NULL; disp <- NULL
  if (!is.null(ora)) {
    gsc  <- oraModules(ora)
    disp <- setNames(ora$Description, ora$ID)   # 图上显示通路名
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
# 整体差异基因模块（可选）：所有 diff genes 作为一个基因集，只报告方向、不出图
gscAll <- NULL
if (file.exists(sigDiffFile)) {
  sig <- read.delim(sigDiffFile, stringsAsFactors = FALSE)
  sg  <- unique(as.vector(sig[, 1]))
  sg  <- intersect(sg, rownames(pseudo))
  if (length(sg) >= minModuleSize)
    gscAll <- GeneSet(sg, setName = paste0(targetGene, "_diff_all"))
  cat("整体 diff genes 模块基因数:", length(sg), "\n")
}
# ===================== 3. 每个集合：GSVA + limma 方向 + 出图 =====================
runSet <- function(nm, gsc, display = NULL) {
  cat("\n########## 正在处理:", nm, "##########\n")
  if (is.null(gsc)) { message("  ", nm, " 无可用基因集，跳过"); return(NULL) }
  if (!is.null(gscAll)) gsc <- GeneSetCollection(c(gsc, gscAll))
  # ---- GSVA 打分（表达矩阵始终全基因） ----
  gp <- gsvaParam(exprData = pseudo, geneSets = gsc,
                  kcdf = "Gaussian", maxDiff = TRUE, tau = 1)
  es <- gsva(gp)
  write.csv(es, file.path(outDir, paste0("GSVA_scores_", nm, ".csv")))
  # ---- limma 比较 HLA-DOA-high vs HLA-DOA-low（校正 sizeFactor + Type） ----
  grp  <- factor(gsub(".*\\.", "", colnames(es)), levels = c("HLA-DOA-low", "HLA-DOA-high"))
  type <- gsub("\\..*$", "", colnames(es))
  sf   <- sizeFactor[colnames(es)]            # 伪样本测序深度协变量
  if (length(unique(type)) > 1) {
    type <- factor(type)
    design <- model.matrix(~ sf + type + grp)
    cat("  设计矩阵: ~ sizeFactor + Type + HLA-DOA_group（校正总表达与",
        paste(levels(type), collapse = "/"), "后比较 HLA-DOA 高低）\n")
  } else {
    design <- model.matrix(~ sf + grp)
    cat("  设计矩阵: ~ sizeFactor + HLA-DOA_group\n")
  }
  # ----【修正 C】动态提取 grp 系数名（model.matrix 会把连字符"转成点"） ----
  coefName <- grep("^grp", colnames(design), value = TRUE)[1]
  cat("  比较系数:", coefName, "\n")
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
  # 整体 diff 模块方向（只打印，不出图）
  sigRow <- grep(paste0("^", targetGene, "_diff_all$"), rownames(res))
  if (length(sigRow))
    cat("  [整体模块] diff genes 全部 logFC =", round(res$logFC[sigRow], 3),
        " P =", round(res$P.Value[sigRow], 4), " ->", res$direction[sigRow], "\n")
  # ---- 方向条形图（topN，清理名称 + 星号） ----
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
  # ---- 热图（topN，列注释含 Type + HLA-DOA） ----
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
# 结果说明
# - GSVA_direction_*.csv：各集合方向表（logFC>0=正向/高表达组活性更高；
#   已校正 sizeFactor 测序深度与 Type 协变量，基于"表达 HLA-DOA 细胞内部"的分组）
# - GSVA_direction_barplot_*.pdf/png：方向条形图（带星号）
# - GSVA_heatmap_*.pdf/png：各通路在各伪样本的 GSVA 分数热图
# - GSVA_scores_*.csv：原始 GSVA 打分矩阵
# - 本脚本回答"diff genes 富集到的通路在内源性 HLA-DOA 高/低表达组之间的评分差异"；
#   该观察性比较不等同于虚拟网络扰动后的直接转录响应。
#   "富集在哪"由 ORA_diffgenes_all.R 回答，两个脚本配套使用。
# - 论文方法学表述建议：在表达 HLA-DOA 的巨噬细胞中按表达中位数分为
#   high/low 两组，伪样本化后以 log 总表达量（sizeFactor）与样本组（Type）
#   为协变量，limma 比较两组 GSVA 通路分数（duplicateCorrelation 配对样本）。
#######################################################################
