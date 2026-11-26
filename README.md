# Enrichment and GSVA direction analysis for HLA-DOA in cervical squamous cell carcinoma macrophages

R code accompanying the manuscript "Exploring HLA-DOA and Candidate Therapeutic Compounds in Cervical Squamous Cell Carcinoma Using Single-Cell Mendelian Randomization and Multi-Omics Integration".

This repository contains the two scripts behind the pathway part of the study: over-representation analysis (ORA) of the differentially regulated genes identified after a network-level gain-of-function perturbation of HLA-DOA regulatory output in tumour macrophages, and GSVA-based scoring that evaluates whether these pathways differ between endogenous HLA-DOA-high and HLA-DOA-low macrophage pseudobulk profiles. The GSVA step is observational and does not directly measure transcriptional responses to the network perturbation.

## 1. What the two scripts do

| Script | Question answered | Method |
|---|---|---|
| `ORA_diffgenes_all.R` | Which pathways contain the HLA-DOA perturbation-associated differentially regulated genes? | Hypergeometric ORA with clusterProfiler on five collections: GO (BP/CC/MF), KEGG, MSigDB Hallmark (H), C7 immune signatures, and C2-Reactome |
| `GSVA_direction_from_ORA.R` | Do these pathway scores differ between endogenous HLA-DOA-high and HLA-DOA-low macrophages? | Gene modules built from the differentially regulated gene members of each nominal ORA term, followed by GSVA scores (kcdf = "Gaussian", maxDiff = TRUE) and limma comparison of HLA-DOA-high versus HLA-DOA-low pseudobulk samples with sizeFactor adjustment and sample pairing |

Run the two scripts in this order: the second one reads the tables written by the first.

## 2. Repository structure

```text
.
├── ORA_diffgenes_all.R          # step 1: ORA (GO / KEGG / H / C7 / C2-Reactome)
├── GSVA_direction_from_ORA.R    # step 2: GSVA + limma direction per collection
├── README.md
└── (generated at run time)
    ├── ORA_out/                 # enrichment tables, bar plots, dot plots
    └── GSVA_out/                # GSVA scores, direction tables, bar plots, heatmaps
```

## 3. Requirements

R 4.4.2 and Seurat 5.4.0 were used for the analysis reported in the manuscript. The principal analysis packages were clusterProfiler 4.14.6, org.Hs.eg.db 3.20.0, msigdbr 26.1.1, GSVA 2.0.7, limma 3.66.0, edgeR 4.4.2, and GSEABase 1.68.0. Annotation resources were accessed on 26 September 2026. Install the remaining packages from CRAN or Bioconductor and check the console log if a function signature has changed.

```r
install.packages(c("dplyr", "ggplot2", "pheatmap", "msigdbr"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot",
                       "GSVA", "GSEABase", "limma", "edgeR"))
install.packages("Seurat")   # Seurat 5.x
```

`msigdbr` is used through its collection/subcollection interface (msigdbr >= 7.5), i.e. `msigdbr(species, collection = "H")` rather than the deprecated `category` argument.

## 4. Inputs

- `sigDiff.txt` — table of the significant differentially regulated genes obtained after the network-level HLA-DOA perturbation (tab-separated, header required, first column = gene symbol). Only the first column is used.
- `seurat.Rdata` — a Seurat object that must be named `pbmc` in the workspace (`stopifnot(exists("pbmc"))`). It is expected to contain:
  - the macrophage subset only (`cellType == "macrophages"` is applied again inside the script);
  - a `Type` metadata column with `Control` / `Treat`;
  - an `orig.ident` column identifying the original sample, used as the block for `duplicateCorrelation`;
  - normalized expression in the RNA assay `data` layer.

Both inputs and the `ORA_out/` directory are referenced by the paths set in the Settings block at the top of each script; edit them there rather than inside the code body. For the Treat-only analysis reported in the manuscript, set `keepTypes <- c("Treat")`.

## 5. Running

```r
# step 1
Rscript ORA_diffgenes_all.R

# step 2 (reads ORA_out/, writes GSVA_out/)
Rscript GSVA_direction_from_ORA.R
```

Both scripts are also runnable interactively in RStudio; `set.seed(123)` is fixed in step 2.

## 6. Outputs

### `ORA_out/`

- `ORA_GO.txt`, `ORA_KEGG.txt`, `ORA_H.txt`, `ORA_C7.txt`, `ORA_C2Reactome.txt` — enrichment tables filtered at nominal P < 0.05; inspect `p.adjust` for multiple-testing-adjusted significance.
- `ORA_<collection>.barplot.pdf`, `ORA_<collection>.bubble.pdf` — bar and dot plots.

### `GSVA_out/`

- `GSVA_scores_<collection>.csv` — raw GSVA score matrix.
- `GSVA_direction_<collection>.csv` — logFC, nominal P, Benjamini–Hochberg-adjusted P, and direction for every gene set (logFC > 0 = higher score in the HLA-DOA-high group).
- `GSVA_direction_barplot_<collection>.pdf/png` — direction bar plot with significance marks.
- `GSVA_heatmap_<collection>.pdf/png` — GSVA scores across pseudobulk samples (row-scaled).

## 7. Statistical design of the direction analysis

For the Treat-only analysis, the CSCC macrophage subset contained 2,930 cells, of which 1,228 (41.9%) had detectable HLA-DOA expression. Zero-expression cells were excluded to avoid a contrast dominated by HLA-DOA presence versus absence. The remaining HLA-DOA-positive cells were split at the median into high and low groups (614 cells each).

Cells were aggregated by Type, original sample, and HLA-DOA group. Pseudobulk profiles containing fewer than 10 cells were excluded, yielding eight pseudobulk samples from four CSCC samples (paired high/low profiles per sample). Genes detected in at least two pseudobulk samples were retained, and counts were converted to log2-CPM over all genes using edgeR.

Each nominal ORA term (P < 0.05) contributed a gene module formed by its differentially regulated gene members; modules with at least five genes detectable in the pseudobulk matrix were retained. A module comprising all significant differentially regulated genes was also evaluated. GSVA scores were compared with limma using log-transformed total aggregated RNA counts (`sizeFactor`) as a covariate. For the Treat-only analysis, the design was `~ sizeFactor + HLA-DOA_group`. When both Control and Treat profiles are present, the script uses `~ sizeFactor + Type + HLA-DOA_group`. `duplicateCorrelation` accounts for pairing by original sample, and `topTable` reports Benjamini–Hochberg-adjusted P values within each collection.

Two points should be kept in mind when reading the direction results:

- ORA tables are filtered on nominal P values. Nominal enrichment and Benjamini–Hochberg-adjusted significance must be distinguished using the `pvalue` and `p.adjust` columns. The GSVA analysis evaluates endogenous expression-associated differences and should not be interpreted as a causal response to the network perturbation.
- C7 immune signatures were tested by ORA but yielded no terms at nominal P < 0.05. Because `ORA_C7.txt` is only written when a nominal term is retained, no C7 ORA table was available for module construction. The GSVA script therefore fell back to the corresponding predefined MSigDB sets and printed the switch to the console. The fallback is flagged in the log, not in the output file name.

## 8. Data

All data analysed with these scripts are public:

- single-cell RNA-seq: GEO GSE208653, GSE197461;
- spatial transcriptomics: CROST VISDS000672, VISDS000674;
- immune-cell eQTL (exposure): OneK1K (Yazar et al., Science 2022);
- cervical cancer GWAS (outcome): BioBank Japan `ebi-a-GCST90018597` (Sakaue et al., Nature Genetics 2021), accessed through IEU OpenGWAS.

## 9. Citation

If you use these scripts, please cite the manuscript above (full reference to be added upon publication) and the methods they depend on: clusterProfiler (Wu et al., Innovation 2021), GSVA (Hänzelmann et al., BMC Bioinformatics 2013), limma (Ritchie et al., Nucleic Acids Research 2015), Seurat (Hao et al., Cell 2021), and MSigDB (Subramanian et al., PNAS 2005).

## 10. Contact

Miao Shengdi, Hunan Normal University — miaoshengdi@hunnu.edu.cn
