# Enrichment and GSVA direction analysis for HLA-DOA virtual overexpression in cervical squamous cell carcinoma macrophages

R code accompanying the manuscript *"Single-cell Mendelian randomization and multi-omics integration reveal the role of HLA-DOA and potential drug targets in cervical squamous cell carcinoma"*.

This repository contains the two scripts behind the pathway part of the study: over-representation analysis (ORA) of the genes differentially expressed after *in silico* overexpression of **HLA-DOA** in tumour macrophages, and GSVA-based scoring that assigns an up/down direction to each enriched pathway.

---

## 1. What the two scripts do

| Script | Question answered | Method |
|---|---|---|
| `ORA_diffgenes_all.R` | *Which pathways do the HLA-DOA perturbation genes fall into?* | Hypergeometric ORA (`clusterProfiler`) on five collections: GO (BP/CC/MF), KEGG, MSigDB Hallmark (H), C7 immune signatures, C2-Reactome |
| `GSVA_direction_from_ORA.R` | *Are those pathways more or less active when HLA-DOA is high?* | Gene modules built from the diff-gene members of each ORA pathway → GSVA scores (`kcdf = "Gaussian"`, `maxDiff = TRUE`) → `limma` comparison of HLA-DOA-high vs HLA-DOA-low pseudobulk samples |

Run the two scripts **in this order**: the second one reads the tables written by the first.

## 2. Repository structure

```
.
├── ORA_diffgenes_all.R          # step 1: ORA (GO / KEGG / H / C7 / C2-Reactome)
├── GSVA_direction_from_ORA.R    # step 2: GSVA + limma direction per collection
├── README.md
└── (generated at run time)
    ├── ORA_out/                 # enrichment tables, bar plots, dot plots
    └── GSVA_out/                # GSVA scores, direction tables, bar plots, heatmaps
```

## 3. Requirements

R 4.4.2 with Seurat 5.4.0 was used for the analysis reported in the manuscript. The remaining packages were installed from CRAN / Bioconductor at the time of the analysis; install the current versions and check the console log if a function signature has changed.

```r
install.packages(c("dplyr", "ggplot2", "pheatmap", "msigdbr"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot",
                       "GSVA", "GSEABase", "limma", "edgeR"))

install.packages("Seurat")   # Seurat 5.x
```

`msigdbr` is used through its `collection` / `subcollection` interface (msigdbr ≥ 7.5), i.e. `msigdbr(species, collection = "H")` rather than the deprecated `category` argument.

## 4. Inputs

**`sigDiff.txt`** — table of the significant genes obtained after virtual overexpression of HLA-DOA (tab-separated, header required, **first column = gene symbol**). Only the first column is used.

**`seurat.Rdata`** — a Seurat object that must be named `pbmc` in the workspace (the script calls `stopifnot(exists("pbmc"))`). It is expected to contain:

- the macrophage subset only (`cellType == "macrophages"` is applied again inside the script),
- a `Type` metadata column with `Control` / `Treat`,
- an `orig.ident` column identifying the original sample (used as the block for `duplicateCorrelation`),
- normalized expression in the `RNA` assay `data` layer.

Both inputs and the `ORA_out/` directory are referenced by the paths set in the *Settings* block at the top of each script; edit them there rather than inside the code body.

## 5. Running

```r
# step 1
Rscript ORA_diffgenes_all.R

# step 2 (reads ORA_out/, writes GSVA_out/)
Rscript GSVA_direction_from_ORA.R
```

Both scripts are also runnable interactively in RStudio; `set.seed(123)` is fixed in step 2.

## 6. Outputs

`ORA_out/`

- `ORA_GO.txt`, `ORA_KEGG.txt`, `ORA_H.txt`, `ORA_C7.txt`, `ORA_C2Reactome.txt` — enrichment tables (nominal *P* < 0.05)
- `ORA_<collection>.barplot.pdf`, `ORA_<collection>.bubble.pdf` — bar and dot plots

`GSVA_out/`

- `GSVA_scores_<collection>.csv` — raw GSVA score matrix
- `GSVA_direction_<collection>.csv` — logFC, *P*, adjusted *P* and direction for every gene set (`logFC > 0` = more active in the HLA-DOA-high group)
- `GSVA_direction_barplot_<collection>.pdf/png` — direction bar plot with significance marks
- `GSVA_heatmap_<collection>.pdf/png` — GSVA scores across pseudobulk samples (row-scaled)

## 7. Statistical design of the direction analysis

Cells were split at the median HLA-DOA expression within macrophages and aggregated into pseudobulk profiles by `Type × sample × HLA-DOA group`, giving 8 pseudobulk samples (4 disease samples × high/low). Profiles were converted to log2-CPM over all genes. Each significant ORA pathway contributed a gene module formed by its diff-gene members (≥ 5 genes retained); GSVA scores were compared between groups with `limma`, using `duplicateCorrelation` blocked by original sample, and adjusted for the sample group (`~ Type + HLA-DOA_group`) when both groups were present.

Two points should be kept in mind when reading the direction results:

1. ORA is filtered on the nominal *P* value. **None of the pathways passed multiple-testing correction**, so the direction results are reported as exploratory trends, not as confirmed activation.
2. When an ORA collection yields no usable pathway (C7 in this dataset), the script falls back to the corresponding pre-defined MSigDB sets and prints the switch to the console; the fallback is flagged in the log, not in the output file name.

## 8. Data

All data analysed with these scripts are public:

- single-cell RNA-seq: GEO **GSE208653**, **GSE197461**
- spatial transcriptomics: CROST **VISDS000672**, **VISDS000674**
- immune-cell eQTL (exposure): OneK1K (Yazar et al., *Science* 2022)
- cervical cancer GWAS (outcome): BioBank Japan **bbj-a-98** (Ishigaki et al., *Nat Genet* 2020)

## 9. Citation

If you use these scripts, please cite the manuscript above (full reference to be added upon publication) and the methods they depend on: `clusterProfiler` (Wu et al., *Innovation* 2021), `GSVA` (Hänzelmann et al., *BMC Bioinformatics* 2013), `limma` (Ritchie et al., *Nucleic Acids Res* 2015), `Seurat` (Hao et al., *Cell* 2021), MSigDB (Subramanian et al., *PNAS* 2005).

## 10. Contact

Miao Shengdi, Hunan Normal University — miaoshengdi@hunnu.edu.cn
