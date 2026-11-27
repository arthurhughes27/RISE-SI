# RISE-SI
Repository for files for a project applying surrogate index methods prior to using nonparametric surrogate evaluation approaches on transcriptomic data.

## Data
The data are not in the repository due to data protection. Put the raw files here (folders `data-raw/` and `data-raw/is2/`):

```
data-raw/
  BTM_for_GSEA_20131008.gmt                   blood transcription modules (Li et al. 2014)
  is2/
    all_noNorm_eset.rds                       Immune Signatures 2 gene expression (ExpressionSet)
    hai_2025-01-10_01-13-41.xlsx              haemagglutination inhibition titres
    neut_ab_titer_2025-01-10_01-13-22.xlsx    neutralising antibody titres
```

Then, open the Rproj. project file in RStudio and build the processed data by running:

```r
source("analysis/preprocessing/make_analysis_data.R")
```

This writes `data/analysis_set.rds` (containing clinical endpoint + study info) and `data/expr.rds` (expression data). All analysis scripts read these files through the
helpers in `R/load_data.R`.

## Analysis
Run every script from the project root. To run the whole main analysis (preprocessing, data description and
`analysis/application/main/`) in order:

```r
source("analysis/run_main_analysis.R")
```

- `analysis/descriptive/descriptive.R`: data description figures.
- `analysis/application/main/`: main analysis (gene expression at day 7).
  `01_screen.R` screens all gene sets; `02`–`07` study the top gene set, M156.0.
- `analysis/application/supplementary/`: other days (`other_days/`), other learners (`learners/`) and
  experiments (`experiments/`).

Figures are written to `output/figures/`, and those used in the latex plan document are also written to `manuscript/latex/figures/`.

## R packages
CRAN: tidyverse, fs, readxl, janitor, glmnet, pls, spls, ranger, xgboost, SuperLearner, scales.
Bioconductor: Biobase (preprocessing only). SurrogateRank (RISE-Meta).
