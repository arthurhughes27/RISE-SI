# Helpers to load the processed data (data/, built by analysis/preprocessing/make_analysis_data.R)
# Source after loading the tidyverse.

# Endpoint at one GE day (e.g. "Day 7"): one row per individual.
# y_pre, y_post: mean of the log2 titres against H1, H3 and B, before vaccination and at day 28.
load_endpoint <- function(ge_day) {
  readRDS(fs::path("data", "analysis_set.rds"))$endpoint %>%
    filter(ge_day == !!ge_day) %>%
    group_by(participant_id, study, age, sex, ge_pre_day, ge_post_day) %>%
    summarise(y_pre = mean(log_value_pre), y_post = mean(log_value_post), .groups = "drop")
}

# Expression of the individuals in `ind` (from load_endpoint): rows are all pre samples, then all post samples.
# Columns: the genes measured in every study (all genes, or only those in `genes`).
load_expr <- function(ind, genes = NULL) {
  expr <- readRDS(fs::path("data", "expr.rds"))
  if (!is.null(genes)) expr <- expr[intersect(genes, rownames(expr)), , drop = FALSE]
  E <- t(expr[, c(paste(ind$participant_id, ind$ge_pre_day), paste(ind$participant_id, ind$ge_post_day)), drop = FALSE])
  E[, colSums(is.na(E)) == 0, drop = FALSE]
}

# Baseline covariates of the individuals in `ind`, rows as in load_expr (identical at both time points):
# age (years) and sexMale (1 = male, 0 = female; female is the reference). Every SI includes them as main effects.
load_covariates <- function(ind) {
  cbind(age = rep(ind$age, 2), sexMale = rep(as.numeric(ind$sex == "Male"), 2))
}

# BTM gene sets: list of gene symbols, named by module (e.g. "plasma cells & B cells, immunoglobulins (M156.0)")
load_btm <- function() {
  gmt <- strsplit(readLines(fs::path("data-raw", "BTM_for_GSEA_20131008.gmt")), "\t")
  setNames(lapply(gmt, `[`, -(1:2)), make.unique(sapply(gmt, `[`, 1)))
}

# Genes of one BTM, by its identifier (e.g. "M156.0")
btm_genes <- function(id) {
  btm <- load_btm()
  btm[[grep(paste0("(", id, ")"), names(btm), fixed = TRUE)]]
}
