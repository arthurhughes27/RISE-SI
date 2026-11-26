# Screen of all gene sets (BTMs) at day 7: mean expression vs surrogate index (ridge) as trial-level surrogates
# For each gene set: mean of the genes, and ridge SI cross-fitted leave-one-study-out (as in the M156.0 example),
# then paired RISE-Meta for each surrogate type over all gene sets.
# The surrogates are slow to compute: they are saved on the first run (output/application/screen_surrogates_day7.rds)
# and loaded on later runs, so the RISE-Meta settings can be changed and re-run quickly.
# Endpoint Y: mean of the log2 titres against H1, H3 and B, before vaccination and at day 28.
# Input: data/analysis/analysis_set_endpoint.rds (from analysis/descriptive/descriptive.R). Run from the project root.
# Output: output/application/screen_all_genesets_day7.rds (RISE-Meta results, read by the figure scripts)

library(tidyverse)
library(Biobase)
library(glmnet)
library(SurrogateRank)

ge_day <- "Day 7"
min_genes <- 5          # gene sets need at least this many genes measured in every study

# RISE-Meta settings (change and re-run; the surrogates are not recomputed)
epsilon <- 0.2          # equivalence margin (pooled and within-study)
test_target <- "pi"     # "pi": prediction interval, "ci": confidence interval
p_correction <- "BH"    # multiplicity correction across gene sets ("none", "BH", ...)
test <- "knha"          # variance estimation of the pooled effect

surr_file <- fs::path("output", "application", "screen_surrogates_day7.rds")

# ---- Surrogates: computed once, then loaded ---------------------------------------------------
if (!file.exists(surr_file)) {
  set.seed(1)

  # ---- Endpoint: one row per individual -----------------------------------------------------------
  ind <- readRDS(fs::path("data", "analysis", "analysis_set_endpoint.rds")) %>%
    filter(ge_day == !!ge_day) %>%
    group_by(participant_id, study, ge_pre_day, ge_post_day) %>%
    summarise(y_pre = mean(log_value_pre), y_post = mean(log_value_post), .groups = "drop")
  n <- nrow(ind)
  pre <- 1:n
  post <- n + 1:n

  # ---- Expression of all genes: rows are all pre samples, then all post samples -------------------
  eset <- readRDS(fs::path("data-raw", "is2", "all_noNorm_eset.rds"))
  key <- paste(eset$participant_id, round(eset$study_time_collected, 2))
  cols <- match(c(paste(ind$participant_id, ind$ge_pre_day), paste(ind$participant_id, ind$ge_post_day)), key)
  E <- t(exprs(eset)[, cols])
  rm(eset)
  E <- E[, colSums(is.na(E)) == 0]     # genes measured in every study

  # Centre and scale genes and Y within study, both time points pooled (removes platform / assay offsets;
  # does not change anyone's direction of change)
  study <- rep(ind$study, 2)
  y <- c(ind$y_pre, ind$y_post)
  for (k in unique(study)) {
    r <- study == k
    E[r, ] <- scale(E[r, ])
    y[r] <- as.numeric(scale(y[r]))
  }
  E[is.na(E)] <- 0                      # genes constant within a study

  # ---- Gene sets: BTMs (excluding unannotated "TBA" modules) with at least min_genes measured genes --
  gmt <- strsplit(readLines(fs::path("data-raw", "BTM_for_GSEA_20131008.gmt")), "\t")
  gene_sets <- setNames(lapply(gmt, function(g) intersect(g[-(1:2)], colnames(E))), make.unique(sapply(gmt, `[`, 1)))
  gene_sets <- gene_sets[!grepl("^TBA", names(gene_sets)) & lengths(gene_sets) >= min_genes]
  length(gene_sets)

  # ---- Surrogates: mean and ridge SI for each gene set ------------------------------------------
  # SI: ridge regression of Y on the genes, cross-fitted leave-one-study-out (inner CV folds by individual)
  fold <- rep(sample(rep(1:5, length.out = n)), 2)
  ridge_si <- function(S) {
    si <- rep(NA, 2 * n)
    for (k in unique(study)) {
      tr <- study != k
      fit <- cv.glmnet(S[tr, ], y[tr], alpha = 0, foldid = fold[tr])
      si[!tr] <- predict(fit, S[!tr, ], s = "lambda.min")
    }
    si
  }
  surr_mean <- sapply(gene_sets, function(g) rowMeans(E[, g]))
  surr_si <- sapply(seq_along(gene_sets), function(i) {
    if (i %% 25 == 0) cat(i, "/", length(gene_sets), "gene sets\n")
    ridge_si(E[, gene_sets[[i]]])
  })
  colnames(surr_si) <- names(gene_sets)

  fs::dir_create(fs::path("output", "application"))
  saveRDS(list(ind = ind, surr = list(Mean = surr_mean, `SI (ridge)` = surr_si), n_genes = lengths(gene_sets)),
          surr_file)
}
saved <- readRDS(surr_file)   # endpoint, surrogates (rows: pre samples, then post samples), gene set sizes
ind <- saved$ind
n <- nrow(ind)
pre <- 1:n
post <- n + 1:n

# ---- Paired RISE-Meta, one call per surrogate type ---------------------------------------------
# Studies with a zero within-study SE for a gene set are dropped by RISE-Meta for that gene set only.
rise_args <- list(yone = ind$y_post, yzero = ind$y_pre, studyone = ind$study, studyzero = ind$study,
                  paired.all = TRUE, epsilon.meta.mode = "user", epsilon.meta = epsilon, epsilon.study = epsilon,
                  test = test, test.target = test_target, meta.analysis.method = "RE", p.correction = p_correction)
run_rise <- function(surr, screen_plot = FALSE) {
  capture.output(res <- do.call(rise.screen.meta, c(rise_args, list(
    sone = surr[post, ], szero = surr[pre, ],
    return.screen.plot = screen_plot, screen.plot.topN = 20, return.forest.plot = FALSE, return.fit.plot = FALSE,
    return.evaluate.results = screen_plot,   # if FALSE, RISE-Meta drops the screening plot when a gene set is significant
    return.study.similarity.plot = FALSE))))
  res
}
rise <- list(Mean = run_rise(saved$surr$Mean), `SI (ridge)` = run_rise(saved$surr$`SI (ridge)`, screen_plot = TRUE))

# Fit and forest plots for the top gene set of the SI (smallest unadjusted p-value), on the same (training) data
top <- rise$`SI (ridge)`$screening.metrics.meta %>% slice_min(p.unadjusted, n = 1, with_ties = FALSE) %>% pull(marker)
top
capture.output(top_eval <- do.call(rise.evaluate.meta, c(rise_args, list(
  sone = saved$surr$`SI (ridge)`[post, top, drop = FALSE], szero = saved$surr$`SI (ridge)`[pre, top, drop = FALSE],
  markers = top, evaluate.weights = FALSE, return.fit.plot = TRUE, return.forest.plot = TRUE))))

# ---- RISE-Meta figures for the SI: printed, then saved ------------------------------------------
rise_plots <- list(screen.plot = rise$`SI (ridge)`$gamma.s.plot$screen.plot,
                   fit.plot = top_eval$gamma.s.plot$fit.plot, forest.plot = top_eval$gamma.s.plot$forest.plot)
# Adapt RISE-Meta's default plots to the document page. The screening plot uses a base font size of 20, which
# inflates the legends and squeezes the panel: reset it to a normal theme (titles described in the captions).
rise_plots$screen.plot <- rise_plots$screen.plot + labs(title = NULL, subtitle = NULL) +
  scale_y_discrete(labels = scales::label_wrap(40)) +   # wrap long gene set names
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(colour = "black"))
rise_plots$fit.plot <- rise_plots$fit.plot + labs(title = NULL, x = "Treatment effect on SI", y = "Treatment effect on endpoint")
rise_plots$screen.plot
rise_plots$fit.plot
rise_plots$forest.plot

plot_sizes <- list(screen.plot = c(6.8, 4.8), fit.plot = c(8, 6), forest.plot = c(12.5, 6))   # width, height (in); the forest plot layout needs a wide canvas
for (d in c(fs::path("output", "figures", "application"), fs::path("manuscript", "latex", "figures"))) {
  fs::dir_create(d)
  for (pl in names(plot_sizes)) {
    if (!is.null(rise_plots[[pl]])) {
      ggsave(fs::path(d, paste0("screen_si_", sub(".plot", "", pl, fixed = TRUE), ".pdf")), rise_plots[[pl]],
             width = plot_sizes[[pl]][1], height = plot_sizes[[pl]][2])
    }
  }
}

results <- map_dfr(rise, ~ as_tibble(.x$screening.metrics.meta), .id = "surrogate") %>%
  transmute(geneset = marker, surrogate, n_genes = saved$n_genes[geneset], n_studies = n.studies,
            delta = mu.delta, lower = ci.delta.lower, upper = ci.delta.upper,
            pi_lower = pi.delta.lower, pi_upper = pi.delta.upper, tau2, p_adjusted = p.adjusted)

results %>% group_by(surrogate) %>% summarise(n_genesets = n(), n_pass = sum(p_adjusted < 0.05))

saveRDS(results, fs::path("output", "application", "screen_all_genesets_day7.rds"))
