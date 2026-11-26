# Screen of all gene sets (BTMs) at one gene expression (GE) day, with the ridge surrogate index (SI).
# Sourced by screen_day1.R, screen_day2.R, screen_day3.R and screen_day1_2.R, which set `ge_day` and `tag`.
#   1. SI (ridge, cross-fitted leave-one-study-out) for every gene set; saved on the first run, loaded afterwards
#   2. Paired RISE-Meta on the SIs: screening plot; fit and forest plots for the top gene set
#   3. Top gene set: individual predictions (SI vs observed) and pooled delta of mean, median, max and SI
# Endpoint Y: mean of the log2 titres against H1, H3 and B, before vaccination and at day 28.
# Input: data/analysis/analysis_set_endpoint.rds (from analysis/descriptive/descriptive.R). Run from the project root.

library(tidyverse)
library(Biobase)
library(glmnet)
library(SurrogateRank)

min_genes <- 5          # gene sets need at least this many genes measured in every study
epsilon <- 0.2          # equivalence margin (pooled and within-study)
test_target <- "pi"     # "pi": prediction interval, "ci": confidence interval
p_correction <- "BH"    # multiplicity correction across gene sets in the screen

surr_file <- fs::path("output", "application", paste0("screen_surrogates_", tag, ".rds"))
fig_dir <- fs::path("output", "figures", "application", tag)
fs::dir_create(c(fs::path("output", "application"), fig_dir))

# ---- Endpoint: one row per individual -----------------------------------------------------------
ind <- readRDS(fs::path("data", "analysis", "analysis_set_endpoint.rds")) %>%
  filter(ge_day == !!ge_day) %>%
  group_by(participant_id, study, ge_pre_day, ge_post_day) %>%
  summarise(y_pre = mean(log_value_pre), y_post = mean(log_value_post), .groups = "drop")
n <- nrow(ind)
if (n == 0) stop(ge_day, " is not in the analysis set: re-run analysis/descriptive/descriptive.R")
pre <- 1:n
post <- n + 1:n
cat(ge_day, ":", n, "individuals in", n_distinct(ind$study), "studies\n")

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

# ---- 1. SI for every gene set (computed once, then loaded) --------------------------------------
# Ridge regression of Y on the genes, cross-fitted leave-one-study-out (inner CV folds by individual)
set.seed(1)
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
if (!file.exists(surr_file) || nrow(readRDS(surr_file)) != 2 * n) {   # recompute if missing or out of date
  si_all <- sapply(gene_sets, function(g) ridge_si(E[, g]))
  saveRDS(si_all, surr_file)
}
si_all <- readRDS(surr_file)

# ---- 2. Paired RISE-Meta on the SIs -------------------------------------------------------------
rise_args <- list(yone = ind$y_post, yzero = ind$y_pre, studyone = ind$study, studyzero = ind$study,
                  paired.all = TRUE, epsilon.meta.mode = "user", epsilon.meta = epsilon, epsilon.study = epsilon,
                  test = "knha", test.target = test_target, meta.analysis.method = "RE")
capture.output(screen <- do.call(rise.screen.meta, c(rise_args, list(
  sone = si_all[post, ], szero = si_all[pre, ], p.correction = p_correction,
  return.screen.plot = TRUE, screen.plot.topN = 20, return.forest.plot = FALSE, return.fit.plot = FALSE,
  return.evaluate.results = TRUE,   # if FALSE, RISE-Meta drops the screening plot whenever a gene set is significant
  return.study.similarity.plot = FALSE))))
screen$screening.metrics.meta %>% arrange(p.unadjusted) %>%
  select(marker, n.studies, mu.delta, pi.delta.lower, pi.delta.upper, p.unadjusted, p.adjusted) %>% head(10)

# Top gene set (smallest unadjusted p-value): fit and forest plots on the same data
top <- screen$screening.metrics.meta %>% slice_min(p.unadjusted, n = 1, with_ties = FALSE) %>% pull(marker)
top
capture.output(top_eval <- do.call(rise.evaluate.meta, c(rise_args, list(
  sone = si_all[post, top, drop = FALSE], szero = si_all[pre, top, drop = FALSE], p.correction = "none",
  markers = top, evaluate.weights = FALSE, return.fit.plot = TRUE, return.forest.plot = TRUE))))

# RISE-Meta's own plots, adapted to page size (its screening plot uses a base font size of 20)
p_screen <- screen$gamma.s.plot$screen.plot + labs(title = NULL, subtitle = NULL) +
  scale_y_discrete(labels = scales::label_wrap(40)) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(colour = "black"))
p_fit <- top_eval$gamma.s.plot$fit.plot + labs(title = NULL, x = "Treatment effect on SI", y = "Treatment effect on endpoint")
p_forest <- top_eval$gamma.s.plot$forest.plot
p_screen
p_fit
p_forest
ggsave(fs::path(fig_dir, "screen_si_screen.pdf"), p_screen, width = 6.8, height = 4.8)
ggsave(fs::path(fig_dir, "screen_si_fit.pdf"), p_fit, width = 8, height = 6)
ggsave(fs::path(fig_dir, "screen_si_forest.pdf"), p_forest, width = 12.5, height = 6)

# ---- 3a. Top gene set: individual predictions -------------------------------------------------
pred <- tibble(study, time = rep(c("Day 0", "Day 28"), each = n), y, si = si_all[, top])
r2 <- pred %>% group_by(study) %>% summarise(label = sprintf("R^2 == %.2f", cor(si, y)^2))   # within-study R2
lims <- range(c(pred$y, pred$si))
p_pred <- ggplot(pred, aes(y, si, colour = time)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +   # perfect prediction
  geom_point(size = 1, alpha = 0.6) +
  geom_text(data = r2, aes(x = lims[1], y = lims[2], label = label), parse = TRUE,
            hjust = 0, vjust = 1, size = 2.8, inherit.aes = FALSE) +
  scale_colour_manual(values = c(`Day 0` = "#257de0", `Day 28` = "#f47230"), name = NULL) +
  facet_wrap(~study, ncol = 5) +
  coord_equal(xlim = lims, ylim = lims) +
  labs(x = "Observed antibody response", y = "Predicted antibody response (SI)",
       title = paste0("Cross-fitted SI predictions in each held-out study (", top, ", ", tolower(ge_day), ")")) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
        axis.text = element_text(colour = "black"), legend.position = "right")
p_pred
ggsave(fs::path(fig_dir, "top_geneset_prediction.pdf"), p_pred, width = 6.3, height = 4.4)

# ---- 3b. Top gene set: mean, median, max vs SI ------------------------------------------------
S <- E[, gene_sets[[top]]]
surr <- cbind(Mean = rowMeans(S), Median = apply(S, 1, median), Max = apply(S, 1, max), `SI (ridge)` = si_all[, top])

# RISE-Meta drops a study for a surrogate when its within-study SE is zero; drop such studies for all four
y_up <- ind$y_post > ind$y_pre
zero_se <- unique(unlist(lapply(colnames(surr), function(s) {
  d <- y_up - (surr[post, s] > surr[pre, s])
  names(which(tapply(d, ind$study, var) == 0))
})))
zero_se
keep <- !ind$study %in% zero_se

args_keep <- modifyList(rise_args, list(yone = ind$y_post[keep], yzero = ind$y_pre[keep],
                                        studyone = ind$study[keep], studyzero = ind$study[keep]))
capture.output(res <- do.call(rise.screen.meta, c(args_keep, list(
  sone = surr[post, ][keep, ], szero = surr[pre, ][keep, ], p.correction = "none",
  return.screen.plot = FALSE, return.forest.plot = FALSE, return.fit.plot = FALSE,
  return.study.similarity.plot = FALSE, return.evaluate.results = FALSE))))

pooled <- as_tibble(res$screening.metrics.meta) %>%
  transmute(surrogate = factor(marker, colnames(surr)), delta = mu.delta, lower = ci.delta.lower,
            upper = ci.delta.upper, pi_lower = pi.delta.lower, pi_upper = pi.delta.upper, tau2)
pooled

p_pooled <- ggplot(pooled, aes(delta, fct_rev(surrogate))) +     # thick = CI, thin = prediction interval
  annotate("rect", xmin = -epsilon, xmax = epsilon, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_linerange(aes(xmin = pi_lower, xmax = pi_upper)) +
  geom_linerange(aes(xmin = lower, xmax = upper), linewidth = 1.2) +
  geom_point(shape = 18, size = 2.5) +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(x = expression("Pooled " * delta), y = NULL, title = paste0(top, ", ", tolower(ge_day))) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text = element_text(colour = "black"))
p_pooled
ggsave(fs::path(fig_dir, "top_geneset_aggregation_pooled.pdf"), p_pooled, width = 4.8, height = 2.2)
