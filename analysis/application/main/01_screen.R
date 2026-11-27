# Main analysis: screen of all gene sets (BTMs) at day 7 with the ridge surrogate index (SI), and the mean of the
# genes as a reference.
#   1. SI (cross-fitted leave-one-study-out) for every gene set; saved on the first run, loaded afterwards
#   2. Paired RISE-Meta on the SIs: screening plot; fit and forest plots for the top gene set
#   3. Pooled delta of the mean vs the SI, one point per gene set
#   4. Top gene set: individual predictions (SI vs observed) and pooled delta of mean, median, max and SI
# The supplementary screens (analysis/application/supplementary/) source this script after setting
#   `ge_day`     GE day, as in data/analysis_set.rds (default "Day 7")
#   `tag`        name of the outputs: output/figures/application/<tag>/, output/application/screen_surrogates_<tag>.rds
#   `learner`    "ridge" (default), "pls" (partial least squares, 1-3 components chosen by CV) or
#                "spls" (sparse PLS, components and sparsity chosen by CV)
#   `covariates` TRUE to add age and sex as baseline predictors (default FALSE)
# The figures of the main analysis (day 7, ridge, no covariates) are also saved in manuscript/latex/figures/.
# Endpoint Y: mean of the log2 titres against H1, H3 and B, before vaccination and at day 28.
# Input: data/ (from analysis/preprocessing/make_analysis_data.R). Run from the project root.

library(tidyverse)
library(glmnet)
library(pls)
library(SurrogateRank)
source(fs::path("R", "load_data.R"))

if (!exists("ge_day")) ge_day <- "Day 7"
if (!exists("tag")) tag <- "day7"
if (!exists("learner")) learner <- "ridge"
if (!exists("covariates")) covariates <- FALSE
si_label <- paste0("SI (", c(ridge = "ridge", pls = "PLS", spls = "sparse PLS")[[learner]], if (covariates) " + age, sex", ")")

min_genes <- 5          # gene sets need at least this many genes measured in every study
epsilon <- 0.2          # equivalence margin (pooled and within-study)
test_target <- "pi"     # "pi": prediction interval, "ci": confidence interval
p_correction <- "BH"    # multiplicity correction across gene sets in the screen
alpha <- 0.05

surr_file <- fs::path("output", "application", paste0("screen_surrogates_", tag, ".rds"))
fig_dir <- fs::path("output", "figures", "application", tag)
fs::dir_create(c(fs::path("output", "application"), fig_dir))
main <- ge_day == "Day 7" && learner == "ridge" && !covariates
save_fig <- function(p, file, w, h, doc = FALSE) {   # doc: figure of the plan document
  for (d in c(fig_dir, if (doc && main) fs::path("manuscript", "latex", "figures")))
    ggsave(fs::path(d, file), p, width = w, height = h)
}

# ---- Endpoint: one row per individual -----------------------------------------------------------
ind <- load_endpoint(ge_day)
n <- nrow(ind)
if (n == 0) stop(ge_day, " is not in the analysis set: re-run analysis/preprocessing/make_analysis_data.R")
pre <- 1:n
post <- n + 1:n
cat(ge_day, ":", n, "individuals in", n_distinct(ind$study), "studies\n")

# ---- Expression of all genes: rows are all pre samples, then all post samples -------------------
E <- load_expr(ind)     # genes measured in every study

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

# Baseline covariates (identical at both time points): age (centred and scaled) and male sex
X_cov <- if (covariates) cbind(age = as.numeric(scale(rep(ind$age, 2))), male = rep(as.numeric(ind$sex == "Male"), 2))

# ---- Gene sets: BTMs (excluding unannotated "TBA" modules) with at least min_genes measured genes --
gene_sets <- lapply(load_btm(), intersect, colnames(E))
gene_sets <- gene_sets[!grepl("^TBA", names(gene_sets)) & lengths(gene_sets) >= min_genes]
length(gene_sets)

# ---- 1. SI for every gene set (computed once, then loaded) --------------------------------------
# Regression of Y on the genes (and covariates), cross-fitted leave-one-study-out. Inner CV folds (ridge penalty,
# number of PLS components) are by individual. Additive covariates cancel in the paired comparison; they can
# only act by changing the gene weights.
set.seed(1)
fold <- rep(sample(rep(1:5, length.out = n)), 2)
# Sparse PLS (spls package): number of components K and sparsity eta chosen by CV over the given folds
spls_tuned <- function(x, y, newx, fold, grid = expand.grid(K = 1:2, eta = c(0.3, 0.6, 0.9))) {
  grid$K <- pmin(grid$K, ncol(x))
  cv_err <- sapply(seq_len(nrow(grid)), function(g) {
    sum(sapply(unique(fold), function(f) {
      fit <- spls::spls(x[fold != f, ], y[fold != f], K = grid$K[g], eta = grid$eta[g], trace = FALSE)
      sum((y[fold == f] - predict(fit, x[fold == f, , drop = FALSE]))^2)
    }))
  })
  best <- grid[which.min(cv_err), ]
  fit <- spls::spls(x, y, K = best$K, eta = best$eta, trace = FALSE)
  as.numeric(predict(fit, newx))
}
fit_predict <- function(X, tr) {
  if (learner == "ridge") {
    pf <- as.numeric(!colnames(X) %in% c("age", "male"))          # covariates unpenalised
    fit <- cv.glmnet(X[tr, ], y[tr], alpha = 0, foldid = fold[tr], penalty.factor = pf)
    as.numeric(predict(fit, X[!tr, ], s = "lambda.min"))
  } else if (learner == "spls") {
    spls_tuned(X[tr, ], y[tr], X[!tr, ], fold[tr])
  } else {
    fit <- plsr(y ~ ., data = data.frame(y = y[tr], X[tr, ]), ncomp = min(3, ncol(X)), validation = "CV",
                segments = split(seq_len(sum(tr)), fold[tr]))
    ncomp <- which.min(RMSEP(fit, estimate = "CV", intercept = FALSE)$val[1, 1, ])
    as.numeric(predict(fit, newdata = data.frame(X[!tr, ]), ncomp = ncomp))
  }
}
cross_fit_si <- function(S) {
  X <- cbind(S, X_cov)
  si <- rep(NA, 2 * n)
  for (k in unique(study)) si[study == k] <- fit_predict(X, study != k)
  si
}
saved <- if (file.exists(surr_file)) readRDS(surr_file)
if (!is.list(saved) || !identical(saved$ind, ind)) {   # recompute if missing or if the analysis set has changed
  saved <- list(ind = ind, si = sapply(gene_sets, function(g) cross_fit_si(E[, g])))
  saveRDS(saved, surr_file)
}
si_all <- saved$si
surr_mean <- sapply(gene_sets, function(g) rowMeans(E[, g]))

# ---- 2. Paired RISE-Meta on the SIs (and on the means, for 3) -----------------------------------
# Studies with a zero within-study SE for a gene set are dropped by RISE-Meta for that gene set only.
rise_args <- list(yone = ind$y_post, yzero = ind$y_pre, studyone = ind$study, studyzero = ind$study,
                  paired.all = TRUE, epsilon.meta.mode = "user", epsilon.meta = epsilon, epsilon.study = epsilon,
                  test = "knha", test.target = test_target, meta.analysis.method = "RE")
run_screen <- function(surr, screen_plot = FALSE) {
  capture.output(res <- do.call(rise.screen.meta, c(rise_args, list(
    sone = surr[post, ], szero = surr[pre, ], p.correction = p_correction,
    return.screen.plot = screen_plot, screen.plot.topN = 20, return.forest.plot = FALSE, return.fit.plot = FALSE,
    return.evaluate.results = screen_plot,   # if FALSE, RISE-Meta drops the screening plot whenever a gene set is significant
    return.study.similarity.plot = FALSE))))
  res
}
screen <- run_screen(si_all, screen_plot = TRUE)
screen_mean <- run_screen(surr_mean)
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
p_fit <- top_eval$gamma.s.plot$fit.plot + labs(title = NULL, x = paste("Treatment effect on", si_label), y = "Treatment effect on endpoint")
p_forest <- top_eval$gamma.s.plot$forest.plot
p_screen
p_fit
p_forest
save_fig(p_screen, "screen_si_screen.pdf", 6.8, 4.8, doc = TRUE)
save_fig(p_fit, "screen_si_fit.pdf", 8, 6, doc = TRUE)
save_fig(p_forest, "screen_si_forest.pdf", 12.5, 6, doc = TRUE)   # the forest plot layout needs a wide canvas

# ---- 3. Pooled delta of the mean vs the SI, one point per gene set ------------------------------
wide <- bind_rows(mean = as_tibble(screen_mean$screening.metrics.meta),
                  si = as_tibble(screen$screening.metrics.meta), .id = "surrogate") %>%
  select(geneset = marker, surrogate, delta = mu.delta, p_adjusted = p.adjusted) %>%
  pivot_wider(names_from = surrogate, values_from = c(delta, p_adjusted)) %>%
  mutate(pass = case_when(p_adjusted_mean < alpha & p_adjusted_si < alpha ~ "Both",
                          p_adjusted_mean < alpha ~ "Mean only",
                          p_adjusted_si < alpha ~ "SI only",
                          TRUE ~ "Neither"),
         pass = factor(pass, c("Neither", "Mean only", "SI only", "Both")))
count(wide, pass)

p_mean_si <- ggplot(wide, aes(delta_mean, delta_si)) +
  annotate("rect", xmin = -epsilon, xmax = epsilon, ymin = -Inf, ymax = Inf, fill = "grey93") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -epsilon, ymax = epsilon, fill = "grey93") +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") +
  geom_point(aes(colour = pass), size = 1.3) +
  scale_colour_manual(values = c(Neither = "grey65", `Mean only` = "#f47230", `SI only` = "#257de0", Both = "black"),
                      name = "Passes (PI test, BH)", drop = FALSE) +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(x = expression("Pooled " * delta * ", mean"), y = bquote("Pooled " * delta * ", " * .(si_label))) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "black"),
        legend.position = "right")
p_mean_si
save_fig(p_mean_si, "screen_mean_vs_si.pdf", 4.8, 3.6, doc = TRUE)

# ---- 4a. Top gene set: individual predictions -------------------------------------------------
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
  labs(x = "Observed antibody response", y = paste("Predicted antibody response,", si_label),
       title = paste0("Cross-fitted SI predictions in each held-out study (", top, ", ", tolower(ge_day), ")")) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
        axis.text = element_text(colour = "black"), legend.position = "right")
p_pred
save_fig(p_pred, "top_geneset_prediction.pdf", 6.3, 4.4)

# ---- 4b. Top gene set: mean, median, max vs SI ------------------------------------------------
S <- E[, gene_sets[[top]]]
surr <- cbind(Mean = rowMeans(S), Median = apply(S, 1, median), Max = apply(S, 1, max), si = si_all[, top])
colnames(surr)[4] <- si_label

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
save_fig(p_pooled, "top_geneset_aggregation_pooled.pdf", 4.8, 2.2)

rm(list = ls())
