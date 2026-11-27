# Experiment: the M156.0 example without cross-fitting
# The ridge SI is fitted once on all studies and applied to the same studies (in-sample), instead of being
# cross-fitted leave-one-study-out. Shows (1) the individual-level predictions and (2) the comparison of simple
# aggregations with the SI in RISE-Meta, with the cross-fitted SI kept as a reference.
# Gene set M156.0, gene expression at day 7. Endpoint Y: mean of the log2 titres against H1, H3 and B.
# Input: data/ (from analysis/preprocessing/make_analysis_data.R). Run from the project root.

library(tidyverse)
library(glmnet)
library(SurrogateRank)
source(fs::path("R", "load_data.R"))

geneset_id <- "M156.0"
ge_day <- "Day 7"
epsilon <- 0.2     # equivalence margin (pooled and within-study)
set.seed(1)

# ---- Endpoint: one row per individual -----------------------------------------------------------
ind <- load_endpoint(ge_day)
n <- nrow(ind)
pre <- 1:n
post <- n + 1:n

# ---- Gene set expression: rows are all pre samples, then all post samples -----------------------
genes <- btm_genes(geneset_id)
S <- load_expr(ind, genes)     # genes measured in every study

# Centre and scale genes and Y within study, both time points pooled
study <- rep(ind$study, 2)
y <- c(ind$y_pre, ind$y_post)
for (k in unique(study)) {
  r <- study == k
  S[r, ] <- scale(S[r, ])
  y[r] <- as.numeric(scale(y[r]))
}
S[is.na(S)] <- 0                      # genes constant within a study
X <- cbind(S, load_covariates(ind))     # predictors of the SI: genes and baseline covariates (age, sexMale)

# ---- Surrogates ---------------------------------------------------------------------------------
# Inner CV folds (for the ridge penalty) are by individual, so an individual's two samples stay together.
fold <- rep(sample(rep(1:5, length.out = n)), 2)
surr <- cbind(Mean = rowMeans(S), Median = apply(S, 1, median), Max = apply(S, 1, max),
              `SI (ridge, cross-fitted)` = NA, `SI (ridge, no cross-fitting)` = NA)

# Cross-fitted SI (as in the example): the SI of study k comes from a model trained on the other studies
for (k in unique(study)) {
  tr <- study != k
  fit <- cv.glmnet(X[tr, ], y[tr], alpha = 0, foldid = fold[tr])
  surr[!tr, "SI (ridge, cross-fitted)"] <- predict(fit, X[!tr, ], s = "lambda.min")
}

# No cross-fitting: one model trained on all studies, applied to the same studies
fit_all <- cv.glmnet(X, y, alpha = 0, foldid = fold)
surr[, "SI (ridge, no cross-fitting)"] <- predict(fit_all, X, s = "lambda.min")

# ---- Figure 1: in-sample predictions vs observed, one panel per study --------------------------
pred <- tibble(study, time = rep(c("Day 0", "Day 28"), each = n), y, si = surr[, "SI (ridge, no cross-fitting)"])
r2 <- pred %>% group_by(study) %>% summarise(label = sprintf("R^2 == %.2f", cor(si, y)^2))   # within-study R2
r2

lims <- range(c(pred$y, pred$si))   # same range on both axes and in every panel
p_pred <- ggplot(pred, aes(y, si, colour = time)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +   # perfect prediction
  geom_point(size = 1, alpha = 0.6) +
  geom_text(data = r2, aes(x = lims[1], y = lims[2], label = label), parse = TRUE,
            hjust = 0, vjust = 1, size = 2.8, inherit.aes = FALSE) +
  scale_colour_manual(values = c(`Day 0` = "#257de0", `Day 28` = "#f47230"), name = NULL) +
  facet_wrap(~study, ncol = 5) +
  coord_equal(xlim = lims, ylim = lims) +
  labs(x = "Observed antibody response", y = "Predicted antibody response (SI)",
       title = paste0("SI predictions without cross-fitting (", geneset_id, ", ", tolower(ge_day), ")")) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
        axis.text = element_text(colour = "black"), legend.position = "right")
p_pred

# ---- Paired RISE-Meta ---------------------------------------------------------------------------
# RISE-Meta drops a study for a surrogate when its within-study SE is zero; drop such studies for all
# surrogates, so that every surrogate is pooled over the same studies
y_up <- ind$y_post > ind$y_pre
zero_se <- unique(unlist(lapply(colnames(surr), function(s) {
  d <- y_up - (surr[post, s] > surr[pre, s])
  names(which(tapply(d, ind$study, var) == 0))
})))
zero_se
keep <- !ind$study %in% zero_se

capture.output(res <- rise.screen.meta(
  yone = ind$y_post[keep], yzero = ind$y_pre[keep],
  sone = surr[post, ][keep, ], szero = surr[pre, ][keep, ],
  studyone = ind$study[keep], studyzero = ind$study[keep], paired.all = TRUE,
  epsilon.meta.mode = "user", epsilon.meta = epsilon, epsilon.study = epsilon,
  test = "knha", test.target = "pi", meta.analysis.method = "RE", p.correction = "none",
  return.screen.plot = FALSE, return.forest.plot = FALSE, return.fit.plot = FALSE,
  return.study.similarity.plot = FALSE, return.evaluate.results = FALSE))

pooled <- as_tibble(res$screening.metrics.meta) %>%
  transmute(surrogate = factor(marker, colnames(surr)), delta = mu.delta, lower = ci.delta.lower,
            upper = ci.delta.upper, pi_lower = pi.delta.lower, pi_upper = pi.delta.upper, tau2)
pooled

# ---- Figure 2: pooled delta, one row per summary (as in the example, plus the SI without cross-fitting) --
# Thick = 95% CI, thin = 95% prediction interval
p_pooled <- ggplot(pooled, aes(delta, fct_rev(surrogate))) +
  annotate("rect", xmin = -epsilon, xmax = epsilon, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_linerange(aes(xmin = pi_lower, xmax = pi_upper)) +
  geom_linerange(aes(xmin = lower, xmax = upper), linewidth = 1.2) +
  geom_point(shape = 18, size = 2.5) +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(x = expression("Pooled " * delta), y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text = element_text(colour = "black"))
p_pooled

for (d in c(fs::path("output", "figures", "application"), fs::path("manuscript", "latex", "figures"))) {
  fs::dir_create(d)
  ggsave(fs::path(d, "example_geneset_no_crossfit_prediction.pdf"), p_pred, width = 6.3, height = 4.4)
  ggsave(fs::path(d, "example_geneset_no_crossfit_pooled.pdf"), p_pooled, width = 4.8, height = 2.2)
}

rm(list = ls())
