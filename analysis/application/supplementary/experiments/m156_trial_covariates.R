# Experiment: trial-level covariates in the ridge surrogate index (SI) of M156.0, day 7
# A trial-level covariate is constant within a study, so it cancels in the paired (within-person) comparison and
# in the within-study centring. It can only act through gene x covariate interactions, i.e. by letting the
# weight of each gene depend on the type of study. Three trial-level covariates are tried one at a time, then together:
#   cell type     whole blood (vs PBMC)
#   platform      Illumina (vs Affymetrix or custom)
#   season        vaccination season (year, centred and scaled)
# The mean expression of the gene set and the SI without covariates are references.
# Endpoint Y: mean of the log2 titres against H1, H3 and B, before vaccination and at day 28.
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

# ---- Gene set expression (rows: all pre samples, then all post samples) and study metadata --------
genes <- btm_genes(geneset_id)
S <- load_expr(ind, genes)     # genes measured in every study

# Trial-level covariates (one row per study)
trial <- readRDS(fs::path("data", "analysis_set.rds"))$studies %>%
  filter(study %in% ind$study) %>%
  transmute(study, whole_blood = as.numeric(cells == "WB"),
            illumina = as.numeric(vendor %in% "Illumina"), season = as.numeric(scale(as.numeric(season))))
trial

# Centre and scale genes and Y within study, both time points pooled
study <- rep(ind$study, 2)
y <- c(ind$y_pre, ind$y_post)
for (k in unique(study)) {
  r <- study == k
  S[r, ] <- scale(S[r, ])
  y[r] <- as.numeric(scale(y[r]))
}
S[is.na(S)] <- 0                      # genes constant within a study

# ---- Surrogates: mean, SI without covariates, and SI with gene x trial covariate interactions ----
# Each SI is cross-fitted leave-one-study-out; inner CV folds (for the ridge penalty) are by individual.
covariates <- c(`cell type` = "whole_blood", platform = "illumina", season = "season")
surr <- matrix(NA, 2 * n, 3 + length(covariates),
               dimnames = list(NULL, c("Mean", "SI (ridge)", paste0("SI (ridge) x ", names(covariates)),
                                       "SI (ridge) x all three")))
surr[, "Mean"] <- rowMeans(S)

fold <- rep(sample(rep(1:5, length.out = n)), 2)
ridge_cf <- function(X) {
  si <- rep(NA, 2 * n)
  for (k in unique(study)) {
    tr <- study != k
    fit <- cv.glmnet(X[tr, ], y[tr], alpha = 0, foldid = fold[tr])
    si[!tr] <- predict(fit, X[!tr, ], s = "lambda.min")
  }
  si
}
surr[, "SI (ridge)"] <- ridge_cf(S)
Z <- as.matrix(trial[match(study, trial$study), covariates])   # trial covariates of each row
for (i in seq_along(covariates)) {
  surr[, 2 + i] <- ridge_cf(cbind(S, S * Z[, i]))               # genes, and genes x covariate
}
surr[, "SI (ridge) x all three"] <- ridge_cf(cbind(S, S * Z[, 1], S * Z[, 2], S * Z[, 3]))

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

# ---- Figure: pooled delta by surrogate ----------------------------------------------------------
# Thick = CI, thin = prediction interval
p <- ggplot(pooled, aes(delta, fct_rev(surrogate))) +
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
p

# Experiment: saved for inspection only (not in the document)
fs::dir_create(fs::path("output", "figures", "application"))
ggsave(fs::path("output", "figures", "application", "example_geneset_trial_covariates.pdf"), p, width = 4.8, height = 2.6)

rm(list = ls())
