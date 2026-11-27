# Example application on a single gene set: surrogate index (ridge) vs each gene of the set on its own
# Gene set M156.0 (plasma cells & B cells, immunoglobulins), gene expression at day 7.
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

# ---- Gene set expression: rows are all pre samples, then all post samples -----------------------
genes <- btm_genes(geneset_id)
S <- load_expr(ind, genes)     # genes measured in every study

# Centre and scale genes and Y within study, both time points pooled (removes platform / assay offsets;
# does not change anyone's direction of change)
study <- rep(ind$study, 2)
y <- c(ind$y_pre, ind$y_post)
for (k in unique(study)) {
  r <- study == k
  S[r, ] <- scale(S[r, ])
  y[r] <- as.numeric(scale(y[r]))
}
S[is.na(S)] <- 0                      # genes constant within a study
X <- cbind(S, load_covariates(ind))     # predictors of the SI: genes and baseline covariates (age, sexMale)

# ---- Surrogates: ridge SI and each gene, cross-fitted leave-one-study-out ---------------------
# SI: ridge regression of Y on the genes and covariates (inner CV folds by individual, for the penalty).
# Single genes: oriented so that larger is better, using the sign of their correlation with Y in the
# other studies (the SI learns this orientation; a gene on its own does not).
surr <- matrix(NA, 2 * n, ncol(S) + 1, dimnames = list(NULL, c("SI (ridge)", colnames(S))))
fold <- rep(sample(rep(1:5, length.out = n)), 2)
for (k in unique(study)) {
  tr <- study != k
  fit <- cv.glmnet(X[tr, ], y[tr], alpha = 0, foldid = fold[tr])
  surr[!tr, "SI (ridge)"] <- predict(fit, X[!tr, ], s = "lambda.min")
  surr[!tr, colnames(S)] <- sweep(S[!tr, ], 2, sign(cor(S[tr, ], y[tr])), "*")
}

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
  transmute(surrogate = marker, delta = mu.delta, lower = ci.delta.lower, upper = ci.delta.upper,
            pi_lower = pi.delta.lower, pi_upper = pi.delta.upper, tau2) %>%
  arrange(surrogate != "SI (ridge)", delta)                          # SI first, then genes by pooled delta
pooled

# ---- Figure: pooled delta for the SI (top) and each gene --------------------------------------
# Thick = 95% CI, thin = 95% prediction interval
p <- ggplot(pooled, aes(delta, factor(surrogate, rev(surrogate)))) +
  annotate("rect", xmin = -epsilon, xmax = epsilon, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_hline(yintercept = nrow(pooled) - 0.5, colour = "grey60") +    # separates the SI from the genes
  geom_linerange(aes(xmin = pi_lower, xmax = pi_upper)) +
  geom_linerange(aes(xmin = lower, xmax = upper), linewidth = 1.2) +
  geom_point(shape = 18, size = 2.5) +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(x = expression("Pooled " * delta), y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text = element_text(colour = "black"))
p

for (d in c(fs::path("output", "figures", "application"), fs::path("manuscript", "latex", "figures"))) {
  fs::dir_create(d)
  ggsave(fs::path(d, "example_geneset_single_genes.pdf"), p, width = 4.8, height = 4)
}

rm(list = ls())
