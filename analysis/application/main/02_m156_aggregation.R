# Example application on a single gene set: simple aggregation (mean, median, max) vs surrogate index (ridge)
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

fig_dirs <- c(fs::path("output", "figures", "application"), fs::path("manuscript", "latex", "figures"))
fs::dir_create(fig_dirs)
save_pdf <- function(p, file, w, h) for (d in fig_dirs) ggsave(fs::path(d, file), p, width = w, height = h)

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

# ---- Surrogates ---------------------------------------------------------------------------------
surr <- cbind(Mean = rowMeans(S), Median = apply(S, 1, median), Max = apply(S, 1, max), `SI (ridge)` = NA)

# SI: ridge regression of Y on the genes, cross-fitted leave-one-study-out.
# Inner CV folds (for the penalty) are by individual, so an individual's two samples stay together.
fold <- rep(sample(rep(1:5, length.out = n)), 2)
for (k in unique(study)) {
  tr <- study != k
  fit <- cv.glmnet(S[tr, ], y[tr], alpha = 0, foldid = fold[tr])
  surr[!tr, "SI (ridge)"] <- predict(fit, S[!tr, ], s = "lambda.min")
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

res$screening.metrics.meta %>%
  select(marker, n.studies, mu.delta, ci.delta.lower, ci.delta.upper, pi.delta.lower, pi.delta.upper, tau2)

# ---- Figure: delta per study and pooled, one panel per surrogate ------------------------------
# Study rows: estimate with 95% CI. Pooled row: thick = 95% CI, thin = 95% prediction interval.
by_study <- as_tibble(res$screening.metrics.study) %>%
  transmute(surrogate = marker, study, delta, lower = delta - 1.96 * sd, upper = delta + 1.96 * sd)
pooled <- as_tibble(res$screening.metrics.meta) %>%
  transmute(surrogate = marker, study = "Pooled", delta = mu.delta, lower = ci.delta.lower, upper = ci.delta.upper,
            pi_lower = pi.delta.lower, pi_upper = pi.delta.upper)
fig <- bind_rows(by_study, pooled) %>%
  mutate(surrogate = factor(surrogate, colnames(surr)),
         study = factor(study, rev(c(sort(unique(by_study$study)), "Pooled"))),
         is_pooled = study == "Pooled")

p <- ggplot(fig, aes(delta, study)) +
  annotate("rect", xmin = -epsilon, xmax = epsilon, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_linerange(aes(xmin = pi_lower, xmax = pi_upper), na.rm = TRUE) +
  geom_linerange(aes(xmin = lower, xmax = upper, linewidth = is_pooled)) +
  geom_point(aes(shape = is_pooled), size = 1.8) +
  scale_linewidth_manual(values = c(`FALSE` = 0.4, `TRUE` = 1.2), guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
  facet_wrap(~surrogate, nrow = 1) +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(x = expression(delta[k] == U[list(Y, k)] - U[list(S, k)]), y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "grey92"), axis.text = element_text(colour = "black"))
p

save_pdf(p, "example_geneset_aggregation.pdf", 6.3, 4)

# ---- Figure: pooled delta only, one row per surrogate ---------------------------------------
# Thick = 95% CI, thin = 95% prediction interval
p_pooled <- ggplot(pooled, aes(delta, factor(surrogate, rev(colnames(surr))))) +
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

save_pdf(p_pooled, "example_geneset_aggregation_pooled.pdf", 4.8, 2)

rm(list = ls())
