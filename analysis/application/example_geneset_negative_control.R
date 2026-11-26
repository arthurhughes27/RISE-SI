# Example application: negative control
# Does the ridge surrogate index (SI) turn any gene set into a good trial-level surrogate?
# Same comparison as for M156.0 (mean, median, max vs SI) on M1.1 (integrin cell surface interactions (II)),
# a module chosen at random among those not expected to respond to the vaccine.
# Gene expression at day 7. Endpoint Y: mean of the log2 titres against H1, H3 and B, before vaccination and at day 28.
# Input: data/analysis/analysis_set_endpoint.rds (from analysis/descriptive/descriptive.R). Run from the project root.

library(tidyverse)
library(Biobase)
library(glmnet)
library(SurrogateRank)

geneset_id <- "M1.1"
ge_day <- "Day 7"
epsilon <- 0.2     # equivalence margin (pooled and within-study)
set.seed(1)

# ---- Endpoint: one row per individual -----------------------------------------------------------
ind <- readRDS(fs::path("data", "analysis", "analysis_set_endpoint.rds")) %>%
  filter(ge_day == !!ge_day) %>%
  group_by(participant_id, study, ge_pre_day, ge_post_day) %>%
  summarise(y_pre = mean(log_value_pre), y_post = mean(log_value_post), .groups = "drop")
n <- nrow(ind)
pre <- 1:n
post <- n + 1:n

# ---- Gene set expression: rows are all pre samples, then all post samples -----------------------
gmt <- strsplit(readLines(fs::path("data-raw", "BTM_for_GSEA_20131008.gmt")), "\t")
genes <- gmt[[which(grepl(paste0("(", geneset_id, ")"), sapply(gmt, `[`, 1), fixed = TRUE))]][-(1:2)]

eset <- readRDS(fs::path("data-raw", "is2", "all_noNorm_eset.rds"))
key <- paste(eset$participant_id, round(eset$study_time_collected, 2))
cols <- match(c(paste(ind$participant_id, ind$ge_pre_day), paste(ind$participant_id, ind$ge_post_day)), key)
S <- t(exprs(eset)[intersect(genes, rownames(eset)), cols])
rm(eset)
S <- S[, colSums(is.na(S)) == 0]     # genes measured in every study

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

pooled <- as_tibble(res$screening.metrics.meta) %>%
  transmute(surrogate = factor(marker, colnames(surr)), delta = mu.delta, lower = ci.delta.lower,
            upper = ci.delta.upper, pi_lower = pi.delta.lower, pi_upper = pi.delta.upper, tau2)
pooled

# ---- Figure: pooled delta, one row per summary ------------------------------------------------
# Thick = 95% CI, thin = 95% prediction interval
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

for (d in c(fs::path("output", "figures", "application"), fs::path("manuscript", "latex", "figures"))) {
  fs::dir_create(d)
  ggsave(fs::path(d, "example_geneset_negative_control.pdf"), p, width = 4.8, height = 2)
}
