# Example application on a single gene set: impact of the learner used to build the surrogate index (SI)
# Gene set M156.0 (plasma cells & B cells, immunoglobulins), gene expression at day 7.
# Endpoint Y: mean of the log2 titres against H1, H3 and B, before vaccination and at day 28.
# The mean expression of the gene set is included as a reference (no learning).
# Input: data/analysis/analysis_set_endpoint.rds (from analysis/descriptive/descriptive.R). Run from the project root.

library(tidyverse)
library(Biobase)
library(glmnet)
library(pls)
library(ranger)
library(xgboost)
library(SuperLearner)
library(SurrogateRank)

geneset_id <- "M156.0"
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

# ---- Learners: fit on (x, y), predict newx -----------------------------------------------------
# `fold`: inner CV fold of each training row (by individual), used to tune glmnet penalties and PLS components
glmnet_fit <- function(x, y, newx, fold, alpha) {
  predict(cv.glmnet(x, y, alpha = alpha, foldid = fold), newx, s = "lambda.min")[, 1]
}
learners <- list(
  Ridge = function(x, y, newx, fold) glmnet_fit(x, y, newx, fold, alpha = 0),
  Lasso = function(x, y, newx, fold) glmnet_fit(x, y, newx, fold, alpha = 1),
  `Elastic net` = function(x, y, newx, fold) glmnet_fit(x, y, newx, fold, alpha = 0.5),
  `Linear regression` = function(x, y, newx, fold) {
    b <- lm.fit(cbind(1, x), y)$coefficients
    drop(cbind(1, newx) %*% replace(b, is.na(b), 0))
  },
  PLS = function(x, y, newx, fold) {
    fit <- plsr(y ~ ., data = data.frame(y = y, x), ncomp = 3, validation = "CV", segments = split(seq_along(y), fold))
    ncomp <- which.min(RMSEP(fit, estimate = "CV", intercept = FALSE)$val[1, 1, ])
    drop(predict(fit, newdata = data.frame(newx), ncomp = ncomp))
  },
  `Random forest` = function(x, y, newx, fold) {
    predict(ranger(x = x, y = y, num.trees = 500, min.node.size = 10), newx)$predictions
  },
  XGBoost = function(x, y, newx, fold) {
    fit <- xgb.train(params = list(objective = "reg:squarederror", max_depth = 2, learning_rate = 0.05, subsample = 0.8),
                     data = xgb.DMatrix(x, label = y), nrounds = 300, verbose = 0)
    predict(fit, xgb.DMatrix(newx))
  }
)

# Super learner: non-negative weights on all learners except XGBoost (run time)
sl_wrap <- function(f) function(Y, X, newX, id, ...) list(pred = f(as.matrix(X), Y, as.matrix(newX), fold_of[id]), fit = list())
sl_library <- c("Ridge", "Lasso", "Elastic net", "Linear regression", "PLS", "Random forest")
sl_names <- paste0("SL.", make.names(sl_library))
for (i in seq_along(sl_library)) assign(sl_names[i], sl_wrap(learners[[sl_library[i]]]))

# ---- Surrogates: mean (reference) and cross-fitted SIs -----------------------------------------
# Each SI is cross-fitted leave-one-study-out: the SI of study k comes from a model trained on the other studies.
surr <- matrix(NA, 2 * n, length(learners) + 2, dimnames = list(NULL, c("Mean", names(learners), "Super learner")))
surr[, "Mean"] <- rowMeans(S)

id <- rep(1:n, 2)                            # individual of each row
fold_of <- sample(rep(1:5, length.out = n))  # inner CV fold of each individual
for (k in unique(study)) {
  tr <- study != k
  for (l in names(learners)) surr[!tr, l] <- learners[[l]](S[tr, ], y[tr], S[!tr, ], fold_of[id[tr]])

  # Super learner weights chosen by leave-one-(training)-study-out CV, mimicking use in a new trial
  sl <- SuperLearner(Y = y[tr], X = as.data.frame(S[tr, ]), newX = as.data.frame(S[!tr, ]), id = id[tr],
                     SL.library = sl_names, method = "method.NNLS",
                     cvControl = list(V = n_distinct(study[tr]), validRows = unname(split(seq_len(sum(tr)), study[tr]))))
  surr[!tr, "Super learner"] <- sl$SL.predict
  cat("Held-out study", k, "done\n")
}

# ---- Paired RISE-Meta ---------------------------------------------------------------------------
# RISE-Meta drops a study for a surrogate when its within-study SE is zero; 
# drop such studies for all
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

# ---- Figure: pooled delta by learner, mean as reference ---------------------------------------
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
  ggsave(fs::path(d, "example_geneset_learners.pdf"), p, width = 4.8, height = 3)
}
