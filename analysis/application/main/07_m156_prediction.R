# Example application on a single gene set: how well does the cross-fitted ridge SI predict the antibody response
# of individuals in a new study?
# Gene set M156.0 (plasma cells & B cells, immunoglobulins), gene expression at day 7.
# Endpoint Y: mean of the log2 titres against H1, H3 and B, at baseline and day 28 (centred and scaled within study).
# The SI of study k comes from a ridge model trained on the other studies only, so it is a prediction for a new study.
# Input: data/ (from analysis/preprocessing/make_analysis_data.R). Run from the project root.

library(tidyverse)
library(glmnet)
source(fs::path("R", "load_data.R"))

geneset_id <- "M156.0"
ge_day <- "Day 7"
set.seed(1)

# ---- Endpoint: one row per individual -----------------------------------------------------------
ind <- load_endpoint(ge_day)
n <- nrow(ind)

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

# ---- SI: ridge regression of Y on the genes and covariates, cross-fitted leave-one-study-out ------------------
# Inner CV folds (for the penalty) are by individual, so an individual's two samples stay together.
si <- rep(NA, 2 * n)
fold <- rep(sample(rep(1:5, length.out = n)), 2)
for (k in unique(study)) {
  tr <- study != k
  fit <- cv.glmnet(X[tr, ], y[tr], alpha = 0, foldid = fold[tr])
  si[!tr] <- predict(fit, X[!tr, ], s = "lambda.min")
}

pred <- tibble(participant_id = rep(ind$participant_id, 2), study, time = rep(c("Day 0", "Day 28"), each = n), y, si)

# Within-study R2: squared correlation between predicted (SI) and observed response, both time points pooled
r2 <- pred %>% group_by(study) %>% summarise(label = sprintf("R^2 == %.2f", cor(si, y)^2))
r2

# ---- Figure: predicted vs observed, one panel per held-out study ------------------------------
lims <- range(c(pred$y, pred$si))   # same range on both axes and in every panel
p <- ggplot(pred, aes(y, si, colour = time)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +   # perfect prediction
  geom_point(size = 1, alpha = 0.6) +
  geom_text(data = r2, aes(x = lims[1], y = lims[2], label = label), parse = TRUE,
            hjust = 0, vjust = 1, size = 2.8, inherit.aes = FALSE) +
  scale_colour_manual(values = c(`Day 0` = "#257de0", `Day 28` = "#f47230"), name = NULL) +
  facet_wrap(~study, ncol = 5) +
  coord_equal(xlim = lims, ylim = lims) +
  labs(x = "Observed antibody response", y = "Predicted antibody response (SI)",
       title = paste0("Cross-fitted SI predictions in each held-out study (", geneset_id, ", ", tolower(ge_day), ")")) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
        axis.text = element_text(colour = "black"), legend.position = "right")
p

for (d in c(fs::path("output", "figures", "application"), fs::path("manuscript", "latex", "figures"))) {
  fs::dir_create(d)
  ggsave(fs::path(d, "example_geneset_prediction.pdf"), p, width = 6.3, height = 4.4)
}

rm(list = ls())
