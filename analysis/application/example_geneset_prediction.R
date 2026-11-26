# Example application on a single gene set: how well does the cross-fitted ridge SI predict the antibody response
# of individuals in a new study?
# Gene set M156.0 (plasma cells & B cells, immunoglobulins), gene expression at day 7.
# Endpoint Y: mean of the log2 titres against H1, H3 and B, at baseline and day 28 (centred and scaled within study).
# The SI of study k comes from a ridge model trained on the other studies only, so it is a prediction for a new study.
# Input: data/analysis/analysis_set_endpoint.rds (from analysis/descriptive/descriptive.R). Run from the project root.

library(tidyverse)
library(Biobase)
library(glmnet)

geneset_id <- "M156.0"
ge_day <- "Day 7"
set.seed(1)

# ---- Endpoint: one row per individual -----------------------------------------------------------
ind <- readRDS(fs::path("data", "analysis", "analysis_set_endpoint.rds")) %>%
  filter(ge_day == !!ge_day) %>%
  group_by(participant_id, study, ge_pre_day, ge_post_day) %>%
  summarise(y_pre = mean(log_value_pre), y_post = mean(log_value_post), .groups = "drop")
n <- nrow(ind)

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

# ---- SI: ridge regression of Y on the genes, cross-fitted leave-one-study-out ------------------
# Inner CV folds (for the penalty) are by individual, so an individual's two samples stay together.
si <- rep(NA, 2 * n)
fold <- rep(sample(rep(1:5, length.out = n)), 2)
for (k in unique(study)) {
  tr <- study != k
  fit <- cv.glmnet(S[tr, ], y[tr], alpha = 0, foldid = fold[tr])
  si[!tr] <- predict(fit, S[!tr, ], s = "lambda.min")
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
