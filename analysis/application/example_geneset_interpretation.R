# Example application on a single gene set: interpreting the ridge surrogate index (SI) through its coefficients
# Gene set M156.0 (plasma cells & B cells, immunoglobulins), gene expression at day 7.
# Endpoint Y: mean of the log2 titres against H1, H3 and B, before vaccination and at day 28.
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

# Centre and scale genes and Y within study, both time points pooled, so that coefficients are comparable
# across genes (per within-study standard deviation of expression)
study <- rep(ind$study, 2)
y <- c(ind$y_pre, ind$y_post)
for (k in unique(study)) {
  r <- study == k
  S[r, ] <- scale(S[r, ])
  y[r] <- as.numeric(scale(y[r]))
}
S[is.na(S)] <- 0                      # genes constant within a study

# ---- Ridge coefficients -------------------------------------------------------------------------
# Point: ridge fitted on all studies. Bar: range over the leave-one-study-out fits used to build the SI.
fold <- rep(sample(rep(1:5, length.out = n)), 2)   # inner CV folds by individual
ridge_coef <- function(tr) coef(cv.glmnet(S[tr, ], y[tr], alpha = 0, foldid = fold[tr]), s = "lambda.min")[-1, 1]

coefs <- tibble(gene = colnames(S), all = ridge_coef(rep(TRUE, 2 * n)))
loso <- sapply(unique(study), function(k) ridge_coef(study != k))
coefs <- coefs %>% mutate(lo = apply(loso, 1, min), hi = apply(loso, 1, max), gene = fct_reorder(gene, all))
coefs %>% arrange(desc(all))

# ---- Figure: ridge coefficient of each gene -----------------------------------------------------
p <- ggplot(coefs, aes(all, gene)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_linerange(aes(xmin = lo, xmax = hi)) +
  geom_point(size = 1.8) +
  labs(x = "Ridge coefficient (standardised expression)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text = element_text(colour = "black"))
p

for (d in c(fs::path("output", "figures", "application"), fs::path("manuscript", "latex", "figures"))) {
  fs::dir_create(d)
  ggsave(fs::path(d, "example_geneset_interpretation.pdf"), p, width = 4.8, height = 3.6)
}
