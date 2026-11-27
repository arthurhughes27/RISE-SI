# Example application on a single gene set: interpreting the ridge surrogate index (SI) through its coefficients
# Gene set M156.0 (plasma cells & B cells, immunoglobulins), gene expression at day 7.
# Endpoint Y: mean of the log2 titres against H1, H3 and B, before vaccination and at day 28.
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

rm(list = ls())
