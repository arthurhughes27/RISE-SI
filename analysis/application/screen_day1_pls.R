# Screen of all gene sets at day 1 with the PLS surrogate index, including age and sex as baseline predictors
# (PLS was the best learner for the top day-1 gene set, example_geneset_learners_day1.R).
# All the analysis is in screen_day_template.R. Run from the project root.

ge_day <- "Day 1"       # as in data/analysis/analysis_set_endpoint.rds
tag <- "day1_pls"       # output names: output/figures/application/day1_pls/, output/application/screen_surrogates_day1_pls.rds
learner <- "pls"
covariates <- TRUE
source(fs::path("analysis", "application", "screen_day_template.R"))
