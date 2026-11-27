# Screen of all gene sets at day 7 with the PLS surrogate index, including age and sex as baseline predictors
# (as screen_day1_pls.R, for comparison with the day-7 ridge screen in screen_all_genesets.R).
# All the analysis is in screen_day_template.R. Run from the project root.

ge_day <- "Day 7"       # as in data/analysis/analysis_set_endpoint.rds
tag <- "day7_pls"       # output names: output/figures/application/day7_pls/, output/application/screen_surrogates_day7_pls.rds
learner <- "pls"
covariates <- TRUE
source(fs::path("analysis", "application", "screen_day_template.R"))
