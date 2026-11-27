# Screen of all gene sets with the ridge surrogate index: Day 1
# All the analysis is in screen_day_template.R. Run from the project root.

ge_day <- "Day 1"     # as in data/analysis/analysis_set_endpoint.rds
tag <- "day1"        # output names: output/figures/application/day1/, output/application/screen_surrogates_day1.rds
learner <- "ridge"
covariates <- FALSE
source(fs::path("analysis", "application", "screen_day_template.R"))
