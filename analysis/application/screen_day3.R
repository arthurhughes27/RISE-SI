# Screen of all gene sets with the ridge surrogate index: Day 3
# All the analysis is in screen_day_template.R. Run from the project root.

ge_day <- "Day 3"     # as in data/analysis/analysis_set_endpoint.rds
tag <- "day3"        # output names: output/figures/application/day3/, output/application/screen_surrogates_day3.rds
learner <- "ridge"
covariates <- FALSE
source(fs::path("analysis", "application", "screen_day_template.R"))
