# Screen of all gene sets with the ridge surrogate index: Day 1 and day 2 pooled (day 1 sample, or day 2 when there is no day 1 sample)
# All the analysis is in screen_day_template.R. Run from the project root.

ge_day <- "Day 1/2"     # as in data/analysis/analysis_set_endpoint.rds
tag <- "day1_2"        # output names: output/figures/application/day1_2/, output/application/screen_surrogates_day1_2.rds
learner <- "ridge"
covariates <- FALSE
source(fs::path("analysis", "application", "screen_day_template.R"))
