# Screen of all gene sets with the ridge surrogate index: Day 2
# All the analysis is in screen_day_template.R. Run from the project root.

ge_day <- "Day 2"     # as in data/analysis/analysis_set_endpoint.rds
tag <- "day2"        # output names: output/figures/application/day2/, output/application/screen_surrogates_day2.rds
source(fs::path("analysis", "application", "screen_day_template.R"))
