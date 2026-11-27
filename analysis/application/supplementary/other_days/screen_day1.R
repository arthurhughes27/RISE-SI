# Screen of all gene sets with the ridge surrogate index: Day 1
# All the analysis is in analysis/application/main/01_screen.R. Run from the project root.

ge_day <- "Day 1"     # as in data/analysis_set.rds
tag <- "day1"        # output names: output/figures/application/day1/, output/application/screen_surrogates_day1.rds
learner <- "ridge"
source(fs::path("analysis", "application", "main", "01_screen.R"))

rm(list = ls())
