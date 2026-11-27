# Screen of all gene sets with the ridge surrogate index: Day 2
# All the analysis is in analysis/application/main/01_screen.R. Run from the project root.

ge_day <- "Day 2"     # as in data/analysis_set.rds
tag <- "day2"        # output names: output/figures/application/day2/, output/application/screen_surrogates_day2.rds
learner <- "ridge"
source(fs::path("analysis", "application", "main", "01_screen.R"))

rm(list = ls())
