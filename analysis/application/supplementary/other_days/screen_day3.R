# Screen of all gene sets with the ridge surrogate index: Day 3
# All the analysis is in analysis/application/main/01_screen.R. Run from the project root.

ge_day <- "Day 3"     # as in data/analysis_set.rds
tag <- "day3"        # output names: output/figures/application/day3/, output/application/screen_surrogates_day3.rds
learner <- "ridge"
covariates <- FALSE
source(fs::path("analysis", "application", "main", "01_screen.R"))

rm(list = ls())
