# Screen of all gene sets at day 1 with the sparse PLS surrogate index, including age and sex as baseline predictors
# All the analysis is in analysis/application/main/01_screen.R. Run from the project root.

ge_day <- "Day 1"       # as in data/analysis_set.rds
tag <- "day1_spls"      # output names: output/figures/application/day1_spls/, output/application/screen_surrogates_day1_spls.rds
learner <- "spls"
covariates <- TRUE
source(fs::path("analysis", "application", "main", "01_screen.R"))

rm(list = ls())
