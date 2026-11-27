# Screen of all gene sets at day 1 with the PLS surrogate index, including age and sex as baseline predictors
# (PLS was the best learner for the top day-1 gene set, supplementary/learners/m71_learners_day1.R).
# All the analysis is in analysis/application/main/01_screen.R. Run from the project root.

ge_day <- "Day 1"       # as in data/analysis_set.rds
tag <- "day1_pls"       # output names: output/figures/application/day1_pls/, output/application/screen_surrogates_day1_pls.rds
learner <- "pls"
covariates <- TRUE
source(fs::path("analysis", "application", "main", "01_screen.R"))

rm(list = ls())
