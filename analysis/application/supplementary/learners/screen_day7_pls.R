# Screen of all gene sets at day 7 with the PLS surrogate index
# (as supplementary/learners/screen_day1_pls.R, for comparison with the day-7 ridge screen in main/01_screen.R).
# All the analysis is in analysis/application/main/01_screen.R. Run from the project root.

ge_day <- "Day 7"       # as in data/analysis_set.rds
tag <- "day7_pls"       # output names: output/figures/application/day7_pls/, output/application/screen_surrogates_day7_pls.rds
learner <- "pls"
source(fs::path("analysis", "application", "main", "01_screen.R"))

rm(list = ls())
