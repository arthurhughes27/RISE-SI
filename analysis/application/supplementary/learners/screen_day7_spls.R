# Screen of all gene sets at day 7 with the sparse PLS surrogate index
# All the analysis is in analysis/application/main/01_screen.R. Run from the project root.

ge_day <- "Day 7"       # as in data/analysis_set.rds
tag <- "day7_spls"      # output names: output/figures/application/day7_spls/, output/application/screen_surrogates_day7_spls.rds
learner <- "spls"
source(fs::path("analysis", "application", "main", "01_screen.R"))

rm(list = ls())
