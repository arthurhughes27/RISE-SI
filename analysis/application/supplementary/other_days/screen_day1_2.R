# Screen of all gene sets with the ridge surrogate index: Day 1 and day 2 pooled (day 1 sample, or day 2 when there is no day 1 sample)
# All the analysis is in analysis/application/main/01_screen.R. Run from the project root.

ge_day <- "Day 1/2"     # as in data/analysis_set.rds
tag <- "day1_2"        # output names: output/figures/application/day1_2/, output/application/screen_surrogates_day1_2.rds
learner <- "ridge"
source(fs::path("analysis", "application", "main", "01_screen.R"))

rm(list = ls())
