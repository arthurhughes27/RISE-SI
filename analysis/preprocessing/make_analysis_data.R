# Preprocessing: builds the processed data used by every analysis script, from the raw IS2 files.
# Inactivated influenza vaccine studies of IS2, without the type 2 diabetes cohorts of SDY1119.
# Inputs:
#   data-raw/all_noNorm_eset.rds
#   data-raw/hai_2025-01-10_01-13-41.xlsx
#   data-raw/neut_ab_titer_2025-01-10_01-13-22.xlsx
# Outputs:
#   data/analysis_set.rds  list of
#     endpoint  complete paired data of the analysis set: one row per individual, GE day and subtype
#     counts    complete individuals per study and GE day (all studies); counted = at least min_n individuals
#     studies   one row per study: season, cells (PBMC or WB), platform vendor
#   data/expr.rds          expression (genes x samples) of the samples in `endpoint`, columns named "participant_id day"
# Run from the project root.

library(tidyverse)
library(Biobase)
library(readxl)
library(janitor)

ge_days <- c(1, 2, 3, 7)  # post-vaccination GE days (exact days, no windows)
min_n <- 5                # a study counts only with at least this many complete individuals

# Sample closest to `target` within [lo, hi], per group
closest <- function(d, target, lo, hi, by) {
  d %>% filter(day >= lo, day <= hi) %>%
    group_by(across(all_of(by))) %>% slice_min(abs(day - target), n = 1, with_ties = FALSE) %>% ungroup()
}

# ---- GE samples: inactivated influenza vaccine studies ----------------------------------------
eset <- readRDS(fs::path("data-raw", "all_noNorm_eset.rds"))
sample_key <- paste(eset$participant_id, round(eset$study_time_collected, 2))
samples <- pData(eset) %>%
  as_tibble() %>%
  filter(pathogen == "Influenza", vaccine_type == "Inactivated", !grepl("T2D", matrix)) %>%   # drop the SDY1119 T2D cohorts
  transmute(
    participant_id,
    study = study_accession,
    day = round(study_time_collected, 2),
    age = age_imputed,
    sex = gender,
    season = str_extract(vaccine, "[0-9]{4}"),
    cells = if_else(cell_type == "PBMC", "PBMC", "WB"),
    vendor = featureSetVendor
  )

participants <- samples %>% distinct(participant_id, study, age, sex)
studies <- samples %>% distinct(study, season, cells, vendor) %>% arrange(season, parse_number(study))

# ---- Antibody titres: log2 value per individual, assay, day and subtype ------------------------
subtype_of <- c(
  "A/California" = "H1",
  "A/California/7/2009" = "H1",
  "A/Brisbane/59/2007" = "H1",
  "A/Solomon Islands/3/2006" = "H1",
  "A/Victoria" = "H3",
  "A/Victoria/361/2011" = "H3",
  "A/Brisbane/10/2007" = "H3",
  "A/Uruguay/716/2007" = "H3",
  "A/Perth/16/2009" = "H3",
  "A/Wisconsin/67/2005" = "H3"
)

read_ab <- function(file, assay) {
  read_excel(fs::path("data-raw", file)) %>% clean_names() %>%
    transmute(participant_id, assay = assay, day = study_time_collected, strain = virus, value = value_preferred)
}
ab <- bind_rows(
  read_ab("hai_2025-01-10_01-13-41.xlsx", "HAI"),
  read_ab("neut_ab_titer_2025-01-10_01-13-22.xlsx", "nAb")
) %>%
  inner_join(participants %>% select(participant_id, study), by = "participant_id") %>%
  mutate(subtype = coalesce(subtype_of[strain], if_else(str_starts(strain, "B/"), "B", NA_character_))) %>%
  filter(!is.na(value), !is.na(subtype)) %>%           # strains not in the list above are not used
  group_by(participant_id, study, assay, day, subtype) %>%
  summarise(log_value = mean(log2(value)), .groups = "drop")   # replicates / strains of a subtype

# ---- Complete paired data -----------------------------------------------------------------------
# Pre-vaccination: closest to day 0 within days -7 to 0. Antibody post: closest to day 28 within days 21-35.
ge_pre <- closest(samples, 0, -7, 0, "participant_id") %>% select(participant_id, ge_pre_day = day)
ge_post <- bind_rows(
  samples %>% filter(day %in% ge_days) %>%
    transmute(participant_id, ge_day = paste("Day", day), ge_post_day = day),
  # "Day 1/2": the day 1 sample, or the day 2 sample when there is no day 1 sample
  samples %>% filter(day %in% c(1, 2)) %>%
    group_by(participant_id) %>% slice_min(day, n = 1) %>% ungroup() %>%
    transmute(participant_id, ge_day = "Day 1/2", ge_post_day = day)
)

by_ab <- c("participant_id", "study", "assay", "subtype")
ab_paired <- inner_join(
  closest(ab, 0, -7, 0, by_ab),
  closest(ab, 28, 21, 35, by_ab),
  by = by_ab,
  suffix = c("_pre", "_post")
)

# Endpoint assay: HAI if paired for all three subtypes, otherwise nAb
assay_used <- ab_paired %>%
  count(participant_id, assay) %>% filter(n == 3) %>%
  group_by(participant_id) %>% summarise(assay = if ("HAI" %in% assay) "HAI" else "nAb")

complete_data <- ge_post %>%
  inner_join(ge_pre, by = "participant_id") %>%
  inner_join(participants %>% filter(!is.na(age), sex %in% c("Female", "Male")), by = "participant_id") %>%
  inner_join(assay_used, by = "participant_id") %>%
  inner_join(ab_paired %>% select(-study), by = c("participant_id", "assay"))

counts <- complete_data %>%
  distinct(study, ge_day, participant_id) %>%
  count(study, ge_day) %>%
  complete(study = unique(studies$study), ge_day, fill = list(n = 0L)) %>%
  mutate(counted = n >= min_n)

# Analysis set: studies with at least `min_n` complete individuals at the given GE day
endpoint <- complete_data %>% semi_join(filter(counts, counted), by = c("study", "ge_day"))
endpoint %>% distinct(ge_day, study, participant_id) %>% group_by(ge_day) %>%
  summarise(K = n_distinct(study), N = n())

# ---- Expression of the samples in the analysis set ----------------------------------------------
keys <- unique(c(paste(endpoint$participant_id, endpoint$ge_pre_day), paste(endpoint$participant_id, endpoint$ge_post_day)))
expr <- exprs(eset)[, match(keys, sample_key)]
colnames(expr) <- keys
rm(eset)

# ---- Save ---------------------------------------------------------------------------------------
fs::dir_create("data")
saveRDS(list(endpoint = endpoint, counts = counts, studies = studies), fs::path("data", "analysis_set.rds"))
saveRDS(expr, fs::path("data", "expr.rds"))

rm(list = ls())
