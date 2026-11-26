# Descriptive figures for the plan document (Data description section)
#   Figure 1: individuals with complete paired data, by study and post-vaccination gene expression (GE) day
#   Figure 2: within-study effect of vaccination on the antibody endpoint, U_Y,k
# Also saves the analysis set used by the SI analyses: data/analysis/analysis_set_endpoint.rds
# Run from the project root.

library(tidyverse)
library(Biobase)
library(readxl)
library(janitor)

ge_days <- c(1, 2, 3, 7)  # post-vaccination GE days in the analysis set (exact days, no windows)
fig_days <- paste("Day", c(1, 3, 7))   # days shown in the figures
min_n <- 5                # a study counts only with at least this many complete individuals

fig_dirs <- c(
  fs::path("output", "figures", "descriptive"),
  fs::path("manuscript", "latex", "figures")
)
fs::dir_create(fig_dirs)
save_pdf <- function(p, file, w, h)
  for (d in fig_dirs)
    ggsave(fs::path(d, file), p, width = w, height = h)

# Sample closest to `target` within [lo, hi], per group
closest <- function(d, target, lo, hi, by) {
  d %>% filter(day >= lo, day <= hi) %>%
    group_by(across(all_of(by))) %>% slice_min(abs(day - target), n = 1, with_ties = FALSE) %>% ungroup()
}

# ---- GE samples: inactivated influenza vaccine studies ----------------------------------------
eset <- readRDS(fs::path("data-raw", "is2", "all_noNorm_eset.rds"))
samples <- pData(eset) %>%
  as_tibble() %>%
  filter(pathogen == "Influenza",
         vaccine_type == "Inactivated",!grepl("T2D", matrix)) %>%                                  # drop the SDY1119 type 2 diabetes cohorts
  transmute(
    participant_id,
    study = study_accession,
    day = round(study_time_collected, 2),
    age = age_imputed,
    sex = gender,
    season = str_extract(vaccine, "[0-9]{4}"),
    cells = if_else(cell_type == "PBMC", "PBMC", "WB")
  )
rm(eset)

participants <- samples %>% distinct(participant_id, study, age, sex)
studies <- samples %>% distinct(study, season, cells) %>% arrange(season, parse_number(study))
study_label <- with(studies, setNames(sprintf("%s  (%s, %s)", study, season, cells), study))

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
  read_excel(fs::path("data-raw", "is2", file)) %>% clean_names() %>%
    transmute(
      participant_id,
      assay = assay,
      day = study_time_collected,
      strain = virus,
      value = value_preferred
    )
}
ab <- bind_rows(
  read_ab("hai_2025-01-10_01-13-41.xlsx", "HAI"),
  read_ab("neut_ab_titer_2025-01-10_01-13-22.xlsx", "nAb")
) %>%
  inner_join(participants %>% select(participant_id, study), by = "participant_id") %>%
  filter(!is.na(value)) %>%
  mutate(
    subtype = coalesce(subtype_of[strain], if_else(
      str_starts(strain, "B/"), "B", NA_character_
    )),
    log_value = if_else(study == "SDY1276", value, log2(value))
  ) %>%   # SDY1276 is already on the log scale
  group_by(participant_id, study, assay, day, subtype) %>%
  summarise(log_value = mean(log_value), .groups = "drop")                 # replicates / strains of a subtype

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
  group_by(participant_id) %>% summarise(assay = if ("HAI" %in% assay)
    "HAI"
    else
      "nAb")

complete_data <- ge_post %>%
  inner_join(ge_pre, by = "participant_id") %>%
  inner_join(participants %>% filter(!is.na(age), sex %in% c("Female", "Male")), by = "participant_id") %>%
  inner_join(assay_used, by = "participant_id") %>%
  inner_join(ab_paired %>% select(-study), by = c("participant_id", "assay"))

counts <- complete_data %>%
  distinct(study, ge_day, participant_id) %>%
  count(study, ge_day) %>%
  complete(study = studies$study, ge_day, fill = list(n = 0L)) %>%
  mutate(counted = n >= min_n)
totals <- counts %>% filter(counted) %>% group_by(ge_day) %>% summarise(K = n(), N = sum(n))
totals

# Analysis set for the SI analyses: studies with at least `min_n` complete individuals
analysis_set <- complete_data %>% semi_join(filter(counts, counted), by = c("study", "ge_day"))
fs::dir_create(fs::path("data", "analysis"))
saveRDS(analysis_set,
        fs::path("data", "analysis", "analysis_set_endpoint.rds"))

# ---- Figure 1: complete individuals per study and GE day --------------------------------------
y_levels <- rev(c(studies$study, "Total"))
fig1 <- counts %>% filter(ge_day %in% fig_days) %>% mutate(row = factor(study, y_levels),
                          n_fill = if_else(counted, n, NA_integer_))
fig1_tot <- totals %>% filter(ge_day %in% fig_days) %>% mutate(row = factor("Total", y_levels),
                              label = sprintf("K = %d\nN = %d", K, N))

p1 <- ggplot(fig1, aes(ge_day, row)) +
  geom_tile(aes(fill = n_fill),
            colour = "white",
            linewidth = 1) +
  geom_text(aes(
    label = if_else(counted, as.character(n), ""),
    colour = n > 100
  ), size = 4.4) +
  geom_text(
    data = fig1_tot,
    aes(label = label),
    size = 4.2,
    lineheight = 0.9
  ) +
  scale_fill_gradientn(
    colours = c("#eaf2fc", "#257de0", "#0b3d78"),
    trans = "sqrt",
    na.value = "#f2f2f2",
    name = "Individuals"
  ) +       # grey: fewer than min_n
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"),
                      guide = "none") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(labels = c(study_label, Total = "Total"),
                   drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.text.x.top = element_text(size = 15)
  )
p1

save_pdf(p1, "analysable_by_day_hai_nab_3subtypes.pdf", 6.3, 5.6)

# ---- Figure 2: U_Y,k = P(Y_post > Y_pre) + 0.5 P(Y_post = Y_pre) --------------------------------
# Y: mean of the log2 titres against H1, H3 and B. 95% Wald CI, truncated to [0, 1].
uy <- analysis_set %>%
  filter(ge_day %in% fig_days) %>%
  group_by(study, ge_day, participant_id) %>%
  summarise(
    y_pre = mean(log_value_pre),
    y_post = mean(log_value_post),
    .groups = "drop"
  ) %>%
  mutate(g = (y_post > y_pre) + 0.5 * (y_post == y_pre)) %>%
  group_by(study, ge_day) %>%
  summarise(U = mean(g),
            se = sd(g) / sqrt(n()),
            .groups = "drop") %>%
  mutate(lower = pmax(0, U - 1.96 * se),
         upper = pmin(1, U + 1.96 * se))

p2 <- ggplot(uy, aes(U, factor(study, rev(studies$study)))) +
  geom_linerange(aes(xmin = lower, xmax = upper)) +
  geom_point(size = 1.5) +
  facet_wrap( ~ ge_day, nrow = 1) +
  scale_x_continuous(
    limits = c(0.6, 1.02),
    breaks = c(0.6, 0.8, 1),
    labels = c("0.6", "0.8", "1")
  ) +
  labs(x = expression(U[list(Y, k)]), y = NULL) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    axis.text = element_text(colour = "black"),
    panel.spacing.x = unit(0.8, "lines")
  )
p2

save_pdf(p2, "endpoint_uy_by_study.pdf", 4.8, 3.6)
