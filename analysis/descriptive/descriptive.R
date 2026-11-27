# Descriptive figures for the plan document (Data description section)
#   Figure 1: individuals with complete paired data, by study and post-vaccination gene expression (GE) day
#   Figure 2: within-study effect of vaccination on the antibody endpoint, U_Y,k
# Input: data/analysis_set.rds (from analysis/preprocessing/make_analysis_data.R). Run from the project root.

library(tidyverse)

fig_days <- paste("Day", c(1, 3, 7))   # days shown in the figures

fig_dirs <- c(
  fs::path("output", "figures", "descriptive"),
  fs::path("manuscript", "latex", "figures")
)
fs::dir_create(fig_dirs)
save_pdf <- function(p, file, w, h)
  for (d in fig_dirs)
    ggsave(fs::path(d, file), p, width = w, height = h)

d <- readRDS(fs::path("data", "analysis_set.rds"))
analysis_set <- d$endpoint
counts <- d$counts      # complete individuals per study and GE day; counted: at least 5
studies <- d$studies %>% distinct(study, season, cells)
study_label <- with(studies, setNames(sprintf("%s  (%s, %s)", study, season, cells), study))
totals <- counts %>% filter(counted) %>% group_by(ge_day) %>% summarise(K = n(), N = sum(n))
totals

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
  ) +       # grey: study not counted (too few individuals)
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

rm(list = ls())
