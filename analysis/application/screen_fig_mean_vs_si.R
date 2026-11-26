# Screen of all gene sets at day 7, figure: pooled delta of the mean vs the surrogate index (ridge), one point per gene set
# Input: output/application/screen_all_genesets_day7.rds (from screen_all_genesets.R). Run from the project root.

library(tidyverse)

epsilon <- 0.2
alpha <- 0.05

results <- readRDS(fs::path("output", "application", "screen_all_genesets_day7.rds"))

# One row per gene set: pooled delta for each surrogate, and which surrogate(s) pass after BH correction
wide <- results %>%
  mutate(surrogate = recode(surrogate, Mean = "mean", `SI (ridge)` = "si")) %>%
  select(geneset, surrogate, delta, p_adjusted) %>%
  pivot_wider(names_from = surrogate, values_from = c(delta, p_adjusted)) %>%
  mutate(pass = case_when(p_adjusted_mean < alpha & p_adjusted_si < alpha ~ "Both",
                          p_adjusted_mean < alpha ~ "Mean only",
                          p_adjusted_si < alpha ~ "SI only",
                          TRUE ~ "Neither"),
         pass = factor(pass, c("Neither", "Mean only", "SI only", "Both")))
count(wide, pass)

p <- ggplot(wide, aes(delta_mean, delta_si)) +
  annotate("rect", xmin = -epsilon, xmax = epsilon, ymin = -Inf, ymax = Inf, fill = "grey93") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -epsilon, ymax = epsilon, fill = "grey93") +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") +
  geom_point(aes(colour = pass), size = 1.3) +
  scale_colour_manual(values = c(Neither = "grey65", `Mean only` = "#f47230", `SI only` = "#257de0", Both = "black"),
                      name = "Passes (PI test, BH)", drop = FALSE) +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(x = expression("Pooled " * delta * ", mean"), y = expression("Pooled " * delta * ", SI (ridge)")) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "black"),
        legend.position = "right")
p

for (d in c(fs::path("output", "figures", "application"), fs::path("manuscript", "latex", "figures"))) {
  fs::dir_create(d)
  ggsave(fs::path(d, "screen_mean_vs_si.pdf"), p, width = 4.8, height = 3.6)
}
