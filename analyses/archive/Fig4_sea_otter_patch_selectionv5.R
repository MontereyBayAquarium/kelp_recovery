################################################################################
# Sea otter energetics and prey use — Panels A & B
#   Panel A: SOFA posterior – PUR foraging effort across patch types
#   Panel B: Focal bouts by prey group and patch type
# Author: Joshua G. Smith — UCSC Nearshore Ecology Research Group
################################################################################

rm(list = ls())

require(librarian)
librarian::shelf(
  tidyverse, janitor, readxl, here,
  ggridges, patchwork, sf, lubridate
)

################################################################################
# 1. Load SOFA data ------------------------------------------------------------
################################################################################

sofa_path <- here::here("output", "sofa", "SOFA_Results_2025-10-31.xlsx")

eta_df <- read_xlsx(sofa_path, sheet = "propforagingeffort")
cr_df  <- read_xlsx(sofa_path, sheet = "meanbiomassconsumption")

################################################################################
# 2. Filter to purple urchins and map to patch types ---------------------------
################################################################################

recode_patches <- function(df) {
  df %>%
    filter(prey_name %in% c("urchin_b", "urchin_f", "urchin_i")) %>%
    mutate(
      patch_type = case_when(
        prey_name == "urchin_b" ~ "BAR",
        prey_name == "urchin_f" ~ "FOR",
        prey_name == "urchin_i" ~ "INCIP"
      ),
      # order for ridge y-axis (FOR, BAR, INCIP)
      patch_type = factor(patch_type, levels = c("FOR", "BAR", "INCIP"))
    )
}

eta_focus <- recode_patches(eta_df)
cr_focus  <- recode_patches(cr_df)

################################################################################
# 3. Remove low-value duplicate summary rows (keep true consumption rates) -----
################################################################################

cr_focus <- cr_focus %>%
  group_by(patch_type) %>%
  slice_max(mean, n = 1, with_ties = FALSE) %>%
  ungroup()

################################################################################
# 4. Expand posterior draws from mean + SD -------------------------------------
################################################################################

expand_posterior <- function(tbl, n = 5000) {
  tbl %>%
    rowwise() %>%
    mutate(draws = list(rnorm(n, mean, sd))) %>%
    unnest(draws)
}

eta_draws <- expand_posterior(eta_focus)
cr_draws  <- expand_posterior(cr_focus)

################################################################################
# 5. Summarise 90% credible intervals ------------------------------------------
################################################################################

posterior_summary <- function(df) {
  df %>%
    group_by(patch_type) %>%
    summarise(
      med = median(draws),
      l90 = quantile(draws, 0.05),
      u90 = quantile(draws, 0.95),
      .groups = "drop"
    )
}

eta_sum <- posterior_summary(eta_draws)
cr_sum  <- posterior_summary(cr_draws)

################################################################################
# Fig 4: Panel A (posterior ridgelines) + Panel B (urchin metrics B/C/D stacked)
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(patchwork)
library(purrr)
library(here)

################################################################################
# 1) Shared aesthetics ----------------------------------------------------------
################################################################################

patch_colors <- c(
  "BAR"   = "#7570B3",
  "INCIP" = "#D95F02",
  "FOR"   = "#1B9E77"
)

theme_panel <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid   = element_blank(),
      axis.text    = element_text(color = "black"),
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8),
      plot.tag     = element_text(size = 10, color = "black")
    )
}

################################################################################
# 2) Posterior plot function (Panel A) ------------------------------------------
################################################################################

plot_posterior <- function(draw_df, sum_df, xlab, tag_title, show_legend = FALSE) {
  
  ggplot() +
    geom_density_ridges(
      data = draw_df,
      aes(x = draws, y = patch_type, fill = patch_type),
      alpha = 0.85, scale = 1.1, color = "black"
    ) +
    geom_errorbar(
      data = sum_df,
      aes(y = patch_type, xmin = l90, xmax = u90),
      orientation = "y", height = 0.25,
      linewidth = 1.2, color = "black"
    ) +
    geom_point(
      data = sum_df,
      aes(y = patch_type, x = med),
      color = "white", fill = "black",
      shape = 21, size = 3, stroke = 1
    ) +
    scale_fill_manual(values = patch_colors, name = "Patch type") +
    labs(x = xlab, y = NULL, tag = tag_title) +
    theme_panel(base_size = 11) +
    theme(
      axis.text.y     = element_blank(),
      legend.position = if (show_legend) "right" else "none"
    )
}

################################################################################
# 3) Load processed survey data ------------------------------------------------
################################################################################

dissect_with_pred <- readRDS("output/survey_data/processed/dissect_with_pred.rds")
sizefq_with_pred  <- readRDS("output/survey_data/processed/sizefq_with_pred.rds")

lvl <- c("FOR", "BAR", "INCIP")

dissect_with_pred <- dissect_with_pred %>%
  mutate(site_type_predicted = factor(site_type_predicted, levels = lvl))

sizefq_with_pred <- sizefq_with_pred %>%
  mutate(site_type_predicted = factor(site_type_predicted, levels = lvl))

################################################################################
# 4) Panel A -------------------------------------------------------------------
################################################################################

pA <- plot_posterior(
  eta_draws,
  eta_sum,
  xlab        = "Foraging effort (proportion of sea \nurchins in patch vs. all other prey)",
  tag_title   = "A",
  show_legend = FALSE
)

################################################################################
# 5) Panel B: Urchin metrics (Purple only) --------------------------------------
################################################################################

# Precompute means for colored vertical lines
gi_means <- dissect_with_pred %>%
  filter(species == "purple_urchin") %>%
  drop_na(site_type_predicted, gonad_index) %>%
  group_by(site_type_predicted) %>%
  summarise(mean_val = mean(gonad_index), .groups = "drop")

gm_means <- dissect_with_pred %>%
  filter(species == "purple_urchin") %>%
  drop_na(site_type_predicted, gonad_mass_g) %>%
  group_by(site_type_predicted) %>%
  summarise(mean_val = mean(gonad_mass_g), .groups = "drop")

# B) Gonad index
pB <- ggplot(
  dissect_with_pred %>%
    filter(species == "purple_urchin") %>%
    drop_na(),
  aes(x = gonad_index, y = site_type_predicted, fill = site_type_predicted)
) +
  geom_density_ridges(alpha = 0.85, scale = 1.1,
                      color = "black", linewidth = 0.3) +
  geom_vline(
    data = gi_means,
    aes(xintercept = mean_val, color = site_type_predicted),
    linewidth = 0.9
  ) +
  scale_fill_manual(values = patch_colors) +
  scale_color_manual(values = patch_colors, guide = "none") +
  coord_cartesian(xlim = c(0, 23)) +
  labs(x = "Gonad index", y = NULL, tag = "B") +
  theme_panel(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    legend.position = "none"
  )

# C) Gonad mass
pC <- ggplot(
  dissect_with_pred %>%
    filter(species == "purple_urchin") %>%
    drop_na(),
  aes(x = gonad_mass_g, y = site_type_predicted, fill = site_type_predicted)
) +
  geom_density_ridges(alpha = 0.85, scale = 1.1,
                      color = "black", linewidth = 0.3) +
  geom_vline(
    data = gm_means,
    aes(xintercept = mean_val, color = site_type_predicted),
    linewidth = 0.9
  ) +
  scale_fill_manual(values = patch_colors) +
  scale_color_manual(values = patch_colors, guide = "none") +
  coord_cartesian(xlim = c(0, 6)) +
  labs(x = "Gonad mass (g)", y = NULL, tag = "C") +
  theme_panel(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    legend.position = "none"
  )

# D) Size-frequency (step/line, Purple only)
size_curve <- sizefq_with_pred %>%
  filter(
    species == "Purple",
    !is.na(site_type_predicted),
    !is.na(size_cm),
    !is.na(count)
  ) %>%
  group_by(site_type_predicted, size_cm) %>%
  summarise(n = sum(count), .groups = "drop")

pD <- ggplot(
  size_curve,
  aes(x = size_cm, y = n, color = site_type_predicted)
) +
  geom_step(linewidth = 1.1) +
  scale_color_manual(values = patch_colors, name = "Patch type") +
  scale_x_continuous(limits = c(1, 8), breaks = 1:8) +
  labs(x = "Size (sm)", y = NULL, tag = "D") +
  theme_panel(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    legend.position = "right"
  )

# Stack B/C/D
panelB <- (pB / pC / pD) +
  plot_layout(ncol = 1, heights = c(1, 1, 1), guides = "collect")

################################################################################
# 6) Combine Panel A + Panel B -------------------------------------------------
################################################################################

energetics_combined <- pA + panelB +
  #plot_layout(ncol = 2, widths = c(1, 1)) &
  plot_layout(ncol = 2, widths = c(1.2, 0.8)) &
  theme(plot.margin = unit(c(2, 2, 2, 2), "pt"))

energetics_combined

################################################################################
# 7) Save ----------------------------------------------------------------------
################################################################################

ggsave(
  here::here("figures", "Fig4_effort_allocationv2.png"),
  energetics_combined,
  width = 8, height = 8, dpi = 600, bg = "white"
)
