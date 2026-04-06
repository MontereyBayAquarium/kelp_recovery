################################################################################
# Fig 4: Sea otter energetics and prey use
#   A = SOFA posterior for PUR foraging effort across patch types
#   B = Mean gonad index by predicted patch type
#   C = Purple urchin biomass density by predicted patch type (g m^-2)
#   D = Purple urchin size-frequency by predicted patch type
################################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

require(librarian)
librarian::shelf(
  tidyverse, janitor, readxl, here,
  ggridges, patchwork, sf, lubridate
)

################################################################################
# 1) Load all data --------------------------------------------------------------
################################################################################

# SOFA results
sofa_path <- here::here("output", "sofa", "SOFA_Results_2025-10-31.xlsx")
eta_df    <- read_xlsx(sofa_path, sheet = "propforagingeffort")
cr_df     <- read_xlsx(sofa_path, sheet = "meanbiomassconsumption")

# Zone-level data (.rda file)
loaded_objs <- load(here::here("output", "survey_data", "processed", "zone_level_data4.rda"))

if (!"quad_build3" %in% loaded_objs) {
  stop("quad_build3 was not found in zone_level_data4.rda")
}

# Size-frequency data
sizefq_with_pred <- readRDS(
  here::here("output", "survey_data", "processed", "sizefq_with_pred.rds")
)

################################################################################
# 2) Shared aesthetics ----------------------------------------------------------
################################################################################

patch_colors <- c(
  "BAR"   = "#7570B3",
  "INCIP" = "#D95F02",
  "FOR"   = "#1B9E77"
)

patch_levels <- c("FOR", "BAR", "INCIP")

theme_panel <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid       = element_blank(),
      axis.text        = element_text(color = "black"),
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8),
      plot.tag         = element_text(size = 10, color = "black"),
      strip.background = element_blank(),
      strip.text       = element_text(color = "black")
    )
}

################################################################################
# 3) Prep SOFA data -------------------------------------------------------------
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
      patch_type = factor(patch_type, levels = patch_levels)
    )
}

eta_focus <- recode_patches(eta_df)
cr_focus  <- recode_patches(cr_df)

cr_focus <- cr_focus %>%
  group_by(patch_type) %>%
  slice_max(mean, n = 1, with_ties = FALSE) %>%
  ungroup()

expand_posterior <- function(tbl, n = 5000) {
  tbl %>%
    rowwise() %>%
    mutate(draws = list(rnorm(n, mean, sd))) %>%
    unnest(draws) %>%
    ungroup()
}

eta_draws <- expand_posterior(eta_focus)
cr_draws  <- expand_posterior(cr_focus)

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
# 4) Prep quad data for panels B and C -----------------------------------------
################################################################################

quad_plot_dat <- quad_build3 %>%
  st_drop_geometry() %>%
  mutate(
    patch_type = factor(pred_patch, levels = patch_levels),
    biomass_density_g_m2 = total_biomass_g / 80
  ) %>%
  filter(!is.na(patch_type))

################################################################################
# 5) Prep size-frequency data for panel D --------------------------------------
################################################################################

sizefq_with_pred <- sizefq_with_pred %>%
  mutate(
    site_type_predicted = factor(site_type_predicted, levels = patch_levels)
  )

size_curve <- sizefq_with_pred %>%
  filter(
    species == "Purple",
    !is.na(site_type_predicted),
    !is.na(size_cm),
    !is.na(count)
  ) %>%
  group_by(site_type_predicted, size_cm) %>%
  summarise(n = sum(count), .groups = "drop")

################################################################################
# 6) Plot functions -------------------------------------------------------------
################################################################################

plot_posterior <- function(draw_df, sum_df, xlab, tag_title, show_legend = FALSE) {
  
  ggplot() +
    geom_density_ridges(
      data = draw_df,
      aes(x = draws, y = patch_type, fill = patch_type),
      alpha = 0.85,
      scale = 1.1,
      color = "black",
      linewidth = 0.3
    ) +
    geom_errorbar(
      data = sum_df,
      aes(y = patch_type, xmin = l90, xmax = u90),
      orientation = "y",
      height = 0.25,
      linewidth = 1.0,
      color = "black"
    ) +
    geom_point(
      data = sum_df,
      aes(y = patch_type, x = med),
      color = "white",
      fill = "black",
      shape = 21,
      size = 2.8,
      stroke = 0.9
    ) +
    scale_fill_manual(values = patch_colors, name = "Patch type") +
    labs(x = xlab, y = NULL, tag = tag_title) +
    theme_panel(base_size = 11) +
    theme(
      axis.text.y     = element_blank(),
      legend.position = if (show_legend) "right" else "none"
    )
}

plot_box <- function(df, xvar, xlab, tag_title, xlimits = NULL) {
  
  p <- ggplot(
    df,
    aes(x = .data[[xvar]], y = patch_type, fill = patch_type)
  ) +
    geom_boxplot(
      width = 0.65,
      outlier.shape = 16,
      outlier.size = 1.1,
      linewidth = 0.35,
      color = "black"
    ) +
    scale_fill_manual(values = patch_colors) +
    labs(x = xlab, y = NULL, tag = tag_title) +
    theme_panel(base_size = 11) +
    theme(
      axis.text.y = element_blank(),
      legend.position = "none"
    )
  
  if (!is.null(xlimits)) {
    p <- p + coord_cartesian(xlim = xlimits)
  }
  
  p
}

################################################################################
# 7) Build panels ---------------------------------------------------------------
################################################################################

# Panel A
pA <- plot_posterior(
  draw_df     = eta_draws,
  sum_df      = eta_sum,
  xlab        = "Foraging effort (proportion of sea\nurchins in patch vs. all other prey)",
  tag_title   = "A",
  show_legend = FALSE
)

# Panel B
pB <- quad_plot_dat %>%
  drop_na(mean_gi) %>%
  plot_box(
    xvar      = "mean_gi",
    xlab      = "Mean gonad index \n (ratio of gonad mass to animal mass)",
    tag_title = "B"
  )

# Panel C
pC <- quad_plot_dat %>%
  drop_na(biomass_density_g_m2) %>%
  filter(biomass_density_g_m2 < 1000) %>%
  plot_box(
    xvar      = "biomass_density_g_m2",
    xlab      = expression("Biomass density (g " * m^{-2} * ")"),
    tag_title = "C"
  )

# Panel D
pD <- ggplot(
  size_curve,
  aes(x = size_cm, y = n, color = site_type_predicted)
) +
  geom_step(linewidth = 1.1) +
  scale_color_manual(values = patch_colors, name = "Patch type") +
  scale_x_continuous(limits = c(1, 8), breaks = 1:8) +
  labs(x = "Size frequency (cm)", y = NULL, tag = "D") +
  theme_panel(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    legend.position = "right"
  )

################################################################################
# 8) Combine panels -------------------------------------------------------------
################################################################################

panelBCD <- (pB / pC / pD) +
  plot_layout(ncol = 1, heights = c(1, 1, 1), guides = "collect")

energetics_combined <- pA + panelBCD +
  plot_layout(ncol = 2, widths = c(1.2, 0.8)) &
  theme(plot.margin = unit(c(2, 2, 2, 2), "pt"))

energetics_combined

################################################################################
# 9) Save -----------------------------------------------------------------------
################################################################################

ggsave(
  here::here("figures", "Fig4_effort_allocation_v6.png"),
  energetics_combined,
  width  = 8,
  height = 8,
  dpi    = 600,
  bg     = "white"
)


################################################################################
# Summary stats for manuscript (Fig 4A) ----------------------------------------
################################################################################

eta_sum %>%
  mutate(
    med_pct = med * 100,
    l90_pct = l90 * 100,
    u90_pct = u90 * 100
  ) %>%
  select(patch_type, med_pct, l90_pct, u90_pct) %>%
  arrange(patch_type) %>%
  print()


################################################################################
# Summary stats for gonad index (Fig 4B) ----------------------------------------
################################################################################

gi_summary <- quad_plot_dat %>%
  dplyr::filter(!is.na(mean_gi)) %>%
  dplyr::group_by(patch_type) %>%
  dplyr::summarise(
    n      = sum(is.finite(mean_gi)),
    mean   = mean(mean_gi, na.rm = TRUE),
    sd     = sd(mean_gi, na.rm = TRUE),
    se     = sd / sqrt(n),
    median = median(mean_gi, na.rm = TRUE),
    q25    = quantile(mean_gi, 0.25, na.rm = TRUE),
    q75    = quantile(mean_gi, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(factor(patch_type, levels = c("BAR", "INCIP", "FOR")))

gi_summary

biomass_summary <- quad_plot_dat %>%
  dplyr::filter(!is.na(biomass_density_g_m2)) %>%
  dplyr::group_by(patch_type) %>%
  dplyr::summarise(
    n      = sum(is.finite(biomass_density_g_m2)),
    mean   = mean(biomass_density_g_m2, na.rm = TRUE),
    sd     = sd(biomass_density_g_m2, na.rm = TRUE),
    se     = sd / sqrt(n),
    median = median(biomass_density_g_m2, na.rm = TRUE),
    q25    = quantile(biomass_density_g_m2, 0.25, na.rm = TRUE),
    q75    = quantile(biomass_density_g_m2, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(factor(patch_type, levels = c("BAR", "INCIP", "FOR")))

biomass_summary
