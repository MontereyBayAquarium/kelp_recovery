################################################################################
# 
################################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

require(librarian)
shelf(
  tidyverse, janitor, lubridate, sf, here,
  randomForest, patchwork, ggrepel, tidytext
)

# Reproducibility
set.seed(1985)

################################################################################
# 0. Setup / load data ---------------------------------------------------------
################################################################################

load(here::here("output", "survey_data", "processed", "zone_level_data4.rda"))

years_keep <- c(2024, 2025)

patch_colors <- c(
  "BAR"   = "#7570B3",
  "INCIP" = "#D95F02",
  "FOR"   = "#1B9E77"
)

str(quad_build3)





################################################################################
# Prepare data
################################################################################

plot_dat <- quad_build3 %>%
  filter(year %in% years_keep) %>%
  mutate(
    exposed_urchin_densitym2 = purple_urchin_densitym2 - purple_urchin_conceiledm2,
    prop_exposed = if_else(
      purple_urchin_densitym2 > 0,
      exposed_urchin_densitym2 / purple_urchin_densitym2,
      NA_real_
    )
  )

################################################################################
# Plot
################################################################################

p <- ggplot(
  plot_dat,
  aes(
    x = purple_urchin_densitym2,
    y = prop_exposed,
    size = macro_stipe_density_20m2,
    color = mean_gi
  )
) +
  geom_point(alpha = 0.85) +
  scale_size_continuous(
    name = "Macro stipe density\n(20 m2)",
    range = c(2, 18)
  ) +
  scale_color_viridis_c(name = "Mean GI") +
  labs(
    x = expression("Purple urchin density (m"^"-2" * ")"),
    y = "Proportion exposed purple urchins"
  ) +
  theme_classic()

p

################################################################################
#Save

ggsave(
  "~/Downloads/urchin_behavior_bubble_plot.png",
  plot = p,
  width = 7,
  height = 5,
  dpi = 300
)










################################################################################
# Calculate prop exposed and year-to-year change in kelp density
################################################################################

plot_dat <- quad_build3 %>%
  filter(year %in% c(2024, 2025)) %>%
  mutate(
    exposed_urchin_densitym2 = purple_urchin_densitym2 - purple_urchin_conceiledm2,
    prop_exposed = if_else(
      purple_urchin_densitym2 > 0,
      exposed_urchin_densitym2 / purple_urchin_densitym2,
      NA_real_
    )
  ) %>%
  st_drop_geometry()

delta_dat <- plot_dat %>%
  select(
    patch_id, site, site_type, year,
    prop_exposed, purple_urchin_densitym2, mean_gi,
    macro_stipe_density_20m2, total_biomass_g
  ) %>%
  pivot_wider(
    names_from = year,
    values_from = c(
      prop_exposed, purple_urchin_densitym2, mean_gi,
      macro_stipe_density_20m2, total_biomass_g
    ),
    names_sep = "_"
  ) %>%
  mutate(
    delta_macro_stipe_density_20m2 =
      macro_stipe_density_20m2_2025 - macro_stipe_density_20m2_2024,
    delta_total_biomass_g =
      total_biomass_g_2025 - total_biomass_g_2024,
    delta_prop_exposed =
      prop_exposed_2025 - prop_exposed_2024,
    delta_mean_gi =
      mean_gi_2025 - mean_gi_2024,
    delta_purple_urchin_densitym2 =
      purple_urchin_densitym2_2025 - purple_urchin_densitym2_2024
  )




p_delta <- ggplot(
  delta_dat,
  aes(
    x = prop_exposed_2024,
    y = delta_macro_stipe_density_20m2,
    size = purple_urchin_densitym2_2024,
    color = mean_gi_2024
  )
) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_point(alpha = 0.85) +
  scale_size_continuous(
    name = "Purple urchin density 2024",
    range = c(2, 14)
  ) +
  scale_color_viridis_c(name = "Mean GI 2024") +
  labs(
    x = "Proportion exposed purple urchins in 2024",
    y = expression(Delta * " Macro stipe density (2025 - 2024)")
  ) +
  theme_classic()

p_delta



ggsave(
  "~/Downloads/urchin_bubble_plot_delta.png",
  plot = p_delta,
  width = 7,
  height = 5,
  dpi = 300
)






