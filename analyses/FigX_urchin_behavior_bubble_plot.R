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

