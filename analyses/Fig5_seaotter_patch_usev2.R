################################################################################
# Supplementary Figure S1: Sea otter occupancy and behavior (Panels A & B)
# Joshua G. Smith – UCSC Nearshore Ecology Research Group
################################################################################

# --- 1. Load packages ---------------------------------------------------------
rm(list = ls())
require(librarian)
librarian::shelf(
  tidyverse, janitor, lubridate, sf,
  gridExtra, RColorBrewer, here
)

# --- 2. Load data -------------------------------------------------------------
load(here::here("output", "survey_data", "processed", "zone_level_data4.rda"))
scan_orig <- readr::read_csv(
  file.path(here::here("output", "scans", "scans_data.csv"))
)

years_keep <- c(2024, 2025)
patch_cols <- c("Barren"   = "#7570B3",
                "Incipient"= "#D95F02",
                "Forest"   = "#1B9E77")

################################################################################
# 3. Process behavioral scan data
################################################################################

scan_clean <- scan_orig %>%
  clean_names() %>%
  mutate(
    lat          = as.numeric(lat),
    long         = as.numeric(long),
    year         = lubridate::year(date),
    quarter      = lubridate::quarter(date),
    behavior_group = ifelse(behav == "foraging", "Foraging", "Other")
  ) %>%
  drop_na(lat, long, ind, behavior_group)

scan_sf <- st_as_sf(
  scan_clean,
  coords = c("long", "lat"),
  crs    = st_crs(quad_build3)
)

pts_joined <- st_join(
  scan_sf,
  quad_build3,
  join = st_intersects,
  left = TRUE
) %>%
  mutate(
    year = lubridate::year(date),
    pred_patch = case_when(
      str_detect(tolower(pred_patch), "bar")   ~ "BAR",
      str_detect(tolower(pred_patch), "incip") ~ "INCIP",
      str_detect(tolower(pred_patch), "for")   ~ "FOR",
      TRUE                                     ~ toupper(pred_patch)
    )
  ) %>%
  filter(year %in% years_keep)

pts_behav <- pts_joined %>%
  st_drop_geometry() %>%
  mutate(
    pred_patch = case_when(
      str_detect(tolower(pred_patch), "bar")   ~ "Barren",
      str_detect(tolower(pred_patch), "incip") ~ "Incipient",
      str_detect(tolower(pred_patch), "for")   ~ "Forest",
      TRUE                                     ~ NA_character_
    ),
    pred_patch = factor(pred_patch, levels = c("Barren", "Incipient", "Forest")),
    behav      = stringr::str_to_title(stringr::str_trim(behav))
  ) %>%
  filter(
    !is.na(pred_patch),
    !is.na(behav),
    !is.na(ind),
    ind > 0
  )

################################################################################
# 4. Otter occupancy by patch type (Panel A)
################################################################################

patch_area_df <- quad_build3 %>%
  st_make_valid() %>%
  mutate(area_m2 = as.numeric(st_area(geometry))) %>%
  st_drop_geometry() %>%
  distinct(patch_id, area_m2)

occupancy_df <- pts_behav %>%
  group_by(year, patch_id, pred_patch) %>%
  summarise(
    mean_otters = mean(ind, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  left_join(patch_area_df, by = "patch_id") %>%
  mutate(otters_per_km2 = mean_otters / (area_m2 / 1e6))

occupancy_summary <- occupancy_df %>%
  group_by(pred_patch) %>%
  summarise(
    mean_density = mean(otters_per_km2, na.rm = TRUE),
    se_density   = sd(otters_per_km2, na.rm = TRUE) / sqrt(n()),
    .groups      = "drop"
  )

################################################################################
# 5. Behavioral composition (Panel B)
################################################################################

behavior_summary <- pts_behav %>%
  count(pred_patch, behav) %>%
  group_by(pred_patch) %>%
  mutate(prop = n / sum(n))

################################################################################
# 6. Plot panels A & B
################################################################################

# Custom theme (matches Josh's usual style)
my_theme <- theme(
  axis.text.x      = ggplot2::element_text(size = 8, color = "black"),
  axis.text.y      = ggplot2::element_text(size = 8, color = "black"),
  axis.title       = ggplot2::element_text(size = 10, color = "black"),
  legend.text      = ggplot2::element_text(size = 8, color = "black"),
  legend.title     = ggplot2::element_text(size = 8, color = "black"),
  plot.tag         = ggplot2::element_text(size = 10, color = "black"),
  panel.grid       = ggplot2::element_blank(),
  panel.background = ggplot2::element_blank(),
  axis.line        = ggplot2::element_line(colour = "black"),
  legend.key       = ggplot2::element_blank(),
  plot.title       = ggplot2::element_text(face = "plain", size = 11, hjust = 0),
  legend.background = ggplot2::element_rect(
    fill = scales::alpha("blue", 0)
  ),
  strip.text       = ggplot2::element_text(
    size = 10, face = "bold", color = "black", hjust = 0
  ),
  strip.background = ggplot2::element_blank()
)

# --- Panel A: Otter occupancy -------------------------------------------------
gA <- ggplot(
  occupancy_summary,
  aes(x = pred_patch, y = mean_density, fill = pred_patch)
) +
  geom_col(color = "black", width = 0.6) +
  geom_errorbar(
    aes(
      ymin = mean_density - se_density,
      ymax = mean_density + se_density
    ),
    width = 0.2
  ) +
  scale_fill_manual(values = patch_cols, name = "Patch type") +
  labs(
    x     = "Patch type",
    y     = expression("Mean sea otter density (individuals km"^{-2}*")"),
    title = "A. Sea otter occupancy by patch type"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  ) + my_theme

# --- Panel B: Behavioral composition -----------------------------------------
gB <- ggplot(
  behavior_summary,
  aes(x = pred_patch, y = prop, fill = behav)
) +
  geom_col(color = "black", width = 0.7, position = "stack") +
  scale_fill_manual(
    values = RColorBrewer::brewer.pal(n = 8, "Set2"),
    name   = "Behavior"
  ) +
  labs(
    x     = "Patch type",
    y     = "Proportion of observed behaviors",
    title = "B. Behavioral composition by patch type"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  ) + my_theme

################################################################################
# 7. Combine panels A & B and save
################################################################################

combined_fig <- gridExtra::grid.arrange(
  gA, gB,
  ncol   = 2,
  widths = c(1, 1.1)
)

ggsave(
  here::here("figures", "Fig5_patch_utilization_AB.png"),
  combined_fig,
  width  = 10,
  height = 4,
  dpi    = 600,
  bg     = "white"
)
