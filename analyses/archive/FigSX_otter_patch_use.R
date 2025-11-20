################################################################################
# Supplementary Figure S1: Sea otter occupancy and behavioral composition
# (aggregated across years)
# Joshua G. Smith – UCSC Nearshore Ecology Research Group
################################################################################

# --- 1. Load packages ---------------------------------------------------------
rm(list = ls())
require(librarian)
librarian::shelf(
  tidyverse, janitor, lubridate, sf,
  gridExtra, RColorBrewer, here
)

datdir <- here::here("output")
load(here::here("output", "survey_data", "processed", "zone_level_data4.rda"))
scan_orig <- read_csv(file.path(here::here("output","scans","scans_data.csv")))
dissection_data <- read_csv("/Volumes/enhydra/data/kelp_recovery/MBA_kelp_forest_database/processed/dissection/dissection_data_cleanedv2.csv")
forage_orig <- read_csv("/Volumes/enhydra/data/foraging_data/processed/foraging_data_2024_2025_processed.csv")

years_keep <- c(2024, 2025)

# --- 2. Clean and prepare behavioral data ------------------------------------
scan_clean <- scan_orig %>%
  clean_names() %>%
  mutate(
    lat  = as.numeric(lat),
    long = as.numeric(long),
    year = year(date),
    quarter = lubridate::quarter(date),
    behavior_group = ifelse(behav == "foraging", "Foraging", "Other")
  ) %>%
  drop_na(lat, long, ind, behavior_group)

# Convert to sf
scan_sf <- scan_clean %>%
  st_as_sf(coords = c("long", "lat"), crs = st_crs(quad_build3))

# Join scans to patch polygons
pts_joined <- st_join(scan_sf, quad_build3, join = st_intersects, left = TRUE) %>%
  mutate(
    year = lubridate::year(date),
    pred_patch = case_when(
      str_detect(tolower(pred_patch), "bar")   ~ "BAR",
      str_detect(tolower(pred_patch), "incip") ~ "INCIP",
      str_detect(tolower(pred_patch), "for")   ~ "FOR",
      TRUE ~ toupper(pred_patch)
    )
  ) %>%
  filter(year %in% years_keep)

# Standardize patch and behavior fields
pts_behav <- pts_joined %>%
  st_drop_geometry() %>%
  mutate(
    pred_patch = case_when(
      str_detect(tolower(pred_patch), "bar")   ~ "Barren",
      str_detect(tolower(pred_patch), "incip") ~ "Incipient",
      str_detect(tolower(pred_patch), "for")   ~ "Forest",
      TRUE ~ NA_character_
    ),
    pred_patch = factor(pred_patch, levels = c("Barren", "Incipient", "Forest")),
    behav = str_to_title(str_trim(behav))
  ) %>%
  filter(!is.na(pred_patch), !is.na(behav), !is.na(ind), ind > 0) %>%
  filter(year %in% years_keep)

# --- 3. Compute otter occupancy per patch-year --------------------------------
patch_area_df <- quad_build3 %>%
  st_make_valid() %>%
  mutate(area_m2 = as.numeric(st_area(geometry))) %>%
  st_drop_geometry() %>%
  distinct(patch_id, area_m2)

occupancy_df <- pts_behav %>%
  group_by(year, patch_id, pred_patch) %>%
  summarise(
    mean_otters = mean(ind, na.rm = TRUE),
    n_scans = n(),
    .groups = "drop"
  ) %>%
  left_join(patch_area_df, by = "patch_id") %>%
  mutate(otters_per_km2 = mean_otters / (area_m2 / 1e6))

# --- 4. Summarize otter density across years ----------------------------------
occupancy_summary <- occupancy_df %>%
  group_by(pred_patch) %>%
  summarise(
    mean_density = mean(otters_per_km2, na.rm = TRUE),
    se_density   = sd(otters_per_km2, na.rm = TRUE) / sqrt(n()),
    n_patches    = n(),
    .groups = "drop"
  )

# --- 5. Summarize behavioral composition (pooled across years) ----------------
behavior_summary <- pts_behav %>%
  count(pred_patch, behav) %>%
  group_by(pred_patch) %>%
  mutate(prop = n / sum(n))

# --- 6. Define color schemes --------------------------------------------------
patch_cols <- c("Barren" = "purple", "Incipient" = "orange", "Forest" = "forestgreen")

behav_levels <- sort(unique(behavior_summary$behav))
behavior_cols <- setNames(
  RColorBrewer::brewer.pal(n = max(3, min(length(behav_levels), 8)), "Set2"),
  behav_levels
)

# --- 7. Plot Panel A: Otter occupancy (density) -------------------------------
gA <- ggplot(occupancy_summary, aes(x = pred_patch, y = mean_density, fill = pred_patch)) +
  geom_col(color = "black", width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_density - se_density, ymax = mean_density + se_density),
    width = 0.2
  ) +
  scale_fill_manual(values = patch_cols) +
  labs(
    x = "Patch type",
    y = expression("Mean otter density (individuals km"^{-2}*")"),
    title = "A. Otter occupancy across patch types"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0)
  )

# --- 8. Plot Panel B: Behavioral composition ----------------------------------
gB <- ggplot(behavior_summary, aes(x = pred_patch, y = prop, fill = behav)) +
  geom_col(color = "black", width = 0.7, position = "stack") +
  scale_fill_manual(values = behavior_cols, name = "Behavior") +
  labs(
    x = "Patch type",
    y = "Proportion of observed behaviors",
    title = "B. Behavioral composition by patch type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0)
  )

# --- 9. Combine panels into Supplementary Figure S1 ---------------------------
gridExtra::grid.arrange(gA, gB, ncol = 2, widths = c(0.9, 1.2))

# --- 10. Save figure ----------------------------------------------------------
ggsave(
  here("figures", "Figure_S1_otter_behavior_patchtypes.png"),
  plot = gridExtra::arrangeGrob(gA, gB, ncol = 2, widths = c(0.9, 1.2)),
  width = 10, height = 5, dpi = 600, bg = "white"
)
