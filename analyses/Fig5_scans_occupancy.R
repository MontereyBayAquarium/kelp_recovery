# jogsmith@ucsc.edu

rm(list = ls())

require(librarian)
librarian::shelf(
  tidyverse, lubridate, sf, stringr, purrr, terra, janitor,
  rnaturalearth, rnaturalearthdata, ggspatial, here
)

datadir  <- "/Volumes/enhydra/data/kelp_recovery/"
figdir   <- here::here("figures")

# load benthic survey data
load(file.path(datadir, "MBA_kelp_forest_database/processed/recovery/kelp_recovery_data.rda"))

# load scans
scan_orig <- read_csv(
  file.path(here::here("output", "scans", "scans_data.csv")),
  show_col_types = FALSE
)

# load site patches
site_patches <- st_read(
  here::here("output", "gis_data", "processed", "site_patch_polygons.shp"),
  quiet = TRUE
)

# load LDA-predicted patch types
load(here::here("output", "lda_patch_transitionsv5.rda"))  # transitions_tbl_constrained

################################################################################
# Helpers
################################################################################

get_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

################################################################################
# Step 1: Average to site, zone, site_type for each year
################################################################################

kelp_avg <- kelp_data %>%
  dplyr::select(-macro_stipe_sd_20m2) %>%
  dplyr::group_by(site, site_type, latitude, longitude, zone, survey_date) %>%
  dplyr::summarise(
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::select(-transect)

quad_avg <- quad_data %>%
  dplyr::group_by(site, site_type, latitude, longitude, zone, survey_date) %>%
  dplyr::summarise(
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::select(-quadrat, -transect)

dat_agg <- kelp_avg %>%
  dplyr::inner_join(
    quad_avg,
    by = c("site", "site_type", "latitude", "longitude", "zone", "survey_date"),
    suffix = c("_kelp", "_quad")
  )

quad_zone <- dat_agg %>%
  dplyr::group_by(latitude, longitude, site, site_type, survey_date, zone) %>%
  dplyr::summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

################################################################################
# Step 2: assign model-predicted patch types
################################################################################

quad_zone_with_pred <- quad_zone %>%
  dplyr::left_join(
    transitions_tbl_constrained %>%
      dplyr::select(site, site_type, zone, patch_2024, patch_2025),
    by = c("site", "site_type", "zone")
  ) %>%
  dplyr::mutate(
    pred_patch = dplyr::case_when(
      format(survey_date, "%Y") == "2024" ~ as.character(patch_2024),
      format(survey_date, "%Y") == "2025" ~ as.character(patch_2025),
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(-patch_2024, -patch_2025)

################################################################################
# Step 3: join survey points to patch geometry
################################################################################

quad_zone_sf <- st_as_sf(
  quad_zone_with_pred,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

site_patches_single <- site_patches %>%
  st_cast("POLYGON", warn = FALSE) %>%
  dplyr::mutate(patch_id = dplyr::row_number())

site_patches_with_points <- site_patches_single %>%
  st_join(quad_zone_sf, join = st_intersects, left = TRUE)

quad_build3 <- site_patches_with_points %>%
  dplyr::mutate(
    patch_cat = ifelse(year(survey_date) == 2024, "predicted 2024", "predicted 2025")
  ) %>%
  dplyr::select(-site_type.x) %>%
  dplyr::select(
    patch_id, latitude, longitude, survey_date, site,
    site_type = site_type.y, pred_patch, everything()
  ) %>%
  dplyr::mutate(
    pred_patch = ifelse(is.na(pred_patch), site_type, pred_patch)
  ) %>%
  dplyr::filter(!is.na(pred_patch))

################################################################################
# Step 4: build patch transition table
################################################################################

patch_transition_lookup <- quad_build3 %>%
  st_drop_geometry() %>%
  dplyr::mutate(Year = year(survey_date)) %>%
  dplyr::filter(Year %in% c(2024, 2025)) %>%
  dplyr::group_by(patch_id, Year) %>%
  dplyr::summarise(
    pred_patch = get_mode(pred_patch),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = Year,
    values_from = pred_patch,
    names_prefix = "patch_"
  ) %>%
  dplyr::mutate(
    patch_changes = dplyr::case_when(
      patch_2024 == patch_2025 ~ "Stable",
      patch_2024 == "BAR"   & patch_2025 %in% c("INCIP", "FOR") ~ "Improved",
      patch_2024 == "INCIP" & patch_2025 == "FOR"               ~ "Improved",
      patch_2024 == "FOR"   & patch_2025 %in% c("INCIP", "BAR") ~ "Declined",
      patch_2024 == "INCIP" & patch_2025 == "BAR"               ~ "Declined",
      TRUE ~ NA_character_
    ),
    patch_changes = factor(
      patch_changes,
      levels = c("Declined", "Stable", "Improved")
    ),
    patch_2024_lab = factor(
      dplyr::recode(
        as.character(patch_2024),
        "BAR"   = "Barren",
        "INCIP" = "Incipient",
        "FOR"   = "Forest"
      ),
      levels = c("Barren", "Incipient", "Forest")
    )
  ) %>%
  dplyr::filter(!is.na(patch_changes))

################################################################################
# Step 5: calculate mean otter occupancy per patch per year
################################################################################

scan_build <- scan_orig %>%
  janitor::clean_names() %>%
  dplyr::mutate(
    scan_date = as.Date(date),
    year = year(scan_date),
    ind = as.numeric(ind)
  ) %>%
  dplyr::filter(year %in% c(2024, 2025)) %>%
  dplyr::filter(!is.na(lat), !is.na(long))

scan_sf <- st_as_sf(
  scan_build,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

scan_in_patches <- st_join(
  scan_sf,
  site_patches_single %>% dplyr::select(patch_id),
  join = st_intersects,
  left = FALSE
) %>%
  dplyr::mutate(year = year(scan_date))

scan_patch_summaries <- scan_in_patches %>%
  st_drop_geometry() %>%
  dplyr::group_by(scan_date, year, patch_id) %>%
  dplyr::summarise(
    patch_tot_otters = sum(ind, na.rm = TRUE),
    .groups = "drop"
  )

meanoccupancy <- scan_patch_summaries %>%
  dplyr::group_by(year, patch_id) %>%
  dplyr::summarise(
    avg_otters = mean(patch_tot_otters, na.rm = TRUE),
    .groups = "drop"
  )

################################################################################
# Step 6: build delta occupancy table
################################################################################

df_delta <- meanoccupancy %>%
  tidyr::pivot_wider(
    names_from = year,
    values_from = avg_otters,
    names_prefix = "y"
  ) %>%
  dplyr::left_join(
    patch_transition_lookup %>%
      dplyr::select(patch_id, patch_changes, patch_2024_lab),
    by = "patch_id"
  ) %>%
  dplyr::mutate(
    delta_otters = y2025 - y2024
  ) %>%
  dplyr::filter(
    !is.na(y2024),
    !is.na(y2025),
    !is.na(delta_otters),
    !is.na(patch_changes),
    !is.na(patch_2024_lab)
  )

################################################################################
# Step 7: plot delta occupancy
################################################################################

change_cols <- c(
  "Declined" = "indianred",
  "Stable"   = "grey60",
  "Improved" = "navyblue"
)

p_delta <- ggplot(
  df_delta,
  aes(x = patch_changes, y = delta_otters, fill = patch_changes)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "black"
  ) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    alpha = 0.75,
    color = "black",
    linewidth = 0.3
  ) +
 # geom_jitter(
#    width = 0.12,
#    size = 1.7,
#    alpha = 0.55,
#    color = "black"
#  ) +
  facet_wrap(~ patch_2024_lab, nrow = 1) +
  scale_fill_manual(
    values = change_cols,
    name = "Patch outcome"
  ) +
  labs(
    x = "Patch outcome",
    y = expression(Delta~"mean otter occupancy (2025 - 2024)")
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(angle = 20, hjust = 1),
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9),
    legend.position = "right"
  )

p_delta

################################################################################
# Step 8: export
################################################################################

ggsave(
  here::here("figures", "Fig5_occupancy_scans.png"),
  plot = p_delta,
  bg = "white",
  width = 8,
  height = 4.5,
  units = "in",
  dpi = 600
)
