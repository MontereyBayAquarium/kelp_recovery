# jogsmith@ucsc.edu

rm(list = ls())

require(librarian)
librarian::shelf(
  tidyverse, lubridate, sf, janitor, here
)

datadir <- "/Volumes/enhydra/data/kelp_recovery/"

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
# Step 0: load data
################################################################################

load(file.path(datadir, "MBA_kelp_forest_database/processed/recovery/kelp_recovery_data.rda"))

scan_orig <- readr::read_csv(
  file.path(here::here("output", "scans", "scans_datav2.csv")),
  show_col_types = FALSE
)

site_patches <- sf::st_read(
  here::here("output", "gis_data", "processed", "site_patch_polygons.shp"),
  quiet = TRUE
)

load(here::here("output", "lda_patch_transitionsv5.rda"))  # transitions_tbl_constrained

################################################################################
# Step 1: rebuild patch lookup from benthic data
################################################################################

kelp_avg <- kelp_data %>%
  dplyr::select(-macro_stipe_sd_20m2) %>%
  dplyr::group_by(site, site_type, latitude, longitude, zone, survey_date) %>%
  dplyr::summarise(
    dplyr::across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::select(-transect)

quad_avg <- quad_data %>%
  dplyr::group_by(site, site_type, latitude, longitude, zone, survey_date) %>%
  dplyr::summarise(
    dplyr::across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
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
    dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

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

quad_zone_sf <- sf::st_as_sf(
  quad_zone_with_pred,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

site_patches_single <- site_patches %>%
  sf::st_cast("POLYGON", warn = FALSE) %>%
  dplyr::mutate(patch_id = dplyr::row_number())

site_patches_with_points <- site_patches_single %>%
  sf::st_join(quad_zone_sf, join = sf::st_intersects, left = TRUE)

quad_build3 <- site_patches_with_points %>%
  dplyr::mutate(
    patch_cat = ifelse(lubridate::year(survey_date) == 2024, "predicted 2024", "predicted 2025")
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
# Step 2: patch outcome lookup
################################################################################

patch_transition_lookup <- quad_build3 %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(year = lubridate::year(survey_date)) %>%
  dplyr::filter(year %in% c(2024, 2025)) %>%
  dplyr::group_by(patch_id, year) %>%
  dplyr::summarise(
    pred_patch = get_mode(pred_patch),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = year,
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
  dplyr::filter(!is.na(patch_changes), !is.na(patch_2024_lab))

################################################################################
# Step 3: build full patch-level time series from scans
################################################################################

scan_build <- scan_orig %>%
  janitor::clean_names() %>%
  dplyr::mutate(
    scan_date = as.Date(date),
    year = lubridate::year(scan_date)
  ) %>%
  dplyr::filter(year %in% c(2024, 2025)) %>%
  dplyr::filter(!is.na(lat), !is.na(long))

scan_sf <- sf::st_as_sf(
  scan_build,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

scan_in_patches <- sf::st_join(
  scan_sf,
  site_patches_single %>% dplyr::select(patch_id),
  join = sf::st_intersects,
  left = FALSE
)

# Sum individuals per patch per scan date
scan_patch_daily <- scan_in_patches %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(scan_date, patch_id) %>%
  dplyr::summarise(
    patch_tot_otters = sum(ind, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    month = lubridate::floor_date(scan_date, "month"),
    month_lab = format(month, "%Y-%m")
  )

# Average within month to smooth the twice/month scan schedule
scan_patch_monthly <- scan_patch_daily %>%
  dplyr::group_by(patch_id, month, month_lab) %>%
  dplyr::summarise(
    monthly_occ = mean(patch_tot_otters, na.rm = TRUE),
    n_scans = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    patch_transition_lookup %>%
      dplyr::select(patch_id, patch_changes, patch_2024_lab),
    by = "patch_id"
  ) %>%
  dplyr::filter(!is.na(patch_changes), !is.na(patch_2024_lab))

################################################################################
# Step 4: plotting
################################################################################

change_cols <- c(
  "Declined" = "indianred",
  "Stable"   = "grey60",
  "Improved" = "navyblue"
)

p_ts <- ggplot(
  scan_patch_monthly,
  aes(x = month, y = monthly_occ, group = patch_id, color = patch_changes)
) +
  geom_line(
    alpha = 0.22,
    linewidth = 0.45
  ) +
  geom_smooth(
    aes(group = patch_changes),
    method = "loess",
    se = FALSE,
    linewidth = 1.3,
    span = 0.9
  ) +
  facet_wrap(~ patch_2024_lab, nrow = 1, scales = "free_y") +
  scale_color_manual(
    values = change_cols,
    name = "Patch outcome"
  ) +
  scale_x_date(
    date_breaks = "2 months",
    date_labels = "%b\n%Y"
  ) +
  labs(
    x = NULL,
    y = "Mean monthly otter occupancy per patch"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )

p_ts

################################################################################
# Step 5: export
################################################################################

ggsave(
  filename = "~/Downloads/otter_patch_timeseries_by_outcome.png",
  plot = p_ts,
  bg = "white",
  width = 9,
  height = 4.2,
  units = "in",
  dpi = 600
)



