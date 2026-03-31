rm(list = ls())

################################################################################
# libraries
################################################################################

require(librarian)
librarian::shelf(
  tidyverse, lubridate, sf, janitor, here,
  broom, patchwork, scales
)

################################################################################
# paths + data
################################################################################

datadir <- "/Volumes/enhydra/data/kelp_recovery/"
figdir  <- here::here("figures")

# load benthic survey data
load(file.path(
  datadir,
  "MBA_kelp_forest_database/processed/recovery/kelp_recovery_data.rda"
))

# load scans
scan_orig <- readr::read_csv(
  file.path(here::here("output", "scans", "scans_datav2.csv")),
  show_col_types = FALSE
)

# load site patches
site_patches <- sf::st_read(
  here::here("output", "gis_data", "processed", "site_patch_polygons.shp"),
  quiet = TRUE
)

# load LDA-predicted patch transitions
load(here::here("output", "lda_patch_transitionsv5.rda"))  # transitions_tbl_constrained

################################################################################
# helpers
################################################################################

get_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

clean_patch <- function(x) {
  factor(
    dplyr::recode(
      as.character(x),
      "BAR"   = "Barren",
      "INCIP" = "Incipient",
      "FOR"   = "Forest"
    ),
    levels = c("Barren", "Incipient", "Forest")
  )
}

scale_vec <- function(x) {
  mu  <- mean(x, na.rm = TRUE)
  sdv <- sd(x, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) sdv <- 1
  (x - mu) / sdv
}

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

make_2d_grid <- function(x, y, n_x = 120, n_y = 120) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  x_min <- max(0, min(x, na.rm = TRUE))
  x_max <- min(1, max(x, na.rm = TRUE))
  y_min <- min(y, na.rm = TRUE)
  y_max <- max(y, na.rm = TRUE)
  
  if (!is.finite(x_min) || !is.finite(x_max) ||
      !is.finite(y_min) || !is.finite(y_max)) {
    stop("Predictors contain no finite values.")
  }
  
  if (x_max == x_min) x_max <- x_min + 1e-6
  if (y_max == y_min) y_max <- y_min + 1e-6
  
  x_mu <- mean(x, na.rm = TRUE)
  x_sd <- sd(x, na.rm = TRUE)
  y_mu <- mean(y, na.rm = TRUE)
  y_sd <- sd(y, na.rm = TRUE)
  
  if (!is.finite(x_sd) || x_sd == 0) x_sd <- 1
  if (!is.finite(y_sd) || y_sd == 0) y_sd <- 1
  
  tidyr::crossing(
    occ_freq_2024 = seq(x_min, x_max, length.out = n_x),
    log_total_otters_2024 = seq(y_min, y_max, length.out = n_y)
  ) %>%
    dplyr::mutate(
      occ_freq_2024_sc = (occ_freq_2024 - x_mu) / x_sd,
      log_total_otters_2024_sc = (log_total_otters_2024 - y_mu) / y_sd
    )
}


################################################################################
# STEP 1: aggregate benthic survey data to zone x survey_date
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
# STEP 2: assign annual model-predicted patch types
################################################################################

quad_zone_with_pred <- quad_zone %>%
  dplyr::left_join(
    transitions_tbl_constrained %>%
      dplyr::select(site, site_type, zone, patch_2024, patch_2025),
    by = c("site", "site_type", "zone")
  ) %>%
  dplyr::mutate(
    pred_patch = dplyr::case_when(
      lubridate::year(survey_date) == 2024 ~ as.character(patch_2024),
      lubridate::year(survey_date) == 2025 ~ as.character(patch_2025),
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(-patch_2024, -patch_2025)

################################################################################
# STEP 3: join survey points to patch geometry
################################################################################

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
# STEP 4: build patch transition table
################################################################################

patch_transition_lookup <- quad_build3 %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(Year = lubridate::year(survey_date)) %>%
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
  dplyr::filter(!is.na(patch_2024), !is.na(patch_2025)) %>%
  dplyr::mutate(
    patch_2024_lab = clean_patch(patch_2024),
    patch_2025_lab = clean_patch(patch_2025),
    
    patch_changes = dplyr::case_when(
      patch_2024 == patch_2025 ~ "Stable",
      patch_2024 == "BAR"   & patch_2025 %in% c("INCIP", "FOR") ~ "Improved",
      patch_2024 == "INCIP" & patch_2025 == "FOR"               ~ "Improved",
      patch_2024 == "FOR"   & patch_2025 %in% c("INCIP", "BAR") ~ "Declined",
      patch_2024 == "INCIP" & patch_2025 == "BAR"               ~ "Declined",
      TRUE ~ NA_character_
    ),
    
    barren_outcome = dplyr::case_when(
      patch_2024 == "BAR" & patch_2025 == "BAR"               ~ "Stayed barren",
      patch_2024 == "BAR" & patch_2025 %in% c("INCIP", "FOR") ~ "Improved",
      TRUE ~ NA_character_
    ),
    barren_improved = dplyr::case_when(
      patch_2024 == "BAR" & patch_2025 == "BAR"               ~ 0,
      patch_2024 == "BAR" & patch_2025 %in% c("INCIP", "FOR") ~ 1,
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::mutate(
    patch_changes = factor(
      patch_changes,
      levels = c("Declined", "Stable", "Improved")
    ),
    barren_outcome = factor(
      barren_outcome,
      levels = c("Stayed barren", "Improved")
    )
  ) %>%
  dplyr::filter(!is.na(patch_changes))

################################################################################
# STEP 5: assign scan observations to patches
################################################################################

scan_build <- scan_orig %>%
  janitor::clean_names() %>%
  dplyr::mutate(
    scan_date = as.Date(date),
    year = lubridate::year(scan_date),
    ind = suppressWarnings(as.numeric(ind))
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
) %>%
  dplyr::mutate(year = lubridate::year(scan_date))

################################################################################
# STEP 6: patch-date summaries
################################################################################

scan_patch_summaries <- scan_in_patches %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(scan_date, year, patch_id) %>%
  dplyr::summarise(
    patch_tot_otters = sum(ind, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    occupied = as.integer(patch_tot_otters > 0)
  )

all_scan_dates <- tibble(scan_date = sort(unique(scan_build$scan_date)))
all_patch_ids  <- tibble(patch_id = sort(unique(site_patches_single$patch_id)))

scan_patch_date <- tidyr::crossing(all_scan_dates, all_patch_ids) %>%
  dplyr::left_join(
    scan_patch_summaries %>%
      dplyr::select(scan_date, year, patch_id, patch_tot_otters, occupied),
    by = c("scan_date", "patch_id")
  ) %>%
  dplyr::mutate(
    year = lubridate::year(scan_date),
    patch_tot_otters = tidyr::replace_na(patch_tot_otters, 0),
    occupied = tidyr::replace_na(occupied, 0L)
  )

################################################################################
# STEP 7: 2024 patch-level otter-use metrics
################################################################################

patch_otter_2024 <- scan_patch_date %>%
  dplyr::filter(year == 2024) %>%
  dplyr::group_by(patch_id) %>%
  dplyr::summarise(
    occ_freq_2024 = mean(occupied, na.rm = TRUE),
    total_otters_2024 = sum(patch_tot_otters, na.rm = TRUE),
    log_total_otters_2024 = log1p(total_otters_2024),
    mean_otters_2024 = mean(patch_tot_otters, na.rm = TRUE),
    n_scans_2024 = dplyr::n(),
    .groups = "drop"
  )

################################################################################
# STEP 8: lower panel delta table
################################################################################

meanoccupancy <- scan_patch_summaries %>%
  dplyr::group_by(year, patch_id) %>%
  dplyr::summarise(
    avg_otters = mean(patch_tot_otters, na.rm = TRUE),
    .groups = "drop"
  )

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
# STEP 9: barren-start modeling table
################################################################################

barren_dat <- patch_transition_lookup %>%
  dplyr::left_join(patch_otter_2024, by = "patch_id") %>%
  dplyr::filter(
    patch_2024 == "BAR",
    !is.na(barren_improved),
    !is.na(barren_outcome),
    !is.na(occ_freq_2024),
    !is.na(log_total_otters_2024)
  )

if (nrow(barren_dat) == 0 || length(unique(barren_dat$barren_improved)) < 2) {
  stop("Barren-start improvement model could not be estimated.")
}

################################################################################
# STEP 10: fit barren-start models
################################################################################

barren_fit_surface <- barren_dat %>%
  dplyr::mutate(
    occ_freq_2024_sc = scale_vec(occ_freq_2024),
    log_total_otters_2024_sc = scale_vec(log_total_otters_2024)
  )

mod_barren_surface <- stats::glm(
  barren_improved ~ occ_freq_2024_sc + log_total_otters_2024_sc,
  family = binomial(),
  data = barren_fit_surface
)

barren_fit_curve <- barren_fit_surface %>%
  dplyr::mutate(
    otter_use_index_raw = stats::predict(
      mod_barren_surface,
      newdata = barren_fit_surface,
      type = "link"
    )
  )

index_min <- min(barren_fit_curve$otter_use_index_raw, na.rm = TRUE)
index_max <- max(barren_fit_curve$otter_use_index_raw, na.rm = TRUE)

if (!is.finite(index_min) || !is.finite(index_max) || index_min == index_max) {
  stop("Otter-use index could not be scaled.")
}

barren_fit_curve <- barren_fit_curve %>%
  dplyr::mutate(
    otter_use_index = (otter_use_index_raw - index_min) / (index_max - index_min)
  )

mod_barren_curve <- stats::glm(
  barren_improved ~ otter_use_index,
  family = binomial(),
  data = barren_fit_curve
)

################################################################################
# STEP 11: predictions for panel A
################################################################################

pred_A <- make_2d_grid(
  x = barren_fit_surface$occ_freq_2024,
  y = barren_fit_surface$log_total_otters_2024,
  n_x = 120,
  n_y = 120
) %>%
  dplyr::mutate(
    pred_prob = stats::predict(
      mod_barren_surface,
      newdata = .,
      type = "response"
    )
  ) %>%
  dplyr::filter(
    is.finite(occ_freq_2024),
    is.finite(log_total_otters_2024),
    is.finite(pred_prob)
  )

rug_A <- barren_fit_surface %>%
  dplyr::transmute(
    occ_freq_2024,
    log_total_otters_2024
  )

################################################################################
# STEP 12: predictions for panel B with +/- 1 SE ribbons
################################################################################

pred_B_grid <- tibble(
  otter_use_index = seq(0, 1, length.out = 200)
)

pred_B_link <- stats::predict(
  mod_barren_curve,
  newdata = pred_B_grid,
  type = "link",
  se.fit = TRUE
)

pred_B_grid <- pred_B_grid %>%
  dplyr::mutate(
    fit_link = pred_B_link$fit,
    se_link = pred_B_link$se.fit,
    improved_fit = inv_logit(fit_link),
    improved_se_low = inv_logit(fit_link - se_link),
    improved_se_high = inv_logit(fit_link + se_link),
    stayed_fit = 1 - improved_fit,
    stayed_se_low = 1 - improved_se_high,
    stayed_se_high = 1 - improved_se_low
  )

pred_B <- dplyr::bind_rows(
  pred_B_grid %>%
    dplyr::transmute(
      otter_use_index,
      outcome = "Stayed barren",
      pred_prob = stayed_fit,
      se.low = stayed_se_low,
      se.high = stayed_se_high
    ),
  pred_B_grid %>%
    dplyr::transmute(
      otter_use_index,
      outcome = "Improved",
      pred_prob = improved_fit,
      se.low = improved_se_low,
      se.high = improved_se_high
    )
)

rug_B <- barren_fit_curve %>%
  dplyr::transmute(
    otter_use_index
  )

################################################################################
# STEP 13: panel A
################################################################################

base_theme <- theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    legend.title = element_text(),
    legend.text = element_text(color = "black"),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

change_cols <- c(
  "Declined" = "indianred",
  "Stable"   = "grey60",
  "Improved" = "navyblue"
)

curve_cols <- c(
  "Stayed barren" = "grey60",
  "Improved"      = "navyblue"
)


p_A <- ggplot(
  pred_A,
  aes(x = occ_freq_2024, y = log_total_otters_2024, fill = pred_prob)
) +
  geom_raster(interpolate = TRUE) +
  geom_contour(
    aes(z = pred_prob),
    color = "white",
    linewidth = 0.35,
    alpha = 0.8,
    bins = 6
  ) +
  geom_point(
    data = rug_A,
    aes(x = occ_freq_2024, y = log_total_otters_2024),
    inherit.aes = FALSE,
    color = "black",
    alpha = 0.45,
    size = 1.2
  ) +
  scale_fill_viridis_c(
    name = "Improvement\nprobability",
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  coord_cartesian(xlim = c(0, 0.8)) +
  labs(
    x = "2024 sea otter \noccupancy frequency",
    y = "log(1 + total otters counted \nin patch during 2024)"
  ) +
  base_theme +
  theme(
    legend.position = "right"
  )

################################################################################
# STEP 14: panel B
################################################################################

p_B <- ggplot(
  pred_B,
  aes(
    x = otter_use_index,
    y = pred_prob,
    color = outcome,
    fill = outcome,
    group = outcome
  )
) +
  geom_ribbon(
    aes(ymin = se.low, ymax = se.high),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 1.5, lineend = "round") +
  geom_rug(
    data = rug_B,
    aes(x = otter_use_index),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.24,
    linewidth = 0.35,
    color = "black"
  ) +
  scale_color_manual(values = curve_cols) +
  scale_fill_manual(values = curve_cols) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0.1, 0.3, 0.5, 0.7, 0.9),
    labels = c("10", "30", "50", "70", "90"),
    expand = expansion(mult = c(0.01, 0.02))
  )+
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "2024 occupancy index \n(percentile from model)",
    y = "Predicted probability of \n2025 outcome",
    color = "2025 outcome",
    fill = "2025 outcome"
  ) +
  base_theme +
  theme(
    legend.position = "right"
  )

################################################################################
# STEP 15: panel C
################################################################################

p_C <- ggplot(
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
  facet_wrap(~ patch_2024_lab, nrow = 1) +
  scale_fill_manual(
    values = change_cols,
    name = "Patch outcome"
  ) +
  scale_y_continuous(
    limits = c(-15, NA),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = "Patch outcome",
    y = expression(Delta~"mean otter occupancy (2025 - 2024)")
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 20, hjust = 1),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 11),
    legend.position = "right",
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

################################################################################
# STEP 16: combine figure
################################################################################

top_row <- p_A | p_B

final_fig <- top_row / p_C +
  patchwork::plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.tag = element_text(face = "bold", size = 14)
    )
  ) +
  patchwork::plot_layout(heights = c(0.92, 1.0), widths = c(1.05, 0.95))

final_fig

ggsave(
  filename = file.path(figdir, "Fig5_otter_occupancy_patch_dynamics.png"),
  plot = final_fig,
  width = 10,
  height = 7,
  units = "in",
  dpi = 600
)
