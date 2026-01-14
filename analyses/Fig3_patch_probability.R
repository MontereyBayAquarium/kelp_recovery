################################################################################
# FIGURE 3 (unchanged): Process-level drivers (PDPs)
#
# Depends on:
#   - output/patch_state_lookup_2024_2025.csv (written by Fig 2 script)
#
# Output:
#   - figures/Fig3_patch_process_PDPs.png
################################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

require(librarian)
shelf(
  tidyverse, janitor, lubridate, sf, here,
  randomForest, pdp, patchwork, purrr
)

set.seed(1985)

################################################################################
# 0. Load data + state lookup --------------------------------------------------
################################################################################

load(here::here("output", "survey_data", "processed", "zone_level_data4.rda"))

years_keep    <- c(2024, 2025)
patch_area_m2 <- 80
focal_months  <- c(6, 7, 8, 9)

patch_colors <- c(
  "BAR"   = "#7570B3",
  "INCIP" = "#D95F02",
  "FOR"   = "#1B9E77"
)

state_lookup_all <- readr::read_csv(
  here::here("output", "patch_state_lookup_2024_2025.csv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(state_resp = factor(state_resp, levels = c("BAR","INCIP","FOR")))

################################################################################
# Helpers ----------------------------------------------------------------------
################################################################################

med_impute <- function(x) {
  if (!is.numeric(x)) return(x)
  x[!is.finite(x)] <- NA_real_
  if (all(is.na(x))) return(x)
  x[is.na(x)] <- stats::median(x, na.rm = TRUE)
  x
}

first_present <- function(dat, choices) {
  out <- intersect(choices, names(dat))
  if (length(out) == 0) return(NA_character_)
  out[1]
}

safe_ratio <- function(num, den) {
  out <- ifelse(is.finite(num) & is.finite(den) & den > 0, num / den, 0)
  pmin(pmax(out, 0), 1)
}

################################################################################
# PART B. Build predictors (urchin + foraging) ---------------------------------
################################################################################

forage_path <- "/Volumes/enhydra/data/foraging_data/processed/foraging_data_2024_2025_processed.csv"
forage_orig <- readr::read_csv(forage_path, show_col_types = FALSE)

core_bio_cols <- c(
  "mean_biomass_g",
  "mean_gonad_mass_g",
  "purple_urchin_densitym2",
  "purple_urchin_conceiledm2",
  "purple_urchin_concealedm2",
  "total_biomass_g",
  "total_gonad_mass_g"
)

patch_predictors_all <- quad_build3 %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(year = lubridate::year(survey_date)) %>%
  dplyr::group_by(patch_id, site, zone, year) %>%
  dplyr::summarise(
    dplyr::across(dplyr::any_of(core_bio_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(urchin_biomass_densitym2 = total_biomass_g / patch_area_m2)

conceal_col <- if ("purple_urchin_conceiledm2" %in% names(patch_predictors_all)) {
  "purple_urchin_conceiledm2"
} else if ("purple_urchin_concealedm2" %in% names(patch_predictors_all)) {
  "purple_urchin_concealedm2"
} else {
  NA_character_
}

patch_predictors_all <- patch_predictors_all %>%
  dplyr::mutate(
    behavior_ratio = if (!is.na(conceal_col)) safe_ratio(.data[[conceal_col]], purple_urchin_densitym2) else NA_real_
  )

# Foraging join (PUR only; Jun–Sep)
lon_col  <- first_present(forage_orig, c("long","lon","longitude","Longitude"))
lat_col  <- first_present(forage_orig, c("lat","latitude","Latitude"))
id_col   <- first_present(forage_orig, c("foragdata_id","foragdiv_id","fid","id"))
date_col <- first_present(forage_orig, c("date","Date"))

if (any(is.na(c(lon_col, lat_col, id_col)))) stop("Missing lon/lat/id columns in foraging data.")

forage_dat <- forage_orig %>%
  dplyr::mutate(
    .date_tmp = if (!is.na(date_col)) suppressWarnings(lubridate::as_date(.data[[date_col]])) else as.Date(NA),
    date = dplyr::if_else(
      !is.na(.date_tmp),
      .date_tmp,
      lubridate::make_date(year = .data$year, month = .data$month, day = .data$day)
    )
  ) %>%
  dplyr::mutate(date = lubridate::as_date(date)) %>%
  dplyr::select(dplyr::all_of(c(id_col, lon_col, lat_col, "prey","size","qualifier","number","date","month","year"))) %>%
  dplyr::rename(fid = dplyr::all_of(id_col), long = dplyr::all_of(lon_col), lat = dplyr::all_of(lat_col))

size_key <- tidyr::expand_grid(size = 1:4, qualifier = c("a","b","c")) %>%
  dplyr::mutate(
    size_cm = (size - 1) * 5 + dplyr::case_when(
      qualifier == "a" ~ 1.66,
      qualifier == "b" ~ 3.32,
      TRUE             ~ 4.98
    )
  )

pur_forage <- forage_dat %>%
  dplyr::filter(month %in% focal_months, tolower(prey) == "pur") %>%
  dplyr::left_join(size_key, by = c("size","qualifier")) %>%
  dplyr::mutate(
    test_diameter_mm = size_cm * 10,
    biomass_g        = -14.2 + 7.44 * exp(0.04 * test_diameter_mm),
    biomass_g        = ifelse(biomass_g < 0.5, 0.5, biomass_g),
    number           = tidyr::replace_na(number, 1),
    total_biomass_g  = biomass_g * number
  ) %>%
  dplyr::filter(size_cm < 8) %>%
  dplyr::filter(is.finite(long), is.finite(lat))

pur_sf <- sf::st_as_sf(pur_forage, coords = c("long","lat"), crs = 4326, remove = FALSE)

quad_same_crs <- quad_build3 %>%
  sf::st_transform(sf::st_crs(pur_sf)) %>%
  dplyr::mutate(
    survey_year = lubridate::year(survey_date),
    survey_doy  = lubridate::yday(survey_date)
  ) %>%
  dplyr::select(patch_id, site, zone, survey_date, survey_year, survey_doy)

pur_spatial <- sf::st_join(pur_sf, quad_same_crs, join = sf::st_intersects, left = TRUE)

pur_patch <- pur_spatial %>%
  dplyr::mutate(
    forag_doy = lubridate::yday(date),
    diff_days = abs(forag_doy - survey_doy)
  ) %>%
  dplyr::group_by(fid) %>%
  dplyr::slice_min(diff_days, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  sf::st_drop_geometry() %>%
  dplyr::filter(survey_year %in% years_keep)

foraging_by_patch <- pur_patch %>%
  dplyr::group_by(patch_id, site, zone, survey_year) %>%
  dplyr::summarise(
    total_urchin_biomass_consumed_g = sum(total_biomass_g, na.rm = TRUE),
    n_foraging_obs                  = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    foraging_biomass_kg = total_urchin_biomass_consumed_g / 1000,
    year                = survey_year
  ) %>%
  dplyr::group_by(year) %>%
  dplyr::mutate(rel_foraging_effort = foraging_biomass_kg / sum(foraging_biomass_kg, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(rel_foraging_effort = tidyr::replace_na(rel_foraging_effort, 0))

################################################################################
# Join predictors + habitat-defined states -------------------------------------
################################################################################

driver_df_all <- patch_predictors_all %>%
  dplyr::filter(year %in% years_keep) %>%
  dplyr::left_join(
    foraging_by_patch %>% dplyr::select(patch_id, site, zone, year, rel_foraging_effort),
    by = c("patch_id","site","zone","year")
  ) %>%
  dplyr::left_join(
    state_lookup_all,
    by = c("patch_id","site","zone","year")
  ) %>%
  dplyr::mutate(state_resp = factor(state_resp, levels = c("BAR","INCIP","FOR")))

################################################################################
# Driver RF + PDPs (unchanged) -------------------------------------------------
################################################################################

rf_driver_input <- driver_df_all %>%
  dplyr::select(
    state_resp,
    behavior_ratio,
    mean_biomass_g,
    mean_gonad_mass_g,
    purple_urchin_densitym2,
    urchin_biomass_densitym2,
    rel_foraging_effort
  ) %>%
  dplyr::filter(!is.na(state_resp)) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), med_impute))

predictor_cols_B <- setdiff(names(rf_driver_input), "state_resp")
rf_driver_input$state_resp <- droplevels(rf_driver_input$state_resp)

set.seed(1985)
rf_driver <- randomForest::randomForest(
  x = rf_driver_input[, predictor_cols_B],
  y = rf_driver_input$state_resp,
  ntree      = 1500,
  mtry       = max(2, floor(sqrt(length(predictor_cols_B)))),
  importance = TRUE,
  na.action  = na.omit
)

pdp_plot_vars <- c(
  "mean_biomass_g",
  "urchin_biomass_densitym2",
  "behavior_ratio",
  "mean_gonad_mass_g",
  "rel_foraging_effort"
)

pretty_lab_map <- c(
  "mean_biomass_g"            = "Urchin biomass\n(g / urchin)",
  "urchin_biomass_densitym2"  = "Urchin biomass density\n(g / m²)",
  "behavior_ratio"            = "Urchin concealment ratio\n(concealed / total)",
  "mean_gonad_mass_g"         = "Urchin gonad mass\n(g / urchin)",
  "rel_foraging_effort"       = "Sea otter foraging effort\n(proportion of annual take)"
)

class_levels <- levels(rf_driver_input$state_resp)

pdp_all <- purrr::map_dfr(pdp_plot_vars, function(v) {
  purrr::map_dfr(class_levels, function(cls) {
    pd_tmp <- pdp::partial(
      rf_driver,
      pred.var        = v,
      which.class     = cls,
      prob            = TRUE,
      grid.resolution = 100,
      train           = rf_driver_input
    )
    tibble::tibble(
      x              = pd_tmp[[v]],
      y              = pd_tmp$yhat,
      state          = cls,
      predictor_var  = v,
      predictor_lab  = pretty_lab_map[[v]]
    )
  })
}) %>%
  tidyr::drop_na(x, y)

gini_order <- rf_driver$importance %>%
  as.data.frame() %>%
  tibble::rownames_to_column("variable") %>%
  dplyr::filter(variable %in% pdp_plot_vars) %>%
  dplyr::mutate(predictor_lab = pretty_lab_map[variable]) %>%
  dplyr::arrange(dplyr::desc(MeanDecreaseGini)) %>%
  dplyr::pull(predictor_lab)

pdp_all_ordered <- pdp_all %>%
  dplyr::mutate(
    predictor_lab = factor(predictor_lab, levels = gini_order),
    state         = factor(state, levels = c("BAR","INCIP","FOR"))
  )

note_labels <- rf_driver$importance %>%
  as.data.frame() %>%
  tibble::rownames_to_column("variable") %>%
  dplyr::filter(variable %in% pdp_plot_vars) %>%
  dplyr::mutate(
    predictor_lab = pretty_lab_map[variable],
    gini_note     = paste0("Gini = ", round(MeanDecreaseGini, 2)),
    predictor_lab = factor(predictor_lab, levels = gini_order)
  ) %>%
  dplyr::distinct(predictor_lab, gini_note)

p_final <- ggplot2::ggplot(
  pdp_all_ordered,
  ggplot2::aes(x = x, y = y, color = state, shape = state)
) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::facet_wrap(~ predictor_lab, scales = "free", nrow = 1, strip.position = "bottom") +
  ggplot2::geom_text(
    data = note_labels,
    ggplot2::aes(x = Inf, y = Inf, label = gini_note),
    hjust = 1.05, vjust = 1.5,
    size  = 2.8,
    color = "gray20",
    inherit.aes = FALSE
  ) +
  ggplot2::scale_color_manual(
    values = patch_colors,
    breaks = c("BAR", "INCIP", "FOR"),
    labels = c("Barren", "Incipient", "Forest"),
    name   = "Patch state"
  ) +
  ggplot2::scale_shape_manual(
    values = c(16, 17, 15),
    breaks = c("BAR", "INCIP", "FOR"),
    labels = c("Barren", "Incipient", "Forest"),
    name   = "Patch state"
  ) +
  ggplot2::labs(y = "Probability", x = NULL) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    strip.background   = ggplot2::element_blank(),
    strip.placement    = "outside",
    strip.text.x       = ggplot2::element_text(size = 8),
    axis.text          = ggplot2::element_text(size = 7),
    axis.title.y       = ggplot2::element_text(size = 9),
    legend.position    = "top",
    panel.spacing      = grid::unit(1, "lines"),
    axis.title.x       = ggplot2::element_blank()
  )

p_final

ggsave(
  filename = here::here("figures", "Fig3_patch_process_PDPs.png"),
  plot     = p_final,
  width    = 9,
  height   = 4,
  dpi      = 600,
  bg       = "white"
)
