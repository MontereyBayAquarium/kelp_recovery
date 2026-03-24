################################################################################
# FIGURE: Recovery signature of incipient forests
# Standardized contrasts (Hedges' g with bootstrap 95% CIs)
#
# Positive values  = higher in incipient forests
# Negative values  = lower in incipient forests
#
# Contrasts:
#   1) Incipient vs Barren
#   2) Incipient vs Forest
#
# Variables included:
#   - Sea otter foraging effort
#   - Purple urchin density
#   - Urchin biomass density
#   - Urchin gonad mass
#   - Urchin concealment ratio
#   - Reef relief
#   - Reef rugosity
#   - Total kelp recruit density
#
# Variable order:
#   Ranked by absolute effect size for Incipient vs Barren
#
# Output:
#   figures/Fig_recovery_signature_effectsizes_ranked.png
################################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

require(librarian)
shelf(
  tidyverse, janitor, lubridate, sf, here
)

set.seed(1985)

################################################################################
# 0. Load data -----------------------------------------------------------------
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
  dplyr::mutate(
    state_resp = factor(state_resp, levels = c("BAR", "INCIP", "FOR"))
  )

################################################################################
# 1. Helper functions ----------------------------------------------------------
################################################################################

first_present <- function(dat, choices) {
  out <- intersect(choices, names(dat))
  if (length(out) == 0) return(NA_character_)
  out[1]
}

safe_ratio <- function(num, den) {
  out <- ifelse(is.finite(num) & is.finite(den) & den > 0, num / den, NA_real_)
  pmin(pmax(out, 0), 1)
}

hedges_g <- function(x1, x2) {
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  
  n1 <- length(x1)
  n2 <- length(x2)
  
  if (n1 < 2 || n2 < 2) return(NA_real_)
  
  m1 <- mean(x1)
  m2 <- mean(x2)
  s1 <- stats::sd(x1)
  s2 <- stats::sd(x2)
  
  s_pooled <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  
  if (!is.finite(s_pooled) || s_pooled == 0) return(NA_real_)
  
  d <- (m1 - m2) / s_pooled
  
  # Small-sample correction
  J <- 1 - (3 / (4 * (n1 + n2) - 9))
  
  J * d
}

boot_hedges_g <- function(x1, x2, n_boot = 2000, conf = 0.95) {
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  
  n1 <- length(x1)
  n2 <- length(x2)
  
  if (n1 < 2 || n2 < 2) {
    return(tibble::tibble(
      effect    = NA_real_,
      conf.low  = NA_real_,
      conf.high = NA_real_,
      n1        = n1,
      n2        = n2
    ))
  }
  
  est <- hedges_g(x1, x2)
  
  boot_vals <- replicate(n_boot, {
    x1b <- sample(x1, size = n1, replace = TRUE)
    x2b <- sample(x2, size = n2, replace = TRUE)
    hedges_g(x1b, x2b)
  })
  
  alpha <- 1 - conf
  qs <- stats::quantile(
    boot_vals,
    probs = c(alpha / 2, 1 - alpha / 2),
    na.rm = TRUE
  )
  
  tibble::tibble(
    effect    = est,
    conf.low  = unname(qs[1]),
    conf.high = unname(qs[2]),
    n1        = n1,
    n2        = n2
  )
}

################################################################################
# 2. Build patch-level predictors ----------------------------------------------
################################################################################

core_patch_cols <- c(
  "mean_biomass_g",
  "mean_gonad_mass_g",
  "purple_urchin_densitym2",
  "purple_urchin_conceiledm2",
  "purple_urchin_concealedm2",
  "total_biomass_g",
  "total_gonad_mass_g",
  "relief_cm",
  "risk_index",
  "lamr",
  "macr",
  "nerj",
  "ptej",
  "lsetj",
  "eisj"
)

patch_predictors_all <- quad_build3 %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(year = lubridate::year(survey_date)) %>%
  dplyr::group_by(patch_id, site, zone, year) %>%
  dplyr::summarise(
    dplyr::across(dplyr::any_of(core_patch_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    urchin_biomass_densitym2 = total_biomass_g / patch_area_m2
  )

conceal_col <- if ("purple_urchin_conceiledm2" %in% names(patch_predictors_all)) {
  "purple_urchin_conceiledm2"
} else if ("purple_urchin_concealedm2" %in% names(patch_predictors_all)) {
  "purple_urchin_concealedm2"
} else {
  NA_character_
}

patch_predictors_all <- patch_predictors_all %>%
  dplyr::mutate(
    behavior_ratio = if (!is.na(conceal_col)) {
      safe_ratio(.data[[conceal_col]], purple_urchin_densitym2)
    } else {
      NA_real_
    },
    kelp_recruit_density = rowSums(
      dplyr::across(dplyr::any_of(c("lamr", "macr", "nerj", "ptej", "lsetj", "eisj"))),
      na.rm = TRUE
    )
  )

################################################################################
# 3. Build sea otter foraging predictor ----------------------------------------
################################################################################

forage_path <- "/Volumes/enhydra/data/foraging_data/processed/foraging_data_2024_2025_processed.csv"
forage_orig <- readr::read_csv(forage_path, show_col_types = FALSE)

lon_col  <- first_present(forage_orig, c("long", "lon", "longitude", "Longitude"))
lat_col  <- first_present(forage_orig, c("lat", "latitude", "Latitude"))
id_col   <- first_present(forage_orig, c("foragdata_id", "foragdiv_id", "fid", "id"))
date_col <- first_present(forage_orig, c("date", "Date"))

if (any(is.na(c(lon_col, lat_col, id_col)))) {
  stop("Missing lon/lat/id columns in foraging data.")
}

forage_dat <- forage_orig %>%
  dplyr::mutate(
    .date_tmp = if (!is.na(date_col)) {
      suppressWarnings(lubridate::as_date(.data[[date_col]]))
    } else {
      as.Date(NA)
    },
    date = dplyr::if_else(
      !is.na(.date_tmp),
      .date_tmp,
      lubridate::make_date(year = .data$year, month = .data$month, day = .data$day)
    )
  ) %>%
  dplyr::mutate(date = lubridate::as_date(date)) %>%
  dplyr::select(dplyr::all_of(c(
    id_col, lon_col, lat_col, "prey", "size", "qualifier",
    "number", "date", "month", "year"
  ))) %>%
  dplyr::rename(
    fid  = dplyr::all_of(id_col),
    long = dplyr::all_of(lon_col),
    lat  = dplyr::all_of(lat_col)
  )

size_key <- tidyr::expand_grid(size = 1:4, qualifier = c("a", "b", "c")) %>%
  dplyr::mutate(
    size_cm = (size - 1) * 5 + dplyr::case_when(
      qualifier == "a" ~ 1.66,
      qualifier == "b" ~ 3.32,
      TRUE             ~ 4.98
    )
  )

pur_forage <- forage_dat %>%
  dplyr::filter(month %in% focal_months, tolower(prey) == "pur") %>%
  dplyr::left_join(size_key, by = c("size", "qualifier")) %>%
  dplyr::mutate(
    test_diameter_mm = size_cm * 10,
    biomass_g        = -14.2 + 7.44 * exp(0.04 * test_diameter_mm),
    biomass_g        = ifelse(biomass_g < 0.5, 0.5, biomass_g),
    number           = tidyr::replace_na(number, 1),
    total_biomass_g  = biomass_g * number
  ) %>%
  dplyr::filter(size_cm < 8) %>%
  dplyr::filter(is.finite(long), is.finite(lat))

pur_sf <- sf::st_as_sf(
  pur_forage,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

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
  dplyr::mutate(
    rel_foraging_effort = foraging_biomass_kg / sum(foraging_biomass_kg, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    rel_foraging_effort = tidyr::replace_na(rel_foraging_effort, 0)
  )

################################################################################
# 4. Join predictors + state response ------------------------------------------
################################################################################

driver_df_all <- patch_predictors_all %>%
  dplyr::filter(year %in% years_keep) %>%
  dplyr::left_join(
    foraging_by_patch %>%
      dplyr::select(patch_id, site, zone, year, rel_foraging_effort),
    by = c("patch_id", "site", "zone", "year")
  ) %>%
  dplyr::left_join(
    state_lookup_all,
    by = c("patch_id", "site", "zone", "year")
  ) %>%
  dplyr::mutate(
    state_resp = factor(state_resp, levels = c("BAR", "INCIP", "FOR"))
  )

################################################################################
# 5. Select focal variables ----------------------------------------------------
################################################################################

pretty_lab_map <- c(
  "rel_foraging_effort"      = "Sea otter foraging effort",
  "purple_urchin_densitym2"  = "Purple urchin density",
  "urchin_biomass_densitym2" = "Urchin biomass density",
  "mean_gonad_mass_g"        = "Urchin gonad mass",
  "behavior_ratio"           = "Urchin concealment ratio",
  "relief_cm"                = "Reef relief (cm)",
  "risk_index"               = "Reef rugosity",
  "kelp_recruit_density"     = "Kelp recruit density"
)

vars_focus <- c(
  "rel_foraging_effort",
  "purple_urchin_densitym2",
  "urchin_biomass_densitym2",
  "mean_gonad_mass_g",
  "behavior_ratio",
  "relief_cm",
  "risk_index",
  "kelp_recruit_density"
)

plot_df <- driver_df_all %>%
  dplyr::filter(!is.na(state_resp)) %>%
  dplyr::select(state_resp, dplyr::all_of(vars_focus)) %>%
  tidyr::pivot_longer(
    cols      = -state_resp,
    names_to  = "variable",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    variable_lab = dplyr::recode(variable, !!!pretty_lab_map)
  )

################################################################################
# 6. Compute effect sizes ------------------------------------------------------
################################################################################

get_contrast_effects <- function(dat, focal = "INCIP", ref = "BAR", contrast_lab) {
  dat %>%
    dplyr::group_by(variable, variable_lab) %>%
    dplyr::group_modify(~{
      x1 <- .x$value[.x$state_resp == focal]
      x2 <- .x$value[.x$state_resp == ref]
      boot_hedges_g(x1, x2, n_boot = 2000, conf = 0.95)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      contrast = contrast_lab
    )
}

eff_bar <- get_contrast_effects(
  plot_df,
  focal = "INCIP",
  ref = "BAR",
  contrast_lab = "Incipient vs Barren"
)

eff_for <- get_contrast_effects(
  plot_df,
  focal = "INCIP",
  ref = "FOR",
  contrast_lab = "Incipient vs Forest"
)

effects_df <- dplyr::bind_rows(eff_bar, eff_for) %>%
  dplyr::mutate(
    contrast = factor(
      contrast,
      levels = c("Incipient vs Barren", "Incipient vs Forest")
    )
  )

################################################################################
# 7. Rank variables by Incipient vs Barren effect size -------------------------
################################################################################

var_order <- effects_df %>%
  dplyr::filter(contrast == "Incipient vs Barren") %>%
  dplyr::mutate(abs_effect = abs(effect)) %>%
  dplyr::arrange(dplyr::desc(abs_effect)) %>%
  dplyr::pull(variable_lab)

effects_df_plot <- effects_df %>%
  dplyr::mutate(
    variable_lab = factor(variable_lab, levels = rev(var_order))
  ) %>%
  dplyr::arrange(variable_lab, contrast) %>%
  dplyr::mutate(
    y_base = as.numeric(variable_lab),
    y = dplyr::case_when(
      contrast == "Incipient vs Barren" ~ y_base - 0.16,
      contrast == "Incipient vs Forest" ~ y_base + 0.16,
      TRUE ~ y_base
    )
  )

################################################################################
# 8. Build figure --------------------------------------------------------------
################################################################################

subtitle_text <- paste(
  "Variables are ranked by the absolute effect size for Incipient vs Barren.",
  "Points show Hedges' g and error bars show bootstrap 95% confidence intervals."
)

x_rng <- range(
  c(effects_df_plot$conf.low, effects_df_plot$conf.high),
  na.rm = TRUE
)
x_pad <- diff(x_rng) * 0.08
if (!is.finite(x_pad) || x_pad == 0) x_pad <- 0.25

p_effects <- ggplot(effects_df_plot, aes(x = effect, y = y, color = contrast)) +
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    linewidth  = 0.5,
    color      = "grey40"
  ) +
  geom_segment(
    aes(
      x    = conf.low,
      xend = conf.high,
      y    = y,
      yend = y
    ),
    linewidth = 1.0,
    alpha     = 0.95,
    lineend   = "round"
  ) +
  geom_point(
    size  = 3.8,
    alpha = 0.98
  ) +
  scale_color_manual(
    values = c(
      "Incipient vs Barren" = patch_colors[["BAR"]],
      "Incipient vs Forest" = patch_colors[["FOR"]]
    ),
    name = NULL
  ) +
  scale_y_continuous(
    breaks = seq_along(rev(var_order)),
    labels = rev(var_order),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  coord_cartesian(
    xlim = c(x_rng[1] - x_pad, x_rng[2] + x_pad),
    clip = "off"
  ) +
  labs(
    x        = "Standardized difference (Hedges' g)",
    y        = NULL,
    title    = "Mechanistic correlates of recovery",
    subtitle = subtitle_text,
    caption  = paste(
      "Positive values indicate variables that are higher in incipient forests.",
      "Confidence intervals are bootstrap 95% intervals and may be asymmetric."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position      = "top",
    legend.justification = "left",
    legend.text          = element_text(size = 10),
    
    axis.text.y          = element_text(size = 10, color = "black"),
    axis.text.x          = element_text(size = 10, color = "black"),
    axis.title.x         = element_text(size = 11),
    
    plot.title           = element_text(size = 14, face = "bold"),
    plot.subtitle        = element_text(size = 10, color = "grey20"),
    plot.caption         = element_text(size = 9, color = "grey30"),
    
    panel.grid.major.x   = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor.x   = element_blank(),
    panel.grid.major.y   = element_blank(),
    
    axis.line.y          = element_blank(),
    axis.ticks.y         = element_blank(),
    
    plot.margin          = ggplot2::margin(t = 10, r = 15, b = 10, l = 10)
  )

p_effects

################################################################################
# 9. Save figure ---------------------------------------------------------------
################################################################################

ggsave(
  filename = here::here("figures", "Fig_recovery_signature_effectsizes_ranked.png"),
  plot     = p_effects,
  width    = 8.6,
  height   = 5.2,
  dpi      = 600,
  bg       = "white"
)