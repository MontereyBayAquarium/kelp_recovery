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
#   - Otter occupancy
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

scan_orig <- readr::read_csv(
  here::here("output", "scans", "scans_data.csv"),
  show_col_types = FALSE
)

years_keep    <- c(2024, 2025)
patch_area_m2 <- 80

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
# 2. Build patch-level benthic predictors --------------------------------------
################################################################################

core_patch_cols <- c(
  "mean_gonad_mass_g",
  "purple_urchin_densitym2",
  "purple_urchin_conceiledm2",
  "purple_urchin_concealedm2",
  "total_biomass_g",
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
# 3. Build patch polygons ------------------------------------------------------
################################################################################

patch_geom <- quad_build3 %>%
  dplyr::select(patch_id, site, zone, geometry) %>%
  dplyr::group_by(patch_id, site, zone) %>%
  dplyr::summarise(
    geometry = sf::st_union(geometry),
    .groups = "drop"
  ) %>%
  sf::st_as_sf()

################################################################################
# 4. Scan-based otter occupancy ------------------------------------------------
################################################################################

scan_build <- scan_orig %>%
  janitor::clean_names() %>%
  dplyr::mutate(
    scan_date = as.Date(date),
    year      = lubridate::year(scan_date),
    ind       = as.numeric(ind)
  ) %>%
  dplyr::filter(year %in% years_keep) %>%
  dplyr::filter(!is.na(lat), !is.na(long))

scan_sf <- sf::st_as_sf(
  scan_build,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

scan_in_patches <- sf::st_join(
  scan_sf,
  patch_geom %>% dplyr::select(patch_id, site, zone),
  join = sf::st_intersects,
  left = FALSE
)

scan_patch_date <- scan_in_patches %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(scan_date, year, patch_id, site, zone) %>%
  dplyr::summarise(
    patch_tot_otters = sum(ind, na.rm = TRUE),
    .groups = "drop"
  )

scan_patch_metrics <- scan_patch_date %>%
  dplyr::group_by(year, patch_id, site, zone) %>%
  dplyr::summarise(
    mean_otter_occupancy = mean(patch_tot_otters, na.rm = TRUE),
    .groups = "drop"
  )

################################################################################
# 5. Join all predictors -------------------------------------------------------
################################################################################

driver_df_all <- patch_predictors_all %>%
  dplyr::filter(year %in% years_keep) %>%
  dplyr::left_join(
    scan_patch_metrics,
    by = c("patch_id", "site", "zone", "year")
  ) %>%
  dplyr::left_join(
    state_lookup_all,
    by = c("patch_id", "site", "zone", "year")
  ) %>%
  dplyr::mutate(
    mean_otter_occupancy = tidyr::replace_na(mean_otter_occupancy, 0),
    state_resp = factor(state_resp, levels = c("BAR", "INCIP", "FOR"))
  )

################################################################################
# 6. Variables -----------------------------------------------------------------
################################################################################

pretty_lab_map <- c(
  "mean_otter_occupancy"     = "Sea otter occupancy",
  "purple_urchin_densitym2"  = "Purple urchin density",
  "urchin_biomass_densitym2" = "Urchin biomass density",
  "mean_gonad_mass_g"        = "Urchin gonad mass",
  "behavior_ratio"           = "Urchin concealment ratio",
  "relief_cm"                = "Reef relief",
  "risk_index"               = "Reef rugosity",
  "kelp_recruit_density"     = "Kelp recruit density"
)

vars_focus <- c(
  "mean_otter_occupancy",
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
# 7. Compute effect sizes ------------------------------------------------------
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
# 8. Rank variables ------------------------------------------------------------
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
    y = dplyr::if_else(
      contrast == "Incipient vs Barren",
      y_base - 0.16,
      y_base + 0.16
    )
  )

################################################################################
# 9. Plot ----------------------------------------------------------------------
################################################################################

p_effects <- ggplot(effects_df_plot, aes(x = effect, y = y, color = contrast)) +
  
  # zero reference
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    linewidth  = 0.6,
    color      = "grey40"
  ) +
  
  # CI
  geom_segment(
    aes(
      x    = conf.low,
      xend = conf.high,
      y    = y,
      yend = y
    ),
    linewidth = 1.2,
    alpha     = 0.95,
    lineend   = "round"
  ) +
  
  # points
  geom_point(
    size  = 4,
    alpha = 1
  ) +
  
  # colors
  scale_color_manual(
    values = c(
      "Incipient vs Barren" = patch_colors[["BAR"]],
      "Incipient vs Forest" = patch_colors[["FOR"]]
    ),
    name = NULL
  ) +
  
  # y-axis
  scale_y_continuous(
    breaks = seq_along(rev(var_order)),
    labels = rev(var_order),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Standardized difference (Hedges' g)",
    y = NULL,
    title = ""
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    # legend RIGHT
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 11),
    
    # TEXT
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 12),
    
    plot.title = element_text(size = 16, face = "bold"),
    plot.caption = element_text(size = 10, color = "grey30"),
    
    # GRID
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    
    # LEAN AXES (both x and y, same thickness)
    axis.line = element_line(linewidth = 0.8, color = "black"),
    
    # TICKS ON BOTH AXES
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    
    # Optional: control tick length (looks nice)
    axis.ticks.length = grid::unit(3, "pt"),
    
    plot.margin = ggplot2::margin(t = 10, r = 20, b = 10, l = 10)
  )

p_effects


################################################################################
# 10. Save figure --------------------------------------------------------------
################################################################################

ggsave(
  filename = here::here("figures", "Fig4_recovery_effectsizes_ranked.png"),
  plot     = p_effects,
  width    = 8.8,
  height   = 5.6,
  dpi      = 600,
  bg       = "white"
)


################################################################################
# 11. make table --------------------------------------------------------------
################################################################################

library(gt)

effect_table_pretty <- effects_df %>%
  dplyr::mutate(
    ci_excludes_zero = dplyr::if_else(
      conf.low > 0 | conf.high < 0,
      "Yes",
      "No"
    ),
    
    effect_direction = dplyr::case_when(
      effect > 0 ~ "Higher in incipient",
      effect < 0 ~ "Lower in incipient",
      TRUE ~ "No difference"
    ),
    
    estimate_ci = paste0(
      round(effect, 2),
      " [",
      round(conf.low, 2),
      ", ",
      round(conf.high, 2),
      "]"
    )
  ) %>%
  dplyr::arrange(
    factor(contrast, levels = c("Incipient vs Barren", "Incipient vs Forest")),
    dplyr::desc(abs(effect))
  ) %>%
  dplyr::select(
    contrast,
    variable = variable_lab,
    `Hedges' g [95% CI]` = estimate_ci,
    `Effect direction` = effect_direction,
    `95% CI excludes zero` = ci_excludes_zero
  )

gt_tbl <- effect_table_pretty %>%
  gt::gt(groupname_col = "contrast") %>%
 # gt::tab_header(
  #  title = "Mechanistic correlates of recovery",
   # subtitle = "Standardized effect sizes comparing incipient forests to barrens and forests"
  #) %>%
  gt::cols_label(
    variable = "Variable"
  ) %>%
  gt::cols_align(
    align = "left",
    columns = c(variable, `Effect direction`)
  ) %>%
  gt::cols_align(
    align = "center",
    columns = c(`Hedges' g [95% CI]`, `95% CI excludes zero`)
  ) %>%
  gt::tab_style(
    style = gt::cell_text(color = "red"),
    locations = gt::cells_body(
      rows = `95% CI excludes zero` == "Yes"
    )
  ) %>%
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_row_groups()
  ) %>%
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_column_labels(everything())
  ) %>%
  gt::tab_options(
    table.font.size = 12,
    data_row.padding = gt::px(5)
  )

gt_tbl

gt::gtsave(
  gt_tbl,
  filename = here::here("tables", "TabS1_effect_size_table.png"),
  vwidth = 1000,
  vheight = 600
)

###############################################################################
#Summary stats


behavior_summary <- driver_df_all %>%
  dplyr::filter(!is.na(state_resp)) %>%
  dplyr::group_by(state_resp) %>%
  dplyr::summarise(
    n      = sum(is.finite(behavior_ratio)),
    mean   = mean(behavior_ratio, na.rm = TRUE),
    sd     = sd(behavior_ratio, na.rm = TRUE),
    se     = sd / sqrt(n),
    median = median(behavior_ratio, na.rm = TRUE),
    q25    = quantile(behavior_ratio, 0.25, na.rm = TRUE),
    q75    = quantile(behavior_ratio, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

behavior_summary

# Incipient vs Barren
test_incip_bar <- wilcox.test(
  behavior_ratio ~ state_resp,
  data = driver_df_all %>%
    dplyr::filter(state_resp %in% c("INCIP", "BAR"))
)

# Incipient vs Forest
test_incip_for <- wilcox.test(
  behavior_ratio ~ state_resp,
  data = driver_df_all %>%
    dplyr::filter(state_resp %in% c("INCIP", "FOR"))
)

test_incip_bar
test_incip_for
