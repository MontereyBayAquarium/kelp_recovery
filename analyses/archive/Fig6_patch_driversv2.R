################################################################################
# Patch State Classification + Ecological Drivers + Habitat PCA
#
# Author: Joshua G. Smith — UCSC Nearshore Ecology Research Group
#
# PART A. Habitat-only RF
#   Train on 2024 diver-called patch states using ONLY structure/cover/physical
#   habitat predictors (no urchin / gonad / foraging). Predict 2025 states.
#
# PART B. Ecological driver RF (2025 only)
#   Response  = predicted 2025 state from PART A
#   Predictors = urchin concealment ratio, urchin biomass (g / urchin),
#                urchin gonad mass (g / urchin), urchin density (ind / m²),
#                urchin biomass density (g / m² standing stock),
#                sea otter foraging effort (proportion of annual urchin take).
#
#   Output: variable importance + PDPs faceted by predictor.
#
# PART C. PCA biplot of habitat structure
#   Multivariate habitat space, colored by patch state.
################################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(librarian)
  librarian::shelf(
    tidyverse, janitor, lubridate, sf, here,
    randomForest, vip, pdp, patchwork, ggrepel, purrr
  )
})

set.seed(42)

################################################################################
# 0. Helper functions / constants ---------------------------------------------
################################################################################

clean_state_label <- function(x) {
  x_up <- toupper(x)
  dplyr::case_when(
    x_up %in% c(
      "BAR","BARREN","BARRENS","URCHIN BARREN",
      "BARREN PATCH","BARREN/URCHIN BARREN"
    ) ~ "BAR",
    x_up %in% c(
      "INCIP","INCIPIENT","INCIPIENT FOREST",
      "TRANSITION","TRANSITIONAL","INCIPIENT PATCH"
    ) ~ "INCIP",
    x_up %in% c(
      "FOR","FOREST","FORESTED","KELP FOREST"
    ) ~ "FOR",
    TRUE ~ NA_character_
  )
}

mode_char <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

asin_sqrt <- function(x) {
  # arcsin(sqrt(p)) for percent cover-style vars (0-100)
  p <- pmin(pmax(x, 0), 100) / 100
  asin(sqrt(p))
}

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
  out <- ifelse(
    is.finite(num) & is.finite(den) & den > 0,
    num / den,
    0
  )
  pmin(pmax(out, 0), 1)
}

# plotting colors for patch states
patch_colors <- c(
  "BAR"   = "purple",
  "INCIP" = "orange",
  "FOR"   = "forestgreen"
)

# analysis constants
focal_months  <- c(6, 7, 8, 9)  # Jun–Sep for foraging data
years_keep    <- c(2024, 2025)
patch_area_m2 <- 80            # reef area mapped per patch, m²

################################################################################
# 1. Load data -----------------------------------------------------------------
################################################################################

# Must contain quad_build3
load(here::here("output", "survey_data", "processed", "zone_level_data4.rda"))
if (!exists("quad_build3")) stop("quad_build3 not found in loaded .rda")

# Sea otter foraging data
forage_path <- "/Volumes/enhydra/data/foraging_data/processed/foraging_data_2024_2025_processed.csv"
if (!file.exists(forage_path)) stop("Foraging file not found at forage_path")
forage_orig <- readr::read_csv(forage_path, show_col_types = FALSE)

################################################################################
# PART A. Habitat-only Random Forest (train 2024; predict 2025) ----------------
################################################################################

dat_raw <- quad_build3 %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(year = lubridate::year(survey_date))

# Diver-truth patch state per patch_id-year
truth_all <- dat_raw %>%
  dplyr::mutate(site_type_clean = clean_state_label(site_type)) %>%
  dplyr::filter(year %in% years_keep) %>%
  dplyr::group_by(patch_id, site, zone, year) %>%
  dplyr::summarise(
    state = mode_char(site_type_clean),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    state = factor(state, levels = c("BAR","INCIP","FOR"))
  )

# Drop biological terms (urchin/gonad/etc.) from predictors: habitat only
ban_regex <- "(urchin|gonad|biomass|forag|foraging|\\bgi\\b|_gi$|mean_gi|sd_gi)"

num_cols_all    <- names(dat_raw)[sapply(dat_raw, is.numeric)]
allowed_numeric <- num_cols_all[
  !grepl(ban_regex, num_cols_all, ignore.case = TRUE)
]
allowed_numeric <- setdiff(allowed_numeric, "year")

# Aggregate numeric habitat/cover structure to patch_id-year
agg_hab <- dat_raw %>%
  dplyr::group_by(patch_id, site, zone, year) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::any_of(allowed_numeric),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

meta_hab <- agg_hab %>%
  dplyr::select(patch_id, site, zone, year)

X_raw <- agg_hab %>%
  dplyr::select(-patch_id, -site, -zone, -year)

# Transform cover % with asin(sqrt(p)), others with log1p
cover_cols <- names(X_raw)[grepl("^cov_", names(X_raw))]
other_cols <- setdiff(names(X_raw), cover_cols)

X_tr <- X_raw %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::any_of(cover_cols),
      ~ asin_sqrt(.x)
    ),
    dplyr::across(
      dplyr::any_of(other_cols),
      ~ log1p(pmax(.x, 0))
    )
  )

# Drop predictors with >50% NA
too_na <- sapply(X_tr, function(x) mean(!is.finite(x) | is.na(x)))
X_tr2  <- X_tr %>%
  dplyr::select(dplyr::any_of(names(too_na)[too_na <= 0.5]))

# Median-impute remaining NAs
X_imp <- X_tr2 %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      med_impute
    )
  )

# Drop zero-variance cols
zv_cols <- names(X_imp)[sapply(X_imp, function(x) sd(x, na.rm = TRUE) == 0)]
if (length(zv_cols)) {
  message("[PART A] Dropping zero-variance predictors: ", paste(zv_cols, collapse = ", "))
  X_imp <- X_imp %>%
    dplyr::select(-dplyr::any_of(zv_cols))
}

# Scale predictors
X_scaled_mat <- scale(X_imp)
X_scaled     <- tibble::as_tibble(X_scaled_mat, .name_repair = "minimal")
colnames(X_scaled) <- colnames(X_imp)

hab_scaled_df <- dplyr::bind_cols(meta_hab, X_scaled)

# Merge with diver-called truth
model_df <- hab_scaled_df %>%
  dplyr::left_join(
    truth_all,
    by = c("patch_id","site","zone","year")
  ) %>%
  dplyr::mutate(state = droplevels(state))

# Split train/predict
train_df <- model_df %>%
  dplyr::filter(year == 2024, !is.na(state))
test_df  <- model_df %>%
  dplyr::filter(year == 2025)

cat("\n[PART A] Training rows (2024): ", nrow(train_df), "\n")
cat("[PART A] Prediction rows (2025): ", nrow(test_df), "\n")
cat("[PART A] Training class balance (2024 diver calls):\n")
print(table(train_df$state, useNA = "ifany"))

predictor_cols_A <- setdiff(
  colnames(train_df),
  c("patch_id","site","zone","year","state")
)

# Train RF on 2024 habitat-only predictors
rf_train2024 <- randomForest(
  x = train_df[, predictor_cols_A],
  y = train_df$state,
  ntree      = 1500,
  mtry       = max(2, floor(sqrt(length(predictor_cols_A)))),
  importance = TRUE,
  na.action  = na.omit
)

# OOB diagnostics
oob_acc <- 1 - rf_train2024$err.rate[nrow(rf_train2024$err.rate), "OOB"]
cat(sprintf("\n[PART A] RF (2024 habitat classifier)\nOOB accuracy = %.3f\n", oob_acc))
cat("[PART A] OOB confusion matrix:\n")
print(rf_train2024$confusion)

vip_train2024 <- vip::vip(
  rf_train2024,
  num_features = min(20, length(predictor_cols_A)),
  bar = TRUE
) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(
    title = "Habitat / structural predictors distinguishing BAR / INCIP / FOR (2024)",
    x     = "Predictor (scaled)",
    y     = "Importance"
  )

# Predict 2025 patch states from habitat-only model
test_pred_class <- predict(
  rf_train2024,
  newdata = test_df[, predictor_cols_A],
  type    = "response"
)

pred_2025 <- test_df %>%
  dplyr::select(patch_id, site, zone, year) %>%
  dplyr::mutate(
    predicted_state_2025 = factor(
      test_pred_class,
      levels = c("BAR","INCIP","FOR")
    )
  ) %>%
  dplyr::distinct(patch_id, site, zone, year, .keep_all = TRUE)

cat("\n[PART A] Predicted 2025 state counts (2024-trained habitat model):\n")
print(table(pred_2025$predicted_state_2025, useNA = "ifany"))

# If we actually have 2025 diver calls, evaluate generalization
if (any(!is.na(test_df$state))) {
  eval_2025 <- test_df %>%
    dplyr::select(patch_id, site, zone, year, state) %>%
    dplyr::left_join(
      pred_2025,
      by = c("patch_id","site","zone","year")
    ) %>%
    dplyr::mutate(state = droplevels(state))
  
  conf_2025 <- table(
    True = eval_2025$state,
    Pred = eval_2025$predicted_state_2025
  )
  cat("\n[PART A] 2025 generalization confusion matrix:\n")
  print(conf_2025)
  
  acc_2025 <- mean(
    eval_2025$state == eval_2025$predicted_state_2025,
    na.rm = TRUE
  )
  cat(sprintf("[PART A] 2025 generalization accuracy = %.3f\n", acc_2025))
} else {
  cat("\n[PART A] No 2025 diver calls to score generalization.\n")
}

################################################################################
# PART B. Ecological drivers of 2025 predicted states -------------------------
################################################################################
# Goal: link predicted_state_2025 to prey/urchin/foraging metrics

# 1. Summarize patch-level urchin metrics from quad_build3
core_bio_cols <- c(
  "mean_biomass_g",               # urchin biomass per urchin (g/urchin)
  "mean_gonad_mass_g",            # gonad mass per urchin (g/urchin)
  "purple_urchin_densitym2",      # density (ind / m²)
  "purple_urchin_conceiledm2",    # concealed density? (typo version)
  "purple_urchin_concealedm2",    # concealed density? (alt spelling)
  "total_gonad_mass_g"            # total gonad mass (for gonad per-area calc if needed)
)

patch_predictors_all <- quad_build3 %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(year = lubridate::year(survey_date)) %>%
  dplyr::group_by(patch_id, site, zone, year) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::any_of(core_bio_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    # standing urchin biomass per m² on the reef
    urchin_biomass_densitym2 = mean_biomass_g * purple_urchin_densitym2
  )

# figure out which concealed-density column exists
conceal_col <- if ("purple_urchin_conceiledm2" %in% names(patch_predictors_all)) {
  "purple_urchin_conceiledm2"
} else if ("purple_urchin_concealedm2" %in% names(patch_predictors_all)) {
  "purple_urchin_concealedm2"
} else {
  NA_character_
}

################################################################################
# 2. Sea otter foraging data: biomass of PUR urchins consumed per patch
################################################################################

# robust column ID for lat/long/etc.
lon_col  <- first_present(forage_orig, c("long","lon","longitude","Longitude"))
lat_col  <- first_present(forage_orig, c("lat","latitude","Latitude"))
id_col   <- first_present(forage_orig, c("foragdata_id","foragdiv_id","fid","id"))
date_col <- first_present(forage_orig, c("date","Date"))

if (any(is.na(c(lon_col, lat_col, id_col)))) {
  stop("Missing lon/lat/id columns in foraging data.")
}

# build "forage_dat" with consistent names
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
      lubridate::make_date(
        year  = .data$year,
        month = .data$month,
        day   = .data$day
      )
    )
  ) %>%
  dplyr::mutate(date = lubridate::as_date(date)) %>%
  dplyr::select(
    dplyr::all_of(c(
      id_col, lon_col, lat_col,
      "prey","size","qualifier","number",
      "date","month","year"
    ))
  ) %>%
  dplyr::rename(
    fid  = dplyr::all_of(id_col),
    long = dplyr::all_of(lon_col),
    lat  = dplyr::all_of(lat_col)
  )

# Convert purple urchin prey items into biomass per observation
size_key <- tidyr::expand_grid(size = 1:4, qualifier = c("a","b","c")) %>%
  dplyr::mutate(
    size_cm = (size - 1) * 5 + dplyr::case_when(
      qualifier == "a" ~ 1.66,
      qualifier == "b" ~ 3.32,
      TRUE             ~ 4.98
    )
  )

pur_forage <- forage_dat %>%
  dplyr::filter(
    month %in% focal_months,
    tolower(prey) == "pur"
  ) %>%
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


# Spatially join foraging obs -> patch polygons, match nearest-in-time survey
pur_sf <- sf::st_as_sf(
  pur_forage,
  coords = c("long","lat"),
  crs    = 4326,
  remove = FALSE
)

quad_same_crs <- quad_build3 %>%
  sf::st_transform(sf::st_crs(pur_sf)) %>%
  dplyr::mutate(
    survey_year = lubridate::year(survey_date),
    survey_doy  = lubridate::yday(survey_date)
  ) %>%
  dplyr::select(
    patch_id, site, zone,
    survey_date, survey_year, survey_doy
  )

pur_spatial <- sf::st_join(
  pur_sf,
  quad_same_crs,
  join = sf::st_intersects,
  left = TRUE
)

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

# Summarize foraging -> per (patch_id,year)
foraging_by_patch <- pur_patch %>%
  dplyr::group_by(patch_id, site, zone, survey_year) %>%
  dplyr::summarise(
    total_urchin_biomass_consumed_g = sum(total_biomass_g, na.rm = TRUE),
    n_foraging_obs                  = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    foraging_biomass_kg   = total_urchin_biomass_consumed_g / 1000,
    year                  = survey_year
  ) %>%
  dplyr::group_by(year) %>%
  dplyr::mutate(
    # proportion of all urchin biomass (kg) taken in that patch that year
    rel_foraging_effort = foraging_biomass_kg /
      sum(foraging_biomass_kg, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    rel_foraging_effort = tidyr::replace_na(rel_foraging_effort, 0)
  )

################################################################################
# 3. Join predictors + predicted 2025 states
################################################################################

driver_df_2025 <- patch_predictors_all %>%
  dplyr::filter(year == 2025) %>%
  dplyr::mutate(
    behavior_ratio = if (!is.na(conceal_col)) {
      safe_ratio(.data[[conceal_col]], purple_urchin_densitym2)
    } else {
      NA_real_
    }
  ) %>%
  dplyr::left_join(
    foraging_by_patch %>%
      dplyr::filter(year == 2025) %>%
      dplyr::select(
        patch_id, site, zone, year,
        rel_foraging_effort
      ),
    by = c("patch_id","site","zone","year")
  ) %>%
  dplyr::left_join(
    pred_2025 %>%
      dplyr::select(
        patch_id, site, zone, year,
        predicted_state_2025
      ),
    by = c("patch_id","site","zone","year")
  ) %>%
  dplyr::mutate(
    state_resp = factor(predicted_state_2025, levels = c("BAR","INCIP","FOR"))
  )

cat("\n[PART B] Rows in driver_df_2025 (2025 only): ", nrow(driver_df_2025), "\n")
cat("[PART B] Class balance in predicted_state_2025:\n")
print(table(driver_df_2025$state_resp, useNA = "ifany"))

################################################################################
# 4. Prep RF input for ecological driver model --------------------------------
################################################################################

rf_driver_input <- driver_df_2025 %>%
  dplyr::select(
    state_resp,
    behavior_ratio,               # concealed / total
    mean_biomass_g,               # biomass per urchin (g / urchin)
    mean_gonad_mass_g,            # gonad mass per urchin (g / urchin)
    purple_urchin_densitym2,      # urchin density (ind / m²)
    urchin_biomass_densitym2,     # standing biomass density (g / m²)
    rel_foraging_effort           # otter foraging effort (proportion of take)
  ) %>%
  dplyr::filter(!is.na(state_resp))

# Median-impute numeric predictors
rf_driver_input <- rf_driver_input %>%
  dplyr::mutate(
    dplyr::across(
      .cols = where(is.numeric),
      .fns  = med_impute
    )
  )

# Drop any all-NA or zero-variance predictors
drop_cols_allNA <- names(rf_driver_input)[names(rf_driver_input) != "state_resp"][
  sapply(
    rf_driver_input[names(rf_driver_input) != "state_resp"],
    function(x) all(!is.finite(x) | is.na(x))
  )
]
if (length(drop_cols_allNA)) {
  message("[PART B] Dropping all-NA predictors: ", paste(drop_cols_allNA, collapse = ", "))
  rf_driver_input <- rf_driver_input %>%
    dplyr::select(-dplyr::any_of(drop_cols_allNA))
}

drop_cols_zv <- names(rf_driver_input)[names(rf_driver_input) != "state_resp"][
  sapply(
    rf_driver_input[names(rf_driver_input) != "state_resp"],
    function(x) is.numeric(x) && sd(x, na.rm = TRUE) == 0
  )
]
if (length(drop_cols_zv)) {
  message("[PART B] Dropping zero-variance predictors: ", paste(drop_cols_zv, collapse = ", "))
  rf_driver_input <- rf_driver_input %>%
    dplyr::select(-dplyr::any_of(drop_cols_zv))
}

predictor_cols_B <- setdiff(names(rf_driver_input), "state_resp")
rf_driver_input$state_resp <- droplevels(rf_driver_input$state_resp)

# Fit ecological driver RF
rf_driver <- randomForest(
  x = rf_driver_input[, predictor_cols_B],
  y = rf_driver_input$state_resp,
  ntree      = 1500,
  mtry       = max(2, floor(sqrt(length(predictor_cols_B)))),
  importance = TRUE,
  na.action  = na.omit
)

cat("\n[PART B] Driver RF confusion matrix (trained on predicted 2025 states):\n")
print(rf_driver$confusion)

vip_driver <- vip::vip(
  rf_driver,
  num_features = min(12, length(predictor_cols_B)),
  bar = TRUE
) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(
    title = "Process-level drivers of predicted 2025 patch states",
    x     = "Predictor",
    y     = "Importance"
  )

################################################################################
# 5. Partial Dependence Curves (PDPs) ------------------------------------------
################################################################################

# Variables we'll visualize
pdp_vars <- predictor_cols_B

# Pretty labels for facets (match final RF predictors)
pretty_lab_map <- c(
  "behavior_ratio"            = "Urchin concealment ratio\n(concealed / total)",
  "mean_biomass_g"            = "Urchin biomass\n(g / urchin)",
  "mean_gonad_mass_g"         = "Urchin gonad mass\n(g / urchin)",
  "purple_urchin_densitym2"   = "Urchin density\n(ind / m²)",
  "urchin_biomass_densitym2"  = "Urchin biomass density\n(g / m²)",
  "rel_foraging_effort"       = "Sea otter foraging effort\n(proportion of annual take)"
)

# Gather PDPs for each predictor × state
class_levels <- levels(rf_driver_input$state_resp)

pdp_all <- purrr::map_dfr(pdp_vars, function(v) {
  purrr::map_dfr(class_levels, function(cls) {
    pd_tmp <- pdp::partial(
      rf_driver,
      pred.var       = v,
      which.class    = cls,
      prob           = TRUE,
      grid.resolution = 100,
      train          = rf_driver_input
    )
    tibble::tibble(
      x              = pd_tmp[[v]],
      y              = pd_tmp$yhat,
      state          = cls,                   # "BAR","INCIP","FOR"
      predictor_var  = v,
      predictor_lab  = pretty_lab_map[[v]]
    )
  })
}) %>%
  tidyr::drop_na(x, y)

# Factor order for facets = order of pdp_vars
predictor_lab_order <- unname(pretty_lab_map[pdp_vars])

pdp_all <- pdp_all %>%
  dplyr::mutate(
    predictor_lab = factor(predictor_lab, levels = predictor_lab_order),
    state         = factor(state, levels = c("BAR","INCIP","FOR"))
  )

# Variable importance (MeanDecreaseGini) to annotate each facet
var_imp <- rf_driver$importance %>%
  as.data.frame() %>%
  tibble::rownames_to_column("variable") %>%
  dplyr::select(variable, MeanDecreaseGini) %>%
  dplyr::mutate(
    predictor_lab = pretty_lab_map[variable],
    gini_note     = paste0("Gini = ", round(MeanDecreaseGini, 2))
  ) %>%
  dplyr::filter(!is.na(predictor_lab)) %>%
  dplyr::mutate(
    predictor_lab = factor(predictor_lab, levels = predictor_lab_order)
  ) %>%
  dplyr::distinct(predictor_lab, gini_note)

note_labels <- var_imp

# Final PDP plot across predictors, colored by state
p_final <- ggplot2::ggplot(
  pdp_all,
  ggplot2::aes(x = x, y = y, color = state)
) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::facet_wrap(
    ~ predictor_lab,
    scales         = "free",
    nrow           = 1,
    strip.position = "bottom"
  ) +
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
    breaks = c("BAR","INCIP","FOR"),
    labels = c("Barren","Incipient","Forest"),
    name   = "Patch state"
  ) +
  ggplot2::labs(
    y = "Probability",
    x = NULL
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    strip.background   = ggplot2::element_blank(),
    strip.placement    = "outside",
    strip.text.x       = ggplot2::element_text(size = 9),
    axis.text          = ggplot2::element_text(size = 8),
    axis.title.y       = ggplot2::element_text(size = 10),
    legend.position    = "top",
    panel.spacing      = grid::unit(1, "lines"),
    axis.title.x       = ggplot2::element_blank()
  )

p_final
################################################################################
# PART C. PCA biplot of habitat structure -------------------------------------
################################################################################

# We'll color each patch by:
#   - 2024 diver-called state (for 2024 points)
#   - predicted 2025 state from PART A (for 2025 points)

state_lookup_2024 <- truth_all %>%
  dplyr::rename(state_2024call = state)

state_lookup_2025 <- pred_2025 %>%
  dplyr::mutate(
    predicted_state_2025 = factor(
      predicted_state_2025,
      levels = c("BAR","INCIP","FOR")
    )
  )

habitat_df <- hab_scaled_df %>%
  dplyr::left_join(
    state_lookup_2024,
    by = c("patch_id","site","zone","year")
  ) %>%
  dplyr::left_join(
    state_lookup_2025 %>%
      dplyr::select(patch_id, site, zone, year, predicted_state_2025),
    by = c("patch_id","site","zone","year")
  ) %>%
  dplyr::mutate(
    state_final = dplyr::case_when(
      year == 2024 ~ as.character(state_2024call),
      year == 2025 ~ as.character(predicted_state_2025),
      TRUE         ~ NA_character_
    ),
    state_final = factor(state_final, levels = c("BAR","INCIP","FOR"))
  ) %>%
  dplyr::filter(!is.na(state_final))

# Run PCA on scaled habitat-only predictors
predictor_cols_pca <- setdiff(
  colnames(habitat_df),
  c("patch_id","site","zone","year",
    "state_2024call","predicted_state_2025","state_final")
)

pca_mat <- habitat_df %>%
  dplyr::select(dplyr::any_of(predictor_cols_pca)) %>%
  as.matrix()

pca_obj <- prcomp(pca_mat, center = FALSE, scale. = FALSE)

scores_df <- tibble::as_tibble(pca_obj$x[, 1:2, drop = FALSE]) %>%
  dplyr::rename(PC1 = PC1, PC2 = PC2) %>%
  dplyr::bind_cols(
    habitat_df %>%
      dplyr::select(patch_id, site, zone, year, state_final)
  )

loadings_df <- tibble::as_tibble(
  pca_obj$rotation[, 1:2, drop = FALSE],
  rownames = "variable"
) %>%
  dplyr::rename(PC1 = PC1, PC2 = PC2) %>%
  dplyr::mutate(vec_len = sqrt(PC1^2 + PC2^2))

# pick top habitat/cover gradients (exclude urchin/gonad/etc.)
arrow_ban_regex <- "(urchin|gonad|biomass|forag|foraging|\\bgi\\b|_gi$|mean_gi|sd_gi)"
arrow_df <- loadings_df %>%
  dplyr::filter(!grepl(arrow_ban_regex, variable, ignore.case = TRUE)) %>%
  dplyr::slice_max(order_by = vec_len, n = 8, with_ties = FALSE)

range_x   <- range(scores_df$PC1, na.rm = TRUE)
range_y   <- range(scores_df$PC2, na.rm = TRUE)
max_span  <- max(diff(range_x), diff(range_y))
max_vec   <- max(arrow_df$vec_len)
arrow_scaler <- 0.4 * max_span / max_vec

arrow_df <- arrow_df %>%
  dplyr::mutate(
    x0 = 0,
    y0 = 0,
    x1 = PC1 * arrow_scaler,
    y1 = PC2 * arrow_scaler
  )

p_biplot <- ggplot2::ggplot() +
  ggplot2::stat_ellipse(
    data = scores_df,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = state_final,
      fill  = state_final
    ),
    type = "norm",
    level = 0.95,
    geom  = "polygon",
    alpha = 0.15,
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  ggplot2::stat_ellipse(
    data = scores_df,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = state_final
    ),
    type = "norm",
    level = 0.95,
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(
    data = scores_df,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = state_final,
      shape = factor(year)
    ),
    size  = 3,
    alpha = 0.8
  ) +
  ggplot2::geom_segment(
    data = arrow_df,
    ggplot2::aes(
      x    = x0,
      y    = y0,
      xend = x1,
      yend = y1
    ),
    arrow = ggplot2::arrow(length = grid::unit(0.25, "cm")),
    linewidth = 0.7,
    color     = "black"
  ) +
  ggrepel::geom_label_repel(
    data = arrow_df,
    ggplot2::aes(
      x     = x1,
      y     = y1,
      label = variable
    ),
    size = 3,
    label.size = 0.2,
    min.segment.length = 0
  ) +
  ggplot2::scale_color_manual(
    values = patch_colors,
    breaks = c("BAR","INCIP","FOR"),
    labels = c("Barren","Incipient","Forest"),
    name   = "Patch state"
  ) +
  ggplot2::scale_fill_manual(
    values = patch_colors,
    guide  = "none"
  ) +
  ggplot2::scale_shape_discrete(name = "Year") +
  ggplot2::labs(
    title    = "Patch habitat structure in multivariate space",
    subtitle = "Ellipses = 95% habitat-state envelopes;\nArrows = strongest habitat/cover gradients (no urchin/gonad/foraging vars)",
    x = paste0(
      "PC1 (",
      round(summary(pca_obj)$importance[2,1] * 100, 1),
      "% var)"
    ),
    y = paste0(
      "PC2 (",
      round(summary(pca_obj)$importance[2,2] * 100, 1),
      "% var)"
    )
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid     = ggplot2::element_blank(),
    legend.position = "right",
    legend.title    = ggplot2::element_text(size = 10),
    legend.text     = ggplot2::element_text(size = 9)
  )


################################################################################
# 6. Print / save key outputs --------------------------------------------------
################################################################################

print(vip_train2024)   # habitat RF importance
print(vip_driver)      # ecological driver RF importance
print(p_final)         # PDP panel across predictors
print(p_biplot)        # PCA/biplot of habitat structure

# Example save of PDP figure (edit filename/path as needed)
#ggplot2::ggsave(
#  filename = here::here("figures", "Fig_ecological_drivers_PDP.png"),
#  plot     = p_final,
#  width    = 10,
#  height   = 4,
#  dpi      = 600,
#  bg       = "white"
#)
