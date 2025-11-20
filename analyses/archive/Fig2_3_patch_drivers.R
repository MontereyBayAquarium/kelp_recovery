################################################################################
# Patch State Classification 
# Joshua G. Smith — UCSC Smith Lab / Nearshore Ecology Research Group
#
# PART A. Habitat-only RF
#   Train on 2024 diver-called patch states using ONLY structure/cover/physical
#   habitat predictors (no urchin / gonad / foraging). Exclude predictors
#   to be used in part B. Predict 2025 states.
#
# PART B. Ecological driver RF (2024 + 2025)
#   Build a unified response "state_resp" for BOTH years:
#     - 2024: diver-called state
#     - 2025: predicted state from PART A habitat RF
#   Predictors:
#     - Urchin concealment ratio (concealed / total)
#     - Urchin biomass (g / urchin)
#     - Urchin gonad mass (g / urchin)
#     - Urchin density (ind / m²)
#     - Urchin biomass density (g / m² standing stock; total_biomass_g / 80 m²)
#     - Sea otter foraging effort (proportion of all urchin biomass taken)
#
#   Output: variable importance + PDPs across predictors.
#
# PART C. PCA biplot of habitat structure
#   Multivariate habitat space, colored by patch state.
#   2024 = diver-called, 2025 = predicted from habitat RF.
################################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

require(librarian)
shelf(
  tidyverse, janitor, lubridate, sf, here,
  randomForest, vip, pdp, patchwork, ggrepel, purrr
)

set.seed(1985)

################################################################################
# 0. Setup / load data ---------------------------------------------------------
################################################################################

# Load quad_build3 etc.
load(here::here("output", "survey_data", "processed", "zone_level_data4.rda"))

# Sea otter foraging data
forage_path <- "/Volumes/enhydra/data/foraging_data/processed/foraging_data_2024_2025_processed.csv"
forage_orig <- readr::read_csv(forage_path, show_col_types = FALSE)

# Constants
focal_months  <- c(6, 7, 8, 9)  # Jun–Sep for foraging data
years_keep    <- c(2024, 2025)
patch_area_m2 <- 80            # reef area swath per patch, m²

# Patch state color palette
patch_colors <- c(
  "BAR"   = "#7570B3",
  "INCIP" = "#D95F02",
  "FOR"   = "#1B9E77"
)

################################################################################
# Helper functions -------------------------------------------------------------
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
  # arcsin(sqrt(p)) for % cover vars (0–100)
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

################################################################################
# PART A. Habitat-only RF
# Train on 2024 diver-called states using only habitat/cover/physical vars,
# then predict 2025 patch states.
################################################################################

dat_raw <- quad_build3 %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(year = lubridate::year(survey_date))

# Diver-called patch state per patch_id/year
truth_all <- dat_raw %>%
  dplyr::mutate(site_type_clean = clean_state_label(site_type)) %>%  # site_type is diver-called patch state
  dplyr::filter(year %in% years_keep) %>%
  dplyr::group_by(patch_id, site, zone, year) %>%
  dplyr::summarise(
    state = mode_char(site_type_clean),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    state = factor(state, levels = c("BAR","INCIP","FOR"))
  )

# Identify numeric predictors that are *not* biological
ban_regex <- "(urchin|gonad|biomass|forag|foraging|\\bgi\\b|_gi$|mean_gi|sd_gi)"

num_cols_all    <- names(dat_raw)[sapply(dat_raw, is.numeric)]
allowed_numeric <- num_cols_all[!grepl(ban_regex, num_cols_all, ignore.case = TRUE)]
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

# Merge with diver-called states
model_df <- hab_scaled_df %>%
  dplyr::left_join(
    truth_all,
    by = c("patch_id","site","zone","year")
  ) %>%
  dplyr::mutate(state = droplevels(state))

# Split: train on 2024 (where divers called state), predict 2025
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

# If 2025 diver calls exist, evaluate generalization
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
# PART B. Ecological drivers of patch state 
# Use BOTH YEARS (2024 + 2025).
#
# We build one response = "state_resp":
#   2024 -> diver-called state
#   2025 -> predicted_state_2025 from habitat-only RF
#
# Then relate this response to urchin / gonad / foraging predictors.
################################################################################

# 1. Summarize patch-level urchin metrics from quad_build3 ---------------------

core_bio_cols <- c(
  "mean_biomass_g",               # biomass per urchin (g / urchin)
  "mean_gonad_mass_g",            # gonad mass per urchin (g / urchin)
  "purple_urchin_densitym2",      # density (ind / m²)
  "purple_urchin_conceiledm2",    # concealed density? (typo version)
  "purple_urchin_concealedm2",    # concealed density? (alt spelling)
  "total_biomass_g",              # TOTAL standing urchin biomass in the 80 m² swath
  "total_gonad_mass_g"            # available if needed
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
    # standing urchin biomass density (g / m²) from TOTAL biomass over 80 m²
    urchin_biomass_densitym2 = total_biomass_g / patch_area_m2
  )

# Figure out which concealed-density column exists so we can build behavior_ratio
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
    }
  )

# 2. Sea otter foraging data: biomass of PUR urchins consumed per patch --------

# robust column IDs for lon/lat/etc. in foraging data
lon_col  <- first_present(forage_orig, c("long","lon","longitude","Longitude"))
lat_col  <- first_present(forage_orig, c("lat","latitude","Latitude"))
id_col   <- first_present(forage_orig, c("foragdata_id","foragdiv_id","fid","id"))
date_col <- first_present(forage_orig, c("date","Date"))

if (any(is.na(c(lon_col, lat_col, id_col)))) {
  stop("Missing lon/lat/id columns in foraging data.")
}

# build forage_dat with consistent names
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

# Summarize foraging -> per (patch_id, year)
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
    # proportion of all urchin biomass (kg) taken in that patch that year
    rel_foraging_effort = foraging_biomass_kg /
      sum(foraging_biomass_kg, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    rel_foraging_effort = tidyr::replace_na(rel_foraging_effort, 0)
  )

# 3. Build a unified "state_resp" for BOTH years -------------------------------

# 2024 diver-called state
state_lookup_2024 <- truth_all %>%
  dplyr::filter(year == 2024) %>%
  dplyr::transmute(
    patch_id, site, zone, year,
    state_resp = state  # diver-called
  )

# 2025 predicted state from habitat RF
state_lookup_2025 <- pred_2025 %>%
  dplyr::filter(year == 2025) %>%
  dplyr::transmute(
    patch_id, site, zone, year,
    state_resp = predicted_state_2025  # RF-predicted from habitat-only model
  )

state_lookup_all <- dplyr::bind_rows(
  state_lookup_2024,
  state_lookup_2025
)

# 4. Join predictors + states for both years -----------------------------------

driver_df_all <- patch_predictors_all %>%
  dplyr::filter(year %in% years_keep) %>%
  dplyr::left_join(
    foraging_by_patch %>%
      dplyr::select(
        patch_id, site, zone, year,
        rel_foraging_effort
      ),
    by = c("patch_id","site","zone","year")
  ) %>%
  dplyr::left_join(
    state_lookup_all,
    by = c("patch_id","site","zone","year")
  ) %>%
  dplyr::mutate(
    state_resp = factor(state_resp, levels = c("BAR","INCIP","FOR"))
  )

cat("\n[PART B] Rows in driver_df_all (2024 + 2025): ", nrow(driver_df_all), "\n")
cat("[PART B] Class balance in state_resp:\n")
print(table(driver_df_all$state_resp, useNA = "ifany"))

# 5. Per-patch polygons + states (2024 diver call + 2025 predicted) -----------
#    This is for exporting patch shapes w/ patch_2024 and patch_2025

# 5a. Diver-called per patch (patch_2024)
patch_calls_tbl <- quad_build3 %>%
  dplyr::mutate(
    patch_call_clean = clean_state_label(site_type)
  ) %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(patch_id, site, zone) %>%
  dplyr::summarise(
    patch_2024 = mode_char(patch_call_clean),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    patch_2024 = factor(patch_2024, levels = c("BAR","INCIP","FOR"))
  )

# 5b. Patch geometry dissolved per patch_id
patch_geom_tbl <- quad_build3 %>%
  dplyr::select(patch_id, site, zone, geometry) %>%
  dplyr::group_by(patch_id, site, zone) %>%
  dplyr::summarise(
    geometry = sf::st_union(geometry),
    .groups  = "drop"
  )

# 5c. Predicted 2025 state
patch_pred2025_tbl <- pred_2025 %>%
  dplyr::transmute(
    patch_id,
    patch_2025 = predicted_state_2025
  ) %>%
  dplyr::mutate(
    patch_2025 = factor(patch_2025, levels = c("BAR","INCIP","FOR"))
  )

# 5d. Join together for export
final_patch_sf <- patch_geom_tbl %>%
  dplyr::left_join(
    patch_calls_tbl,
    by = c("patch_id","site","zone")
  ) %>%
  dplyr::left_join(
    patch_pred2025_tbl,
    by = "patch_id"
  ) %>%
  dplyr::distinct(patch_id, .keep_all = TRUE) %>%
  dplyr::select(
    patch_id,
    site,
    zone,
    patch_2024,
    patch_2025,
    geometry
  ) %>%
  sf::st_as_sf()

transitions_tbl_constrained <- final_patch_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(
    # site_type is the diver-called 2024 patch state
    site_type = factor(
      as.character(patch_2024),
      levels = c("BAR","FOR","INCIP")
    ),
    # patch_2024 stored as *character* (not factor)
    patch_2024 = as.character(patch_2024),
    # patch_2025 stored as factor with same levels/order
    patch_2025 = factor(
      as.character(patch_2025),
      levels = c("BAR","FOR","INCIP")
    )
  ) %>%
  dplyr::select(site, zone, site_type, patch_2024, patch_2025)

# sanity check
str(transitions_tbl_constrained)

# save in the same way as before
#save(
#  transitions_tbl_constrained,
#  file = here::here("output", "lda_patch_transitionsv5.rda")
#)

# 5e. (optional) Write shapefile to disk
#out_dir <- here::here("output", "gis_data","processed","patch_state_RFsummary_2024_2025_shp")
#if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

#out_path <- file.path(out_dir, "patch_state_RFsummary_2024_2025.shp")

# sf::st_write(
#   obj          = final_patch_sf,
#   dsn          = out_path,
#   driver       = "ESRI Shapefile",
#   delete_layer = TRUE
# )

################################################################################
# 6. Prep RF input for ecological driver model --------------------------------
################################################################################

rf_driver_input <- driver_df_all %>%
  dplyr::select(
    state_resp,
    behavior_ratio,               # concealed / total
    mean_biomass_g,               # urchin biomass per urchin (g / urchin)
    mean_gonad_mass_g,            # gonad mass per urchin (g / urchin)
    purple_urchin_densitym2,      # urchin density (ind / m²)
    urchin_biomass_densitym2,     # standing biomass density (g / m²; total_biomass_g / 80 m²)
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

# Fit ecological driver RF on BOTH YEARS
rf_driver <- randomForest(
  x = rf_driver_input[, predictor_cols_B],
  y = rf_driver_input$state_resp,
  ntree      = 1500,
  mtry       = max(2, floor(sqrt(length(predictor_cols_B)))),
  importance = TRUE,
  na.action  = na.omit
)

cat("\n[PART B] Driver RF confusion matrix (trained on 2024+2025 habitat-defined states):\n")
print(rf_driver$confusion)

vip_driver <- vip::vip(
  rf_driver,
  num_features = min(12, length(predictor_cols_B)),
  bar = TRUE
) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(
    title = "Process-level drivers of patch state (2024 + 2025)",
    x     = "Predictor",
    y     = "Importance"
  )

################################################################################
# 6. Partial Dependence Curves (PDPs) ------------------------------------------
# Custom ordering and filtering for plotting only
################################################################################

# 1. Full set of predictors actually used in the RF:
#    (we keep them all in the model — including purple_urchin_densitym2 —
#     but we won't plot density)
all_predictors <- predictor_cols_B

# 2. Subset + order JUST the variables we want to show in p_final
#    Order requested:
#    1) urchin biomass (g / urchin)
#    2) urchin biomass density (g / m²)
#    3) concealment ratio
#    4) gonad mass (g / urchin)
#    5) sea otter foraging effort
pdp_plot_vars <- c(
  "mean_biomass_g",
  "urchin_biomass_densitym2",
  "behavior_ratio",
  "mean_gonad_mass_g",
  "rel_foraging_effort"
)

# 3. Pretty facet labels for ONLY those five
pretty_lab_map <- c(
  "mean_biomass_g"            = "Urchin biomass\n(g / urchin)",
  "urchin_biomass_densitym2"  = "Urchin biomass density\n(g / m²)",
  "behavior_ratio"            = "Urchin concealment ratio\n(concealed / total)",
  "mean_gonad_mass_g"         = "Urchin gonad mass\n(g / urchin)",
  "rel_foraging_effort"       = "Sea otter foraging effort\n(proportion of annual take)"
)

class_levels <- levels(rf_driver_input$state_resp)

# 4. Build long PDP data frame across predictors × states
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
      state          = cls,                   # "BAR","INCIP","FOR"
      predictor_var  = v,
      predictor_lab  = pretty_lab_map[[v]]
    )
  })
}) %>%
  tidyr::drop_na(x, y)

# 5. Make sure facet order matches pdp_plot_vars
predictor_lab_order <- unname(pretty_lab_map[pdp_plot_vars])

pdp_all <- pdp_all %>%
  dplyr::mutate(
    predictor_lab = factor(predictor_lab, levels = predictor_lab_order),
    state         = factor(state, levels = c("BAR","INCIP","FOR"))
  )

# 6. Pull variable importance to annotate facets (Gini from rf_driver)
note_labels <- rf_driver$importance %>%
  as.data.frame() %>%
  tibble::rownames_to_column("variable") %>%
  dplyr::filter(variable %in% pdp_plot_vars) %>%             # <- only keep plotted vars
  dplyr::mutate(
    predictor_lab = pretty_lab_map[variable],
    gini_note     = paste0("Gini = ", round(MeanDecreaseGini, 2))
  ) %>%
  dplyr::mutate(
    predictor_lab = factor(predictor_lab, levels = predictor_lab_order)
  ) %>%
  dplyr::distinct(predictor_lab, gini_note)

# 7. Final PDP figure (no purple_urchin_densitym2 panel, custom order)
p_final <- ggplot2::ggplot(
  pdp_all,
  ggplot2::aes(x = x, y = y, color = state, shape = state)
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
  ggplot2::scale_shape_manual(
    values = c(16, 17, 15),  # circles / triangles / squares
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

################################################################################
# PART C. PCA biplot of habitat structure
# Show multivariate habitat/cover space for each patch.
#
# Definition of "state_final":
#   - For 2024 points: diver-called patch state
#   - For 2025 points: RF-predicted patch state from Part A
#
# We generate THREE plots with identical axis ranges:
#   (1) Barrens + Forests only (no Incipient), NO vectors
#   (2) Barrens + Forests + Incipient, NO vectors
#   (3) Barrens + Forests + Incipient, WITH habitat/cover vectors
#
# Outputs:
#   p_biplot_BAR_FOR_novecs
#   p_biplot_ALL_novecs
#   p_biplot_ALL_withvecs
################################################################################


################################################################################
# 1. Build state lookup and habitat_df for PCA
################################################################################

# truth_all: from Part A/B (diver-called states by patch_id/year)
# pred_2025: from Part A (predicted_state_2025 for 2025 patches)
# hab_scaled_df: from Part A (scaled habitat/cover/physical vars per patch_id-year)

state_lookup_2024_for_biplot <- truth_all %>%
  dplyr::rename(state_2024call = state)

state_lookup_2025_for_biplot <- pred_2025 %>%
  dplyr::mutate(
    predicted_state_2025 = factor(
      predicted_state_2025,
      levels = c("BAR","INCIP","FOR")
    )
  )

habitat_df <- hab_scaled_df %>%
  dplyr::left_join(
    state_lookup_2024_for_biplot,
    by = c("patch_id","site","zone","year")
  ) %>%
  dplyr::left_join(
    state_lookup_2025_for_biplot %>%
      dplyr::select(patch_id, site, zone, year, predicted_state_2025),
    by = c("patch_id","site","zone","year")
  ) %>%
  dplyr::mutate(
    state_final = dplyr::case_when(
      year == 2024 ~ as.character(state_2024call),
      year == 2025 ~ as.character(predicted_state_2025),
      TRUE         ~ NA_character_
    ),
    # Relevel so BAR, FOR, INCIP have stable order for shapes/legend
    state_final = factor(state_final, levels = c("BAR","FOR","INCIP"))
  ) %>%
  dplyr::filter(!is.na(state_final))

################################################################################
# 2. Run PCA on habitat-only predictors (no urchin/gonad/foraging vars)
################################################################################

# We'll PCA the scaled habitat predictors from Part A.
# These were already scaled/cleaned in hab_scaled_df.

predictor_cols_pca <- setdiff(
  colnames(habitat_df),
  c("patch_id","site","zone","year",
    "state_2024call","predicted_state_2025","state_final")
)

pca_mat <- habitat_df %>%
  dplyr::select(dplyr::any_of(predictor_cols_pca)) %>%
  as.matrix()

pca_obj <- prcomp(
  pca_mat,
  center = FALSE,  # already scaled in Part A
  scale. = FALSE
)

# Scores (patch locations in PCA space)
scores_df <- tibble::as_tibble(pca_obj$x[, 1:2, drop = FALSE]) %>%
  dplyr::rename(PC1 = PC1, PC2 = PC2) %>%
  dplyr::bind_cols(
    habitat_df %>%
      dplyr::select(patch_id, site, zone, year, state_final)
  )

# Loadings (direction of each predictor in PCA space)
loadings_df <- tibble::as_tibble(
  pca_obj$rotation[, 1:2, drop = FALSE],
  rownames = "variable"
) %>%
  dplyr::rename(PC1 = PC1, PC2 = PC2) %>%
  dplyr::mutate(
    vec_len = sqrt(PC1^2 + PC2^2)
  )

################################################################################
# 3. Build arrow_df for habitat/cover gradients
#    We exclude biological/process vars like urchins, gonads, biomass, foraging.
#    Then we recode variable names to pretty labels for vectors.
################################################################################

arrow_ban_regex <- "(urchin|gonad|biomass|forag|foraging|\\bgi\\b|_gi$|mean_gi|sd_gi)"

arrow_df <- loadings_df %>%
  dplyr::filter(!grepl(arrow_ban_regex, variable, ignore.case = TRUE)) %>%
  dplyr::slice_max(order_by = vec_len, n = 8, with_ties = FALSE) %>%
  dplyr::mutate(
    variable_pretty = dplyr::recode(
      variable,
      "n_macro_plants_20m2"       = "Kelp density",
      "macro_stipe_density_20m2"  = "Kelp stipe density",
      "cov_crustose_coralline"    = "Crustose coralline cover",
      "tegula_densitym2"          = "Tegula spp density",
      "cov_desmarestia_spp"       = "Desmarestia spp cover",
      "density20m2_nerlue"        = "Bull kelp density",
      "density20m2_ptecal"        = "Pterygophora density",
      "cov_lam_holdfast_live"     = "Laminaria spp. holdfast",
      .default = variable
    )
  )

# We scale arrow vectors so they sit nicely inside the point cloud
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

################################################################################
# 4. Define shared axis limits with padding so nothing gets clipped
#    We'll reuse these exact bounds for all three plots so they match in talks.
################################################################################

lims_x_raw <- range(scores_df$PC1, na.rm = TRUE)
lims_y_raw <- range(scores_df$PC2, na.rm = TRUE)

pad_x <- diff(lims_x_raw) * 0.3  # 30% padding each side in x
pad_y <- diff(lims_y_raw) * 0.3  # 30% padding each side in y

lims_x_pad <- c(lims_x_raw[1] - pad_x, lims_x_raw[2] + pad_x)
lims_y_pad <- c(lims_y_raw[1] - pad_y, lims_y_raw[2] + pad_y)

# We want color AND shape to encode patch state.
# Lock stable shapes so slides are consistent.
shape_vals <- c(
  "BAR"   = 16,  # filled circle
  "FOR"   = 15,  # filled square
  "INCIP" = 17   # filled triangle
)

################################################################################
# 5. Helper plotting functions
#    - make_biplot_base(): no vectors
#    - make_biplot_arrows(): with vectors
#
# Both use coord_cartesian(xlim=..., ylim=...) to zoom but DO NOT drop data,
# so ellipses and arrow labels won't get clipped.
################################################################################

make_biplot_base <- function(score_data,
                             lims_x_pad,
                             lims_y_pad,
                             patch_colors) {
  
  ggplot2::ggplot() +
    # Filled 95% ellipse per state
    ggplot2::stat_ellipse(
      data = score_data,
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
    # Ellipse outline
    ggplot2::stat_ellipse(
      data = score_data,
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
    # Patch points
    ggplot2::geom_point(
      data = score_data,
      ggplot2::aes(
        x = PC1,
        y = PC2,
        color = state_final,
        shape = state_final
      ),
      size  = 3,
      alpha = 0.8
    ) +
    ggplot2::scale_color_manual(
      values = patch_colors,
      breaks = c("BAR","FOR","INCIP"),
      labels = c("Barren","Forest","Incipient"),
      name   = "Patch state"
    ) +
    ggplot2::scale_fill_manual(
      values = patch_colors,
      breaks = c("BAR","FOR","INCIP"),
      guide  = "none"
    ) +
    ggplot2::scale_shape_manual(
      values = shape_vals,
      breaks = c("BAR","FOR","INCIP"),
      labels = c("Barren","Forest","Incipient"),
      name   = "Patch state"
    ) +
    ggplot2::coord_cartesian(
      xlim   = lims_x_pad,
      ylim   = lims_y_pad,
      expand = FALSE
    ) +
    ggplot2::labs(
      x = "PC1",
      y = "PC2"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid      = ggplot2::element_blank(),
      legend.position = "right",
      legend.title    = ggplot2::element_text(size = 10),
      legend.text     = ggplot2::element_text(size = 9)
    )
}

make_biplot_arrows <- function(score_data,
                               arrow_data,
                               lims_x_pad,
                               lims_y_pad,
                               patch_colors) {
  
  ggplot2::ggplot() +
    # Filled 95% ellipse per state
    ggplot2::stat_ellipse(
      data = score_data,
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
    # Ellipse outline
    ggplot2::stat_ellipse(
      data = score_data,
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
    # Patch points
    ggplot2::geom_point(
      data = score_data,
      ggplot2::aes(
        x = PC1,
        y = PC2,
        color = state_final,
        shape = state_final
      ),
      size  = 3,
      alpha = 0.8
    ) +
    # Habitat/cover vectors
    ggplot2::geom_segment(
      data = arrow_data,
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
    # Vector labels
    ggrepel::geom_label_repel(
      data = arrow_data,
      ggplot2::aes(
        x     = x1,
        y     = y1,
        label = variable_pretty
      ),
      size = 3,
      label.size = 0.2,
      min.segment.length = 0
    ) +
    ggplot2::scale_color_manual(
      values = patch_colors,
      breaks = c("BAR","FOR","INCIP"),
      labels = c("Barren","Forest","Incipient"),
      name   = "Patch state"
    ) +
    ggplot2::scale_fill_manual(
      values = patch_colors,
      breaks = c("BAR","FOR","INCIP"),
      guide  = "none"
    ) +
    ggplot2::scale_shape_manual(
      values = shape_vals,
      breaks = c("BAR","FOR","INCIP"),
      labels = c("Barren","Forest","Incipient"),
      name   = "Patch state"
    ) +
    ggplot2::coord_cartesian(
      xlim   = lims_x_pad,
      ylim   = lims_y_pad,
      expand = FALSE
    ) +
    ggplot2::labs(
      x = "PC1",
      y = "PC2"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid      = ggplot2::element_blank(),
      legend.position = "right",
      legend.title    = ggplot2::element_text(size = 10),
      legend.text     = ggplot2::element_text(size = 9)
    )
}

################################################################################
# 6. Build the three final biplots
################################################################################

# Plot 1: Barrens + Forests only (no Incipient), NO vectors
scores_BAR_FOR <- scores_df %>%
  dplyr::filter(state_final %in% c("BAR","FOR"))

p_biplot_BAR_FOR_novecs <- make_biplot_base(
  score_data   = scores_BAR_FOR,
  lims_x_pad   = lims_x_pad,
  lims_y_pad   = lims_y_pad,
  patch_colors = patch_colors
)

# Plot 2: Barrens + Forests + Incipient, NO vectors
scores_ALL <- scores_df %>%
  dplyr::filter(state_final %in% c("BAR","FOR","INCIP"))

p_biplot_ALL_novecs <- make_biplot_base(
  score_data   = scores_ALL,
  lims_x_pad   = lims_x_pad,
  lims_y_pad   = lims_y_pad,
  patch_colors = patch_colors
)

# Plot 3: Barrens + Forests + Incipient, WITH habitat/cover vectors
p_biplot_ALL_withvecs <- make_biplot_arrows(
  score_data   = scores_ALL,
  arrow_data   = arrow_df,
  lims_x_pad   = lims_x_pad,
  lims_y_pad   = lims_y_pad,
  patch_colors = patch_colors
)

################################################################################
# 7. Print and save
################################################################################

print(p_biplot_BAR_FOR_novecs)
print(p_biplot_ALL_novecs)
print(p_biplot_ALL_withvecs)

# Save high-res PNGs for slides / manuscript figures
#ggplot2::ggsave(
#  filename = here::here("figures", "biplot_BAR_FOR_novecs.png"),
#  plot     = p_biplot_BAR_FOR_novecs,
#  width    = 6,
#  height   = 5,
#  dpi      = 600,
#  bg       = "white"
#)

#ggplot2::ggsave(
#  filename = here::here("figures", "biplot_ALL_novecs.png"),
#  plot     = p_biplot_ALL_novecs,
#  width    = 6,
#  height   = 5,
#  dpi      = 600,
#  bg       = "white"
#)

#ggplot2::ggsave(
#  filename = here::here("figures", "biplot_ALL_withvecs.png"),
#  plot     = p_biplot_ALL_withvecs,
#  width    = 6,
#  height   = 5,
#  dpi      = 600,
#  bg       = "white"
#)


#ggsave(
#  filename = here::here("figures", "Fig6_incipient_correlatesv2.png"),
#  plot = p_final,
#  width = 9,        # in inches
#  height = 4,        # adjust as needed
#  dpi = 600,         # high-res for publication
#  bg = "white"       # ensures white background if saving to PNG
#)




