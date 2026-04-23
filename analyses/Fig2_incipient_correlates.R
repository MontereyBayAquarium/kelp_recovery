################################################################################
# FIGURE 2: PCA habitat space + contrast-specific RF importance + boxplots
# PLUS: Patch transitions objects used downstream:
#   - final_patch_sf
#   - transitions_tbl_constrained   (EXACT NAME as your downstream code expects)
# and saves transitions_tbl_constrained to:
#   - output/lda_patch_transitionsv5.rda
#
# Joshua G. Smith — UCSC Smith Lab / Nearshore Ecology Research Group
################################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

require(librarian)
shelf(
  tidyverse, janitor, lubridate, sf, here,
  randomForest, patchwork, ggrepel, tidytext
)

# Reproducibility
set.seed(1985)

################################################################################
# 0. Setup / load data ---------------------------------------------------------
################################################################################

load(here::here("output", "survey_data", "processed", "zone_level_data4.rda"))

years_keep <- c(2024, 2025)

patch_colors <- c(
  "BAR"   = "#7570B3",
  "INCIP" = "#D95F02",
  "FOR"   = "#1B9E77"
)

################################################################################
# Helper functions -------------------------------------------------------------
################################################################################

mode_char <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

asin_sqrt <- function(x) {
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

################################################################################
# PART A. Habitat-only RF: train 2024 diver calls, predict 2025 ----------------
################################################################################

dat_raw <- quad_build3 %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(year = lubridate::year(survey_date))

# Diver-called patch state per patch_id/year (2024 truth used for RF training)
truth_all <- dat_raw %>%
  dplyr::filter(year %in% years_keep) %>%
  dplyr::group_by(patch_id, site, zone, year) %>%
  dplyr::summarise(
    state = mode_char(site_type),
    .groups = "drop"
  ) %>%
  dplyr::mutate(state = factor(state, levels = c("BAR","INCIP","FOR")))

# Identify numeric habitat predictors (exclude bio drivers & GI)
ban_regex <- "(urchin|gonad|biomass|forag|foraging|\\bgi\\b|_gi$|mean_gi|sd_gi)"
num_cols_all    <- names(dat_raw)[sapply(dat_raw, is.numeric)]
allowed_numeric <- num_cols_all[!grepl(ban_regex, num_cols_all, ignore.case = TRUE)]
allowed_numeric <- setdiff(allowed_numeric, "year")

# Aggregate numeric habitat/cover structure to patch_id-year
agg_hab <- dat_raw %>%
  dplyr::group_by(patch_id, site, zone, year) %>%
  dplyr::summarise(
    dplyr::across(dplyr::any_of(allowed_numeric), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

meta_hab <- agg_hab %>% dplyr::select(patch_id, site, zone, year)
X_raw <- agg_hab %>% dplyr::select(-patch_id, -site, -zone, -year)

# Transform cover % with asin(sqrt(p)), others with log1p
cover_cols <- names(X_raw)[grepl("^cov_", names(X_raw))]
other_cols <- setdiff(names(X_raw), cover_cols)

X_tr <- X_raw %>%
  dplyr::mutate(
    dplyr::across(dplyr::any_of(cover_cols), ~ asin_sqrt(.x)),
    dplyr::across(dplyr::any_of(other_cols), ~ log1p(pmax(.x, 0)))
  )

# Drop predictors with >50% NA
too_na <- sapply(X_tr, function(x) mean(!is.finite(x) | is.na(x)))
X_tr2  <- X_tr %>% dplyr::select(dplyr::any_of(names(too_na)[too_na <= 0.5]))

# Median impute
X_imp <- X_tr2 %>% dplyr::mutate(dplyr::across(dplyr::everything(), med_impute))

# Drop zero-variance cols
zv_cols <- names(X_imp)[sapply(X_imp, function(x) sd(x, na.rm = TRUE) == 0)]
if (length(zv_cols)) {
  message("[FIG 2] Dropping zero-variance predictors: ", paste(zv_cols, collapse = ", "))
  X_imp <- X_imp %>% dplyr::select(-dplyr::any_of(zv_cols))
}

# Scale predictors
X_scaled_mat <- scale(X_imp)
X_scaled     <- tibble::as_tibble(X_scaled_mat, .name_repair = "minimal")
colnames(X_scaled) <- colnames(X_imp)
hab_scaled_df <- dplyr::bind_cols(meta_hab, X_scaled)

# Merge with diver-called states
model_df <- hab_scaled_df %>%
  dplyr::left_join(truth_all, by = c("patch_id","site","zone","year")) %>%
  dplyr::mutate(state = droplevels(state))

# Split
train_df <- model_df %>% dplyr::filter(year == 2024, !is.na(state))
test_df  <- model_df %>% dplyr::filter(year == 2025)

predictor_cols_A <- setdiff(colnames(train_df), c("patch_id","site","zone","year","state"))

# --- Drop spatial coordinates from habitat predictors (for contrast RFs + importance) ---
spatial_regex <- "(^lat$|^lon$|^long$|latitude|longitude)"

predictor_cols_A_nospatial <- predictor_cols_A[
  !grepl(spatial_regex, predictor_cols_A, ignore.case = TRUE)
]


# Fit habitat RF (2024)
set.seed(1985)
rf_train2024 <- randomForest::randomForest(
  x = train_df[, predictor_cols_A],
  y = train_df$state,
  ntree      = 1500,
  mtry       = max(2, floor(sqrt(length(predictor_cols_A)))),
  importance = TRUE,
  na.action  = na.omit
)

# Predict 2025
test_pred_class <- predict(rf_train2024, newdata = test_df[, predictor_cols_A], type = "response")

pred_2025 <- test_df %>%
  dplyr::select(patch_id, site, zone, year) %>%
  dplyr::mutate(predicted_state_2025 = factor(test_pred_class, levels = c("BAR","INCIP","FOR"))) %>%
  dplyr::distinct(patch_id, site, zone, year, .keep_all = TRUE)

################################################################################
# Patch-level transitions table for LDA (2024 diver calls → 2025 RF states) -----
# (This is your original block, with the same object names)
################################################################################

# 1. Diver-called per patch (patch_2024)
patch_calls_tbl <- quad_build3 %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(year = lubridate::year(survey_date)) %>%
  dplyr::filter(year == 2024) %>%  # make explicit: diver calls are from 2024
  dplyr::group_by(patch_id, site, zone) %>%
  dplyr::summarise(
    patch_2024 = mode_char(site_type),
    .groups    = "drop"
  ) %>%
  dplyr::mutate(
    patch_2024 = factor(patch_2024, levels = c("BAR","INCIP","FOR"))
  )

# 2. Patch geometry dissolved per patch_id
patch_geom_tbl <- quad_build3 %>%
  dplyr::select(patch_id, site, zone, geometry) %>%
  dplyr::group_by(patch_id, site, zone) %>%
  dplyr::summarise(
    geometry = sf::st_union(geometry),
    .groups  = "drop"
  )

# 3. Predicted 2025 state from habitat RF
patch_pred2025_tbl <- pred_2025 %>%
  dplyr::transmute(
    patch_id,
    patch_2025 = predicted_state_2025
  ) %>%
  dplyr::mutate(
    patch_2025 = factor(patch_2025, levels = c("BAR","INCIP","FOR"))
  )

# 4. Join together: patch_2024 (diver) + patch_2025 (RF) + geometry
final_patch_sf <- patch_geom_tbl %>%
  dplyr::left_join(patch_calls_tbl, by = c("patch_id","site","zone")) %>%
  dplyr::left_join(patch_pred2025_tbl, by = "patch_id") %>%
  dplyr::distinct(patch_id, .keep_all = TRUE) %>%
  dplyr::select(patch_id, site, zone, patch_2024, patch_2025, geometry) %>%
  sf::st_as_sf()

# 5. LDA-friendly transitions table (no geometry)
transitions_tbl_constrained <- final_patch_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(
    site_type = factor(as.character(patch_2024), levels = c("BAR","FOR","INCIP")),
    patch_2024 = as.character(patch_2024),
    patch_2025 = factor(as.character(patch_2025), levels = c("BAR","FOR","INCIP"))
  ) %>%
  dplyr::select(site, zone, site_type, patch_2024, patch_2025)

str(transitions_tbl_constrained)

# SAVE with the exact same object name as downstream expects
#save(
#  transitions_tbl_constrained,
#  file = here::here("output", "lda_patch_transitionsv6.rda")
#)

################################################################################
# Export: state lookup for Fig 3 (year-indexed, lightweight) -------------------
################################################################################

state_lookup_2024 <- truth_all %>%
  dplyr::filter(year == 2024) %>%
  dplyr::transmute(patch_id, site, zone, year, state_resp = state)

state_lookup_2025 <- pred_2025 %>%
  dplyr::filter(year == 2025) %>%
  dplyr::transmute(patch_id, site, zone, year, state_resp = predicted_state_2025)

state_lookup_all <- dplyr::bind_rows(state_lookup_2024, state_lookup_2025)

#readr::write_csv(
#  state_lookup_all,
#  here::here("output", "patch_state_lookup_2024_2025.csv")
#)

################################################################################
# PART C (Fig 2A): PCA biplot of habitat structure -----------------------------
################################################################################

state_lookup_2024_for_biplot <- truth_all %>% dplyr::rename(state_2024call = state)
state_lookup_2025_for_biplot <- pred_2025 %>%
  dplyr::mutate(predicted_state_2025 = factor(predicted_state_2025, levels = c("BAR","INCIP","FOR")))

habitat_df <- hab_scaled_df %>%
  dplyr::left_join(state_lookup_2024_for_biplot, by = c("patch_id","site","zone","year")) %>%
  dplyr::left_join(
    state_lookup_2025_for_biplot %>% dplyr::select(patch_id, site, zone, year, predicted_state_2025),
    by = c("patch_id","site","zone","year")
  ) %>%
  dplyr::mutate(
    state_final = dplyr::case_when(
      year == 2024 ~ as.character(state_2024call),
      year == 2025 ~ as.character(predicted_state_2025),
      TRUE ~ NA_character_
    ),
    state_final = factor(state_final, levels = c("BAR","FOR","INCIP"))
  ) %>%
  dplyr::filter(!is.na(state_final))

predictor_cols_pca <- setdiff(
  colnames(habitat_df),
  c("patch_id","site","zone","year","state_2024call","predicted_state_2025","state_final")
)

pca_mat <- habitat_df %>%
  dplyr::select(dplyr::any_of(predictor_cols_pca)) %>%
  as.matrix()

pca_obj <- prcomp(pca_mat, center = FALSE, scale. = FALSE)

scores_df <- tibble::as_tibble(pca_obj$x[, 1:2, drop = FALSE]) %>%
  dplyr::rename(PC1 = PC1, PC2 = PC2) %>%
  dplyr::bind_cols(habitat_df %>% dplyr::select(patch_id, site, zone, year, state_final))

loadings_df <- tibble::as_tibble(pca_obj$rotation[, 1:2, drop = FALSE], rownames = "variable") %>%
  dplyr::rename(PC1 = PC1, PC2 = PC2) %>%
  dplyr::mutate(vec_len = sqrt(PC1^2 + PC2^2))

arrow_ban_regex <- "(urchin|gonad|biomass|forag|foraging|\\bgi\\b|_gi\\b|mean_gi|sd_gi|lat|lon|long)"
arrow_df <- loadings_df %>%
  dplyr::filter(!grepl(arrow_ban_regex, variable, ignore.case = TRUE)) %>%
  dplyr::slice_max(order_by = vec_len, n = 8, with_ties = FALSE) %>%
  dplyr::mutate(
    variable_pretty = dplyr::recode(
      variable,
      "n_macro_plants_20m2"       = "Macrocystis density",
      "macro_stipe_density_20m2"  = "Macrocystis stipe density",
      "cov_crustose_coralline"    = "Crustose coralline cover",
      "tegula_densitym2"          = "Tegula spp density",
      "cov_desmarestia_spp"       = "Desmarestia spp cover",
      "density20m2_nerlue"        = "Bull kelp density",
      "density20m2_ptecal"        = "Pterygophora density",
      "cov_lam_holdfast_live"     = "Laminaria spp. holdfast",
      .default = variable
    )
  )

range_x <- range(scores_df$PC1, na.rm = TRUE)
range_y <- range(scores_df$PC2, na.rm = TRUE)
max_span <- max(diff(range_x), diff(range_y))
max_vec  <- max(arrow_df$vec_len)
arrow_scaler <- 0.4 * max_span / max_vec

arrow_df <- arrow_df %>%
  dplyr::mutate(x0 = 0, y0 = 0, x1 = PC1 * arrow_scaler, y1 = PC2 * arrow_scaler)

shape_vals <- c("BAR" = 16, "FOR" = 15, "INCIP" = 17)

make_biplot_arrows <- function(score_data, arrow_data, patch_colors) {
  ggplot2::ggplot() +
    ggplot2::stat_ellipse(
      data = score_data,
      ggplot2::aes(x = PC1, y = PC2, color = state_final, fill = state_final),
      type = "norm", level = 0.95, geom = "polygon",
      alpha = 0.15, linewidth = 0.6, show.legend = FALSE
    ) +
    ggplot2::stat_ellipse(
      data = score_data,
      ggplot2::aes(x = PC1, y = PC2, color = state_final),
      type = "norm", level = 0.95,
      linewidth = 0.6, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = score_data,
      ggplot2::aes(x = PC1, y = PC2, color = state_final, shape = state_final),
      size = 3, alpha = 0.8
    ) +
    ggplot2::geom_segment(
      data = arrow_data,
      ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1),
      arrow = ggplot2::arrow(length = grid::unit(0.25, "cm")),
      linewidth = 0.7, color = "black"
    ) +
    ggrepel::geom_label_repel(
      data = arrow_data,
      ggplot2::aes(x = x1, y = y1, label = variable_pretty),
      size = 3, label.size = 0.2, min.segment.length = 0
    ) +
    ggplot2::scale_color_manual(
      values = patch_colors,
      breaks = c("BAR","FOR","INCIP"),
      labels = c("Barren","Forest","Incipient"),
      name   = "Patch state"
    ) +
    ggplot2::scale_fill_manual(values = patch_colors, guide = "none") +
    ggplot2::scale_shape_manual(
      values = shape_vals,
      breaks = c("BAR","FOR","INCIP"),
      labels = c("Barren","Forest","Incipient"),
      name   = "Patch state"
    ) +
    ggplot2::coord_equal(xlim = c(-8.8, 8), ylim = c(-8.8, 8), expand = FALSE) +
    ggplot2::labs(x = "PC1", y = "PC2") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = "right")
}

scores_ALL <- scores_df %>% dplyr::filter(state_final %in% c("BAR","FOR","INCIP"))
p_biplot_ALL_withvecs <- make_biplot_arrows(scores_ALL, arrow_df, patch_colors)

################################################################################
# Fig 2B/C: contrast-specific RF importance + boxplots -------------------------
################################################################################

fit_binary_rf_hab <- function(df, pos = "INCIP", neg = "BAR", predictors) {
  df2 <- df %>%
    dplyr::filter(state %in% c(pos, neg)) %>%
    dplyr::mutate(y = factor(ifelse(state == pos, pos, neg), levels = c(neg, pos)))
  
  randomForest::randomForest(
    x = df2[, predictors],
    y = df2$y,
    ntree = 1500,
    mtry  = max(2, floor(sqrt(length(predictors)))),
    importance = TRUE,
    na.action = na.omit
  )
}

hab_train_bin <- train_df %>%
  dplyr::select(state, dplyr::all_of(predictor_cols_A_nospatial)) %>%
  dplyr::filter(!is.na(state))


set.seed(1985)
rf_incip_vs_bar <- fit_binary_rf_hab(
  hab_train_bin,
  pos = "INCIP", neg = "BAR",
  predictors = predictor_cols_A_nospatial
)

set.seed(1985)
rf_incip_vs_for <- fit_binary_rf_hab(
  hab_train_bin,
  pos = "INCIP", neg = "FOR",
  predictors = predictor_cols_A_nospatial
)


get_perm_imp <- function(rf_obj) {
  imp <- randomForest::importance(rf_obj, type = 1)
  if ("MeanDecreaseAccuracy" %in% colnames(imp)) imp[, "MeanDecreaseAccuracy"] else imp[, 1]
}

imp_incip_bar <- get_perm_imp(rf_incip_vs_bar)
imp_incip_for <- get_perm_imp(rf_incip_vs_for)

pretty_labs <- c(
  cov_fleshy_red            = "Fleshy red algae cov.",
  cov_mac_holdfast_live     = "Macrocystis cov.",
  cov_phragmatopoma         = "Phragmatopoma spp. cov.",
  density20m2_macstump      = "Dead holdfast cov.",
  macj                      = "Juv. Macrocystis den.",
  ptej                      = "Juv. Pterygophora den.",
  cov_diopatra_chaetopterus = "Diopatra spp. cov.",
  density20m2_lamset        = "L. setchelli den.",
  cov_crustose_coralline    = "Crustose coralline cov.",
  tegula_densitym2          = "Tegula spp. den.",
  n_macro_plants_20m2       = "Macrocystis den.",
  density20m2_ptecal        = "Pterygophora den.",
  macro_stipe_density_20m2  = "Macrocystis stipe den.",
  cov_demarestia_spp        = "Desmarestia cov.",
  density20m2_nerlue        = "Nereocystis luetkeana den.",
  cov_lam_holdfast_live     = "Lam. holdfast cov.",
  cov_articulated_coralline = "Articulated coralline cov.",
  pomaulax_densitym2        = "Pomaulax spp. den.",
  cov_encrusting_red        = "Encrusting red algae cov.",
  relief_cm                 = "Vertical relief",
  lamr                      = "Lam. recruit den.",
  cov_bare_sand             = "Bare sand cov.",
  risk_index                = "Reef rugosity",
  cov_desmarestia_spp       = "Desmarestia spp. cov.",
  cov_dodecaceria_spp       = "Dodecaceria spp den.",
  cov_barnacle              = "Barnacle cov.",
  cov_dictyoneurum_spp      = "Dictyoneurum spp cov.",
  density20m2_eisarb        = "Eisenia arborea den.",
  cov_hydroids              = "Hydroids cov."
)

topN_B <- 12

imp_tbl <- tibble::tibble(
  variable   = c(names(imp_incip_bar), names(imp_incip_for)),
  importance = c(unname(imp_incip_bar), unname(imp_incip_for)),
  contrast   = c(rep("Incipient vs Barren", length(imp_incip_bar)),
                 rep("Incipient vs Forest", length(imp_incip_for)))
) %>%
  dplyr::mutate(variable_pretty = dplyr::recode(variable, !!!pretty_labs, .default = variable)) %>%
  dplyr::group_by(contrast) %>%
  dplyr::arrange(dplyr::desc(importance), .by_group = TRUE) %>%
  dplyr::slice_head(n = topN_B) %>%
  dplyr::ungroup()

p_B <- imp_tbl %>%
  dplyr::mutate(
    # reorder *within each facet*
    variable_pretty = tidytext::reorder_within(variable_pretty, importance, contrast)
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = variable_pretty, y = importance)) +
  ggplot2::geom_col(fill = "grey40") +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(~ contrast, scales = "free_y") +
  tidytext::scale_x_reordered() +
  ggplot2::theme_bw() +
  ggplot2::labs(x = NULL, y = "Permutation importance (Mean decrease accuracy)") +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 10)
  )


vars_C <- imp_tbl %>% dplyr::pull(variable) %>% unique()

patch_colors_named <- c(
  "Barren"    = patch_colors[["BAR"]],
  "Forest"    = patch_colors[["FOR"]],
  "Incipient" = patch_colors[["INCIP"]]
)

my_theme <- ggplot2::theme(
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
  legend.background = ggplot2::element_rect(fill = scales::alpha("blue", 0)),
  strip.text       = ggplot2::element_text(size = 10, face = "bold", color = "black", hjust = 0),
  strip.background = ggplot2::element_blank()
)

trim_vars <- c("n_macro_plants_20m2","density20m2_ptecal","lamr","density20m2_nerlue","cov_dodecaceria_spp")

box_df <- dat_raw %>%
  dplyr::filter(density20m2_lamset < 3 | is.na(density20m2_lamset)) %>%
  dplyr::left_join(truth_all, by = c("patch_id","site","zone","year")) %>%
  dplyr::select(state, dplyr::all_of(vars_C)) %>%
  tidyr::pivot_longer(cols = -state, names_to = "variable", values_to = "value") %>%
  dplyr::group_by(variable) %>%
  dplyr::mutate(
    q_low  = if (first(variable) %in% trim_vars) stats::quantile(value, 0.01, na.rm = TRUE) else -Inf,
    q_high = if (first(variable) %in% trim_vars) stats::quantile(value, 0.90, na.rm = TRUE) else  Inf
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(value >= q_low, value <= q_high | is.na(value)) %>%
  dplyr::mutate(
    variable_pretty = dplyr::recode(variable, !!!pretty_labs, .default = variable),
    state_label = factor(state, levels = c("BAR","FOR","INCIP"), labels = c("Barren","Forest","Incipient"))
  )

facet_order <- imp_tbl %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(max_imp = max(importance, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(max_imp)) %>%
  dplyr::mutate(variable_pretty = dplyr::recode(variable, !!!pretty_labs, .default = variable)) %>%
  dplyr::pull(variable_pretty)

box_df <- box_df %>% dplyr::mutate(variable_pretty = factor(variable_pretty, levels = facet_order))

p_C <- ggplot2::ggplot(box_df, ggplot2::aes(x = state_label, y = value, fill = state_label, color = state_label)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.75, linewidth = 0.4) +
  ggplot2::facet_wrap(~ variable_pretty, scales = "free_y", ncol = 5) +
  ggplot2::scale_fill_manual(values = patch_colors_named, name = "Patch state") +
  ggplot2::scale_color_manual(values = patch_colors_named, guide = "none") +
  ggplot2::theme_bw() + my_theme +
  ggplot2::theme(
    legend.position  = "bottom",
    strip.text       = ggplot2::element_text(size = 8, face = "bold"),
    axis.title.x     = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(size = 7, angle = 20, hjust = 1),
    axis.text.y      = ggplot2::element_text(size = 7),
    panel.spacing    = grid::unit(0.4, "lines"),
    panel.border     = ggplot2::element_rect(color = "black", linewidth = 0.5)
  ) +
  ggplot2::labs(y = "Observed value", x = NULL)

################################################################################
# Assemble Fig 2 ---------------------------------------------------------------
################################################################################

p_A <- p_biplot_ALL_withvecs + ggplot2::theme(legend.position = "none")

#top_row <- p_A + p_B + patchwork::plot_spacer() +
 # patchwork::plot_layout(widths = c(1.15, 0.85, 0.05))

top_row <- p_A + p_B + patchwork::plot_spacer() +
  patchwork::plot_layout(widths = c(1.05, 1.05, 0.02))


Fig2 <- (top_row) / p_C +
  patchwork::plot_layout(heights = c(1, 1.25)) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(
    plot.tag    = ggplot2::element_text(size = 10),
    plot.margin = ggplot2::margin(t = 5, r = 10, b = 5, l = 10)
  )

Fig2

#ggsave(
#  filename = here::here("figures", "Fig2_patch_habitat_correlatesv3.png"),
#  plot     = Fig2,
#  width    = 8.5,
#  height   = 8.5,
#  dpi      = 600,
#  bg       = "white"
#)



################################################################################
#Alternate plot with just Panels A and B


################################################################################
# Assemble Fig 2: A above B (drop C)
################################################################################

# If you want wrapped facet labels in B (prevents cutoff)
imp_tbl <- imp_tbl %>%
  dplyr::mutate(
    contrast = dplyr::recode(
      contrast,
      "Incipient vs Barren" = "Incipient vs\nBarren",
      "Incipient vs Forest" = "Incipient vs\nForest"
    )
  )

# Rebuild p_B AFTER changing imp_tbl if your p_B uses it
p_B <- imp_tbl %>%
  dplyr::mutate(variable_pretty = tidytext::reorder_within(variable_pretty, importance, contrast)) %>%
  ggplot2::ggplot(ggplot2::aes(x = variable_pretty, y = importance)) +
  ggplot2::geom_col(fill = "grey40") +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(~ contrast, scales = "free_y") +
  tidytext::scale_x_reordered() +
  ggplot2::theme_bw() +
  ggplot2::labs(x = NULL, y = "Permutation importance (Mean decrease accuracy)") +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 10, margin = ggplot2::margin(t = 3, b = 3)),
    plot.margin = ggplot2::margin(t = 6, r = 6, b = 6, l = 6)
  )

p_A <- p_biplot_ALL_withvecs + ggplot2::theme(legend.position = "right")  # or "none"

Fig2 <- (p_A / p_B) +
  patchwork::plot_layout(heights = c(1.2, 0.8)) +   # tweak as desired
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(size = 10),
    plot.margin = ggplot2::margin(t = 5, r = 10, b = 5, l = 10)
  )

Fig2

ggsave(
  filename = here::here("figures", "Fig2_patch_habitat_correlates_ABv2.png"),
  plot     = Fig2,
  width    = 8.5,
  height   = 8.5,
  dpi      = 600,
  bg       = "white"
)

###############################################################################
#Move old panel C to supp


# Variables to show = union of top variables from both contrasts
vars_C <- imp_tbl %>% dplyr::pull(variable) %>% unique()

# Color mapping (labels used in the plot)
patch_colors_named <- c(
  "Barren"    = patch_colors[["BAR"]],
  "Forest"    = patch_colors[["FOR"]],
  "Incipient" = patch_colors[["INCIP"]]
)


trim_vars <- c("n_macro_plants_20m2","density20m2_ptecal","lamr","density20m2_nerlue","cov_dodecaceria_spp")

# Facet order = by max importance across contrasts 
facet_order <- imp_tbl %>%
 group_by(variable) %>%
  summarise(max_imp = max(importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(dplyr::desc(max_imp)) %>%
 mutate(variable_pretty = dplyr::recode(variable, !!!pretty_labs, .default = variable)) %>%
  pull(variable_pretty)

# Build long data for boxplots (RAW values), 
box_df <- dat_raw %>%
  filter(density20m2_lamset < 3 | is.na(density20m2_lamset)) %>%
  left_join(truth_all, by = c("patch_id","site","zone","year")) %>%
  select(state, dplyr::all_of(vars_C)) %>%
pivot_longer(cols = -state, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  mutate(
    q_low  = if (first(variable) %in% trim_vars) stats::quantile(value, 0.01, na.rm = TRUE) else -Inf,
    q_high = if (first(variable) %in% trim_vars) stats::quantile(value, 0.90, na.rm = TRUE) else  Inf
  ) %>%
  ungroup() %>%
  filter(value >= q_low, value <= q_high | is.na(value)) %>%
  mutate(
    variable_pretty = dplyr::recode(variable, !!!pretty_labs, .default = variable),
    state_label = factor(state, levels = c("BAR","FOR","INCIP"),
                         labels = c("Barren","Forest","Incipient"))
  ) %>%
  mutate(variable_pretty = factor(variable_pretty, levels = facet_order))

# Plot Panel C only
p_C <- ggplot(
  box_df,
  aes(x = state_label, y = value, fill = state_label, color = state_label)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, linewidth = 0.4) +
  facet_wrap(~ variable_pretty, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = patch_colors_named, name = "Patch state") +
  scale_color_manual(values = patch_colors_named, guide = "none") +
  theme_bw() + my_theme +
  theme(
    legend.position  = "bottom",
    strip.text       = ggplot2::element_text(size = 8, face = "bold"),
    axis.title.x     = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(size = 7, angle = 20, hjust = 1),
    axis.text.y      = ggplot2::element_text(size = 7),
    panel.spacing    = grid::unit(0.4, "lines"),
    panel.border     = ggplot2::element_rect(color = "black", linewidth = 0.5)
  ) +
  labs(y = "Observed value", x = NULL)

p_C

# Save
ggsave(
  filename = here::here("figures", "FigS_panelC_boxplots.png"),
  plot     = p_C,
  width    = 8.5,
  height   = 8.5,
  dpi      = 600,
  bg       = "white"
)







################################################################################
# Summary statistics for Fig. 2 ------------------------------------------------
################################################################################



library(vegan)

################################################################################
# 1. Sample sizes by patch state ------------------------------------------------
################################################################################

patch_state_ns <- habitat_df %>%
  dplyr::count(state_final, name = "n") %>%
  dplyr::mutate(
    state_final = factor(state_final,
                         levels = c("BAR", "FOR", "INCIP"),
                         labels = c("Barren", "Forest", "Incipient"))
  )

print(patch_state_ns)

################################################################################
# 2. PCA variance explained -----------------------------------------------------
################################################################################

pca_var <- (pca_obj$sdev)^2
pca_var_exp <- pca_var / sum(pca_var)
pca_var_tbl <- tibble::tibble(
  PC = paste0("PC", seq_along(pca_var_exp)),
  eigenvalue = pca_var,
  prop_var = pca_var_exp,
  cum_var = cumsum(pca_var_exp)
)

print(pca_var_tbl %>% dplyr::slice(1:10))

cat("\nPCA variance explained:\n")
cat("PC1 =", round(100 * pca_var_exp[1], 1), "%\n")
cat("PC2 =", round(100 * pca_var_exp[2], 1), "%\n")
cat("PC1 + PC2 =", round(100 * sum(pca_var_exp[1:2]), 1), "%\n")

################################################################################
# 3. Top PCA loadings -----------------------------------------------------------
################################################################################

loading_tbl <- tibble::as_tibble(
  pca_obj$rotation[, 1:2, drop = FALSE],
  rownames = "variable"
) %>%
  dplyr::rename(PC1 = PC1, PC2 = PC2) %>%
  dplyr::mutate(
    abs_PC1 = abs(PC1),
    abs_PC2 = abs(PC2),
    vec_len = sqrt(PC1^2 + PC2^2)
  ) %>%
  dplyr::arrange(dplyr::desc(vec_len))

cat("\nTop loadings by overall vector length:\n")
print(loading_tbl %>% dplyr::slice(1:15))

cat("\nTop positive loadings on PC1:\n")
print(loading_tbl %>% dplyr::arrange(dplyr::desc(PC1)) %>% dplyr::slice(1:10))

cat("\nTop negative loadings on PC1:\n")
print(loading_tbl %>% dplyr::arrange(PC1) %>% dplyr::slice(1:10))

cat("\nTop positive loadings on PC2:\n")
print(loading_tbl %>% dplyr::arrange(dplyr::desc(PC2)) %>% dplyr::slice(1:10))

cat("\nTop negative loadings on PC2:\n")
print(loading_tbl %>% dplyr::arrange(PC2) %>% dplyr::slice(1:10))

################################################################################
# 4. PERMANOVA -----------------------------------------------------------------
################################################################################

# This uses the same transformed/scaled habitat matrix that went into the PCA.
# Euclidean distance is appropriate here because PCA is based on Euclidean space.

permanova_res <- vegan::adonis2(
  pca_mat ~ state_final,
  data = habitat_df,
  method = "euclidean",
  permutations = 9999
)

cat("\nOverall PERMANOVA:\n")
print(permanova_res)

################################################################################
# 5. Dispersion test (recommended alongside PERMANOVA) --------------------------
################################################################################

# Important: PERMANOVA can be influenced by differences in within-group dispersion.
# This checks whether spread differs among groups.

dist_mat <- dist(pca_mat, method = "euclidean")

disp_mod <- vegan::betadisper(dist_mat, habitat_df$state_final)
disp_anova <- anova(disp_mod)
disp_perm  <- permutest(disp_mod, permutations = 9999)

cat("\nHomogeneity of multivariate dispersion:\n")
print(disp_anova)
print(disp_perm)

################################################################################
# 6. Optional pairwise PERMANOVA ------------------------------------------------
################################################################################

# Uncomment if you want pairwise comparisons among states
# library(pairwiseAdonis)
#
# pairwise_perm <- pairwiseAdonis::pairwise.adonis2(
#   x = pca_mat,
#   factors = habitat_df$state_final,
#   sim.method = "euclidean",
#   p.adjust.m = "BH"
# )
#
# cat("\nPairwise PERMANOVA:\n")
# print(pairwise_perm)

################################################################################
# 7. Nicely formatted values for manuscript text --------------------------------
################################################################################

# Pull overall PERMANOVA values
perm_F  <- permanova_res$F[1]
perm_R2 <- permanova_res$R2[1]
perm_p  <- permanova_res$`Pr(>F)`[1]

cat("\nManuscript-ready summary:\n")
cat(
  paste0(
    "PC1 explained ", round(100 * pca_var_exp[1], 1),
    "% and PC2 explained ", round(100 * pca_var_exp[2], 1),
    "% of the variation (", round(100 * sum(pca_var_exp[1:2]), 1),
    "% cumulative). PERMANOVA indicated that habitat structure differed among patch states ",
    "(F = ", round(perm_F, 2),
    ", R2 = ", round(perm_R2, 3),
    ", p = ", format.pval(perm_p, digits = 3), ")."
  ),
  "\n"
)





