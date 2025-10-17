# ==============================================================================
# Sea otter foraging: Prey-choice GLMM + manual prediction plots
# Focus: probability of urchin-focused foraging vs. density × gonad index + patch
# Joshua G. Smith – UCSC Nearshore Ecology Research Group
# ==============================================================================

rm(list = ls())

require(librarian)
librarian::shelf(
  tidyverse, janitor, lubridate, sf,
  glmmTMB, here, patchwork, ggplot2
)

# ------------------------------------------------------------------------------
# 1) Load data
# ------------------------------------------------------------------------------
load(here::here("output","survey_data","processed","zone_level_data3.rda")) # quad_build3
forage_orig     <- read_csv("/Volumes/enhydra/data/foraging_data/processed/foraging_data_2024_2025_processed.csv")
dissection_data <- read_csv("/Volumes/enhydra/data/kelp_recovery/MBA_kelp_forest_database/processed/dissection/dissection_data_recovery.csv")

years_keep   <- c(2024, 2025)
patch_colors <- c("BAR"="purple", "INCIP"="orange", "FOR"="forestgreen")

# ------------------------------------------------------------------------------
# 2) Dissection summaries → join to patch data (quad_build4)
# ------------------------------------------------------------------------------
dissect_build1 <- dissection_data %>%
  mutate(
    year = year(survey_date)
  ) %>%
  group_by(year, site_official, site_type_official, zone, species) %>%
  summarise(
    #mean_gonad_mass_g = mean(gonad_mass_g, na.rm = TRUE),
    mean_gonad_index  = mean(gonad_index,  na.rm = TRUE),
    sd_gonad_index = sd(gonad_index,  na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = species,
    values_from = c(sd_gonad_index, mean_gonad_index, n),
    names_sep   = "_"
  ) %>% select(-mean_gonad_index_red_urchin, -n_red_urchin, -sd_gonad_index_red_urchin)

quad_build4 <- quad_build3 %>%
  mutate(year = lubridate::year(survey_date)) %>%
  left_join(dissect_build1,
            by = c("year","site"="site_number","site_type","zone")) %>%
  mutate(pred_patch = case_when(
    str_detect(tolower(pred_patch), "bar")   ~ "BAR",
    str_detect(tolower(pred_patch), "incip") ~ "INCIP",
    str_detect(tolower(pred_patch), "for")   ~ "FOR",
    TRUE ~ toupper(pred_patch)
  ))

# ------------------------------------------------------------------------------
# 3) Patch-year predictor table (keep density & gonad index; others optional)
# ------------------------------------------------------------------------------
predictors_focus <- c(
  "purple_urchin_densitym2",
  "mean_gonad_index_purple_urchin",
  "n_macro_plants_20m2",
  "density20m2_nerlue"
)

quad_year <- quad_build4 %>%
  mutate(year = year(survey_date)) %>%
  filter(year %in% years_keep) %>%
  st_drop_geometry() %>%
  group_by(patch_id, year, pred_patch) %>%
  summarise(across(all_of(predictors_focus), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

# ------------------------------------------------------------------------------
# 4) Build bouts with prey composition (counterfactual: urchin vs other prey)
# ------------------------------------------------------------------------------
# Assumes forage_orig has columns: long, lat, year, month, bout, prey (strings)
# Keep summer/fall months if desired:
forage_sf <- forage_orig %>%
  #filter(month %in% c(5,6,7,8,9)) %>%
  st_as_sf(coords = c("long","lat"), crs = 4326, remove = FALSE)

# Helper: flag urchin dives
is_urchin <- function(x) {
  # Expand patterns as needed (scientific names, etc.)
  str_detect(x, regex("urch|red|pur|uni|Strongylocentrotus|Mesocentrotus", ignore_case = TRUE))
}

# Summarise each bout by prey composition and centroid
#bout_summary <- forage_sf %>%
#  group_by(bout) %>%
#  summarise(
#    year         = first(year),
#    total_dives  = n(),
#    n_urchin     = sum(is_urchin(prey), na.rm = TRUE),
#    n_other      = sum(!is_urchin(prey), na.rm = TRUE),
#    prop_urchin  = if_else(total_dives > 0, n_urchin / total_dives, NA_real_),
#    prey_types   = paste(sort(unique(prey)), collapse = ","),
#    geometry     = st_centroid(st_union(geometry)),
#    .groups      = "drop"
#  ) %>%
#  mutate(
#    # Dominant-prey classification (you can adjust threshold)
#    prey_focus = case_when(
#      is.na(prop_urchin) | total_dives == 0 ~ NA_character_,
#      prop_urchin > 0.5                     ~ "urchin_focus",
#      prop_urchin > 0 & prop_urchin <= 0.5  ~ "other_focus",
#      TRUE                                  ~ NA_character_
#    )
#  )

# Summarise each bout by prey composition and centroid
bout_summary <- forage_sf %>%
  group_by(bout) %>%
  summarise(
    year         = first(year),
    total_dives  = n(),
    n_urchin     = sum(is_urchin(prey), na.rm = TRUE),
    n_other      = sum(!is_urchin(prey), na.rm = TRUE),
    prop_urchin  = if_else(total_dives > 0, n_urchin / total_dives, NA_real_),
    prey_types   = paste(sort(unique(prey)), collapse = ","),
    geometry     = st_centroid(st_union(geometry)),
    .groups      = "drop"
  ) %>%
  mutate(
    # Revised classifier:
    # - urchin_focus if >50% of dives are urchin AND n_urchin > 3
    # - other_focus otherwise
    prey_focus = case_when(
      is.na(total_dives) | total_dives == 0 ~ NA_character_,
      n_urchin > 2 & prop_urchin > 0.5      ~ "urchin_focus",
      TRUE                                  ~ "other_focus"
    )
  )


# ------------------------------------------------------------------------------
# 5) Join bouts to polygons (patch_id, pred_patch) by spatial intersection & year
# ------------------------------------------------------------------------------
quad4_same_crs <- st_transform(quad_build4, st_crs(bout_summary))

bouts_with_patch <- st_join(
  bout_summary,
  quad4_same_crs %>% dplyr::select(patch_id, pred_patch, survey_date),
  join = st_intersects,
  left = TRUE
) %>%
  mutate(patch_year = lubridate::year(survey_date)) %>%
  filter(year %in% years_keep, !is.na(patch_id), patch_year == year) %>%
  st_drop_geometry() %>%
  distinct(bout, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# 6) Create "used patch" table with predictors (restrict to bouts with a prey_focus)
# ------------------------------------------------------------------------------
rf_choice <- bouts_with_patch %>%
  filter(!is.na(prey_focus)) %>%
  left_join(quad_year, by = c("patch_id","year","pred_patch")) %>%
  # Drop rows lacking key predictors
  drop_na(purple_urchin_densitym2, mean_gonad_index_purple_urchin, pred_patch) %>%
  mutate(
    pred_patch = factor(pred_patch, levels = c("BAR","INCIP","FOR")),
    urchin_focus_bin = if_else(prey_focus == "urchin_focus", 1L, 0L)
  )

# Quick QC: composition by patch
rf_choice %>%
  count(pred_patch, prey_focus) %>%
  arrange(pred_patch, desc(n)) %>%
  print(n = 100)

# ------------------------------------------------------------------------------
# 7) Fit GLMM: P(urchin-focus | used patch) ~ density * GI + patch
#    (Optionally weight by total_dives to upweight better-observed bouts)
# ------------------------------------------------------------------------------
# Toggle weighting (set weights = NULL to turn off)
use_weights <- FALSE

if (use_weights) {
  m_choice <- glmmTMB(
    urchin_focus_bin ~ purple_urchin_densitym2 * mean_gonad_index_purple_urchin + pred_patch +
      (1 | patch_id),
    family = binomial,
    data   = rf_choice
  )
} else {
  m_choice <- glmmTMB(
    urchin_focus_bin ~ purple_urchin_densitym2 * mean_gonad_index_purple_urchin + pred_patch,
    family = binomial,
    data   = rf_choice
  )
}

summary(m_choice)

# ------------------------------------------------------------------------------
# 8) Manual prediction grids (prey-choice scale) + plots
# ------------------------------------------------------------------------------
linkfun <- family(m_choice)$linkinv

# --- 8a) Vary density (hold GI at observed mean within used patches) ----------
x_density <- seq(
  min(rf_choice$purple_urchin_densitym2, na.rm = TRUE),
  max(rf_choice$purple_urchin_densitym2, na.rm = TRUE),
  length.out = 100
)

new_density <- expand.grid(
  purple_urchin_densitym2 = x_density,
  mean_gonad_index_purple_urchin = mean(rf_choice$mean_gonad_index_purple_urchin, na.rm = TRUE),
  pred_patch = levels(rf_choice$pred_patch)
)

pred_density <- predict(m_choice, newdata = new_density, type = "link", se.fit = TRUE)

df_density <- new_density %>%
  mutate(
    fit_link = pred_density$fit,
    se_link  = pred_density$se.fit,
    fit_prob = linkfun(fit_link),
    lo = linkfun(fit_link - 1.0 * se_link),  # ~68% interval
    hi = linkfun(fit_link + 1.0 * se_link)
  )

# --- 8b) Vary gonad index (hold density at observed mean within used patches) --
x_gonad <- seq(
  min(rf_choice$mean_gonad_index_purple_urchin, na.rm = TRUE),
  max(rf_choice$mean_gonad_index_purple_urchin, na.rm = TRUE),
  length.out = 100
)

new_gonad <- expand.grid(
  purple_urchin_densitym2 = mean(rf_choice$purple_urchin_densitym2, na.rm = TRUE),
  mean_gonad_index_purple_urchin = x_gonad,
  pred_patch = levels(rf_choice$pred_patch)
)

pred_gonad <- predict(m_choice, newdata = new_gonad, type = "link", se.fit = TRUE)

df_gonad <- new_gonad %>%
  mutate(
    fit_link = pred_gonad$fit,
    se_link  = pred_gonad$se.fit,
    fit_prob = linkfun(fit_link),
    lo = linkfun(fit_link - 1.0 * se_link),  # ~68% interval
    hi = linkfun(fit_link + 1.0 * se_link)
  )

# ------------------------------------------------------------------------------
# 9) Plots (probability of urchin-focused foraging)
# ------------------------------------------------------------------------------
g_density <- ggplot(df_density,
                    aes(x = purple_urchin_densitym2, y = fit_prob, color = pred_patch)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = pred_patch), alpha = 0.25, color = NA) +
  scale_color_manual(values = patch_colors, name = "Patch type") +
  scale_fill_manual(values = patch_colors,  name = "Patch type") +
  labs(
    x = expression("Purple urchin density (ind. m"^-2*")"),
    y = "P(urchin-focused foraging | used patch)",
    title = "Prey choice vs. urchin density (within used patches)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

g_gonad <- ggplot(df_gonad,
                  aes(x = mean_gonad_index_purple_urchin, y = fit_prob, color = pred_patch)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = pred_patch), alpha = 0.25, color = NA) +
  scale_color_manual(values = patch_colors, name = "Patch type") +
  scale_fill_manual(values = patch_colors,  name = "Patch type") +
  labs(
    x = "Purple urchin gonad index (mean per patch)",
    y = "P(urchin-focused foraging | used patch)",
    title = "Prey choice vs. urchin gonad condition (within used patches)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

library(patchwork)
g_density + g_gonad

# ------------------------------------------------------------------------------
# 10) Diagnostics & quick checks (optional but recommended)
# ------------------------------------------------------------------------------
# Correlation among predictors within used patches
cor(
  rf_choice$purple_urchin_densitym2,
  rf_choice$mean_gonad_index_purple_urchin,
  use = "complete.obs"
)

# Joint distribution of density vs GI (used patches)
ggplot(rf_choice, aes(purple_urchin_densitym2, mean_gonad_index_purple_urchin, color = pred_patch)) +
  geom_point(alpha = 0.5) +
  geom_density_2d() +
  scale_color_manual(values = patch_colors) +
  theme_bw() +
  labs(x = expression("Density (ind. m"^-2*")"),
       y = "Gonad index",
       title = "Observed combinations of density and GI (used patches)")

# Raw prey-focus vs GI by patch type
ggplot(rf_choice, aes(mean_gonad_index_purple_urchin, urchin_focus_bin, color = pred_patch)) +
  geom_jitter(height = 0.05, alpha = 0.4) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE) +
  scale_color_manual(values = patch_colors) +
  theme_bw() +
  labs(x = "Gonad index", y = "Urchin-focused (0/1)",
       title = "Bivariate prey choice vs. gonad index")

# ------------------------------------------------------------------------------
# 11) Optional export
# ------------------------------------------------------------------------------
# ggsave("output/figures/foraging_choice_glmm_manual_predictions.png",
#        plot = g_density + g_gonad, width = 10, height = 4.5, dpi = 600)
