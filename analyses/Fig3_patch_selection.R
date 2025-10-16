# ==============================================================================
# Sea otter foraging: focal-bout GLMM + manual prediction plots
# Focus: density × gonad index + patch type
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
dissection_data <- read_csv("/Volumes/enhydra/data/kelp_recovery/MBA_kelp_forest_database/processed/dissection/dissection_data_cleanedv2.csv")

years_keep   <- c(2024, 2025)
patch_colors <- c("BAR"="purple", "INCIP"="orange", "FOR"="forestgreen")

# ------------------------------------------------------------------------------
# 2) Dissection summaries → join to patch data (quad_build4)
# ------------------------------------------------------------------------------
dissect_build1 <- dissection_data %>%
  mutate(
    species = str_to_lower(str_trim(species)),
    species = case_when(
      str_detect(species, "red") ~ "red_urchin",
      str_detect(species, "pur") ~ "purple_urchin",
      TRUE ~ species
    ),
    year = year(date_collected)
  ) %>%
  group_by(year, site_number, site_type, zone, species) %>%
  summarise(
    mean_gonad_mass_g = mean(gonad_mass_g, na.rm = TRUE),
    mean_gonad_index  = mean(gonad_index,  na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = species,
    values_from = c(mean_gonad_mass_g, mean_gonad_index, n),
    names_sep   = "_"
  )

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
# 4) Foraging bouts → centroids → join to polygons → rf_bout
# ------------------------------------------------------------------------------
# Note: forage_orig is expected to have columns `long` and `lat`
forage_sf <- forage_orig %>%
  st_as_sf(coords = c("long","lat"), crs = 4326, remove = FALSE)

focal_bouts <- forage_sf %>%
  group_by(bout) %>%
  summarise(
    year       = first(year),
    n_dives    = n(),
    prey_types = paste(unique(prey), collapse = ","),
    geometry   = st_centroid(st_union(geometry)),
    .groups    = "drop"
  ) %>%
  mutate(
    focal3 = if_else(n_dives > 3 & grepl("urc|red|pur", prey_types, ignore.case = TRUE), "yes", "no"),
    focal5 = if_else(n_dives > 5, "yes", "no")
  )

# Ensure same CRS before spatial join
quad4_same_crs <- st_transform(quad_build4, st_crs(focal_bouts))

bouts_with_patch <- st_join(
  focal_bouts,
  quad4_same_crs %>% dplyr::select(patch_id, pred_patch, survey_date),
  join = st_intersects,
  left = TRUE
) %>%
  mutate(patch_year = lubridate::year(survey_date)) %>%
  filter(year %in% years_keep, !is.na(patch_id), patch_year == year) %>%
  st_drop_geometry() %>%
  distinct(bout, .keep_all = TRUE)

rf_bout <- bouts_with_patch %>%
  inner_join(quad_year, by = c("patch_id","year","pred_patch")) %>%
  mutate(
    focal3_bin = factor(focal3, levels = c("no","yes")),
    focal5_bin = factor(focal5, levels = c("no","yes"))
  )

# ------------------------------------------------------------------------------
# 5) GLMM data prep (unscaled predictors)
# ------------------------------------------------------------------------------
rf_bout_glmm <- rf_bout %>%
  drop_na(purple_urchin_densitym2, mean_gonad_index_purple_urchin, pred_patch) %>%
  mutate(
    focal3_num = if_else(focal3 == "yes", 1, 0),
    pred_patch = factor(pred_patch, levels = c("BAR", "INCIP", "FOR"))
  )

# ------------------------------------------------------------------------------
# 6) Fit GLMM (unscaled) and manual prediction plots
# ------------------------------------------------------------------------------
m_glmm_raw <- glmmTMB(
  focal3_num ~ purple_urchin_densitym2 * mean_gonad_index_purple_urchin + pred_patch,
  family = binomial,
  data = rf_bout_glmm
)

# Inverse-logit
linkfun <- family(m_glmm_raw)$linkinv

# --- 6a) Vary density (hold GI at mean) ---------------------------------------
x_density <- seq(
  min(rf_bout_glmm$purple_urchin_densitym2, na.rm = TRUE),
  max(rf_bout_glmm$purple_urchin_densitym2, na.rm = TRUE),
  length.out = 100
)

new_density <- expand.grid(
  purple_urchin_densitym2 = x_density,
  mean_gonad_index_purple_urchin = mean(rf_bout_glmm$mean_gonad_index_purple_urchin, na.rm = TRUE),
  pred_patch = levels(rf_bout_glmm$pred_patch)
)

pred_density <- predict(m_glmm_raw, newdata = new_density, type = "link", se.fit = TRUE)

df_density <- new_density %>%
  mutate(
    fit_link = pred_density$fit,
    se_link  = pred_density$se.fit,
    fit_prob = linkfun(fit_link),
    lo = linkfun(fit_link - 1.0 * se_link),  # ~68% interval
    hi = linkfun(fit_link + 1.0 * se_link)
  )

# --- 6b) Vary gonad index (hold density at mean) ------------------------------
x_gonad <- seq(
  min(rf_bout_glmm$mean_gonad_index_purple_urchin, na.rm = TRUE),
  max(rf_bout_glmm$mean_gonad_index_purple_urchin, na.rm = TRUE),
  length.out = 100
)

new_gonad <- expand.grid(
  purple_urchin_densitym2 = mean(rf_bout_glmm$purple_urchin_densitym2, na.rm = TRUE),
  mean_gonad_index_purple_urchin = x_gonad,
  pred_patch = levels(rf_bout_glmm$pred_patch)
)

pred_gonad <- predict(m_glmm_raw, newdata = new_gonad, type = "link", se.fit = TRUE)

df_gonad <- new_gonad %>%
  mutate(
    fit_link = pred_gonad$fit,
    se_link  = pred_gonad$se.fit,
    fit_prob = linkfun(fit_link),
    lo = linkfun(fit_link - 1.0 * se_link),  # ~68% interval
    hi = linkfun(fit_link + 1.0 * se_link)
  )

# ------------------------------------------------------------------------------
# 7) Plots (your preferred style)
# ------------------------------------------------------------------------------
# Density effect
g_density <- ggplot(df_density,
                    aes(x = purple_urchin_densitym2, y = fit_prob, color = pred_patch)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = pred_patch), alpha = 0.25, color = NA) +
  scale_color_manual(values = patch_colors, name = "Patch type") +
  scale_fill_manual(values = patch_colors,  name = "Patch type") +
  labs(
    x = expression("Purple urchin density (ind. m"^-2*")"),
    y = "Predicted probability of focal foraging",
    title = "Effect of urchin density across patch types"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

# Gonad index effect
g_gonad <- ggplot(df_gonad,
                  aes(x = mean_gonad_index_purple_urchin, y = fit_prob, color = pred_patch)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = pred_patch), alpha = 0.25, color = NA) +
  scale_color_manual(values = patch_colors, name = "Patch type") +
  scale_fill_manual(values = patch_colors,  name = "Patch type") +
  labs(
    x = "Purple urchin gonad index (mean per patch)",
    y = "Predicted probability of focal foraging",
    title = "Effect of urchin gonad condition across patch types"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

# Side-by-side
library(patchwork)
g_density + g_gonad

# ------------------------------------------------------------------------------
# Optional export
# ------------------------------------------------------------------------------
# ggsave("output/figures/foraging_glmm_manual_predictions.png",
#        plot = g_density + g_gonad, width = 10, height = 4.5, dpi = 600)
