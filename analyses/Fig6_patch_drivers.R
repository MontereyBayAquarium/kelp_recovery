################################################################################
# Multiclass Sensitivity Analysis:
# What defines Barrens, Incipient Forests, and Forests?
#
# Model:
#   Random Forest classifier of patch state (BAR / INCIP / FOR)
#
# Outputs:
#   1. Variable importance (which predictors matter most overall)
#   2. Per-state importance (what best predicts INCIP vs BAR vs FOR)
#   3. Directionality (does more X increase or decrease P(INCIP)? etc.)
#
# Author: Joshua G. Smith — UCSC Nearshore Ecology Research Group
################################################################################

rm(list = ls())

require(librarian)
librarian::shelf(
  tidyverse, janitor, lubridate, sf, here,
  randomForest, vip, pdp, patchwork
)

################################################################################
# 1. Load data -----------------------------------------------------------------
################################################################################

# quad_build3 must include patch polygons + benthic survey data
load(here::here("output", "survey_data", "processed", "zone_level_data4.rda"))

# foraging observations (dives etc.)
forage_orig <- read_csv(
  "/Volumes/enhydra/data/foraging_data/processed/foraging_data_2024_2025_processed.csv"
)

years_keep   <- c(2024, 2025)
focal_months <- c(6, 7, 8, 9)  # Jun–Sep
patch_colors <- c("BAR"="purple","INCIP"="orange","FOR"="forestgreen")

################################################################################
# 2. Convert purple urchin prey to biomass -------------------------------------
################################################################################

# map size bins to estimated test diameter
size_key <- expand_grid(size = 1:4, qualifier = c("a","b","c")) %>%
  mutate(
    size_class = paste0(size, qualifier),
    size_cm = (size - 1) * 5 + case_when(
      qualifier == "a" ~ 1.66,
      qualifier == "b" ~ 3.32,
      qualifier == "c" ~ 4.98
    )
  )

# keep purple urchin prey in focal months, estimate biomass per prey item
pur_forage <- forage_orig %>%
  filter(
    month %in% focal_months,
    str_detect(prey, regex("^pur$", ignore_case = TRUE))  # purple urchin only
  ) %>%
  left_join(size_key, by = c("size","qualifier")) %>%
  mutate(
    test_diameter_mm = size_cm * 10,
    biomass_g        = -14.2 + 7.44 * exp(0.04 * test_diameter_mm),
    biomass_g        = if_else(biomass_g < 0.5, 0.5, biomass_g),  # floor for sanity
    total_biomass_g  = biomass_g * number
  ) %>%
  filter(size_cm < 8) # discard impossible monsters

################################################################################
# 3. Match foraging records to habitat patches ---------------------------------
################################################################################

# foraging locations as sf
pur_sf <- st_as_sf(
  pur_forage,
  coords = c("long","lat"),
  crs = 4326,
  remove = FALSE
)

# survey polygons in same CRS, with survey metadata
quad_same_crs <- st_transform(quad_build3, st_crs(pur_sf)) %>%
  mutate(
    survey_year = year(survey_date),
    survey_doy  = yday(survey_date)
  )

# spatial join (which patch did the otter forage in?)
pur_spatial <- st_join(
  pur_sf,
  quad_same_crs %>%
    select(patch_id, pred_patch, survey_date, survey_year, survey_doy),
  join = st_intersects,
  left = TRUE
)

# match each dive/bout to the *closest-in-time* survey of that patch
pur_patch <- pur_spatial %>%
  mutate(
    forag_doy = yday(date),
    diff_days = abs(forag_doy - survey_doy)
  ) %>%
  group_by(foragdata_id) %>%
  slice_min(diff_days, with_ties = FALSE) %>%   # best temporal match
  ungroup() %>%
  st_drop_geometry() %>%
  filter(survey_year %in% years_keep)

################################################################################
# 4. Summarize to patch-year and compute *relative* otter foraging effort -------
################################################################################

# collapse to patch_id × year × patch_type
# - total purple urchin biomass consumed
# - how "important" this patch was to otter foraging that year
biomass_patch <- pur_patch %>%
  group_by(survey_year, patch_id, pred_patch) %>%
  summarise(
    total_biomass_g = sum(total_biomass_g, na.rm = TRUE),
    n_foraging_obs  = n(),  # how many feeding obs in that patch-year
    .groups = "drop"
  ) %>%
  mutate(
    foraging_biomass_kg = total_biomass_g / 1000,
    year                = survey_year
  ) %>%
  group_by(year) %>%
  mutate(
    rel_foraging_effort =
      foraging_biomass_kg / sum(foraging_biomass_kg, na.rm = TRUE)
  ) %>%
  ungroup()

# rel_foraging_effort ~ "how much of all otter urchin take this year
# happened in this specific patch?"
# High value = otters are targeting this patch heavily (selection pressure).

################################################################################
# 5. Join with habitat / prey / behavior / physical attributes -----------------
################################################################################

# predictors of interest (benthic survey level)
predictors_focus <- c(
  "mean_gonad_mass_g",          # gonad mass per urchin (food quality)
  "mean_biomass_g",             # urchin biomass per individual (size/condition)
  "purple_urchin_densitym2",    # prey availability
  "purple_urchin_conceiledm2",  # refuge-seeking behavior
  "relief_cm",                  # rugosity / structure
  "risk_index"                  # physical risk / exposure
)

patch_predictors <- quad_same_crs %>%
  mutate(year = year(survey_date)) %>%
  st_drop_geometry() %>%
  group_by(patch_id, year, pred_patch) %>%
  summarise(
    across(all_of(predictors_focus), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# merge foraging summary with benthic predictors
model_data <- biomass_patch %>%
  left_join(
    patch_predictors,
    by = c("patch_id","year","pred_patch")
  ) %>%
  mutate(
    pred_patch = factor(pred_patch, levels = c("BAR","INCIP","FOR")),
    behavior_ratio = purple_urchin_conceiledm2 / purple_urchin_densitym2
  ) %>%
  filter(
    !is.na(rel_foraging_effort),
    !is.na(mean_gonad_mass_g),
    !is.na(mean_biomass_g),
    !is.na(purple_urchin_densitym2),
    !is.na(purple_urchin_conceiledm2),
    !is.na(relief_cm),
    !is.na(risk_index),
    is.finite(behavior_ratio) # excludes div-by-zero weirdness
  )

# rf_dat is now "one row per patch-year" with all predictors + class label
rf_dat <- model_data %>%
  select(
    pred_patch,
    rel_foraging_effort,
    mean_biomass_g,
    mean_gonad_mass_g,
    behavior_ratio,
    relief_cm,
    risk_index
  )

################################################################################
# 6. Fit a single multiclass Random Forest -------------------------------------
################################################################################

set.seed(42)
rf_multi <- randomForest(
  pred_patch ~ rel_foraging_effort +
    mean_biomass_g +
    mean_gonad_mass_g +
    behavior_ratio +
    relief_cm +
    risk_index,
  data = rf_dat,
  ntree = 2000,
  mtry  = 3,
  importance = TRUE,
  probability = TRUE
)

print(rf_multi)

# Notes:
# - rf_multi$confusion gives class-by-class accuracy (BAR vs INCIP vs FOR)
# - rf_multi$votes gives predicted probability for each class

################################################################################
# 7. Variable importance (overall and per-state) ------------------------------
################################################################################

# Overall permutation importance (MeanDecreaseAccuracy)
varImpPlot(
  rf_multi,
  type = 1,
  main = "Predictors of Patch State (overall importance)"
)

# Publication-style variable importance per class (BAR, INCIP, FOR)
# vip::vip(..., target="INCIP") asks:
#   "Which predictors most affect P(pred_patch == INCIP) ?"
vip_incip <- vip(
  rf_multi,
  target = "INCIP",
  geom = "col",
  aesthetics = list(fill = "orange", alpha = 0.8)
) +
  labs(
    title    = "Predictors of Incipient Forest",
    subtitle = "Permutation importance on P(incipient)",
    x = "Importance",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank())

vip_bar <- vip(
  rf_multi,
  target = "BAR",
  geom = "col",
  aesthetics = list(fill = "purple", alpha = 0.8)
) +
  labs(
    title    = "Predictors of Barren",
    subtitle = "Permutation importance on P(barren)",
    x = "Importance",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank())

vip_for <- vip(
  rf_multi,
  target = "FOR",
  geom = "col",
  aesthetics = list(fill = "forestgreen", alpha = 0.8)
) +
  labs(
    title    = "Predictors of Forest",
    subtitle = "Permutation importance on P(forest)",
    x = "Importance",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank())

# You can view them individually:
# print(vip_incip); print(vip_bar); print(vip_for)

# Or combine in a 1x3 panel:
vip_panel <- vip_bar + vip_incip + vip_for +
  plot_layout(ncol = 3)

print(vip_panel)

################################################################################
# 8. Directionality via partial dependence ------------------------------------
################################################################################
# PDP tells you:
#  as predictor X increases, does P(class=k) go up or down?

# helper to generate PDP for one predictor/one class
make_pdp <- function(var, class_label, line_col) {
  pd <- partial(
    rf_multi,
    pred.var    = var,
    which.class = class_label,
    prob        = TRUE,
    grid.resolution = 100
  )
  
  pretty_x <- c(
    rel_foraging_effort = "Relative otter urchin foraging effort\n(proportion of annual take)",
    mean_biomass_g      = "Urchin biomass (g / urchin)",
    mean_gonad_mass_g   = "Urchin gonad mass (g)",
    behavior_ratio      = "Urchin concealment ratio\n(concealed / total)",
    relief_cm           = "Reef relief (cm)",
    risk_index          = "Risk index"
  )
  
  ggplot(pd, aes_string(x = var, y = "yhat")) +
    geom_line(linewidth = 1.2, color = line_col) +
    labs(
      x = pretty_x[[var]],
      y = paste0("P(", class_label, ")"),
      title = paste0(class_label, ": effect of ", var)
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text  = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = 10)
    )
}


################################################################################
# 5. Generate PDPs for each patch type -----------------------------------------
################################################################################
################################################################################
# Multi-panel PDP Figure — Directionality of Key Predictors Across Patch States
# Clean version with "Probability" y-axis and right-side row labels
################################################################################

################################################################################
# Multi-panel PDP Figure — Directionality of Key Predictors Across Patch States
# Clean, publication-style version (no titles, just axes + tags)
################################################################################

require(patchwork)

# --- Generate PDPs for each state ---------------------------------------------

# --- Barren ---
p_bar_behavior <- make_pdp("behavior_ratio",      "BAR", "purple") +
  labs(title = NULL)
p_bar_biomass  <- make_pdp("mean_biomass_g",     "BAR", "purple") +
  labs(title = NULL)
p_bar_foraging <- make_pdp("rel_foraging_effort","BAR", "purple") +
  labs(title = NULL)

# --- Incipient ---
p_incip_behavior <- make_pdp("behavior_ratio",      "INCIP", "orange") +
  labs(title = NULL)
p_incip_biomass  <- make_pdp("mean_biomass_g",     "INCIP", "orange") +
  labs(title = NULL)
p_incip_foraging <- make_pdp("rel_foraging_effort","INCIP", "orange") +
  labs(title = NULL)

# --- Forest ---
p_for_behavior <- make_pdp("behavior_ratio",      "FOR", "forestgreen") +
  labs(title = NULL)
p_for_biomass  <- make_pdp("mean_biomass_g",     "FOR", "forestgreen") +
  labs(title = NULL)
p_for_foraging <- make_pdp("rel_foraging_effort","FOR", "forestgreen") +
  labs(title = NULL)


library(grid)
library(patchwork)
library(cowplot)

library(grid)
library(patchwork)
library(cowplot)

# --- Base grid with tags -----------------------------------------------------
pdp_grid_clean <- (
  (p_bar_behavior | p_bar_biomass | p_bar_foraging) /
    (p_incip_behavior | p_incip_biomass | p_incip_foraging) /
    (p_for_behavior | p_for_biomass | p_for_foraging)
) +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.title   = element_blank(),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_blank(),  # remove individual y labels
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    plot.margin  = unit(c(4, 4, 4, 4), "mm"),
    panel.grid   = element_blank()
  )

# --- Draw figure with slightly tighter spacing -------------------------------
pdp_grid_labeled <- ggdraw() +
  # bring panels closer to left margin (was x = 0.045)
  draw_plot(pdp_grid_clean, x = 0.035, y = 0, width = 0.92, height = 1) +
  
  # shared y-axis label, snug to tick labels
  draw_label("Probability",
             x = 0.022, y = 0.5,
             angle = 90, vjust = 0.5, size = 11) +
  
  # right-side labels (moved up slightly and inward)
  draw_label("Barren",
             x = 0.952, y = 0.87,
             angle = 270, fontface = "bold", size = 11) +
  draw_label("Incipient",
             x = 0.952, y = 0.54,
             angle = 270, fontface = "bold", size = 11) +
  draw_label("Forest",
             x = 0.952, y = 0.21,
             angle = 270, fontface = "bold", size = 11)

print(pdp_grid_labeled)
