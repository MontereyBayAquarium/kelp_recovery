################################################################################
# Multiclass Sensitivity Analysis:
# What defines Barrens, Incipient Forests, and Forests?
#
# Model:
#   Random Forest classifier of patch state (BAR / INCIP / FOR)
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

load(here::here("output", "survey_data", "processed", "zone_level_data4.rda"))

forage_orig <- read_csv(
  "/Volumes/enhydra/data/foraging_data/processed/foraging_data_2024_2025_processed.csv"
)

years_keep   <- c(2024, 2025)
focal_months <- c(6, 7, 8, 9)  # Jun–Sep
patch_colors <- c("BAR"="purple","INCIP"="orange","FOR"="forestgreen")

################################################################################
# 2. Convert purple urchin prey to biomass -------------------------------------
################################################################################

size_key <- expand_grid(size = 1:4, qualifier = c("a","b","c")) %>%
  mutate(
    size_class = paste0(size, qualifier),
    size_cm = (size - 1) * 5 + case_when(
      qualifier == "a" ~ 1.66,
      qualifier == "b" ~ 3.32,
      qualifier == "c" ~ 4.98
    )
  )

pur_forage <- forage_orig %>%
  filter(
    month %in% focal_months,
    str_detect(prey, regex("^pur$", ignore_case = TRUE))  # purple urchin only
  ) %>%
  left_join(size_key, by = c("size","qualifier")) %>%
  mutate(
    test_diameter_mm = size_cm * 10,
    biomass_g        = -14.2 + 7.44 * exp(0.04 * test_diameter_mm),
    biomass_g        = if_else(biomass_g < 0.5, 0.5, biomass_g),
    total_biomass_g  = biomass_g * number
  ) %>%
  filter(size_cm < 8) # discard impossible monsters >8 cm

################################################################################
# 3. Match foraging records to habitat patches ---------------------------------
################################################################################

pur_sf <- st_as_sf(
  pur_forage,
  coords = c("long","lat"),
  crs = 4326,
  remove = FALSE
)

quad_same_crs <- st_transform(quad_build3, st_crs(pur_sf)) %>%
  mutate(
    survey_year = year(survey_date),
    survey_doy  = yday(survey_date)
  )

pur_spatial <- st_join(
  pur_sf,
  quad_same_crs %>%
    select(patch_id, pred_patch, survey_date, survey_year, survey_doy),
  join = st_intersects,
  left = TRUE
)

pur_patch <- pur_spatial %>%
  mutate(
    forag_doy = yday(date),
    diff_days = abs(forag_doy - survey_doy)
  ) %>%
  group_by(foragdata_id) %>%
  slice_min(diff_days, with_ties = FALSE) %>%
  ungroup() %>%
  st_drop_geometry() %>%
  filter(survey_year %in% years_keep)

################################################################################
# 4. Summarize to patch-year and compute foraging metrics ----------------------
################################################################################

biomass_patch <- pur_patch %>%
  group_by(survey_year, patch_id, pred_patch) %>%
  summarise(
    total_biomass_g = sum(total_biomass_g, na.rm = TRUE),
    n_foraging_obs  = n(),
    .groups = "drop"
  ) %>%
  mutate(
    foraging_biomass_kg = total_biomass_g / 1000,
    biomass_density_gm2 = total_biomass_g / 80,
    year                = survey_year
  ) %>%
  group_by(year) %>%
  mutate(
    rel_foraging_effort = foraging_biomass_kg / sum(foraging_biomass_kg, na.rm = TRUE)
  ) %>%
  ungroup()

################################################################################
# 5. Join with habitat / prey / physical attributes ----------------------------
################################################################################

predictors_focus <- c(
  "mean_gonad_mass_g",
  "total_gonad_mass_g",         # for gonad mass per m2
  "mean_biomass_g",
  "purple_urchin_densitym2",
  "purple_urchin_conceiledm2",
  "relief_cm",
  "risk_index"
)

patch_predictors <- quad_same_crs %>%
  mutate(year = year(survey_date)) %>%
  st_drop_geometry() %>%
  group_by(patch_id, year, pred_patch) %>%
  summarise(across(all_of(predictors_focus), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(gonad_mass_gm2 = total_gonad_mass_g / 80)   # NEW predictor

biomass_patch_full <- patch_predictors %>%
  select(patch_id, year, pred_patch) %>%
  distinct() %>%
  left_join(
    biomass_patch %>%
      select(patch_id, year, rel_foraging_effort, biomass_density_gm2),
    by = c("patch_id", "year")
  ) %>%
  mutate(
    rel_foraging_effort  = replace_na(rel_foraging_effort, 0),
    biomass_density_gm2  = replace_na(biomass_density_gm2, 0)
  )

model_data <- biomass_patch_full %>%
  left_join(patch_predictors,
            by = c("patch_id", "year", "pred_patch")) %>%
  mutate(
    pred_patch = factor(pred_patch, levels = c("BAR", "INCIP", "FOR")),
    behavior_ratio = purple_urchin_conceiledm2 / purple_urchin_densitym2
  ) %>%
  filter(
    !is.na(mean_gonad_mass_g),
    !is.na(mean_biomass_g),
    !is.na(gonad_mass_gm2),
    !is.na(purple_urchin_densitym2),
    !is.na(purple_urchin_conceiledm2),
    !is.na(relief_cm),
    !is.na(risk_index),
    is.finite(behavior_ratio)
  )

################################################################################
# 6. Fit multiclass Random Forest ----------------------------------------------
################################################################################

rf_dat <- model_data %>%
  select(
    pred_patch,
    rel_foraging_effort,
    mean_biomass_g,
    mean_gonad_mass_g,
    gonad_mass_gm2,          # NEW predictor
    biomass_density_gm2,
    behavior_ratio,
    relief_cm,
    risk_index
  )

set.seed(42)
rf_multi <- randomForest(
  pred_patch ~ rel_foraging_effort +
    mean_biomass_g +
    mean_gonad_mass_g +
    gonad_mass_gm2 +          # NEW predictor
    biomass_density_gm2 +
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

################################################################################
# 7. Variable importance (overall and per-state) -------------------------------
################################################################################

varImpPlot(
  rf_multi,
  type = 1,
  main = "Predictors of Patch State (overall importance)"
)

################################################################################
# 8. Partial dependence helper -------------------------------------------------
################################################################################

make_pdp <- function(var, class_label, line_col) {
  pd <- partial(
    rf_multi,
    pred.var    = var,
    which.class = class_label,
    prob        = TRUE,
    grid.resolution = 100
  )
  
  pretty_x <- c(
    behavior_ratio        = "Urchin concealment ratio\n(concealed / total)",
    mean_biomass_g        = "Urchin biomass (g / urchin)",
    biomass_density_gm2   = "Urchin biomass density (g / m²)",
    rel_foraging_effort   = "Relative otter urchin foraging effort\n(proportion of annual take)",
    mean_gonad_mass_g     = "Urchin gonad mass (g)",
    gonad_mass_gm2        = "Gonad mass density (g / m²)",   # NEW label
    relief_cm             = "Reef relief (cm)",
    risk_index            = "Risk index"
  )
  
  ggplot(pd, aes_string(x = var, y = "yhat")) +
    geom_line(linewidth = 1.2, color = line_col) +
    labs(
      x = pretty_x[[var]],
      y = paste0("P(", class_label, ")")
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text  = element_text(color = "black")
    )
}

################################################################################
# 9. Generate PDPs for each patch type -----------------------------------------
################################################################################

# BAR
p_bar_behavior  <- make_pdp("behavior_ratio",        "BAR",   "purple")
p_bar_biomass   <- make_pdp("mean_biomass_g",       "BAR",   "purple")
p_bar_biomassD  <- make_pdp("biomass_density_gm2",  "BAR",   "purple")
p_bar_foraging  <- make_pdp("rel_foraging_effort",  "BAR",   "purple")
p_bar_gonadm2   <- make_pdp("gonad_mass_gm2",       "BAR",   "purple")  # NEW PDP

# INCIP
p_incip_behavior <- make_pdp("behavior_ratio",        "INCIP", "orange")
p_incip_biomass  <- make_pdp("mean_biomass_g",       "INCIP", "orange")
p_incip_biomassD <- make_pdp("biomass_density_gm2",  "INCIP", "orange")
p_incip_foraging <- make_pdp("rel_foraging_effort",  "INCIP", "orange")
p_incip_gonadm2  <- make_pdp("gonad_mass_gm2",       "INCIP", "orange") # NEW PDP

# FOR
p_for_behavior  <- make_pdp("behavior_ratio",        "FOR",   "forestgreen")
p_for_biomass   <- make_pdp("mean_biomass_g",       "FOR",   "forestgreen")
p_for_biomassD  <- make_pdp("biomass_density_gm2",  "FOR",   "forestgreen")
p_for_foraging  <- make_pdp("rel_foraging_effort",  "FOR",   "forestgreen")
p_for_gonadm2   <- make_pdp("gonad_mass_gm2",       "FOR",   "forestgreen") # NEW PDP

################################################################################
# 10. Stack PDPs for facet_grid ------------------------------------------------
################################################################################

extract_pdp <- function(p, predictor_label, state_label) {
  df <- p$data
  names(df) <- c("x", "y")
  df %>%
    mutate(
      state     = state_label,
      predictor = predictor_label
    )
}

pdp_all <- bind_rows(
  # Barren
  extract_pdp(p_bar_behavior,  "Behavior ratio",                     "Barren"),
  extract_pdp(p_bar_biomass,   "Mean biomass (g)",                   "Barren"),
  extract_pdp(p_bar_biomassD,  "Biomass density (g / m²)",           "Barren"),
  extract_pdp(p_bar_foraging,  "Relative foraging effort",           "Barren"),
  extract_pdp(p_bar_gonadm2,   "Gonad mass density (g / m²)",        "Barren"), # NEW
  
  # Incipient
  extract_pdp(p_incip_behavior,"Behavior ratio",                     "Incipient"),
  extract_pdp(p_incip_biomass, "Mean biomass (g)",                   "Incipient"),
  extract_pdp(p_incip_biomassD,"Biomass density (g / m²)",           "Incipient"),
  extract_pdp(p_incip_foraging,"Relative foraging effort",           "Incipient"),
  extract_pdp(p_incip_gonadm2, "Gonad mass density (g / m²)",        "Incipient"), # NEW
  
  # Forest
  extract_pdp(p_for_behavior,  "Behavior ratio",                     "Forest"),
  extract_pdp(p_for_biomass,   "Mean biomass (g)",                   "Forest"),
  extract_pdp(p_for_biomassD,  "Biomass density (g / m²)",           "Forest"),
  extract_pdp(p_for_foraging,  "Relative foraging effort",           "Forest"),
  extract_pdp(p_for_gonadm2,   "Gonad mass density (g / m²)",        "Forest") # NEW
) %>%
  drop_na(x, y)

predictor_labels <- c(
  "Behavior ratio"             = "Urchin concealment \nratio (concealed / total)",
  "Mean biomass (g)"           = "Urchin biomass \n(g / urchin)",
  "Biomass density (g / m²)"   = "Urchin biomass \ndensity (g / m²)",
  "Relative foraging effort"   = "Sea otter \nforaging effort \n(proportion take in patch)",
  "Gonad mass density (g / m²)"= "Urchin gonad mass \ndensity (g / m²)" # NEW label
)

ggplot(pdp_all, aes(x = x, y = y, color = state)) +
  geom_line(linewidth = 1) +
  facet_wrap(
    ~ predictor,
    labeller = labeller(predictor = predictor_labels),
    scales = "free",
    strip.position = "bottom",
    nrow=1
  ) +
  scale_color_manual(values = c(
    "Barren"    = "purple",
    "Incipient" = "orange",
    "Forest"    = "forestgreen"
  )) +
  labs(
    y = "Probability",
    x = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background   = element_blank(),
    strip.placement    = "outside",
    strip.text.x       = element_text(size = 10, face = "plain"),
    axis.text          = element_text(size = 8),
    axis.title.y       = element_text(size = 10),
    legend.position    = "none",
    panel.spacing      = unit(1, "lines"),
    axis.title.x       = element_blank()
  )

