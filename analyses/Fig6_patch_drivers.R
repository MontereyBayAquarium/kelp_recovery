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


# ---------------------------------------------------------------------------
# Variable importance (MeanDecreaseGini) to annotate each predictor
# ---------------------------------------------------------------------------
var_imp <- importance(rf_multi, type = 2) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename(MeanDecreaseGini = MeanDecreaseGini)

# map RF variable names -> the human-readable predictor names we used in pdp_all
var_label_map <- tibble(
  variable = c(
    "behavior_ratio",
    "mean_biomass_g",
    "mean_gonad_mass_g",
    "biomass_density_gm2",
    "rel_foraging_effort"
  ),
  predictor = c(
    "Behavior ratio",
    "Mean biomass (g)",
    "Gonad mass (g / urchin)",
    "Biomass density (g / m²)",
    "Relative foraging effort"
  )
)

# join to get Gini per pretty predictor label
var_imp_pretty <- var_imp %>%
  inner_join(var_label_map, by = "variable") %>%
  mutate(
    gini_note = paste0("Gini = ", round(MeanDecreaseGini, 2))
  ) %>%
  select(predictor, gini_note)

# ---------------------------------------------------------------------------
# Set the predictor facet order MANUALLY the way you want:
# behavior → urchin biomass → gonad mass → biomass density → foraging effort
# ---------------------------------------------------------------------------
desired_order <- c(
  "Behavior ratio",
  "Mean biomass (g)",
  "Gonad mass (g / urchin)",
  "Biomass density (g / m²)",
  "Relative foraging effort"
)

# also set nice multiline facet labels for each predictor
predictor_labels <- c(
  "Behavior ratio"              = "Urchin concealment \nratio (concealed / total)",
  "Mean biomass (g)"            = "Urchin biomass \n(g / urchin)",
  "Gonad mass (g / urchin)"     = "Urchin gonad mass \n(g / urchin)",
  "Biomass density (g / m²)"    = "Urchin biomass \ndensity (g / m²)",
  "Relative foraging effort"    = "Sea otter \nforaging effort \n(proportion take in patch)"
)

# ---------------------------------------------------------------------------
# Make sure pdp_all has the variables that match these labels
# NOTE: we need to rebuild pdp_all with "Gonad mass (g / urchin)" instead of
# "Gonad mass density (g / m²)" since you want gonad mass per urchin
# in slot #3 of the order.
# We'll regenerate pdp_all below to guarantee consistency.
# ---------------------------------------------------------------------------

extract_pdp_df <- function(var, class_label, line_col, pretty_name) {
  pd_tmp <- partial(
    rf_multi,
    pred.var    = var,
    which.class = class_label,
    prob        = TRUE,
    grid.resolution = 100
  )
  
  # pd_tmp has columns var and yhat.
  # Standardize to x/y plus state/predictor for stacking.
  out <- pd_tmp %>%
    rename(x = !!sym(var), y = yhat) %>%
    mutate(
      state     = case_when(
        class_label == "BAR"   ~ "Barren",
        class_label == "INCIP" ~ "Incipient",
        class_label == "FOR"   ~ "Forest",
        TRUE ~ class_label
      ),
      predictor = pretty_name,
      color     = line_col
    )
  out
}

# Build a long df for all predictors × states in the order we care about
pdp_all <- bind_rows(
  # Behavior ratio
  extract_pdp_df("behavior_ratio",       "BAR",   patch_colors["BAR"],   "Behavior ratio"),
  extract_pdp_df("behavior_ratio",       "INCIP", patch_colors["INCIP"], "Behavior ratio"),
  extract_pdp_df("behavior_ratio",       "FOR",   patch_colors["FOR"],   "Behavior ratio"),
  
  # Mean biomass per urchin
  extract_pdp_df("mean_biomass_g",       "BAR",   patch_colors["BAR"],   "Mean biomass (g)"),
  extract_pdp_df("mean_biomass_g",       "INCIP", patch_colors["INCIP"], "Mean biomass (g)"),
  extract_pdp_df("mean_biomass_g",       "FOR",   patch_colors["FOR"],   "Mean biomass (g)"),
  
  # Gonad mass per urchin
  extract_pdp_df("mean_gonad_mass_g",    "BAR",   patch_colors["BAR"],   "Gonad mass (g / urchin)"),
  extract_pdp_df("mean_gonad_mass_g",    "INCIP", patch_colors["INCIP"], "Gonad mass (g / urchin)"),
  extract_pdp_df("mean_gonad_mass_g",    "FOR",   patch_colors["FOR"],   "Gonad mass (g / urchin)"),
  
  # Biomass density m^-2
  extract_pdp_df("biomass_density_gm2",  "BAR",   patch_colors["BAR"],   "Biomass density (g / m²)"),
  extract_pdp_df("biomass_density_gm2",  "INCIP", patch_colors["INCIP"], "Biomass density (g / m²)"),
  extract_pdp_df("biomass_density_gm2",  "FOR",   patch_colors["FOR"],   "Biomass density (g / m²)"),
  
  # Relative otter foraging effort
  extract_pdp_df("rel_foraging_effort",  "BAR",   patch_colors["BAR"],   "Relative foraging effort"),
  extract_pdp_df("rel_foraging_effort",  "INCIP", patch_colors["INCIP"], "Relative foraging effort"),
  extract_pdp_df("rel_foraging_effort",  "FOR",   patch_colors["FOR"],   "Relative foraging effort")
) %>%
  drop_na(x, y)

# apply manual facet order
pdp_all <- pdp_all %>%
  mutate(
    predictor = factor(predictor, levels = desired_order),
    state     = factor(state, levels = c("Barren","Incipient","Forest"))
  ) 

# prep panel-wise Gini notes, in that same order
note_labels <- var_imp_pretty %>%
  filter(predictor %in% desired_order) %>%
  mutate(predictor = factor(predictor, levels = desired_order))

# ---------------------------------------------------------------------------
# Final multi-panel PDP plot
# ---------------------------------------------------------------------------

p_final <- ggplot(pdp_all, aes(x = x, y = y, color = state)) +
  geom_line(linewidth = 1.5) +
  facet_wrap(
    ~ predictor,
    labeller = labeller(predictor = predictor_labels),
    scales   = "free",
    strip.position = "bottom",
    nrow = 1
  ) +
  # add MeanDecreaseGini text to each facet, upper right
  geom_text(
    data = note_labels,
    aes(x = Inf, y = Inf, label = gini_note),
    hjust = 1.05, vjust = 1.5,
    size = 2.8,
    color = "gray20",
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Barren"    = unname(patch_colors["BAR"]),
      "Incipient" = unname(patch_colors["INCIP"]),
      "Forest"    = unname(patch_colors["FOR"])
    ),
    breaks = c("Barren","Incipient","Forest"),
    name   = "Patch state"
  )+
  labs(
    y = "Probability",
    x = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background   = element_blank(),
    strip.placement    = "outside",
    strip.text.x       = element_text(size = 10),
    axis.text          = element_text(size = 8),
    axis.title.y       = element_text(size = 10),
    legend.position    = "top",
    panel.spacing      = unit(1, "lines"),
    axis.title.x       = element_blank()
  ) 

# Show it
p_final



ggsave(
  filename = here::here("figures", "Fig6_incipient_correlates.png"),
  plot = p_final,
  width = 9,        # in inches
  height = 4,        # adjust as needed
  dpi = 600,         # high-res for publication
 bg = "white"       # ensures white background if saving to PNG
)

