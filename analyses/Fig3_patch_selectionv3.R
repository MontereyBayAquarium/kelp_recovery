################################################################################
# Sea otter foraging vs. patch type (absolute foraging biomass + smooth PDPs)
# Author: Joshua G. Smith — UCSC Nearshore Ecology Research Group
# Focus: Predict patch type (BAR / INCIP / FOR) from purple urchin density,
#        mean gonad mass, and absolute urchin biomass consumed by sea otters.
################################################################################

rm(list = ls())

require(librarian)
librarian::shelf(
  tidyverse, janitor, lubridate, sf, here,
  randomForest, pdp, patchwork
)

################################################################################
#Load data

load(here::here("output", "survey_data", "processed", "zone_level_data4.rda"))  
forage_orig <- read_csv("/Volumes/enhydra/data/foraging_data/processed/foraging_data_2024_2025_processed.csv")

years_keep   <- c(2024, 2025)
focal_months <- c(6, 7, 8, 9)  # June–September
patch_colors <- c("BAR"="purple", "INCIP"="orange", "FOR"="forestgreen")

################################################################################
#Load data
# 2. Convert purple urchin prey to biomass in foraging data

# Helper table converting categorical size to numeric diameter (cm)
size_key <- expand_grid(size = 1:4, qualifier = c("a","b","c")) %>%
  mutate(
    size_class = paste0(size, qualifier),
    size_cm = (size - 1) * 5 + case_when(
      qualifier == "a" ~ 1.66,
      qualifier == "b" ~ 3.32,
      qualifier == "c" ~ 4.98
    )
  )

# Filter for purple urchin prey, focal months, convert to biomass
pur_forage <- forage_orig %>%
  filter(
    month %in% focal_months,
    str_detect(prey, regex("^pur$", ignore_case = TRUE))
  ) %>%
  left_join(size_key, by = c("size", "qualifier")) %>%
  mutate(
    test_diameter_mm = size_cm * 10,
    biomass_g = -14.2 + 7.44 * exp(0.04 * test_diameter_mm),
    biomass_g = if_else(biomass_g < 0.5, 0.5, biomass_g),
    total_biomass_g = biomass_g * number
  ) %>%
  filter(size_cm < 8)   # realistic prey cutoff (<8 cm)

################################################################################
#Join with benthic patches

pur_sf <- pur_forage %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326, remove = FALSE)

quad4_same_crs <- st_transform(quad_build3, st_crs(pur_sf))

#pur_patch <- st_join(
#  pur_sf,
#  quad4_same_crs %>% select(patch_id, pred_patch, survey_date),
#  join = st_intersects
#) %>%
#  mutate(year = year(survey_date)) %>%
#  filter(year %in% years_keep) %>%
#  st_drop_geometry()

# Create one polygon per patch_id–year–pred_patch
# Create a unique patch-year layer (keep geometry)
# Ensure patch polygons are unique by patch_id and year



# Ensure patch polygons have numeric survey dates
quad4_same_crs <- quad4_same_crs %>%
  mutate(
    survey_year = year(survey_date),
    survey_doy  = yday(survey_date)
  )

# Convert foraging data to sf (already done)
pur_sf <- st_as_sf(pur_forage, coords = c("long", "lat"), crs = 4326, remove = FALSE)

# Spatial join to all polygons (no year filtering yet)
pur_spatial <- st_join(pur_sf, quad4_same_crs %>%
                         select(patch_id, pred_patch, survey_date, survey_year, survey_doy),
                       join = st_intersects,
                       left = TRUE)

# Compute difference in days between foraging date and survey date
pur_patch <- pur_spatial %>%
  mutate(
    forag_doy = yday(date),
    diff_days = abs(forag_doy - survey_doy)
  ) %>%
  group_by(foragdata_id) %>%
  slice_min(diff_days, with_ties = FALSE) %>%  # keep nearest-in-time match
  ungroup() %>%
  st_drop_geometry() %>%
  filter(survey_year %in% years_keep)


################################################################################
#Aggregate biomass per patch-year (absolute kg)

#biomass_patch <- pur_patch %>%
#  group_by(year, patch_id, pred_patch) %>%
#  summarise(total_biomass_g = sum(total_biomass_g, na.rm = TRUE), .groups = "drop") %>%
#  mutate(foraging_biomass_kg = total_biomass_g / 1000)  # convert to kg

# Aggregate biomass per patch-year
biomass_patch <- pur_patch %>%
  group_by(survey_year, patch_id, pred_patch) %>%
  summarise(
    total_biomass_g = sum(total_biomass_g, na.rm = TRUE),
    n_foraging_obs = n(),
    .groups = "drop"
  ) %>%
  mutate(foraging_biomass_kg = total_biomass_g / 1000) %>%
  rename(year = survey_year) 

################################################################################
#Add benthic patch-level predictors

predictors_focus <- c("purple_urchin_densitym2", "mean_gonad_mass_g")

patch_predictors <- quad_build3 %>%
  mutate(year = year(survey_date)) %>%
  st_drop_geometry() %>%
  group_by(patch_id, year, pred_patch) %>%
  summarise(across(all_of(predictors_focus), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

model_data <- left_join(biomass_patch, patch_predictors,
                        by = c("patch_id", "year", "pred_patch")) %>%
  filter(!is.na(foraging_biomass_kg)) %>%
  mutate(pred_patch = factor(pred_patch, levels = c("BAR","INCIP","FOR")))

################################################################################
#Build random forest model (classification)

rf_dat <- model_data %>%
  select(pred_patch, purple_urchin_densitym2, mean_gonad_mass_g, foraging_biomass_kg) %>%
  drop_na()

set.seed(42)
rf_patch <- randomForest(
  pred_patch ~ purple_urchin_densitym2 + mean_gonad_mass_g + foraging_biomass_kg,
  data = rf_dat,
  ntree = 1000,
  mtry = 2,
  importance = TRUE
)

print(rf_patch)
importance(rf_patch)

################################################################################
#Partial dependence plots (smoothed)
# ------------------------------------------------------------------------------
patch_cols <- c("BAR"="purple","INCIP"="orange","FOR"="forestgreen")

my_theme <- theme(
  axis.text.x = element_text(size=10, color = "black"),
  axis.text.y = element_text(size=10, color = "black"),
  axis.title  = element_text(size=12,color = "black"),
  legend.text = element_text(size=8,color = "black"),
  legend.title= element_text(size=8,color = "black"),
  plot.tag = element_text(size = 10, color = "black"),
  # Gridlines
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_blank(), 
  axis.line = element_line(colour = "black"),
  # Legend
  legend.key = element_blank(),
  legend.background = element_rect(fill=alpha('blue', 0)),
  # Facets
  strip.text = element_text(size=10, face = "bold",color = "black", hjust=0),
  strip.background = element_blank()
)

make_pdp <- function(var) {
  pd_all <- lapply(levels(rf_dat$pred_patch), function(cl) {
    pdp::partial(
      rf_patch,
      pred.var = var,
      which.class = cl,
      prob = TRUE,
      grid.resolution = 150
    ) %>% mutate(patch_type = cl)
  }) %>% bind_rows()
  
  pretty_names <- c(
    "purple_urchin_densitym2" = "Purple urchin density (No. per m²)",
    "mean_gonad_mass_g"       = "Mean gonad mass (g)",
    "foraging_biomass_kg"     = "Sea otter foraging effort \n(kg urchin biomass per patch-year)"
  )
  
  ggplot(pd_all, aes_string(x = var, y = "yhat", color = "patch_type")) +
    geom_smooth(se = FALSE, method = "loess", span = 0.6, linewidth = 1.4) +  # smoother lines
    scale_color_manual(values = patch_cols) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    labs(x = pretty_names[[var]], y = NULL) + my_theme
}

p1 <- make_pdp("purple_urchin_densitym2")
p2 <- make_pdp("mean_gonad_mass_g")
p3 <- make_pdp("foraging_biomass_kg")

pdp_grid <- (p1 | p2 | p3) +
  plot_annotation(
    title = "Partial dependence of patch type probability",
    #subtitle = "Predictors: sea otter foraging effort, urchin density, and gonad mass (Jun–Sep)",
    tag_levels = "A"
  ) &
  theme(plot.tag = element_text(size = 10))

pdp_grid

################################################################################
#Summarize range of foraging effort (for reporting)
# ------------------------------------------------------------------------------
summary_stats <- rf_dat %>%
  summarise(
    min_kg = min(foraging_biomass_kg, na.rm = TRUE),
    median_kg = median(foraging_biomass_kg, na.rm = TRUE),
    max_kg = max(foraging_biomass_kg, na.rm = TRUE)
  )
print(summary_stats)

################################################################################
#Export figure

 ggsave(
   here::here("figures","Fig3_foraging_patchtype_rf_absolute_smooth.png"),
   pdp_grid,
   width = 9, height = 4, dpi = 600
 )







