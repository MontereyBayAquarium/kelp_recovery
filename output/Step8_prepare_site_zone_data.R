
#jogsmith@ucsc.edu

rm(list = ls())

#data required for merge
#1. Benthic survey data averaged to zone level
#2. Dissection data averaged to zone level
#3. GIS isobath layers to split patch type
#4. Sea otter scan data 
#5. Patch area as determined by landsat clustering

################################################################################
#Step 0: set paths and load data
require(librarian)
librarian::shelf(tidyverse, lubridate, sf, stringr, purrr, terra)

datadir <- "/Volumes/enhydra/data/kelp_recovery/"
localdir <- here::here("output")

#load benthic survey data
load(file.path(datadir, "MBA_kelp_forest_database/processed/recovery/kelp_recovery_data.rda"))

#load dissection data
dissection_orig <- read_csv(file.path(datadir, "MBA_kelp_forest_database/processed/dissection/dissection_data_recovery.csv"))

#GIS layers
#bathy_5m_raw <- st_read(file.path(datadir, "gis_data/raw/bathymetry/contours_5m/contours_5m.shp"))
bathy_2m_raw <- rast("/Users/jossmith/Downloads/bat_ccsr_n_2m_bathy.tif")

#laod scans
scan_orig <- read_csv(file.path(here::here("output","scans","scans_data.csv")))

#load site patches
site_patches <- st_read(here::here("output","gis_data","processed","site_patch_polygons.shp"))

#load LDA-predicted patch types
lda_patch <- load(here::here("output","lda_patch_transitions.rda"))


################################################################################
#Step 1: process benthic survey data to obtain zone-level averages
quad_zone <- quad_data %>%
  dplyr::select(-substrate)%>%
  group_by(survey_type, region, latitude, longitude, site, site_type,
           survey_date, zone) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

#intermediate step: prepare site metadata table

site_table <- quad_zone %>% 
  group_by(site, site_type, zone) %>%
  distinct(latitude, longitude)

#write_csv(site_table, here::here("output","site_meta_data","site_table.csv"))


################################################################################
#Step 2: assign model-predicted patch types

str(quad_zone)
str(transitions_tbl_constrained)

quad_zone_with_pred <- quad_zone %>%
  left_join(
    transitions_tbl_constrained %>%
      dplyr::select(site, site_type, zone, patch_2024, patch_2025),
    by = c("site","site_type", "zone")
  ) %>%
  mutate(
    pred_patch = case_when(
      format(survey_date, "%Y") == "2024" ~ as.character(patch_2024),
      format(survey_date, "%Y") == "2025" ~ as.character(patch_2025),
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(-patch_2024, -patch_2025)


################################################################################
#Step 5: join with patch geometry

# Convert to sf using lat/lon
quad_zone_sf <- st_as_sf(
  quad_zone_with_pred,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)


site_patches_single <- site_patches %>%
  st_cast("POLYGON") %>%              # split MULTIPOLYGON into indiv POLYGON
  mutate(patch_id = row_number())     # assign unique polygon ID

plot(site_patches_single)

#join points to polygons
site_patches_with_points <- site_patches_single %>%
  st_join(quad_zone_sf, join = st_intersects, left = TRUE) %>%
  filter(!(is.na(survey_type)))


#inspect
ggplot(site_patches_with_points %>% filter(year(survey_date) == 2024)) +
  geom_sf(aes(fill = pred_patch), color = "black") +
  theme_minimal() +
  labs(
    title = "Predicted Patch Type by Independent Polygon (2024)",
    fill = "Predicted Patch"
  )

str(site_patches_with_points)



quad_build3 <- site_patches_with_points %>%
                mutate(patch_cat = ifelse(year(survey_date) == 2024,"predicted 2024","predicted 2025"))

ggplot(quad_build3) +
  geom_sf(aes(fill = pred_patch), color = "black") +
  facet_wrap(~patch_cat, nrow=1)+
  theme_minimal() +
  labs(
    title = "Predicted Patch Type by Independent Polygon (2024)",
    fill = "Predicted Patch"
  )


