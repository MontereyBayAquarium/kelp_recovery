
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
librarian::shelf(tidyverse, lubridate, sf, stringr, purrr)

datadir <- "/Volumes/enhydra/data/kelp_recovery/"
localdir <- here::here("output")

#load benthic survey data
load(file.path(datadir, "MBA_kelp_forest_database/processed/recovery/kelp_recovery_data.rda"))

#load dissection data
dissection_orig <- read_csv(file.path(datadir, "MBA_kelp_forest_database/processed/dissection/dissection_data_recovery.csv"))

#GIS layers
bathy_10m <- st_read(file.path(datadir, "gis_data/raw/bathymetry/contours_5m/contours_5m.shp")) %>% filter(CONTOUR == "-10")

#laod scans
scan_orig <- read_csv(file.path(here::here("output","scans","scans_data.csv")))

#set patch for patch polygons
poly_dir <- file.path(localdir, "recovery_polygons")


################################################################################
#Step 1: process benthic survey data to obtain zone-level averages
quad_zone <- quad_data %>%
              dplyr::select(-substrate)%>%
              group_by(survey_type, region, latitude, longitude, site, site_type,
                       survey_date, zone) %>%
              summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
                        .groups = "drop")



################################################################################
# Step 2: Read and prepare individual polygon patches (KMZ to KML to separate polygons)
kmz_files <- list.files(poly_dir, pattern = "\\.kmz$", full.names = TRUE, ignore.case = TRUE)

read_kmz <- function(kmz) {
  tmpdir <- tempfile(); dir.create(tmpdir)
  unzip(kmz, exdir = tmpdir)
  kmls <- list.files(tmpdir, pattern = "\\.kml$", full.names = TRUE)
  if (length(kmls) == 0) stop(paste("No KML found in KMZ:", kmz))
  sf::st_read(kmls[1], quiet = TRUE)
}

poly_list <- lapply(kmz_files, function(f) {
  obj <- read_kmz(f)
  obj$source_file <- basename(f)
  obj
})
polys_raw <- do.call(rbind, poly_list)
polys_individual <- polys_raw %>%
  filter(st_geometry_type(.) %in% c("POLYGON","MULTIPOLYGON")) %>%
  st_cast("POLYGON") %>%
  mutate(
    patch_type = case_when(
      str_detect(source_file, "incipient_forests")  ~ "Incipient Forests",
      str_detect(source_file, "persistent_barrens") ~ "Persistent Barrens",
      str_detect(source_file, "persistent_forests") ~ "Persistent Forests",
      TRUE                                            ~ source_file
    ),
    patch_id = row_number()
  ) %>%
  dplyr::select(patch_id, patch_type, geometry)





