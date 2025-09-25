
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
bathy_5m_raw <- st_read(file.path(datadir, "gis_data/raw/bathymetry/contours_5m/contours_5m.shp"))


# Download to local working directory (change path if you prefer)
drive_download(file, path = "bat_nccsr_2m_bathy.tif", overwrite = TRUE)

# Load into R
bathy_raster <- rast("bat_nccsr_2m_bathy.tif")
plot(bathy_raster)


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

plot(polys_individual)


################################################################################
#Step 3: intersect quad data with patches

#convert quad_data to sf points
quad_sf <- st_as_sf(
  quad_data,
  coords = c("longitude", "latitude"),
  crs = 4326,     
  remove = FALSE  
)

#Check CRS and align polygons
#transform quad points to match polygon CRS if needed
if (st_crs(quad_sf) != st_crs(polys_individual)) {
  quad_sf <- st_transform(quad_sf, st_crs(polys_individual))
}

#Spatial join: assign patch_type & geometry
quad_with_patch <- st_join(
  quad_sf,
  polys_individual %>% dplyr::select(patch_id, patch_type),
  join = st_intersects,
  left = TRUE
)

#Check what didn' match
# Plot polygons with fill = patch_type, then overlay quadrat points
ggplot() +
  # polygons
  geom_sf(data = polys_individual, aes(fill = patch_type), alpha = 0.4, color = "black") +
  # quadrats with patch type
  geom_sf(data = filter(quad_with_patch, !is.na(patch_type)),
          aes(color = patch_type), size = 2) +
  # quadrats with no patch assignment
  geom_sf(data = filter(quad_with_patch, is.na(patch_type)),
          color = "black", fill = "black", shape = 21, size = 2, stroke = 0.5) +
  scale_fill_viridis_d(option = "plasma") +
  scale_color_viridis_d(option = "plasma", na.value = "black") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_blank()
  ) 


# Find which points had no intersection
quad_na <- quad_with_patch %>% filter(is.na(patch_type))

# For each NA point, find nearest polygon
nearest_ids <- st_nearest_feature(quad_na, polys_individual)

# Join patch attributes back
quad_na_fixed <- quad_na %>%
  mutate(
    patch_id   = polys_individual$patch_id[nearest_ids],
    patch_type = polys_individual$patch_type[nearest_ids]
  )

# Replace NA rows with fixed ones
quad_with_patch_snapped <- quad_with_patch %>%
  filter(!is.na(patch_type)) %>%   # keep valid matches
  bind_rows(quad_na_fixed)         # add snapped points back in


#Reexamine
ggplot() +
  # polygons
  geom_sf(data = polys_individual, aes(fill = patch_type), alpha = 0.4, color = "black") +
  # quadrats with patch type
  geom_sf(data = filter(quad_with_patch_snapped, !is.na(patch_type)),
          aes(color = patch_type), size = 2) +
  # quadrats with no patch assignment
  geom_sf(data = filter(quad_with_patch_snapped, is.na(patch_type)),
          color = "black", fill = "black", shape = 21, size = 2, stroke = 0.5) +
  scale_fill_viridis_d(option = "plasma") +
  scale_color_viridis_d(option = "plasma", na.value = "black") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_blank()
  ) 

#tidy up
quad_build1 <- quad_with_patch_snapped %>%
                  dplyr::select(-patch_id) %>%
                  rename(landsat_patch = patch_type)

rm(quad_with_patch)
rm(quad_with_patch_snapped)
rm(quad_na_fixed)
rm(dup_counts)

View(quad_build1)


################################################################################
#Step 4: process bathy

monterey_bbox <- st_as_sfc(st_bbox(c(
  xmin = -122.0,
  xmax = -121.85,
  ymin = 36.5,
  ymax = 36.65
), crs = 4326))

monterey_bbox_proj <- st_transform(monterey_bbox, st_crs(bathy_raw))

bathy_mpen <- bathy_raw %>%
  filter(CONTOUR %in% c(-5, -10, -15, -20)) %>%
  st_intersection(monterey_bbox_proj)

plot(bathy_mpen)
rm(bath_raw)

