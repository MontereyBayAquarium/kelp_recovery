
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

#intermediate step: prepare site metadata table

site_table <- quad_zone %>% 
                group_by(site, site_type, zone) %>%
                distinct(latitude, longitude)

write_csv(site_table, here::here("output","site_meta_data","site_table.csv"))

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

#--------------------------------------------------
# 2. Define Monterey Peninsula bounding box (WGS84)
#--------------------------------------------------
monterey_bbox <- st_as_sfc(st_bbox(c(
  xmin = -122.0,
  xmax = -121.85,
  ymin = 36.5,
  ymax = 36.65
), crs = 4326))

# Reproject bbox to raster CRS
monterey_bbox_proj <- st_transform(monterey_bbox, st_crs(bathy_2m_raw))

# Convert sf polygon to terra vector
monterey_bbox_spat <- vect(monterey_bbox_proj)

#--------------------------------------------------
# 3. Crop and mask raster
#--------------------------------------------------
bathy_mpen_raster <- crop(bathy_2m_raw, monterey_bbox_spat)
bathy_mpen_raster <- mask(bathy_mpen_raster, monterey_bbox_spat)

#--------------------------------------------------
# 4. Define 2-m breakpoints and build reclass matrix
#--------------------------------------------------
levels <- seq(-20, 0, by = 2)

# Class values = 1, 2, 3, … instead of depths
rcl <- cbind(head(levels, -1), tail(levels, -1), seq_along(head(levels, -1)))

# Add NA bins for outside range
min_val <- min(values(bathy_mpen_raster), na.rm = TRUE)
max_val <- max(values(bathy_mpen_raster), na.rm = TRUE)

rcl <- rbind(c(min_val, -20, NA), rcl, c(0, max_val, NA))

#--------------------------------------------------
# 5. Reclassify raster into 2-m bands
#--------------------------------------------------
bathy_class <- classify(bathy_mpen_raster, rcl, include.lowest = TRUE)

# Convert raster to polygons
bathy_polys <- as.polygons(bathy_class, dissolve = TRUE, values = TRUE)
bathy_polys_sf <- st_as_sf(bathy_polys)

#--------------------------------------------------
# 6. Add human-readable depth labels
#--------------------------------------------------
depth_labels <- paste(head(levels, -1), "to", tail(levels, -1))
bathy_polys_sf$depth_band <- factor(
  depth_labels[bathy_polys_sf$bat_ccsr_n_2m_bathy],
  levels = depth_labels
)

#--------------------------------------------------
# 7. Plot polygons
#--------------------------------------------------
ggplot() +
  geom_sf(data = bathy_polys_sf, aes(fill = depth_band), color = NA) +
  geom_sf(data = monterey_bbox_proj, fill = NA, color = "red", linetype = "dashed") +
  #scale_fill_viridis_d(option = "plasma", name = "Depth (m)") +
  theme_minimal() +
  labs(
    title = "Bathymetry of Monterey Peninsula",
    subtitle = "2-m depth bands (-20 m to 0 m)",
    x = "Longitude", y = "Latitude"
  )

#export this as a shapefile for use in QGIS
st_write(
  bathy_polys_sf,
  here::here("output","gis_data","processed","mpen_2m_bathy.shp"),     #
  driver = "ESRI Shapefile",
  delete_dsn = TRUE         # overwrite if it already exists
)


################################################################################
#Step 4: split polygons by depth band

# --- libraries ---
library(sf)
library(dplyr)
library(ggplot2)
library(lwgeom)   # for st_split
sf_use_s2(FALSE)  # planar ops, fewer surprises

# --- inputs assumed in memory ---
# quad_build1: sf POINTs with columns site, site_type, zone (Shallow/Deep), etc.
# polys_individual: sf (MULTI)POLYGON with patch_id, patch_type
coast_path <- "/Volumes/enhydra/data/kelp_recovery/gis_data/raw/Coastn83/coastn83.shp"

# ===============================
# 0) CRS & read / prepare layers
# ===============================
target_crs <- 3310
coastline  <- st_read(coast_path, quiet = TRUE) |> st_transform(target_crs)
polys_individual <- st_transform(polys_individual, target_crs) |> st_make_valid()
quad_build1      <- st_transform(quad_build1,      target_crs)

# ============================================================
# 1) Join quads to patches & keep only quads on a patch
# ============================================================
quad_with_patch <- st_join(
  quad_build1,
  polys_individual |> dplyr::select(patch_id, patch_type),
  join = st_intersects
)
quad_valid <- quad_with_patch |> filter(!is.na(patch_id))

# =====================================================================
# 2) One representative point per site/site_type/patch_id/zone (unique)
# =====================================================================
zone_points <- quad_valid |>
  group_by(site, site_type, patch_id, zone) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  st_centroid()

# ==============================================================
# 3) Simplify & COMBINE coastline (single MULTILINESTRING)
# ==============================================================
coast_tol  <- 150  # meters; tune (e.g., 100–300)
coast_simpl <- st_simplify(coastline, dTolerance = coast_tol)
coast_ml    <- st_cast(st_combine(coast_simpl), "MULTILINESTRING") |> st_make_valid()

# --------------------------------------------------------------
# 3a) Monterey viewport (in 3310) + quick diagnostic plot
# --------------------------------------------------------------
monterey_wgs  <- st_as_sfc(st_bbox(c(xmin = -122, xmax = -121.88,
                                     ymin = 36.52,  ymax = 36.65), crs = 4326))
monterey_3310 <- st_transform(monterey_wgs, target_crs)
mbb <- st_bbox(monterey_3310)

ggplot() +
  geom_sf(data = polys_individual, fill = "grey85", color = "black", alpha = 0.3) +
  geom_sf(data = st_as_sf(coast_simpl), color = "red", linewidth = 0.25) +
  geom_sf(data = zone_points, aes(color = zone), size = 1.2) +
  scale_color_manual(values = c(Shallow = "#1f78b4", Deep = "#33a02c")) +
  coord_sf(xlim = c(mbb["xmin"], mbb["xmax"]), ylim = c(mbb["ymin"], mbb["ymax"])) +
  theme_minimal() +
  labs(title = sprintf("Diagnostic: Simplified Coastline (tol = %dm)", coast_tol),
       subtitle = "Grey = patches; Red = simplified coast; points = Shallow/Deep",
       color = "Zone")

# ==========
# helpers
# ==========
# min scalar distance from point to MULTILINESTRING
dist_min <- function(pt, lines_multi) as.numeric(min(st_distance(pt, lines_multi)))

# keep the longest line feature (after clipping etc.)
longest_line <- function(lines_sf) {
  if (!inherits(lines_sf, "sf") || nrow(lines_sf) == 0) return(lines_sf)
  lines_sf[which.max(st_length(lines_sf)), , drop = FALSE]
}

# ===========================================================
# 4) Split one patch using coastline-parallel offset at
#    midpoint distance between Shallow/Deep to coastline
# ===========================================================
split_patch_by_offset <- function(g_one) {
  if (!all(c("Shallow","Deep") %in% g_one$zone)) return(NULL)
  
  pid   <- unique(g_one$patch_id)
  site  <- unique(g_one$site)
  stype <- unique(g_one$site_type)
  
  patch_poly <- polys_individual |> dplyr::filter(patch_id == pid)
  if (nrow(patch_poly) == 0 || st_is_empty(patch_poly)) return(NULL)
  patch_poly <- st_make_valid(patch_poly)
  
  shallow_pt <- g_one |> dplyr::filter(zone=="Shallow") |> st_geometry() |> st_union() |> st_centroid()
  deep_pt    <- g_one |> dplyr::filter(zone=="Deep")    |> st_geometry() |> st_union() |> st_centroid()
  
  ds <- as.numeric(min(st_distance(shallow_pt, coast_ml)))
  dd <- as.numeric(min(st_distance(deep_pt,    coast_ml)))
  m  <- (ds + dd) / 2
  if (!is.finite(m) || m <= 0) return(NULL)
  
  win <- st_buffer(st_geometry(patch_poly), dist = (2*m + 1500))
  local_coast <- tryCatch(st_intersection(coast_ml, win), error = function(e) NULL)
  if (is.null(local_coast) || length(local_coast) == 0) local_coast <- coast_ml
  
  iso_poly <- st_buffer(local_coast, dist = m)
  iso_line <- st_boundary(iso_poly) |> st_collection_extract("LINESTRING")
  if (length(iso_line) == 0) return(NULL)
  
  iso_local <- tryCatch(
    st_intersection(st_as_sf(iso_line), st_buffer(st_geometry(patch_poly), m + 200)),
    error = function(e) st_as_sf(iso_line)
  )
  if (nrow(iso_local) == 0) return(NULL)
  
  cut_line <- iso_local[which.max(st_length(iso_local)), , drop = FALSE]
  
  if (!any(st_intersects(st_geometry(cut_line), st_geometry(patch_poly), sparse = FALSE))) {
    iso_local2 <- tryCatch(
      st_intersection(st_as_sf(iso_line), st_buffer(st_geometry(patch_poly), 3*m + 200)),
      error = function(e) NULL
    )
    if (!is.null(iso_local2) && nrow(iso_local2) > 0) {
      cut_line <- iso_local2[which.max(st_length(iso_local2)), , drop = FALSE]
    }
  }
  if (!any(st_intersects(st_geometry(cut_line), st_geometry(patch_poly), sparse = FALSE))) return(NULL)
  
  pieces <- tryCatch({
    st_collection_extract(st_split(st_make_valid(patch_poly), st_geometry(cut_line)), "POLYGON")
  }, error = function(e) NULL)
  if (is.null(pieces) || nrow(pieces) == 0) return(NULL)
  
  pieces <- st_make_valid(pieces)
  pieces <- st_collection_extract(pieces, "POLYGON")
  pieces <- pieces[!st_is_empty(pieces), ]
  if (nrow(pieces) == 0) return(NULL)
  
  A_all  <- as.numeric(st_area(st_union(patch_poly)))
  A_part <- as.numeric(st_area(pieces))
  keep   <- A_part > max(A_all * 0.01, A_all * 0.005)
  pieces <- pieces[keep, ]
  if (nrow(pieces) == 0) return(NULL)
  
  # --- fixed labeling ---
  pc <- st_point_on_surface(pieces)
  d_sh <- as.numeric(st_distance(pc, shallow_pt))
  d_dp <- as.numeric(st_distance(pc, deep_pt))
  bad <- is.na(d_sh) | is.na(d_dp)
  if (any(bad)) {
    reps <- st_sfc(st_geometry(shallow_pt)[[1]], st_geometry(deep_pt)[[1]], crs = st_crs(pieces))
    nn   <- st_nearest_feature(pc[bad], st_as_sf(data.frame(id=1:2), geometry=reps))
    d_sh[bad] <- ifelse(nn == 1, 0, Inf)
    d_dp[bad] <- ifelse(nn == 2, 0, Inf)
  }
  pieces$zone <- ifelse(d_sh < d_dp, "Shallow", "Deep")
  
  pieces <- pieces |>
    mutate(area_m2 = as.numeric(st_area(geometry))) |>
    group_by(zone) |>
    slice_max(order_by = area_m2, n = 1, with_ties = FALSE) |>
    ungroup()
  
  pieces$site       <- site
  pieces$site_type  <- stype
  pieces$patch_id   <- pid
  pieces$patch_type <- unique(patch_poly$patch_type)
  pieces$mid_dist_m <- m
  pieces$coast_tol  <- coast_tol
  
  pieces |> dplyr::select(patch_id, patch_type, site, site_type, zone, mid_dist_m, coast_tol, geometry)
}

# =====================================================
# 5) Apply per (site, site_type, patch_id)
# =====================================================
groups <- split(
  zone_points,
  interaction(zone_points$site, zone_points$site_type, zone_points$patch_id, drop = TRUE)
)

split_list <- lapply(groups, split_patch_by_offset)
split_list <- split_list[!sapply(split_list, is.null)]
patches_split <- if (length(split_list)) dplyr::bind_rows(split_list) else NULL

message("✅ Groups attempted: ", length(groups),
        " | successful splits: ", ifelse(is.null(patches_split), 0, nrow(patches_split)))

# =========================================
# 6) Plot Monterey results
# =========================================
if (!is.null(patches_split)) {
  ggplot() +
    geom_sf(data = polys_individual, aes(fill = patch_type),
            color = "grey80", alpha = 0.25, show.legend = FALSE) +
    geom_sf(data = patches_split, aes(fill = zone),
            color = "black", alpha = 0.6) +
    geom_sf(data = zone_points, aes(color = zone), size = 1.2) +
    geom_sf(data = st_as_sf(coast_simpl), color = "red", linewidth = 0.25) +
    scale_fill_manual(values = c(Shallow = "#1f78b4", Deep = "#33a02c")) +
    scale_color_manual(values = c(Shallow = "#1f78b4", Deep = "#33a02c")) +
    coord_sf(xlim = c(mbb["xmin"], mbb["xmax"]), ylim = c(mbb["ymin"], mbb["ymax"])) +
    theme_minimal() +
    labs(title = "Patch splits by coastline-parallel midpoint offset",
         subtitle = sprintf("Coastline simplified (tol = %dm). Slivers dropped; robust zone labeling.", coast_tol),
         fill = "Zone", color = "Zone")
} else {
  message("⚠️ No split patches were created")
}













# --- libraries ---
library(sf)
library(dplyr)
library(ggplot2)
sf_use_s2(FALSE)  # planar ops, fewer surprises

# --- inputs in memory assumed ---
# quad_build1: sf POINTs with site, site_type, zone ("Shallow"/"Deep"), etc.
# polys_individual: sf (MULTI)POLYGON with patch_id, patch_type
coast_path <- "/Volumes/enhydra/data/kelp_recovery/gis_data/raw/Coastn83/coastn83.shp"

# --- 0) CRS & read / prepare ---
target_crs <- 3310
coastline  <- st_read(coast_path, quiet = TRUE) |> st_transform(target_crs)
polys_individual <- st_transform(polys_individual, target_crs) |> st_make_valid()
quad_build1      <- st_transform(quad_build1,      target_crs)

# --- 1) Join quads to patches & keep those that hit a patch ---
quad_with_patch <- st_join(
  quad_build1,
  polys_individual |> dplyr::select(patch_id, patch_type),
  join = st_intersects
)
quad_valid <- quad_with_patch |> filter(!is.na(patch_id))

# --- 2) One representative point per site/site_type/patch_id/zone ---
zone_points <- quad_valid |>
  group_by(site, site_type, patch_id, zone) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  st_centroid()

# --- 3) Simplify & COMBINE coastline once (important) ---
coast_tol  <- 100   # meters; adjust if needed (100–300)
coast_simpl <- st_simplify(coastline, dTolerance = coast_tol)
coast_ml     <- st_cast(st_combine(coast_simpl), "MULTILINESTRING")
coast_ml_sf  <- st_sf(geometry = coast_ml)  # handy sf wrapper

# --- (optional) Monterey viewport & diagnostic plot ---
monterey_wgs  <- st_as_sfc(st_bbox(c(xmin = -122, xmax = -121.88,
                                     ymin = 36.52,  ymax = 36.65), crs = 4326))
monterey_3310 <- st_transform(monterey_wgs, target_crs)
mbb <- st_bbox(monterey_3310)

ggplot() +
  geom_sf(data = polys_individual, fill = "grey85", color = "black", alpha = 0.3) +
  geom_sf(data = st_as_sf(coast_simpl), color = "red", linewidth = 0.25) +
  geom_sf(data = zone_points, aes(color = zone), size = 1.2) +
  scale_color_manual(values = c(Shallow = "#1f78b4", Deep = "#33a02c")) +
  coord_sf(xlim = c(mbb["xmin"], mbb["xmax"]), ylim = c(mbb["ymin"], mbb["ymax"])) +
  theme_minimal() +
  labs(title = sprintf("Diagnostic: Simplified Coastline (tol = %dm)", coast_tol),
       subtitle = "Grey = patches; Red = simplified coast; points = Shallow/Deep",
       color = "Zone")

# --- helpers ---
has_nonempty <- function(x) {
  if (is.null(x)) return(FALSE)
  if (inherits(x, "sf"))   return(nrow(x) > 0 && any(!st_is_empty(x)))
  if (inherits(x, "sfc"))  return(length(x) > 0 && any(!st_is_empty(x)))
  FALSE
}

# scalar min distance point -> (multi)lines
dist_min <- function(pt, lines_multi) as.numeric(min(st_distance(pt, lines_multi)))

# dissolve to clean sf MULTIPOLYGON, ensure geometry column exists (no hard-coded name)
dissolve_clean <- function(x) {
  if (!has_nonempty(x)) return(NULL)
  g <- st_make_valid(st_union(x))
  g <- st_buffer(g, 0)
  g <- st_collection_extract(g, "POLYGON")
  out <- st_as_sf(g)
  if (!has_nonempty(out)) return(NULL)
  out
}

# split function: area-preserving using coastline buffer
split_patch_coastbuffer <- function(g_one,
                                    scale_tries   = c(1.00, 1.15, 0.90, 1.30, 0.75),
                                    gap_tol_frac  = 1e-4) {
  # need both zones
  if (!all(c("Shallow","Deep") %in% g_one$zone)) return(NULL)
  
  pid   <- unique(g_one$patch_id)
  site  <- unique(g_one$site)
  stype <- unique(g_one$site_type)
  
  patch_poly <- polys_individual |> filter(patch_id == pid)
  if (nrow(patch_poly) == 0 || st_is_empty(patch_poly)) return(NULL)
  patch_poly <- st_make_valid(patch_poly)
  
  shallow_pt <- g_one |> filter(zone == "Shallow") |> st_geometry() |> st_union() |> st_centroid()
  deep_pt    <- g_one |> filter(zone == "Deep")    |> st_geometry() |> st_union() |> st_centroid()
  
  ds <- dist_min(shallow_pt, coast_ml)
  dd <- dist_min(deep_pt,    coast_ml)
  m0 <- (ds + dd) / 2
  if (!is.finite(m0) || m0 <= 0) return(NULL)
  
  base_win <- st_buffer(st_geometry(patch_poly), dist = (2*m0 + 1500))
  A_patch  <- as.numeric(st_area(st_union(patch_poly)))
  
  best <- NULL
  
  for (f in scale_tries) {
    m <- m0 * f
    if (!is.finite(m) || m <= 0) next
    
    # local coast lines (crop for speed/stability)
    local_raw <- tryCatch(st_intersection(coast_ml_sf, base_win), error = function(e) NULL)
    local_lines <- if (has_nonempty(local_raw)) local_raw else coast_ml_sf
    
    # coastline buffer at m -> nearshore mask
    coast_buf <- tryCatch(st_buffer(st_union(st_geometry(local_lines)), dist = m),
                          error = function(e) NULL)
    if (is.null(coast_buf) || all(st_is_empty(coast_buf))) next
    coast_buf <- st_make_valid(coast_buf)
    
    # nearshore = patch ∩ buffer; offshore = patch \ nearshore
    shallow_raw <- tryCatch(st_intersection(patch_poly, coast_buf), error = function(e) NULL)
    if (!has_nonempty(shallow_raw)) next
    shallow_raw <- st_collection_extract(shallow_raw, "POLYGON", warn = FALSE)
    shallow_raw <- shallow_raw[!st_is_empty(shallow_raw), ]
    
    shallow_mask <- st_make_valid(st_union(shallow_raw))
    deep_raw     <- tryCatch(st_difference(patch_poly, shallow_mask), error = function(e) NULL)
    if (!has_nonempty(deep_raw)) next
    deep_raw <- st_collection_extract(deep_raw, "POLYGON", warn = FALSE)
    deep_raw <- deep_raw[!st_is_empty(deep_raw), ]
    
    # dissolve & clean
    shallow_1 <- dissolve_clean(st_as_sf(shallow_raw))
    deep_1    <- dissolve_clean(st_as_sf(deep_raw))
    if (is.null(shallow_1) || is.null(deep_1)) next
    
    # repair tiny leftover gaps (donate to nearest side)
    union_sd <- st_make_valid(st_union(shallow_1, deep_1))
    gap <- tryCatch(st_difference(st_geometry(patch_poly), union_sd), error = function(e) NULL)
    if (!is.null(gap)) {
      gap <- st_collection_extract(gap, "POLYGON", warn = FALSE)
      if (inherits(gap, "sfc")) gap <- st_as_sf(gap)
      if (has_nonempty(gap)) {
        A_gap <- sum(as.numeric(st_area(gap)))
        if (A_gap / A_patch > gap_tol_frac) {
          # nearest by centroid to shallow or deep rep points
          gc <- st_point_on_surface(gap)
          dS <- as.numeric(st_distance(gc, shallow_pt))
          dD <- as.numeric(st_distance(gc, deep_pt))
          to_S <- gap[dS <= dD, , drop = FALSE]
          to_D <- gap[dD  <  dS, , drop = FALSE]
          if (has_nonempty(to_S)) shallow_1 <- dissolve_clean(rbind(shallow_1, to_S))
          if (has_nonempty(to_D)) deep_1    <- dissolve_clean(rbind(deep_1,    to_D))
        }
      }
    }
    
    # final union sanity; if some tiny leftovers persist, accept (below tol)
    union_sd2 <- st_make_valid(st_union(shallow_1, deep_1))
    leftover <- tryCatch(st_difference(st_geometry(patch_poly), union_sd2), error = function(e) NULL)
    if (!is.null(leftover)) {
      leftover <- st_collection_extract(leftover, "POLYGON", warn = FALSE)
      if (inherits(leftover, "sfc")) leftover <- st_as_sf(leftover)
      if (has_nonempty(leftover)) {
        A_left <- sum(as.numeric(st_area(leftover)))
        if (A_left / A_patch > gap_tol_frac) next  # try next scale
      }
    }
    
    # attach labels & attributes (DON'T select geometry by name)
    shallow_1$zone <- "Shallow"
    deep_1$zone    <- "Deep"
    
    out <- dplyr::bind_rows(shallow_1, deep_1)
    out$site       <- site
    out$site_type  <- stype
    out$patch_id   <- pid
    out$patch_type <- unique(patch_poly$patch_type)
    out$mid_dist_m <- m0
    out$coast_tol  <- coast_tol
    
    best <- out
    break
  }
  
  best
}

# --- 5) Apply per (site, site_type, patch_id) ---
groups <- split(zone_points, interaction(zone_points$site,
                                         zone_points$site_type,
                                         zone_points$patch_id,
                                         drop = TRUE))
split_list <- lapply(groups, split_patch_coastbuffer)
split_list <- split_list[!sapply(split_list, is.null)]
patches_split <- if (length(split_list)) dplyr::bind_rows(split_list) else NULL

message("✅ Groups attempted: ", length(groups),
        " | successful splits: ", ifelse(is.null(patches_split), 0, nrow(patches_split)))

# --- 6) Plot Monterey results ---
if (!is.null(patches_split)) {
  ggplot() +
    geom_sf(data = polys_individual, aes(fill = patch_type),
            color = "grey80", alpha = 0.25, show.legend = FALSE) +
    geom_sf(data = patches_split, aes(fill = zone),
            color = "black", alpha = 0.6) +
    geom_sf(data = zone_points, aes(color = zone), size = 1.2) +
    geom_sf(data = st_as_sf(coast_simpl), color = "red", linewidth = 0.25) +
    scale_fill_manual(values = c(Shallow="#1f78b4", Deep="#33a02c")) +
    scale_color_manual(values = c(Shallow="#1f78b4", Deep="#33a02c")) +
    coord_sf(xlim = c(mbb["xmin"], mbb["xmax"]), ylim = c(mbb["ymin"], mbb["ymax"])) +
    theme_minimal() +
    labs(title = "Area-preserving onshore/offshore patch split",
         subtitle = sprintf("Nearshore = patch ∩ buffer(coast, m̄). coast_tol=%dm; scales=%s",
                            coast_tol, paste0(c(1.00,1.15,0.90,1.30,0.75), collapse=", ")),
         fill = "Zone", color = "Zone")
} else {
  message("⚠️ No split patches were created")

}

























# ======================================================
# Split each patch by coastline-parallel midpoint offset
# (shallow/deep defined by distance to coastline)
# ======================================================

# --- libraries ---
library(sf)
library(dplyr)
library(ggplot2)
library(lwgeom)        # st_split backend
sf_use_s2(FALSE)       # planar ops for stability

# --- PARAMETERS ---
target_crs  <- 3310    # California Albers (meters)
coast_path  <- "/Volumes/enhydra/data/kelp_recovery/gis_data/raw/Coastn83/coastn83.shp"
coast_tol   <- 150     # coastline simplification tolerance (m); try 100–300
m_nudge     <- 1.00    # 0.98 or 1.02 to bias split slightly shoreward/offshore
sliver_frac <- 0.01    # drop pieces < 1% of original patch area
tiny_abs_m2 <- 5       # also drop absolute tinies (< 5 m^2)

# Monterey viewport for plotting (in WGS84; will transform below)
monterey_wgs  <- st_as_sfc(st_bbox(c(xmin = -122.00, xmax = -121.88,
                                     ymin = 36.52,  ymax = 36.65), crs = 4326))

# --- INPUTS IN MEMORY REQUIRED ---
# quad_build1: sf POINTS with at least: site, site_type, zone (Shallow/Deep)
# polys_individual: sf (MULTI)POLYGON with: patch_id, patch_type

stopifnot(inherits(quad_build1, "sf"), inherits(polys_individual, "sf"))

# --- LOAD & PREPARE COASTLINE ---
coastline <- st_read(coast_path, quiet = TRUE) |>
  st_transform(target_crs)

# --- REPROJECT OTHER LAYERS & MAKE VALID ---
quad_build1      <- st_transform(quad_build1,      target_crs)
polys_individual <- st_transform(polys_individual, target_crs) |> st_make_valid()

# --- JOIN QUADS → PATCH (keep only those on a patch) ---
quad_with_patch <- st_join(
  quad_build1,
  polys_individual |> dplyr::select(patch_id, patch_type),
  join = st_intersects
)
quad_valid <- quad_with_patch |> filter(!is.na(patch_id))

# One representative point per (site, site_type, patch_id, zone)
zone_points <- quad_valid |>
  group_by(site, site_type, patch_id, zone) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  st_centroid()

# --- SIMPLIFY & COMBINE COASTLINE (single MULTILINESTRING) ---
coast_simpl <- st_simplify(coastline, dTolerance = coast_tol)
coast_ml    <- st_cast(st_combine(coast_simpl), "MULTILINESTRING") |> st_make_valid()

# --- DIAGNOSTIC PLOT (before splitting) ---
monterey_3310 <- st_transform(monterey_wgs, target_crs)
mbb <- st_bbox(monterey_3310)

ggplot() +
  geom_sf(data = polys_individual, fill = "grey85", color = "black", alpha = 0.3) +
  geom_sf(data = st_as_sf(coast_simpl), color = "red", linewidth = 0.25) +
  geom_sf(data = zone_points, aes(color = zone), size = 1.2) +
  scale_color_manual(values = c(Shallow = "#1f78b4", Deep = "#33a02c")) +
  coord_sf(xlim = c(mbb["xmin"], mbb["xmax"]), ylim = c(mbb["ymin"], mbb["ymax"])) +
  theme_minimal() +
  labs(title = sprintf("Diagnostic: Simplified Coastline (tol = %dm)", coast_tol),
       subtitle = "Grey = patches; Red = simplified coast; points = Shallow/Deep",
       color = "Zone")

# --- HELPERS ---
has_nonempty <- function(x) {
  if (is.null(x)) return(FALSE)
  if (inherits(x, "sf"))  return(nrow(x) > 0 && any(!st_is_empty(x)))
  if (inherits(x, "sfc")) return(length(x) > 0 && any(!st_is_empty(x)))
  FALSE
}

dissolve_clean <- function(x) {
  if (!has_nonempty(x)) return(NULL)
  g <- st_make_valid(st_union(x))
  g <- st_buffer(g, 0)
  g <- st_collection_extract(g, "POLYGON")
  out <- st_as_sf(g)
  if (!has_nonempty(out)) return(NULL)
  out
}

# --- constants used as defaults (avoid same-name defaults) ---
COAST_TOL_DEFAULT   <- coast_tol    # you already set coast_tol earlier
M_NUDGE_DEFAULT     <- 1.00         # try 0.98..1.02 to bias shore/offshore
SLIVER_FRAC_DEFAULT <- 0.01
TINY_ABS_M2_DEFAULT <- 5

# --- REPLACE your split_patch_by_offset() with this version ---
split_patch_by_offset <- function(
    g_one,
    coast_ml,
    coast_tol_used  = COAST_TOL_DEFAULT,
    m_nudge_val     = M_NUDGE_DEFAULT,
    sliver_frac_val = SLIVER_FRAC_DEFAULT,
    tiny_abs_m2_val = TINY_ABS_M2_DEFAULT
) {
  # must have both zones present
  if (!all(c("Shallow","Deep") %in% g_one$zone)) return(NULL)
  
  pid   <- unique(g_one$patch_id)
  site  <- unique(g_one$site)
  stype <- unique(g_one$site_type)
  
  patch_poly <- polys_individual |>
    dplyr::filter(patch_id == pid) |>
    st_make_valid()
  if (nrow(patch_poly) == 0 || st_is_empty(patch_poly)) return(NULL)
  
  # representative points (only for distance calc)
  shallow_pt <- g_one |> dplyr::filter(zone=="Shallow") |> st_geometry() |> st_union() |> st_centroid()
  deep_pt    <- g_one |> dplyr::filter(zone=="Deep")    |> st_geometry() |> st_union() |> st_centroid()
  
  # midpoint distance to coastline
  ds <- as.numeric(min(st_distance(shallow_pt, coast_ml)))
  dd <- as.numeric(min(st_distance(deep_pt,    coast_ml)))
  if (!is.finite(ds) || !is.finite(dd)) return(NULL)
  m  <- ((ds + dd) / 2) * m_nudge_val
  if (!is.finite(m) || m <= 0) return(NULL)
  
  # big enough window so isoline is long enough
  win <- st_buffer(st_geometry(patch_poly), dist = (2*m + 1500))
  local_coast <- tryCatch(st_intersection(coast_ml, win), error = function(e) coast_ml)
  if (is.null(local_coast) || length(local_coast) == 0) local_coast <- coast_ml
  
  # isodistance isoline at distance m
  iso_poly <- st_buffer(local_coast, dist = m)
  iso_line <- st_boundary(iso_poly) |> st_collection_extract("LINESTRING")
  if (length(iso_line) == 0) return(NULL)
  
  # clip isoline near patch & pick longest segment
  iso_local <- tryCatch(st_intersection(st_as_sf(iso_line), st_buffer(st_geometry(patch_poly), m + 200)),
                        error = function(e) st_as_sf(iso_line))
  if (nrow(iso_local) == 0) return(NULL)
  cut_line <- iso_local[which.max(st_length(iso_local)), , drop = FALSE]
  
  # ensure it intersects the patch; one retry with bigger clip
  if (!any(st_intersects(st_geometry(cut_line), st_geometry(patch_poly), sparse = FALSE))) {
    iso_local2 <- tryCatch(st_intersection(st_as_sf(iso_line), st_buffer(st_geometry(patch_poly), 3*m + 200)),
                           error = function(e) NULL)
    if (!is.null(iso_local2) && nrow(iso_local2) > 0)
      cut_line <- iso_local2[which.max(st_length(iso_local2)), , drop = FALSE]
  }
  if (!any(st_intersects(st_geometry(cut_line), st_geometry(patch_poly), sparse = FALSE))) return(NULL)
  
  # split
  pieces <- tryCatch({
    st_collection_extract(st_split(st_make_valid(patch_poly), st_geometry(cut_line)), "POLYGON")
  }, error = function(e) NULL)
  if (is.null(pieces) || nrow(pieces) == 0) return(NULL)
  
  # drop slivers
  A_all <- as.numeric(st_area(st_union(patch_poly)))
  pieces$area_m2 <- as.numeric(st_area(pieces))
  pieces <- pieces |> dplyr::filter(area_m2 > max(A_all * sliver_frac_val, tiny_abs_m2_val))
  if (nrow(pieces) == 0) return(NULL)
  
  # label by coast distance vs m (more robust than point-centroid proximity)
  pc <- st_centroid(pieces)
  d_piece <- as.numeric(st_distance(pc, coast_ml))
  pieces$zone <- ifelse(d_piece <= m, "Shallow", "Deep")
  
  # dissolve by zone
  dissolve_clean <- function(x) {
    if (is.null(x) || nrow(x) == 0) return(NULL)
    g <- st_make_valid(st_union(x))
    g <- st_buffer(g, 0)
    g <- st_collection_extract(g, "POLYGON")
    out <- st_as_sf(g)
    if (is.null(out) || nrow(out) == 0) return(NULL)
    out
  }
  shallow_sf <- pieces |> dplyr::filter(zone=="Shallow") |> dissolve_clean()
  deep_sf    <- pieces |> dplyr::filter(zone=="Deep")    |> dissolve_clean()
  
  # if one missing, nudge m a bit
  if (is.null(shallow_sf) || is.null(deep_sf)) {
    for (f in c(0.95, 1.05, 0.90, 1.10)) {
      m_try <- m * f
      pieces$zone <- ifelse(d_piece <= m_try, "Shallow", "Deep")
      shallow_sf <- pieces |> dplyr::filter(zone=="Shallow") |> dissolve_clean()
      deep_sf    <- pieces |> dplyr::filter(zone=="Deep")    |> dissolve_clean()
      if (!is.null(shallow_sf) && !is.null(deep_sf)) break
    }
  }
  if (is.null(shallow_sf) || is.null(deep_sf)) return(NULL)
  
  # fill any gap: patch - (shallow ∪ deep)
  union_sd <- st_union(shallow_sf, deep_sf)
  gap <- suppressWarnings(st_difference(st_geometry(patch_poly), union_sd))
  if (!is.null(gap) && length(gap) > 0 && any(!st_is_empty(gap))) {
    gap_sf <- st_as_sf(st_collection_extract(gap, "POLYGON", warn = FALSE))
    if (!is.null(gap_sf) && nrow(gap_sf) > 0) {
      gc <- st_centroid(gap_sf)
      dg <- as.numeric(st_distance(gc, coast_ml))
      gap_sf$zone <- ifelse(dg <= m, "Shallow", "Deep")
      shallow_sf <- dissolve_clean(rbind(shallow_sf, gap_sf |> dplyr::filter(zone=="Shallow")))
      deep_sf    <- dissolve_clean(rbind(deep_sf,    gap_sf |> dplyr::filter(zone=="Deep")))
    }
  }
  
  # remove overlaps (assign overlap to nearer side)
  ov <- suppressWarnings(st_intersection(st_geometry(shallow_sf), st_geometry(deep_sf)))
  if (!is.null(ov) && length(ov) > 0 && any(!st_is_empty(ov))) {
    oc  <- st_centroid(st_as_sf(ov))
    d_o <- as.numeric(st_distance(oc, coast_ml))
    if (mean(d_o, na.rm = TRUE) <= m) {
      deep_sf <- suppressWarnings(st_difference(deep_sf, ov)) |> dissolve_clean()
    } else {
      shallow_sf <- suppressWarnings(st_difference(shallow_sf, ov)) |> dissolve_clean()
    }
  }
  
  shallow_sf$zone <- "Shallow"
  deep_sf$zone    <- "Deep"
  out <- dplyr::bind_rows(shallow_sf, deep_sf)
  if (is.null(out) || nrow(out) == 0) return(NULL)
  
  out$site       <- site
  out$site_type  <- stype
  out$patch_id   <- pid
  out$patch_type <- unique(patch_poly$patch_type)
  out$mid_dist_m <- m
  out$coast_tol  <- coast_tol_used
  
  out |> dplyr::select(patch_id, patch_type, site, site_type, zone, mid_dist_m, coast_tol, geometry)
}


# --- IMPORTANT: pass arguments explicitly in lapply (no defaults evaluated) ---
split_list <- lapply(
  groups,
  function(g) split_patch_by_offset(
    g_one          = g,
    coast_ml       = coast_ml,
    coast_tol_used = COAST_TOL_DEFAULT,
    m_nudge_val    = M_NUDGE_DEFAULT,
    sliver_frac_val= SLIVER_FRAC_DEFAULT,
    tiny_abs_m2_val= TINY_ABS_M2_DEFAULT
  )
)
split_list   <- split_list[!sapply(split_list, is.null)]
patches_split<- if (length(split_list)) dplyr::bind_rows(split_list) else NULL


# --- APPLY PER (site, site_type, patch_id) ---
groups <- split(zone_points, interaction(zone_points$site,
                                         zone_points$site_type,
                                         zone_points$patch_id,
                                         drop = TRUE))

split_list <- lapply(groups, function(g) split_patch_by_offset(g, coast_ml = coast_ml))
split_list <- split_list[!sapply(split_list, is.null)]
patches_split <- if (length(split_list)) dplyr::bind_rows(split_list) else NULL

message("✅ Groups attempted: ", length(groups),
        " | successful splits: ", ifelse(is.null(patches_split), 0, nrow(patches_split)))

# --- PLOT (Monterey) ---
if (!is.null(patches_split)) {
  ggplot() +
    geom_sf(data = polys_individual, aes(fill = patch_type),
            color = "grey80", alpha = 0.25, show.legend = FALSE) +
    geom_sf(data = patches_split, aes(fill = zone),
            color = "black", alpha = 0.6) +
    geom_sf(data = zone_points, aes(color = zone), size = 1.2) +
    geom_sf(data = st_as_sf(coast_simpl), color = "red", linewidth = 0.25) +
    scale_fill_manual(values = c(Shallow = "#1f78b4", Deep = "#33a02c")) +
    scale_color_manual(values = c(Shallow = "#1f78b4", Deep = "#33a02c")) +
    coord_sf(xlim = c(mbb["xmin"], mbb["xmax"]), ylim = c(mbb["ymin"], mbb["ymax"])) +
    theme_minimal() +
    labs(title = "Patch splits by coastline-parallel midpoint offset",
         subtitle = sprintf("Coastline simplified (tol = %dm). Coverage repaired; slivers dropped.", coast_tol),
         fill = "Zone", color = "Zone")
} else {
  message("⚠️ No split patches were created")
}
