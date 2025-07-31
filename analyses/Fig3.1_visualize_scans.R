rm(list = ls())

################################################################################
#Load packages and set directories

require(librarian)
librarian::shelf(
  tidyverse, janitor, lubridate, sf,
  geodata, viridis, scales, here
)

datdir   <- here::here("output")
poly_dir <- file.path(datdir, "recovery_polygons")

################################################################################
#Read and prepare individual polygon patches (KMZ to KML to separate polygons)
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

################################################################################
#Read and prepare scan data
scan_dat <- read_csv(file.path(datdir, "scans", "scans_data.csv"))
scan_clean <- scan_dat %>%
  clean_names() %>%
  mutate(
    lat            = as.numeric(lat),
    long           = as.numeric(long),
    year           = year(date),
    quarter        = lubridate::quarter(date),
    period         = paste0(year, " Q", quarter),
    behavior_group = ifelse(behav == "foraging", "Foraging", "Other")
  ) %>%
  drop_na(lat, long, ind, behavior_group)

################################################################################
#Intersect points with individual polygons
tmp_sf      <- st_as_sf(scan_clean, coords = c("long","lat"), crs = 4326)
pts_in_poly <- st_join(tmp_sf, polys_individual, join = st_within)

################################################################################
#Summarize total individuals per polygon patch × behavior × year & quarter
tpatch_counts <- pts_in_poly %>%
  st_drop_geometry() %>%
  group_by(patch_id, patch_type, behavior_group, year, quarter) %>%
  summarise(total_ind = sum(ind, na.rm = TRUE), .groups = "drop")

#Ensure all polygons appear for all years, quarters, and behaviors
years      <- unique(pts_in_poly$year)
quarters   <- unique(pts_in_poly$quarter)
behaviors  <- unique(pts_in_poly$behavior_group)
patch_expand <- expand_grid(
  patch_id       = polys_individual$patch_id,
  year           = years,
  quarter        = quarters,
  behavior_group = behaviors
)

patch_sf <- patch_expand %>%
  left_join(tpatch_counts, by = c("patch_id", "behavior_group", "year", "quarter")) %>%
  mutate(total_ind = replace_na(total_ind, 0)) %>%
  left_join(polys_individual, by = "patch_id") %>%
  st_as_sf() %>%
  # compute relative proportion per year, quarter, and behavior
  group_by(year, quarter, behavior_group) %>%
  mutate(rel_ind = total_ind / sum(total_ind)) %>%
  ungroup()

#Fetch high-resolution land mask clipped to expanded bbox
usa_gadm    <- geodata::gadm(country = "USA", level = 0, path = tempdir())
usa_sf      <- st_as_sf(usa_gadm) %>% st_make_valid()

#Compute and expand bounding box of polygons
tmp_bbox    <- st_bbox(polys_individual)
xpad        <- (tmp_bbox["xmax"] - tmp_bbox["xmin"]) * 0.05
ypad        <- (tmp_bbox["ymax"] - tmp_bbox["ymin"]) * 0.05
bbox_expanded <- tmp_bbox + c(-xpad, -ypad, xpad, ypad)
#Clip land to that expanded bbox
bbox_poly   <- st_as_sfc(bbox_expanded) %>% st_set_crs(st_crs(polys_individual))
land_mask   <- st_intersection(usa_sf, bbox_poly)

#Prepare Foraging-only data and plotting
plot_sf <- patch_sf %>% filter(behavior_group == "Foraging") %>%
  mutate(period = paste(year, quarter, sep = " Q"))

#replicate land for each facet year and quarter
years    <- sort(unique(plot_sf$year))
quarters <- sort(unique(plot_sf$quarter))
# Use tidyr::crossing to expand land for each year/quarter combination
land_plot <- land_mask %>%
  st_sf() %>%
  tidyr::crossing(year = years, quarter = quarters)

################################################################################
#Plot 

p <- ggplot() +
  # land beneath, has year and quarter
  geom_sf(
    data        = land_plot,
    aes(geometry = geometry),
    fill        = "gray90",
    color       = NA,
    inherit.aes = FALSE
  ) +
  # patch polygons colored by relative proportion
  geom_sf(
    data = plot_sf,
    aes(fill = rel_ind, geometry = geometry),
    color = "black",
    size  = 0.3,
    inherit.aes = FALSE
  ) +
  # custom fill: 0 gray, positives viridis
  scale_fill_gradientn(
    colors = c("gray80", viridis::viridis(100, option = "A")),
    values = scales::rescale(c(0, 1e-6, 1)),
    name   = "Proportion",
    breaks  = c(0.01, 0.25, 0.50, 0.75, 1),
    labels  = scales::percent_format(accuracy = 1)
  ) +
  facet_grid(year ~ quarter, switch = "y") +
  coord_sf(
    xlim   = c(bbox_expanded["xmin"], bbox_expanded["xmax"]),
    ylim   = c(bbox_expanded["ymin"], bbox_expanded["ymax"]),
    expand = FALSE,
    clip   = "on"
  ) +
  labs(
    title    = "Relative foraging activity across patches",
    subtitle = "By year and quarter",
    x        = NULL,
    y        = NULL
  ) +
  theme_void(base_size = 14) +
  theme(
    strip.placement     = "outside",
    strip.text.y.left   = element_text(angle = 90, hjust = 0.5, vjust = 0.5,
                                       face = "bold"),
    strip.text.x        = element_text(margin = margin(t = 5), face = "bold"),
    panel.grid          = element_blank(),
    axis.text           = element_blank(),
    axis.ticks          = element_blank(),
    legend.position     = "right"
  )

################################################################################
#Export

ggsave(
  filename = here::here("figures", "Fig14_foraging_scan_relative_activity.png"),
  plot     = p,
  width    = 8,
  height   = 10,
  dpi      = 600,
  bg = "white"
)

