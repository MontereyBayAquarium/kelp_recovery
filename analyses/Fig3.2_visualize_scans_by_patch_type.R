

rm(list = ls())

################################################################################
# Load packages and set directories
require(librarian)
librarian::shelf(
  tidyverse, janitor, lubridate, sf,
  geodata, viridis, scales, here
)

datdir   <- here::here("output")
poly_dir <- file.path(datdir, "recovery_polygons")

################################################################################
# Read and prepare individual polygon patches (KMZ to KML to separate polygons)
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
# Read and prepare scan data
scan_dat <- read_csv(file.path(datdir, "scans", "scans_data.csv"))
scan_clean <- scan_dat %>%
  clean_names() %>%
  mutate(
    lat            = as.numeric(lat),
    long           = as.numeric(long),
    year           = year(date),
    quarter        = lubridate::quarter(date),
    behavior_group = ifelse(behav == "foraging", "Foraging", "Other")
  ) %>%
  drop_na(lat, long, ind, behavior_group)

################################################################################
# Intersect points with individual polygons
tmp_sf      <- st_as_sf(scan_clean, coords = c("long","lat"), crs = 4326)
pts_in_poly <- st_join(tmp_sf, polys_individual, join = st_within)

################################################################################
# Summarize total individuals per patch × behavior × year × quarter
tpatch_counts <- pts_in_poly %>%
  st_drop_geometry() %>%
  group_by(patch_id, patch_type, behavior_group, year, quarter) %>%
  summarise(total_ind = sum(ind, na.rm = TRUE), .groups = "drop")

################################################################################
# Ensure all polygons appear for all years, quarters, and behaviors
years     <- sort(unique(pts_in_poly$year))
quarters  <- sort(unique(pts_in_poly$quarter))
behaviors <- unique(pts_in_poly$behavior_group)

# xpand grid then bring in patch_type mapping
patch_expand <- expand_grid(
  patch_id       = polys_individual$patch_id,
  year           = years,
  quarter        = quarters,
  behavior_group = behaviors
) %>%
  # add patch_type from polys_individual
  left_join(
    st_drop_geometry(polys_individual) %>% dplyr::select(patch_id, patch_type),
    by = "patch_id"
  )

patch_sf <- patch_expand %>% 
  left_join(tpatch_counts, by = c("patch_id","patch_type","behavior_group","year","quarter")) %>%
  mutate(total_ind = replace_na(total_ind, 0)) %>%
  left_join(polys_individual, by = "patch_id") %>%
  st_as_sf() %>%
  group_by(year, quarter, behavior_group) %>%
  mutate(rel_ind = total_ind / sum(total_ind)) %>%
  ungroup()


################################################################################
# Barplot: percent distribution of foraging counts across patch types per year × quarter
dist_data <- pts_in_poly %>%
  st_drop_geometry() %>%
  filter(!is.na(patch_type), behavior_group == "Foraging") %>%
  group_by(patch_type, year, quarter) %>%
  summarise(total_foraging = sum(ind, na.rm = TRUE), .groups = "drop") %>%
  group_by(year, quarter) %>%
  mutate(pct_foraging = total_foraging / sum(total_foraging)) %>%
  ungroup() %>%
  mutate(period = paste0(year, " Q", quarter))

# Plot stacked horizontal bars per period
p <- ggplot(dist_data, aes(x = pct_foraging, y = forcats::fct_rev(period), fill = patch_type)) +
  geom_col(color = "black", size = 0.2) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_viridis_d(option = "A", end = 0.8, name = "Patch Type") +
  labs(
    title = "Distribution of foraging activity by patch type",
    x     = "Percent of foraging",
    y     = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y      = element_text(size = 12),
    panel.grid.major = element_blank(),
    legend.position  = "right"
  )

p

ggsave(
  filename = here::here("figures", "Fig15_foraging_scan_patch_type.png"),
  plot     = p,
  width    = 6,
  height   = 10,
  dpi      = 600,
  bg = "white"
)




