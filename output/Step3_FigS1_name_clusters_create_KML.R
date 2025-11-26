
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap)

#set directories 
basedir <- "/Volumes/enhydra/data/kelp_recovery/"

figdir <- here::here("figures")
output <- here::here("output")

#read landsat 
final_data <- st_read(file.path(output, "/landsat/processed/kelp_area_by_cluster2.geojson")) %>%
  mutate(site_numeric = as.numeric(cluster)) %>%
  st_cast("POINT")

#read state
ca_counties <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")



###############################################################################
#clean up spatial extent 

#set bouding box for study area
xlims <- c(-121.937462, -121.935460, -121.996530, -121.996530)
ylims <- c(36.640507, 36.579086, 36.581153, 36.641310)

box_coords <- tibble(x = xlims, y = ylims) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  st_set_crs(st_crs(final_data))

#get the bounding box of the two x & y coordintates, make sfc
bounding_box <- st_bbox(box_coords) %>% st_as_sfc()
plot(bounding_box)

# Filter out the points falling within the bounding box
landsat_build1 <- st_difference(final_data, bounding_box) 


#clean up
xlims <- c(-121.891489, -122.010802, -121.991380, -121.891489)
ylims <- c(36.651895, 36.658547, 36.519945, 36.519945)


box_coords <- tibble(x = xlims, y = ylims) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  st_set_crs(st_crs(landsat_build1))

#get the bounding box of the two x & y coordintates, make sfc
bounding_box <- st_bbox(box_coords) %>% st_as_sfc()

# Filter out the points falling within the bounding box
landsat_build2 <- st_intersection(landsat_build1, bounding_box) 


plot(landsat_build2 %>% filter(year == 2023))

################################################################################
#rename clusters

#rename sites and create table
landsat_build3 <- landsat_build2 %>%
  mutate(
    #set cluster order
    site_num = case_when(
      cluster == 79 ~ 1,
      cluster == 78 ~ 2,
      cluster == 77 ~ 3,
      cluster == 76 ~ 4,
      cluster == 75 ~ 5,
      cluster == 74 ~ 6,
      cluster == 73 ~ 7,
      cluster == 72 ~ 8, 
      cluster == 71 ~ 9,
      cluster == 70 ~ 10,
      cluster == 69 ~ 11,
      cluster == 68 ~ 12,
      cluster == 67 ~ 13,
      cluster == 66 ~ 14,
      cluster == 64 ~ 15,
      cluster == 101 ~ 16,
      cluster == 100 ~ 17,
      cluster == 45 ~ 18,
      cluster == 44 ~ 19,
      cluster == 1 ~ 20,
      cluster == 2 ~ 21,
      cluster == 4 ~ 22,
      cluster == 85 ~ 23,
      cluster == 86 ~ 24,
      cluster == 6 ~ 25,
      cluster == 7 ~ 26,
      cluster == 103 ~ 27,
      cluster == 14 ~ 28,
      cluster == 11 ~ 29,
      cluster == 16 ~ 30,
      cluster == 17 ~ 31,
      cluster ==22 ~ 32,
      cluster == 21 ~ 33,
      cluster == 9 ~ 34,
      cluster == 88 ~ 35,
      cluster == 23 ~ 36,
      cluster == 24 ~ 37,
      cluster == 119 ~ 38,
      cluster == 28 ~ 39,
      cluster == 96 ~ 40,
      cluster == 97 ~ 41,
      cluster == 31 ~ 42,
      cluster == 33 ~ 43,
      cluster == 37 ~ 44,
      cluster == 39 ~ 45,
      cluster == 40 ~ 46,
      cluster == 43 ~ 47,
      cluster == 51 ~ 48,
      cluster == 56 ~ 49,
      cluster == 47 ~ 50,
      cluster == 46 ~ 51,
      cluster == 55 ~ 52,
      cluster == 59 ~ 53,
      cluster == 99 ~ 54,
      cluster == 53 ~ 55,
      cluster == 60 ~ 56,
      cluster == 50 ~ 57,
      cluster == 54 ~ 58,
      cluster == 62 ~ 59,
      cluster == 49 ~ 60,
      cluster == 52 ~ 61,
      cluster == 58 ~ 62,
      cluster == 61 ~ 63,
      cluster == 117 ~ 64,
      cluster == 65 ~ 65,
      cluster == 63 ~ 66,
      cluster == 57 ~ 67,
      cluster == 48 ~ 68,
      cluster == 41 ~ 69,
      cluster == 38 ~ 70,
      cluster == 98 ~ 71,
      cluster == 123 ~ 72,
      cluster == 110 ~ 73,
      TRUE ~ NA)
  ) %>%
  #filter to data extent
  filter(!(is.na(cluster)))


################################################################################
#plot timeseries for each cluster to inspect

# Theme
base_theme <-  theme(axis.text=element_text(size=7, color = "black"),
                     axis.title=element_text(size=8,color = "black"),
                     legend.text=element_text(size=7,color = "black"),
                     legend.title=element_text(size=8,color = "black"),
                     plot.tag=element_text(size=8,color = "black"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key = element_rect(fill=alpha('blue', 0)),
                     legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=6, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())


# Plot 
ggplot(landsat_build3, aes(x = year, y = perc_of_baseline)) +
  geom_point() + 
  geom_smooth(se = TRUE) +
  facet_wrap(~ site_num, scales = "free_y") + theme_bw() + base_theme


################################################################################
#determine site type

# Updated classification function using a threshold of 5% of the baseline mean.
# Any measurement less than 0.1 is considered equivalent to 0.
#classify_site <- function(values, threshold = 0.1) {
#  if (!any(values < threshold)) { 
#    # If none of the values are below 1, then the site is persistent.
#    return("Forest")
#  } else {
#    # Get the index of the first occurrence where the measurement is below the threshold.
#    first_low <- which(values < threshold)[1]
#    # If all values from that point onward are below the threshold, classify as barren.
#    if (all(values[first_low:length(values)] < threshold)) {
#      return("Barren")
#    } else {
#      # Otherwise, the site went low and then recovered, so classify as incipient.
#      return("Incipient")
#    }
#  }
#}

classify_site <- function(values, threshold = 0.05) {
  if (!any(values < threshold)) {
    return("Forest")
  }
  
  below_thresh <- values < threshold
  
  # Find runs of values above/below threshold
  runs <- rle(below_thresh)
  
  # Count the number of below-threshold periods (TRUE runs)
  n_below_periods <- sum(runs$values)
  
  if (n_below_periods > 1) {
    return("Barren")  # dropped, recovered, dropped again
  } else {
    # One drop and stayed above afterward
    # Check if that drop was followed by a sustained recovery
    drop_idx <- which(below_thresh)
    last_drop <- max(drop_idx)
    
    # If all values after the last drop are above threshold, it's incipient
    if (all(values[(last_drop + 1):length(values)] >= threshold)) {
      return("Incipient")
    } else {
      return("Barren")  # dropped and didn't sustain recovery
    }
  }
}




# Apply the function to data (for years after 2013) to summarize each site.
site_class <- landsat_build3 %>% 
  filter(year > 2013) %>% 
  group_by(site_num) %>% 
  arrange(year) %>% 
  summarize(site_type = classify_site(perc_of_baseline)) %>% 
  ungroup() %>% st_drop_geometry()

# Inspect the classification summary
print(site_class)

# Merge the classification back into the main dataset
landsat_classified <- landsat_build3 %>%
  #filter(year > 2013) %>%
  left_join(site_class, by = "site_num")

################################################################################
#Figure S1

my_theme <- theme(
  axis.text.x      = element_text(size = 8, color = "black"),
  axis.text.y      = element_text(size = 8, color = "black"),
  axis.title       = element_text(size = 12, color = "black"),
  legend.text      = element_text(size = 12, color = "black"),
  legend.title     = element_text(size = 12, color = "black"),
  plot.tag         = element_text(size = 10, color = "black"),
  panel.grid       = element_blank(),
  panel.background = element_blank(),
  axis.line        = element_line(colour = "black"),
  legend.background = element_rect(
    fill = scales::alpha("blue", 0)
  ),
  legend.key = element_rect(fill = "transparent", color = NA),
  legend.position = "top",
  strip.text       = element_text(
    size = 8, face = "plain", color = "black", hjust = 0
  ),
  strip.background = element_blank()
)

# Plot percent of max for each site, colored by site type
# Patch-state colors (same as other figs)
patch_cols <- c(
  "Barren"    = "#7570B3",
  "Incipient" = "#D95F02",
  "Forest"    = "#1B9E77"
)

g <- ggplot(
  landsat_classified,
  aes(x = year, y = perc_of_baseline, color = site_type)
) +
  geom_point(alpha = 0.8, size = 1.1) +
  geom_smooth(se = FALSE, linewidth = 0.6) +
  facet_wrap(~ site_num, scales = "free_y") +
  scale_color_manual(
    values = patch_cols,
    breaks = c("Barren", "Incipient", "Forest"),
    name   = "Patch state"
  ) +
  #add heatwave period
  annotate("rect",
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf,
           fill = "indianred", alpha = 0.3) +
  theme_bw() +
  labs(
    x = "Year",
    y = "Proportion of baseline",
    color = "Patch state"
  ) + my_theme 
g

# Fig 2: Patch state habitat correlates
ggsave(
  filename = here::here("figures", "FigS1_patch_canopy_trends.png"),
  plot     = g,
  width    = 14,
  height   = 8,
  dpi      = 600,
  bg       = "white"
)


################################################################################
#Export named clusters

# Export landsat_build4 as a GeoJSON file
st_write(landsat_classified, file.path(output, "landsat/processed/named_clustersv2.geojson"), delete_dsn = TRUE)


################################################################################
#prep for export to KML to visualize in Google maps

##OLD

# transform landsat data to Teale Albers

landsat_build5 <- landsat_classified %>% st_as_sf(crs = 4326) %>% filter(year == 2023)


#Build barren layer
rast_build1 <- st_transform(landsat_build5, crs = 3310) %>% filter (site_type == "Barren") 
r <- rast(rast_build1, res=30)
landsat_rast_barren <- rasterize(rast_build1, r)
plot(landsat_rast_barren)

#Build forest layer
rast_build1 <- st_transform(landsat_build5, crs = 3310) %>% filter (site_type == "Forest") 
r <- rast(rast_build1, res=30)
landsat_rast_forest <- rasterize(rast_build1, r)
plot(landsat_rast_forest)

#Build Incipient
rast_build1 <- st_transform(landsat_build5, crs = 3310) %>% filter (site_type == "Incipient") 
r <- rast(rast_build1, res=30)
landsat_rast_incipient <- rasterize(rast_build1, r)
plot(landsat_rast_incipient)

#prep KML files
rast_my_spat <- raster(landsat_rast_barren)
poly_rast <- rasterToPolygons(rast_my_spat)
barren_shape <- shapefile(poly_rast, file.path(basedir,"kelp_landsat/processed/monterey_peninsula/ArcGIS_files/persistent_barrens.shp"), overwrite=TRUE)
barren_shp <- st_read(file.path(basedir,"kelp_landsat/processed/monterey_peninsula/ArcGIS_files/persistent_barrens.shp"))
dissolved_barren<- st_union(barren_shp)
st_write(dissolved_barren, file.path(basedir, "/kelp_landsat/processed/monterey_peninsula/ArcGIS_files/persistent_barrens.kml"), driver = "KML", append=TRUE)


#prep KML files
rast_my_spat <- raster(landsat_rast_forest)
poly_rast <- rasterToPolygons(rast_my_spat)
forest_shape <- shapefile(poly_rast, file.path(basedir,"kelp_landsat/processed/monterey_peninsula/ArcGIS_files/persistent_forests.shp"), overwrite=TRUE)
forest_shp <- st_read(file.path(basedir,"kelp_landsat/processed/monterey_peninsula/ArcGIS_files/persistent_forests.shp"))
dissolved_forest<- st_union(forest_shp)
st_write(dissolved_forest, file.path(basedir, "/kelp_landsat/processed/monterey_peninsula/ArcGIS_files/persistent_forests.kml"), driver = "KML", append=TRUE)


#prep KML files
rast_my_spat <- raster(landsat_rast_incipient)
poly_rast <- rasterToPolygons(rast_my_spat)
incipient_shape <- shapefile(poly_rast, file.path(basedir,"kelp_landsat/processed/monterey_peninsula/ArcGIS_files/incipient_forests.shp"), overwrite=TRUE)
incipient_shp <- st_read(file.path(basedir,"kelp_landsat/processed/monterey_peninsula/ArcGIS_files/incipient_forests.shp"))
dissolved_incipient<- st_union(incipient_shp)
st_write(dissolved_incipient, file.path(basedir, "/kelp_landsat/processed/monterey_peninsula/ArcGIS_files/incipient_forests.kml"), driver = "KML", append=TRUE)

