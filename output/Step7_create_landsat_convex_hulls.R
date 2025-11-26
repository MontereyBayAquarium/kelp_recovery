
#Joshua G. Smith; jogsmith@ucsc.edu

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap, scales, readxl, smoothr,
                 purrr, lwgeom)

#set directories 
basedir <- "/Volumes/enhydra/data/kelp_recovery/"
localdir <- "/Users/jossmith/Documents/Data/landsat"
figdir <- here::here("figures")
output <- here::here("output")

#read landsat 
#clusters <- st_read(file.path(output, "/landsat/processed/named_clusters.geojson")) #old clustering
clusters <- st_read(file.path(output, "/landsat/processed/named_clustersv2.geojson"))

#read state polyline
coast_path  <- "/Volumes/enhydra/data/kelp_recovery/gis_data/raw/Coastn83/coastn83.shp"

#read state
ca_counties <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

#multivariate patch types
load(here::here("output","lda_patch_transitionsv3.rda"))

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")

################################################################################
#plot individual points and clusters


cluster_coord <- clusters %>%
  filter(year == 2023) %>%
  mutate(longitude = sf::st_coordinates(.)[,1],
         latitude = sf::st_coordinates(.)[,2])%>%
  st_drop_geometry()%>%
  group_by(site_num)%>%
  dplyr::summarize(lat = mean(latitude),
                   long = mean(longitude))

base_theme <-  theme(axis.text=element_text(size=7, color = "black"),
                     axis.title=element_text(size=8,color = "black"),
                     legend.text=element_text(size=7,color = "black"),
                     legend.title=element_text(size=8,color = "black"),
                     plot.tag=element_text(size=10,color = "black"),
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

# Build inset
g1_inset <- ggplotGrob(
  ggplot() +
    # Plot land
    geom_sf(data = foreign, fill = "grey80", color = "white", lwd = 0.3) +
    geom_sf(data = usa, fill = "grey80", color = "white", lwd = 0.3) +
    # Plot box (if needed; you can remove this if your panel border is enough)
    annotate("rect", xmin = -122.6, xmax = -121, ymin = 36.2, ymax = 37.1, 
             color = "black", fill = NA, lwd = 0.6) +
    labs(x = "", y = "") +
    # Crop to desired extent
    coord_sf(xlim = c(-124.5, -117), ylim = c(32.5, 42)) +
    # Start with an empty theme
    theme_void() +
    # Add back just the panel border so that the inset has an outline.
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.6))
)
#g1_inset

p1 <- ggplot() +
  # Add clusters
  geom_sf(data = clusters %>% filter(year == 2023), aes(color = site_type)) +
  scale_color_manual(values = c("Forest" = "#1B9E77", "Barren" = "#7570B3", "Incipient" = "#D95F02"), name = "Site type") +
  geom_sf(data = ca_counties, fill = "gray", color = "gray80") +
  labs(title = "", tag = "") +
  # Add landmarks
  #geom_text(data = monterey_label, mapping = aes(x = x, y = y, label = label),
  #        size = 3, fontface = "bold") +
  # Add CA inset
  annotation_custom(grob = g1_inset, 
                    xmin = -122.01, 
                    xmax = -121.96,
                    ymin = 36.625) +
  # Add jittered labels for site_name with boxes
  #ggrepel::geom_label_repel(
  #  data = cluster_coord,
  #  aes(x = long, y = lat, label = site_num),
  #  box.padding = 0.3,
  #  point.padding = 0.5,
  #  force = 18,
  #  size = 2,
  #  min.segment.length = 0.1,
  #  segment.color = "black"
  #) +
  #add scale bar
  # ggsn::scalebar(x.min = -121.99, x.max = -121.88, 
  #               y.min = 36.519, y.max = 36.645,
  #              #anchor=c(x=-124.7,y=41),
  #             location="bottomright",
  #            dist = 2, dist_unit = "km",
  #           transform=TRUE, 
  #          model = "WGS84",
  #         st.dist=0.02,
  #        st.size=2,
  #       border.size=.5,
  #      height=.02
  #)+
  #add north arrow
  # ggsn::north(x.min = -121.99, x.max = -121.88, 
  #            y.min = 36.519, y.max = 36.65,
  #           location = "topright", 
  #          scale = 0.05, 
  #         symbol = 10)+
  theme_bw() +  theme(
    plot.tag.position = c(-0.03, 1),
    axis.title = element_blank(),
    legend.position = "none") +
  labs(title = "",
       x="",
       y="",
       tag = "A")+
  #theme(axis.text.x = element_blank(),
  #      axis.text.y = element_blank())+
  #guides(fill = guide_legend(override.aes = list(size = 3))) +
  base_theme+
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326) 

p1


###############################################################################
#plot convex hulls


# Create a convex hull for each cluster while keeping the incipient attribute
cluster_polygons <- clusters %>%
  filter(year == 2023) %>%
  group_by(site_num, site_type) %>%              # group by your cluster ID and type
  summarize(geometry = st_combine(geometry), .groups = "drop") %>%  # combine points in each cluster
  st_convex_hull()                                # compute the convex hull
# Compute the convex hull for each group



p2 <- ggplot() +
  # Plot cluster polygons with the fill based on the incipient category
  geom_sf(data = cluster_polygons, 
          aes(fill = site_type),      # fill mapped to incipient
          color = "black")+            # black border for each cluster            # optional transparency
  scale_fill_manual(values = c("Forest" = "#1B9E77", 
                               "Barren" = "#7570B3", 
                               "Incipient" = "#D95F02"), 
                    name = "Site type") +
  # Plot other spatial layers (e.g., county boundaries)
  geom_sf(data = ca_counties, fill = "gray", color = "gray80") +
  # Include the inset map
 # annotation_custom(grob = g1_inset, 
  #                  xmin = -122.01, xmax = -121.96, ymin = 36.625) +
  # Add cluster labels with repelling to avoid overlaps
  ggrepel::geom_label_repel(
    data = cluster_coord,
    aes(x = long, y = lat, label = site_num),
    box.padding = 0.3,
    point.padding = 0.5,
    force = 18,
    size = 2,
    min.segment.length = 0.1,
    segment.color = "black"
  ) +
  # Add scalebar and north arrow (as in your original code)
  #ggsn::scalebar(
  #  x.min = -121.99, x.max = -121.88,  
  #  y.min = 36.519, y.max = 36.645,
  #  location = "bottomright",
  #  dist = 2, dist_unit = "km",
  #  transform = TRUE, model = "WGS84",
  #  st.dist = 0.02, st.size = 4,
  #  border.size = 0.5, height = 0.02
  #) +
  #ggsn::north(
  #  x.min = -121.99, x.max = -121.88,  
  #  y.min = 36.519, y.max = 36.65,
  #  location = "topright", scale = 0.05, symbol = 10
  #) +
  # Apply theme adjustments
  theme_bw() +
  base_theme +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326)+
  theme(
    #axis.text.x = element_blank(),   
    axis.text.y = element_blank(),   
    #axis.ticks.x = element_blank(),  
    axis.ticks.y = element_blank(),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  ) + labs(tag = "B")

p2



################################################################################
#resolve adjacent and overlapping clusters

# Step 1: Project to UTM Zone 10N
cluster_proj <- st_transform(cluster_polygons, 32610)

# Step 2: Create individual 100 m buffers
indiv_buffers <- st_buffer(cluster_proj, dist = 40)

# Step 3: Subtract overlapping original polygons
non_overlap_buffers <- map_dfr(1:nrow(indiv_buffers), function(i) {
  this_buf <- indiv_buffers[i, ]
  other_polys <- st_union(cluster_proj[-i, ])
  trimmed_geom <- st_difference(this_buf$geometry, other_polys)
  
  # Return only valid, non-empty geometries
  if (!st_is_empty(trimmed_geom) && length(trimmed_geom) > 0) {
    st_sf(site_type = cluster_proj$site_type[i], geometry = trimmed_geom)
  } else {
    NULL
  }
})

# Step 4: Transform back to WGS84
non_overlap_buffers_wgs <- st_transform(non_overlap_buffers, 4326)
ca_counties_wgs <- st_transform(ca_counties, 4326)

# Step 5: Plot
ggplot() +
  geom_sf(data = non_overlap_buffers_wgs, aes(fill = site_type, alpha=0.2), color = "black") +
  geom_sf(data = ca_counties_wgs, fill = "gray", color = "gray80") +
  scale_fill_manual(values = c("Forest" = "#1B9E77",
                               "Barren" = "#7570B3",
                               "Incipient" = "#D95F02")) +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645)) +
  theme_minimal() +
  ggtitle("Non-Overlapping 100m Buffers from Cluster Polygons")


# Step 4: Merge adjacent polygons by site_type
non_overlap_merged <- non_overlap_buffers %>%
  group_by(site_type) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()


# Step 5: Transform to WGS84 for plotting
non_overlap_merged_wgs <- st_transform(non_overlap_merged, 4326)


# Step 6: Plot
ggplot() +
  geom_sf(data = non_overlap_merged_wgs, aes(fill = site_type), color = "black") +
  geom_sf(data = ca_counties_wgs, fill = "gray", color = "gray80") +
  scale_fill_manual(values = c("Forest" = "#1B9E77",
                               "Barren" = "#7570B3",
                               "Incipient" = "#D95F02")) +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645)) +
  theme_minimal() +
  ggtitle("Non-Overlapping, Merged 75m Buffers from Cluster Polygons")



# Remove all interior holes from the merged polygons
non_overlap_merged_clean <- nngeo::st_remove_holes(non_overlap_merged)
non_overlap_merged_clean_wgs <- st_transform(non_overlap_merged_clean, 4326)


# Step 6: Plot
ggplot() +
  geom_sf(data = non_overlap_merged_clean_wgs, aes(fill = site_type), color = "black") +
  geom_sf(data = ca_counties_wgs, fill = "gray", color = "gray80") +
  scale_fill_manual(values = c("Forest" = "#1B9E77",
                               "Barren" = "#7570B3",
                               "Incipient" = "#D95F02")) +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645)) +
  theme_minimal() +
  ggtitle("Non-Overlapping, Merged 75m Buffers from Cluster Polygons")


# Step 5.5: Smooth boundaries (UTM projection)
non_overlap_smoothed <- smooth(non_overlap_merged_clean, method = "ksmooth", smoothness = 15)


non_overlap_smoothed_wgs <- st_transform(non_overlap_smoothed, 4326)


p3 <- ggplot() +
  geom_sf(data = non_overlap_smoothed_wgs, aes(fill = site_type), color = "black") +
  geom_sf(data = ca_counties_wgs, fill = "gray", color = "gray80") +
  scale_fill_manual(values = c("Forest" = "#1B9E77",
                               "Barren" = "#7570B3",
                               "Incipient" = "#D95F02"),
                    name = "Patch type") +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645)) +
  theme_bw() +
  base_theme +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326)+
  theme(
    #axis.text.x = element_blank(),   
    axis.text.y = element_blank(),   
    #axis.ticks.x = element_blank(),  
    axis.ticks.y = element_blank(),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank()
  ) +
  labs(tag = "C")


p3


################################################################################
#Figure S2

p <- ggarrange(p1, p2, p3, nrow = 1, align = "h",
               widths = c(0.5,0.5,0.5)) 

p

ggsave(
  filename = here::here("figures", "FigS2_cluster_assignments.png"),
  plot     = p,
  width    = 14,
  height   = 5,
  dpi      = 600,
  bg       = "white"
)


################################################################################
#INTERMEDIATE STEP

#Export shapwfile and load into QGIS for processing.
st_write(non_overlap_smoothed_wgs, file.path(here::here("output","landsat",
                                                        "processed",
                                                        "smoothed_patch_polygons.shp")), delete_layer = TRUE)

#Files required for QGIS (currently live in output folder):
  #- high res coast polyline - coastn83.shp
  #- Polygons for patches (created in this script) - smoothed_patch_polygons.shp
  #- Isobaths - mpen_2m_bathy.shp
  #- site locations - site_table.csv

#steps conducted in QGIS:
  #- Import patch polygons and extend patches to the shoreline
  #- Import site locations as points
  #- Import isobaths
  #- Use split features tool to clip each polygons at the isobath that is 
      #between the deep and shallow sites. This creates individual shallow vs.
      # deep polygons. 
  #- Save the new polygon layer and load into R. 
  #- The next step (in R) is to reassign each patch type based on the predicted
      #patch type from the LDA model in step 6. 









