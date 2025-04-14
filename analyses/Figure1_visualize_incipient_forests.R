
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
localdir <- "/Users/jossmith/Documents/Data/landsat"
figdir <- here::here("figures")
output <- here::here("output")

#read landsat 
clusters <- st_read(file.path(output, "/landsat/processed/named_clusters.geojson"))

#read state
ca_counties <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")



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


################################################################################
#determine center coords (THIS WILL BE IMPORTANT TO SAVE LATER)

cluster_coord <- clusters %>%
  filter(year == 2023) %>%
  mutate(longitude = sf::st_coordinates(.)[,1],
         latitude = sf::st_coordinates(.)[,2])%>%
  st_drop_geometry()%>%
  group_by(site_num)%>%
  dplyr::summarize(lat = mean(latitude),
                   long = mean(longitude))%>%
  ungroup() %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)


#inspect

ggplot(data = cluster_coord) +
  #geom_sf() +  # Plot the spatial data
  geom_sf(data = ca_counties, fill = "gray", color = "gray80") +
  geom_sf_text(aes(label = site_num), size = 4, color = "black")+
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326) 


ggplot() +
  geom_sf(data = clusters, aes(color = incipient)) +
  scale_color_manual(values = c("Forest" = "forestgreen", "Barren" = "purple", "Incipient" = "orange"), name = "Incipient") +
  coord_sf(crs = 4326) +
  theme_minimal() 



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


# Build inset
g1_inset <-  ggplotGrob(
  ggplot() +
    # Plot land
    geom_sf(data=foreign, fill="grey80", color="white", lwd=0.3) +
    geom_sf(data=usa, fill="grey80", color="white", lwd=0.3) +
    # Plot box
    annotate("rect", xmin=-122.6, xmax=-121, ymin=36.2, ymax=37.1, color="black", fill=NA, lwd=0.6) +
    # Label regions
    #geom_text(data=region_labels, mapping=aes(y=lat_dd, label=region), x= -124.4, hjust=0, size=2) +
    # Labels
    labs(x="", y="") +
    # Crop
    coord_sf(xlim = c(-124.5, -117), ylim = c(32.5, 42)) +
    # Theme
    theme_bw() + base_theme +
    theme( plot.margin = unit(rep(0, 4), "null"),
           panel.margin = unit(rep(0, 4), "null"),
           panel.background = element_rect(fill='transparent'), #transparent panel bg
           # plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
           axis.ticks = element_blank(),
           axis.ticks.length = unit(0, "null"),
           axis.ticks.margin = unit(0, "null"),
           axis.text = element_blank(),
           axis.title=element_blank(),
           axis.text.y = element_blank())
)
#g1_inset

p1 <- ggplot() +
  # Add clusters
  geom_sf(data = clusters %>% filter(year == 2023), aes(color = incipient)) +
  scale_color_manual(values = c("Forest" = "forestgreen", "Barren" = "purple", "Incipient" = "orange"), name = "Site type") +
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
  #add scale bar
  ggsn::scalebar(x.min = -121.99, x.max = -121.88, 
                 y.min = 36.519, y.max = 36.645,
                 #anchor=c(x=-124.7,y=41),
                 location="bottomright",
                 dist = 2, dist_unit = "km",
                 transform=TRUE, 
                 model = "WGS84",
                 st.dist=0.02,
                 st.size=2,
                 border.size=.5,
                 height=.02
  )+
  #add north arrow
  ggsn::north(x.min = -121.99, x.max = -121.88, 
              y.min = 36.519, y.max = 36.65,
              location = "topright", 
              scale = 0.05, 
              symbol = 10)+
  theme_bw() +  theme(
    plot.tag.position = c(-0.03, 1),
    axis.title = element_blank()) +
  labs(title = "",
       x="",
       y="")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())+
  #guides(fill = guide_legend(override.aes = list(size = 3))) +
  base_theme+
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326) 

p1




################################################################################
#plot convex hulls


# Create a convex hull for each cluster while keeping the incipient attribute
cluster_polygons <- clusters %>%
  filter(year == 2023) %>%
  group_by(site_num, incipient) %>%              # group by your cluster ID and type
  summarize(geometry = st_combine(geometry), .groups = "drop") %>%  # combine points in each cluster
  st_convex_hull()                                # compute the convex hull
# Compute the convex hull for each group


p2 <- ggplot() +
  # Plot cluster polygons with the fill based on the incipient category
  geom_sf(data = cluster_polygons, 
          aes(fill = incipient),      # fill mapped to incipient
          color = "black",            # black border for each cluster
          size = 0.5,                 # adjust line width as needed
          alpha = 0.7) +              # optional transparency
  scale_fill_manual(values = c("Forest" = "forestgreen", 
                               "Barren" = "purple", 
                               "Incipient" = "orange"), 
                    name = "Site type") +
  # Plot other spatial layers (e.g., county boundaries)
  geom_sf(data = ca_counties, fill = "gray", color = "gray80") +
  # Include the inset map
  annotation_custom(grob = g1_inset, 
                    xmin = -122.01, xmax = -121.96, ymin = 36.625) +
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
  ggsn::scalebar(
    x.min = -121.99, x.max = -121.88,  
    y.min = 36.519, y.max = 36.645,
    location = "bottomright",
    dist = 2, dist_unit = "km",
    transform = TRUE, model = "WGS84",
    st.dist = 0.02, st.size = 2,
    border.size = 0.5, height = 0.02
  ) +
  ggsn::north(
    x.min = -121.99, x.max = -121.88,  
    y.min = 36.519, y.max = 36.65,
    location = "topright", scale = 0.05, symbol = 10
  ) +
  # Apply theme adjustments
  theme_bw() +
  base_theme +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326)

p2






################################################################################
#plot everything


#ggsave(p1,  filename=file.path(figdir, "Cluster_site_type_map.png"), width = 7.5, height = 9.5, units = "in",
#      bg = "white", dpi = 600)


# Theme
base_theme <-  theme(axis.text.x=element_text(size=10, color = "black"),
                     axis.text.y=element_text(size=9, color = "black"),
                     axis.title=element_text(size=12,color = "black"),
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
                     strip.text = element_text(size=12, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())


# Plot 


# Define the sequence of years to label every 3 years
years_to_label <- seq(2014, max(clusters$year), by = 3)

#Make sure to turn off year filter for landsay_build4
p3 <- ggplot(clusters %>% filter(year > 2013), aes(x = year, y = perc_of_max_3,
                                                         color = incipient)) +
  geom_point() + 
  geom_smooth(se = TRUE) +
  scale_color_manual(values = c("Forest" = "forestgreen", "Barren" = "purple", "Incipient" = "orange"), name = "Site type") +
  labs(y = "Percent of max (relative to 2009-2013)",
       title = "")+
  facet_wrap(~ site_num, scales = "free_y") + 
  theme_bw() + 
  base_theme +
  scale_x_continuous(breaks = years_to_label)  # Setting breaks every 3 years

p3


#ggsave(p2,  filename=file.path(figdir, "Cluster_timeseries.png"), width = 16, height = 10, units = "in",
#      bg = "white", dpi = 600)


################################################################################
#resolve adjacent clustering

# Transform to a projected CRS (e.g., UTM zone 10N)
cluster_polygons_proj <- st_transform(cluster_polygons, 32610)

# Function to join polygons (for one incipient type) using a 50-meter buffer
join_polygons <- function(polys, distance = 50) {
  # Buffer by 50 m
  buffered <- st_buffer(polys, distance)
  # Union the buffered geometries; this dissolves adjacent features
  unioned <- st_union(buffered)
  # If unioned results in disjoint pieces, break them into individual polygons
  unioned_polys <- st_cast(unioned, "POLYGON")
  # Compute convex hulls for neat boundaries
  hulls <- st_convex_hull(unioned_polys)
  return(hulls)
}

# Split the polygons by incipient type and process each group
joined_list <- cluster_polygons_proj %>% 
  split(.$incipient) %>% 
  map(~ join_polygons(.x, distance = 50))

# Add the incipient attribute back for each piece (imap gives you both the object and its name)
joined_list <- imap(joined_list, ~ st_sf(incipient = .y, geometry = .x))

# Combine the resulting polygons into one sf object
joined_polygons <- do.call(rbind, joined_list)

# Optionally, transform back to WGS84 for mapping
joined_polygons <- st_transform(joined_polygons, st_crs(cluster_polygons))


p4 <- ggplot() +
  # Plot the joined patches with fill by incipient type
  geom_sf(data = joined_polygons, aes(fill = incipient), color = "black", size = 0.5
          #, alpha = 0.7
          ) +
  scale_fill_manual(values = c("Forest" = "forestgreen", 
                               "Barren" = "purple", 
                               "Incipient" = "orange"),
                    name = "Site type") +
  # Add additional layers (e.g., county boundaries, inset map, labels)
  geom_sf(data = ca_counties, fill = "gray", color = "gray80") +
  annotation_custom(grob = g1_inset, xmin = -122.01, xmax = -121.96, ymin = 36.625) +
 # ggrepel::geom_label_repel(
#    data = cluster_coord,
#    aes(x = long, y = lat, label = site_num),
#    box.padding = 0.3, point.padding = 0.5, force = 18,
#    size = 2, min.segment.length = 0.1, segment.color = "black"
#  ) +
  ggsn::scalebar(x.min = -121.99, x.max = -121.88,  
                 y.min = 36.519, y.max = 36.645,
                 location = "bottomright",
                 dist = 2, dist_unit = "km",
                 transform = TRUE, model = "WGS84",
                 st.dist = 0.02, st.size = 2,
                 border.size = 0.5, height = 0.02
  ) +
  ggsn::north(x.min = -121.99, x.max = -121.88,  
              y.min = 36.519, y.max = 36.65,
              location = "topright", scale = 0.05, symbol = 10
  ) +
  theme_bw() +
  base_theme +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326)

p4


################################################################################
#reassign clusters and plot new patch timeseries


# ===== STEP 0: Ensure joined_polygons Has the Required Attributes =====

# Check the column names of joined_polygons
print("Columns in joined_polygons:")
print(names(joined_polygons))

# If new_cluster is missing, create it (assigning a unique row number to each polygon)
if (!"new_cluster" %in% names(joined_polygons)) {
  cat("new_cluster not found in joined_polygons. Creating it...\n")
  joined_polygons <- joined_polygons %>% 
    mutate(new_cluster = row_number())
}

# Likewise, ensure there is an incipient column
if (!"incipient" %in% names(joined_polygons)) {
  if ("Incipient" %in% names(joined_polygons)) {
    joined_polygons <- joined_polygons %>% rename(incipient = Incipient)
  } else {
    stop("Error: Neither 'incipient' nor 'Incipient' exist in joined_polygons.")
  }
}

# ===== STEP 1: Prepare Clusters for the Spatial Join =====
# If clusters already has an "incipient" column that conflicts, remove it.
if ("incipient" %in% names(clusters)) {
  clusters_clean <- clusters %>% select(-incipient)
} else {
  clusters_clean <- clusters
}

# ===== STEP 2: Perform the Spatial Join =====
# Use st_intersects to include features that may touch the boundary.
clusters_joined <- st_join(clusters_clean, 
                           joined_polygons[, c("new_cluster", "incipient")],
                           join = st_intersects)

# ===== STEP 3: Handle Duplicate Columns (if present) =====
# Sometimes the join creates duplicate columns (e.g., incipient.x, incipient.y)
if (all(c("incipient.x", "incipient.y") %in% names(clusters_joined))) {
  clusters_joined <- clusters_joined %>%
    mutate(incipient = coalesce(incipient.x, incipient.y)) %>%
    select(-incipient.x, -incipient.y)
}

# (Optional) Inspect the resulting joined dataset:
print("Columns in clusters_joined:")
print(names(clusters_joined))
print(head(st_drop_geometry(clusters_joined)))

# ===== STEP 4: Aggregate the Time Series Data =====
# For each resolved cluster (new_cluster) and each year, calculate the mean of perc_of_max_3.
agg_ts <- clusters_joined %>%
  group_by(new_cluster, year) %>%
  summarise(mean_perc = mean(perc_of_max_3, na.rm = TRUE),
            incipient = first(incipient)) %>%
  ungroup()

# Optionally inspect aggregated data:
print(head(agg_ts))

# ===== STEP 5: Plot the Aggregated Time Series =====
p2_agg <- ggplot(agg_ts %>% filter(year > 2013),
                 aes(x = year, y = mean_perc, color = incipient)) +
  geom_point() +
  geom_smooth(se = TRUE) +
  scale_color_manual(values = c("Forest" = "forestgreen",
                                "Barren" = "purple",
                                "Incipient" = "orange"),
                     name = "Site type") +
  labs(y = "Percent of max (relative to 2009-2013)",
       title = "Aggregated Time Series for Resolved Clusters") +
  facet_wrap(~ new_cluster, scales = "free_y") +
  theme_bw() +
  base_theme +               # Ensure base_theme is defined in your workspace
  scale_x_continuous(breaks = years_to_label)  # Ensure years_to_label is defined

# Display the plot
print(p2_agg)


################################################################################
#plot final cluster map with new labels

### STEP 1: Compute Cluster Centroids and Jitter their Positions

# Assume joined_polygons is your resolved clusters sf object with attributes "new_cluster" and "incipient"
# Transform to a metric CRS (e.g., UTM Zone 10N: EPSG:32610) for accurate distance calculations
joined_polygons_metric <- st_transform(joined_polygons, 32610)

# Compute centroids for each resolved cluster
centroids <- st_centroid(joined_polygons_metric)

# Get the coordinates as a numeric matrix
coords <- st_coordinates(centroids)

# Set a seed for reproducibility and define jitter amount (in meters)
set.seed(123)
jitter_amount <- 20  # Adjust as needed

# Apply random jitter to the coordinates
jittered_coords <- coords + cbind(rnorm(nrow(coords), mean = 0, sd = jitter_amount),
                                  rnorm(nrow(coords), mean = 0, sd = jitter_amount))

# Create a new sf object from the jittered coordinates while preserving the original attributes
jittered_labels <- st_as_sf(as.data.frame(jittered_coords), coords = c("X", "Y"), crs = 32610)
# Copy the new_cluster and incipient attributes from the centroids
jittered_labels$new_cluster <- centroids$new_cluster
jittered_labels$incipient   <- centroids$incipient

# Transform the jittered label positions back to the original CRS (assumed to be the same as joined_polygons)
jittered_labels <- st_transform(jittered_labels, st_crs(joined_polygons))


### STEP 2: Create/Update Your Map with the New Labels

# Suppose p_map is your original map created using joined_polygons. We'll add the jittered labels.
p_map <- ggplot() +
  # Plot your resolved clusters (polygons) colored by incipient type
  geom_sf(data = joined_polygons, aes(fill = incipient), color = "black", size = 0.5, alpha = 0.7) +
  scale_fill_manual(values = c("Forest" = "forestgreen",
                               "Barren" = "purple",
                               "Incipient" = "orange"),
                    name = "Site type") +
  # Add county boundaries for context
  geom_sf(data = ca_counties, fill = "gray", color = "gray80") +
  # Add your inset map (assuming g1_inset is defined)
  annotation_custom(grob = g1_inset, xmin = -122.01, xmax = -121.96, ymin = 36.625) +
  # Use ggrepel to add labels for each new cluster at the jittered locations
  ggrepel::geom_label_repel(data = jittered_labels,
                            aes(label = new_cluster, geometry = geometry),
                            stat = "sf_coordinates",  # This tells ggplot2 to extract x and y from the geometry column
                            size = 3, color = "black") +
  # Define the map extent (adjust as needed)
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = st_crs(joined_polygons)) +
  theme_bw() +
  base_theme   # Ensure your base_theme is defined

# Display the final map
print(p_map)



