
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap, scales)

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

ggplot(clusters %>% filter(year > 2013), aes(x = year, y = perc_of_max_3,
                                             color = incipient)) +
  geom_point() + 
  geom_smooth(se = TRUE) +
  facet_wrap(~ site_num, scales = "free_y") + theme_bw() + base_theme


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
          color = "black")+            # black border for each cluster            # optional transparency
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
    st.dist = 0.02, st.size = 4,
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
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326)+
  theme(
    axis.text.x = element_blank(),   
    axis.text.y = element_blank(),   
    axis.ticks.x = element_blank(),  
    axis.ticks.y = element_blank(),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank()
  )

p2



################################################################################
#plot trends by each patch 


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
                 st.dist = 0.02, st.size = 4,
                 border.size = 0.5, height = 0.02
  ) +
  ggsn::north(x.min = -121.99, x.max = -121.88,  
              y.min = 36.519, y.max = 36.65,
              location = "topright", scale = 0.05, symbol = 10
  ) +
  theme_bw() +
  base_theme +
  theme(axis.title = element_blank())+
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326)+
  theme(
    axis.text.x = element_blank(),   
    axis.text.y = element_blank(),   
    axis.ticks.x = element_blank(),  
    axis.ticks.y = element_blank(),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank()
  )

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
  clusters_clean <- clusters %>% dplyr::select(-incipient)
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
p5 <- ggplot(agg_ts %>% filter(year > 2013),
             aes(x = year, y = mean_perc, color = incipient, fill = incipient)) +
  # Add a red rectangle for the 2014-2016 period:
  annotate(geom = "rect", xmin = 2014, xmax = 2016, ymin = -Inf, ymax = Inf, 
           fill = "indianred", alpha = 0.7) +
  geom_point() +
  geom_smooth(se = TRUE) +
  scale_color_manual(values = c("Forest" = "forestgreen",
                                "Barren" = "purple",
                                "Incipient" = "orange"),
                     name = "Site type") +
  scale_fill_manual(values = c("Forest" = "forestgreen",
                               "Barren" = "purple",
                               "Incipient" = "orange"),
                    name = "Site type") +
  labs(y = "Percent of max (relative to 2009-2013)", title = "",
       x = "Year") +
  facet_wrap(~ new_cluster, scales = "free_y") +
  theme_bw() +
  base_theme +               # Ensure base_theme is defined
  scale_x_continuous(breaks = years_to_label) 

p5




################################################################################
#plot final cluster map with new labels


# Create a mapping vector as character (keys as the original new_cluster values, as characters)
mapping <- c("23" = "1",  "11" = "2",  "24" = "3",  "12" = "4",
             "37" = "5",  "13" = "6",  "38" = "7",  "9"  = "8",
             "30" = "9",  "20" = "10", "29" = "11", "7"  = "12",
             "27" = "13", "26" = "14", "25" = "15", "6"  = "16",
             "4"  = "17", "18" = "18", "8"  = "19", "28" = "20",
             "19" = "21", "5"  = "22", "17" = "23", "3"  = "24",
             "36" = "25", "16" = "26", "10" = "27", "35" = "28",
             "22" = "29", "21" = "30", "34" = "31", "33" = "32",
             "31" = "33", "1"  = "34", "32" = "35", "15" = "36",
             "2"  = "37", "14" = "38")

# Recode the new_cluster field in joined_polygons:
joined_polygons <- joined_polygons %>%
  mutate(new_cluster_relabel = recode(as.character(new_cluster), !!!mapping))


# Assuming you already computed centroids in a metric CRS for jittering:
# (This code is from our previous snippet)
joined_polygons_metric <- st_transform(joined_polygons, 32610)
centroids <- st_centroid(joined_polygons_metric)
coords <- st_coordinates(centroids)

set.seed(123)
jitter_amount <- 20  # jitter amount in meters
jittered_coords <- coords + cbind(rnorm(nrow(coords), mean = 0, sd = jitter_amount),
                                  rnorm(nrow(coords), mean = 0, sd = jitter_amount))

jittered_labels <- st_as_sf(as.data.frame(jittered_coords), coords = c("X", "Y"), crs = 32610)
# Copy the new_cluster attribute from centroids (it should match row order)
jittered_labels$new_cluster <- centroids$new_cluster

# Now, recode to the new labels using the same mapping:
jittered_labels <- jittered_labels %>%
  mutate(new_cluster_relabel = recode(as.character(new_cluster), !!!mapping))

# Transform back to your map's CRS (likely WGS84)
jittered_labels <- st_transform(jittered_labels, st_crs(joined_polygons))



p6 <- ggplot() +
  # Plot your resolved clusters (polygons) with fill based on incipient
  geom_sf(data = joined_polygons, aes(fill = incipient), color = "black", size = 0.5) +
  scale_fill_manual(values = c("Forest" = "forestgreen",
                               "Barren" = "purple",
                               "Incipient" = "orange"),
                    name = "Site type") +
  # Add county boundaries for context (assuming ca_counties is defined)
  geom_sf(data = ca_counties, fill = "gray", color = "gray80") +
  # Add your inset map (g1_inset assumed to be defined)
  annotation_custom(grob = g1_inset, xmin = -122.01, xmax = -121.96, ymin = 36.625) +
  # Add jittered labels for each new cluster with the recoded label:
  ggrepel::geom_label_repel(data = jittered_labels,
                            aes(label = new_cluster_relabel, geometry = geometry),
                            stat = "sf_coordinates",
                            size = 3, color = "black") +
  #add north arrow
  ggsn::north(x.min = -121.99, x.max = -121.88, 
              y.min = 36.519, y.max = 36.65,
              location = "topright", 
              scale = 0.05, 
              symbol = 10)+
  # Define the map extent. Adjust xlim and ylim as needed.
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = st_crs(joined_polygons)) +
  theme_bw() +
  labs(x = "",y="")+
  base_theme  # Ensure base_theme is defined

p6

################################################################################
#Visualize occupiable habitat

p6 <- ggplot() +
  # Add clusters
  geom_sf(data = clusters %>% filter(year == 2023), fill = "forestgreen", color = "forestgreen") +
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

p6


################################################################################
#check out plots and export

p1
p2
p3
p4
p5
p6


#save unresolved clusters
ggsave(p2, filename = file.path(figdir, "Fig2_clusters_unresolved.png"), 
      width = 6, height = 8, units = "in", dpi = 600, bg = "white") #last write 26 Sept 2024


#save clusters resolved with convex hulls
ggsave(p4, filename = file.path(figdir, "Fig3_clusters_convex_hull.png"), 
       width = 6, height = 8, units = "in", dpi = 600, bg = "white") #last write 26 Sept 2024


#save clusters resolved with convex hulls
ggsave(p5, filename = file.path(figdir, "Fig4_cluster_timerseries.png"), 
       width = 13, height = 7.5, units = "in", dpi = 600, bg = "white") #last write 26 Sept 2024




