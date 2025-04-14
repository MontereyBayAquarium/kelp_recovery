
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
  geom_sf(data = clusters, aes(color = incipient)) +
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
  group_by(site_num, incipient) %>%              # group by your cluster ID and type
  summarize(geometry = st_combine(geometry), .groups = "drop") %>%  # combine points in each cluster
  st_convex_hull()                                # compute the convex hull
# Compute the convex hull for each group


p1 <- ggplot() +
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

p1






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
years_to_label <- seq(2014, max(landsat_build4$year), by = 3)

#Make sure to turn off year filter for landsay_build4
p2 <- ggplot(landsat_build4 %>% filter(year > 2013), aes(x = year, y = perc_of_max_3,
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

p2


#ggsave(p2,  filename=file.path(figdir, "Cluster_timeseries.png"), width = 16, height = 10, units = "in",
#      bg = "white", dpi = 600)


################################################################################
#resolve adjacent clustering







