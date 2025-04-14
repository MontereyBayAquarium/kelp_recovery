
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
final_data <- st_read(file.path(output, "/landsat/processed/kelp_area_by_cluster.geojson")) %>%
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
ggplot(landsat_build3 %>% filter(year > 2013), aes(x = year, y = perc_of_max_3)) +
  geom_point() + 
  geom_smooth(se = TRUE) +
  facet_wrap(~ site_num, scales = "free_y") + theme_bw() + base_theme


################################################################################
#determine center coords (THIS WILL BE IMPORTANT TO SAVE LATER)

cluster_coord <- landsat_build4 %>%
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



################################################################################
#label incipient clusters

#rename sites and create table
landsat_build4 <- landsat_build3 %>%
  mutate(
    #set cluster order
    incipient = case_when(
      site_num == 1 ~ "Forest",
      site_num == 2 ~ "Forest",
      site_num == 3 ~ "Forest",
      site_num == 4 ~ "Forest",
      site_num == 5 ~ "Barren",
      site_num == 6 ~ "Barren",
      site_num == 7 ~ "Barren",
      site_num == 8 ~ "Barren",
      site_num == 9 ~ "Barren",
      site_num == 10 ~ "Forest",
      site_num == 11 ~ "Barren",
      site_num == 12 ~ "Barren",
      site_num == 13 ~ "Barren",
      site_num == 14 ~ "Barren",
      site_num == 15 ~ "Incipient",
      site_num == 16 ~ "Incipient",
      site_num == 17 ~ "Barren",
      site_num == 18 ~ "Barren",
      site_num == 19 ~ "Incipient",
      site_num == 20 ~ "Barren",
      site_num == 21 ~ "Barren",
      site_num == 22 ~ "Barren",
      site_num == 23 ~ "Incipient",
      site_num == 24 ~ "Barren",
      site_num == 25 ~ "Forest",
      site_num == 26 ~ "Incipient",
      site_num == 27 ~ "Incipient",
      site_num == 28 ~ "Barren",
      site_num == 29 ~ "Incipient",
      site_num == 30 ~ "Barren",
      site_num == 31 ~ "Barren",
      site_num == 32 ~ "Incipient",
      site_num == 33 ~ "Incipient",
      site_num == 34 ~ "Incipient",
      site_num == 35 ~ "Incipient",
      site_num == 36 ~ "Barren",
      site_num == 37 ~ "Barren",
      site_num == 38 ~ "Barren",
      site_num == 39 ~ "Forest",
      site_num == 40 ~ "Forest",
      site_num == 41 ~ "Forest",
      site_num == 42 ~ "Barren",
      site_num == 43 ~ "Incipient",
      site_num == 44 ~ "Incipient",
      site_num == 45 ~ "Forest",
      site_num == 46 ~ "Forest",
      site_num == 47 ~ "Barren",
      site_num == 48 ~ "Forest",
      site_num == 49 ~ "Forest",
      site_num == 50 ~ "Forest",
      site_num == 51 ~ "Barren",
      site_num == 52 ~ "Forest",
      site_num == 53 ~ "Forest",
      site_num == 54 ~ "Barren",
      site_num == 55 ~ "Forest",
      site_num == 56 ~ "Incipient",
      site_num == 57 ~ "Barren",
      site_num == 58 ~ "Barren",
      site_num == 59 ~ "Forest",
      site_num == 60 ~ "Forest",
      site_num == 61 ~ "Barren",
      site_num == 62 ~ "Incipient",
      site_num == 63 ~ "Forest",
      site_num == 64 ~ "Forest",
      site_num == 65 ~ "Incipient",
      site_num == 66 ~ "Incipient",
      site_num == 67 ~ "Incipient",
      site_num == 68 ~ "Incipient",
      site_num == 69 ~ "Barren",
      site_num == 70 ~ "Incipient",
      site_num == 71 ~ "Forest",
      site_num == 72 ~ "Barren",
      site_num == 73 ~ "Forest",
      TRUE ~ NA)
  ) 

################################################################################
#Export named clusters

# Export landsat_build4 as a GeoJSON file
st_write(landsat_build4, file.path(output, "landsat/processed/named_clusters.geojson"), delete_dsn = TRUE)


################################################################################
#prep for export to KML to visualize in Google maps

# transform landsat data to Teale Albers

landsat_build5 <- landsat_build4 %>% st_as_sf(crs = 4326) %>% filter(year == 2023)


#Build barren layer
rast_build1 <- st_transform(landsat_build5, crs = 3310) %>% filter (incipient == "Barren") 
r <- rast(rast_build1, res=30)
landsat_rast_barren <- rasterize(rast_build1, r)

#Build forest layer
rast_build1 <- st_transform(landsat_build5, crs = 3310) %>% filter (incipient == "Forest") 
r <- rast(rast_build1, res=30)
landsat_rast_forest <- rasterize(rast_build1, r)

#Build Incipient
rast_build1 <- st_transform(landsat_build5, crs = 3310) %>% filter (incipient == "Incipient") 
r <- rast(rast_build1, res=30)
landsat_rast_incipient <- rasterize(rast_build1, r)


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

