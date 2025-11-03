
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
librarian::shelf(tidyverse, lubridate, sf, stringr, purrr, terra, janitor,
                 rnaturalearth, rnaturalearthdata)

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

#load site patches
site_patches <- st_read(here::here("output","gis_data","processed","site_patch_polygons.shp"))

#load LDA-predicted patch types
#lda_patch <- load(here::here("output","lda_patch_transitionsv2.rda")) #old
lda_patch <- load(here::here("output","lda_patch_transitionsv3.rda")) #new


# read CA state
ca_state <- st_read("/Volumes/enhydra/data/kelp_recovery/gis_data/raw/CA_state/ca_boundary_wgs84.shp", quiet=TRUE) |> st_transform(4326)


################################################################################
#Step 1: Average to site, zone, site_type for each year 
kelp_avg <- kelp_data %>%
  dplyr::select(-macro_stipe_sd_20m2) %>%
  dplyr::group_by(site, site_type, latitude, longitude, zone, survey_date) %>%
  dplyr::summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  dplyr::select(-transect)

quad_avg <- quad_data %>%
  dplyr::group_by(site, site_type, latitude, longitude, zone, survey_date) %>%
  dplyr::summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  dplyr::select(-quadrat, -transect)

dat_agg <- kelp_avg %>%
  dplyr::inner_join(
    quad_avg,
    by = c("site","site_type","latitude","longitude","zone","survey_date"),
    suffix = c("_kelp","_quad")
  )


quad_zone <- dat_agg %>%
  group_by(latitude, longitude, site, site_type,
           survey_date, zone) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

#intermediate step: prepare site metadata table

site_table <- quad_zone %>% 
  group_by(site, site_type, zone) %>%
  distinct(latitude, longitude)

#write_csv(site_table, here::here("output","site_meta_data","site_table.csv"))


################################################################################
#Step 2: assign model-predicted patch types

str(quad_zone)
str(transitions_tbl_constrained)

quad_zone_with_pred <- quad_zone %>%
  left_join(
    transitions_tbl_constrained %>%
      dplyr::select(site, site_type, zone, patch_2024, patch_2025),
    by = c("site","site_type", "zone")
  ) %>%
  mutate(
    pred_patch = case_when(
      format(survey_date, "%Y") == "2024" ~ as.character(patch_2024),
      format(survey_date, "%Y") == "2025" ~ as.character(patch_2025),
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(-patch_2024, -patch_2025)


################################################################################
#Step 5: join with patch geometry

# Convert to sf using lat/lon
quad_zone_sf <- st_as_sf(
  quad_zone_with_pred,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)


site_patches_single <- site_patches %>%
  st_cast("POLYGON") %>%              # split MULTIPOLYGON into indiv POLYGON
  mutate(patch_id = row_number())     # assign unique polygon ID

plot(site_patches_single)

#join points to polygons
site_patches_with_points <- site_patches_single %>%
  st_join(quad_zone_sf, join = st_intersects, left = TRUE)


#inspect
ggplot(site_patches_with_points %>% filter(year(survey_date) == 2024)) +
  geom_sf(aes(fill = pred_patch), color = "black") +
  theme_minimal() +
  labs(
    title = "Predicted Patch Type by Independent Polygon (2024)",
    fill = "Predicted Patch"
  )

str(site_patches_with_points)


#tidy up and add classifier
quad_build3 <- site_patches_with_points %>%
                mutate(patch_cat = ifelse(year(survey_date) == 2024,"predicted 2024","predicted 2025")) %>%
                dplyr::select(-site_type.x) %>%
                dplyr::select(patch_id, latitude, longitude, 
                       survey_date,site, site_type = site_type.y, pred_patch, 
                       everything()) %>%
                mutate(pred_patch = ifelse(is.na(pred_patch),site_type,pred_patch)) %>%
                filter(!(is.na(pred_patch)))

ggplot(quad_build3) +
  geom_sf(aes(fill = pred_patch), color = "black") +
  geom_sf(data = st_centroid(quad_build3), color = "black", size = 1) +   # overlay points
  facet_wrap(~patch_cat, nrow=1)+
  theme_minimal() +
  labs(
    title = "Predicted Patch Type by Independent Polygon (2024)",
    fill = "Predicted Patch"
  )


#save(quad_build3, file = here::here("output","survey_data","processed","zone_level_data3.rda"))

################################################################################
#Step 6: prepare scan data for plot

# Convert to sf using lat/lon
scan_sf <- st_as_sf(
  scan_orig,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
) %>% filter(year == 2024 | year == 2025)



################################################################################
#Build Figure1

my_theme <- theme(axis.text=element_text(size=6),
                  axis.text.y = element_text(angle = 90, hjust = 0.5),
                  axis.title=element_text(size=8),
                  plot.tag=element_text(size=8, face="bold"),
                  plot.title=element_text(size=7, face="bold"),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(colour = "black"),
                  legend.key = element_blank(),
                  legend.background = element_rect(fill=alpha('blue', 0)),
                  legend.key.height = unit(1, "lines"), 
                  legend.text = element_text(size = 6),
                  legend.title = element_text(size = 7),
                  strip.background = element_blank(),
                  strip.text = element_text(size = 6 ,face="bold"),
                  panel.background = element_rect(fill="white"))


# Inset map (California with MB bounding box)
usa     <- ne_states(country="United States of America", returnclass="sf")
foreign <- ne_countries(country=c("Canada","Mexico"), returnclass="sf")
#read state
ca_counties <- st_read(file.path(datadir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 


# Build inset
g1_inset <-  ggplotGrob(
  ggplot() +
    # Plot land
    geom_sf(data=foreign, fill="grey80", color="white", lwd=0.3) +
    geom_sf(data=usa, fill="grey80", color="white", lwd=0.3) +
    # Plot box
    annotate("rect", xmin=-122.6, xmax=-121, ymin=36.2, ymax=37.1, color="black", fill=NA, lwd=0.8) +
    # Label regions
    #geom_text(data=region_labels, mapping=aes(y=lat_dd, label=region), x= -124.4, hjust=0, size=2) +
    # Labels
    labs(x="", y="") +
    # Crop
    coord_sf(xlim = c(-124.5, -117), ylim = c(32.5, 42)) +
    # Theme
    theme_bw() + my_theme +
    theme( plot.margin = unit(rep(0, 4), "null"),
           panel.margin = unit(rep(0, 4), "null"),
           panel.background = element_rect(fill='transparent'), #transparent panel bg
           # plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
           axis.ticks = element_blank(),
           axis.ticks.length = unit(0, "null"),
           axis.ticks.margin = unit(0, "null"),
           axis.text = element_blank(),
           axis.title=element_blank(),
           axis.text.y = element_blank(),
           plot.title =  element_text(hjust = 0, size=14, face = "bold"))
)






g1 <-ggplot(quad_build3 %>% filter(year(survey_date) == 2024)) +
  # plot patches
  geom_sf(aes(fill = pred_patch), color = "black") +
  # add foraging locations (sea otter observations)
  geom_sf(
    data = scan_sf %>% filter(year == 2024),
    aes(shape = "Sea otter \nobservation"),
    color = "black",
    size = 0.5,
    alpha = 0.3
  ) +
  # add survey sites (subtidal sites)
  geom_sf(
    data = st_centroid(quad_build3 %>% filter(year(survey_date) == 2024)),
    aes(shape = "Subtidal survey \nsite"),
    fill = "yellow",
    #alpha=0.7,
    color = "black",
    size = 1.5
  ) +
  # plot land
  geom_sf(data = ca_counties, color = "grey70", fill = "grey85") +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.53, 36.64)) +
  theme_bw() +
  # Add plot inset
  annotation_custom(
    grob = g1_inset,
    xmin = -122.001,
    xmax = -121.97,
    ymin = 36.62
  ) +
  labs(
    title = "2024",
    fill = "Patch type",
    shape = ""  # legend title for shapes
  ) +
  scale_fill_manual(
    values = c(
      "BAR"   = "#7570B3",
      "FOR"   = "#1B9E77",
      "INCIP" = "#D95F02"
    )
  ) +
  # custom shape legend
  scale_shape_manual(
    values = c(
      "Subtidal survey \nsite"        = 24,  # filled triangle
      "Sea otter \nobservation" = 21  # filled circle
    )
  ) +
  guides(
    fill = guide_legend(order = 1),
    shape = guide_legend(order = 2, override.aes = list(
      fill = c("black", "yellow"),
      color = c("black", "black"),
      alpha = c(1, 0.6),
      size = c(3, 2)
    ))
  ) +
  my_theme 
#g1


g2 <-ggplot(quad_build3 %>% filter(year(survey_date) == 2025)) +
  # plot patches
  geom_sf(aes(fill = pred_patch), color = "black") +
  # add foraging locations (sea otter observations)
  geom_sf(
    data = scan_sf %>% filter(year == 2025),
    aes(shape = "Sea otter \nobservation"),
    color = "black",
    size = 0.5,
    alpha = 0.3
  ) +
  # add survey sites (subtidal sites)
  geom_sf(
    data = st_centroid(quad_build3 %>% filter(year(survey_date) == 2025)),
    aes(shape = "Subtidal survey \nsite"),
    fill = "yellow",
    color = "black",
    size = 1.5,
   # alpha=0.7
  ) +
  # plot land
  geom_sf(data = ca_counties, color = "grey70", fill = "grey85") +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.53, 36.64)) +
  theme_bw() +
  labs(
    title = "2025",
    fill = "Patch type",
    shape = ""  # legend title for shapes
  ) +
  scale_fill_manual(
    values = c(
      "BAR"   = "#7570B3",
      "FOR"   = "#1B9E77",
      "INCIP" = "#D95F02"
    )
  ) +
  # custom shape legend
  scale_shape_manual(
    values = c(
      "Subtidal survey \nsite"        = 24,  # filled triangle
      "Sea otter \nobservation" = 21  # filled circle
    )
  ) +
  guides(
    fill = guide_legend(order = 1),
    shape = guide_legend(order = 2, override.aes = list(
      fill = c("black", "yellow"),
      color = c("black", "black"),
      alpha = c(1, 0.6),
      size = c(3, 2)
    ))
  ) +
  my_theme
#g2

g <- ggarrange(g1, g2, common.legend = TRUE, legend = "right")

ggsave(
   here::here("figures", "Fig1_map_figure.png"),
   g,
   width = 7, height = 4, dpi = 600, bg = "white"
 )




