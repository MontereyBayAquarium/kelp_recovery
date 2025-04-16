

rm(list=ls())

################################################################################
#prep workspace and load data

require(librarian)

librarian::shelf(tidyverse, sf, ggplot2, gganimate, gifski)

#set directories
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
localdir <- "/Users/jossmith/Documents/Data/landsat"

# Read in the Landsat data (using your provided shapefile)
landsat_orig <- st_read(file.path(localdir, "processed/monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))
#Note: landsat data were obstained from the EDI Data Portal 
# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sbc.54.7

# Read in the California counties data
ca_counties <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp"))
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")

################################################################################
#prep data

# Filter to focal years and create presence/absence
landsat_filtered <- landsat_orig %>% 
  filter(year >= 1990, year <= 2023, quarter == 3) %>%
  mutate(kelp = factor(ifelse(!is.na(area) & area != 0, "present", "absent"),
                       levels = c("present", "absent")))

################################################################################
#buildanimation

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

# Build inset
g1_inset <- ggplotGrob(
  ggplot() +
    # Plot land
    geom_sf(data = foreign, fill = "grey80", color = "white", lwd = 0.3) +
    geom_sf(data = usa, fill = "grey80", color = "white", lwd = 0.3) +
    annotate("rect", xmin = -122.6, xmax = -121, ymin = 36.2, ymax = 37.1, 
             color = "black", fill = NA, lwd = 0.6) +
    labs(x = "", y = "") +
    coord_sf(xlim = c(-124.5, -117), ylim = c(32.5, 42)) +
    theme_void() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.6))
)



p <- ggplot() +
  #add basemap
  geom_sf(data = ca_counties, fill = "grey80", color = "white", size = 0.3) +
  #add kelp presence
  geom_sf(data = landsat_filtered, aes(fill = kelp), shape = 22, size = 1, alpha = 0.2,
          color = "transparent") +
  # Inset map 
  annotation_custom(grob = g1_inset, 
                    xmin = -122.024, xmax = -121.96,
                    ymin = 36.62) +
  scale_fill_manual(values = c("present" = "forestgreen", "absent" = "transparent")) +
  labs(title = "Kelp cover (presence only) - year: {current_frame}", x = "", y = "") +
  # North arrow
  ggsn::north(x.min = -122, x.max = -121.87, 
              y.min = 36.5, y.max = 36.65,
              location = "topright", scale = 0.07, symbol = 10) +
  # select focal area (Monterey Peninsula)
  coord_sf(xlim = c(-122, -121.87), ylim = c(36.50, 36.65), crs = 4326) +
  guides(fill = "none") +
  transition_manual(year) +
  ease_aes('linear') +
  enter_fade() +
  exit_fade() +
  annotate("text",
           x = -121.87,    
           y = 36.50,     
           label = paste("Code for reproducibility:",
                         "https://github.com/MontereyBayAquarium/kelp_recovery/",
                         "Data accessed from:",
                         "Santa Barbara Coastal LTER,",
                         "K.C. Cavanaugh, D.A. Siegel, D.C. Reed, and T.W. Bell. 2020.",
                         "SBC LTER: Time series of kelp biomass in the canopy from Landsat 5,",
                         "1984-2011 ver 7. Environmental Data Initiative.",
                         "https://doi.org/10.6073/pasta/4ae0597b9f0c6763215d10d4102a6067",
                         "(Accessed 2025-04-15).", sep = "\n"),
           hjust = 1,
           vjust = 0,
           size = 3,      
           color = "black",
           lineheight = 0.9) +
  theme_bw() +
  base_theme +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),  
    axis.text.x = element_blank(),   
    axis.text.y = element_blank(),   
    axis.ticks.x = element_blank(),  
    axis.ticks.y = element_blank(),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank()
  )

p


################################################################################
#render animation

anim <- animate(p,
                nframes = length(unique(landsat_filtered$year)),
                fps = 3,
                width = 600,
                height = 800,
                renderer = gifski_renderer())

anim 

# Save 
anim_save(file.path("/Users/jossmith/Downloads", "kelp_presence.gif"), animation = anim)





