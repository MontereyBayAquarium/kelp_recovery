

rm(list=ls())


################################################################################
#Prep workspace

#required packages
require(librarian)
librarian::shelf(tidyverse, readxl, gstat, sp, raster, lubridate, shiny, leaflet,
                 rnaturalearth, rnaturalearthdata, sf, viridis)

#set directories 
datdir <- here::here("output")

#read scan data
scan_dat <- read_csv(file.path(datdir,"scans","scans_data.csv"))

################################################################################
#Process for visualization

#prep for shiny
scan_prepped <- scan_dat %>%
  mutate(
    year = year(date),
    quarter = quarter(date),
    behav = as.factor(behav)
  )

#get coastline and buffer
coastline <- ne_download(scale = 10, type = "coastline", category = "physical", returnclass = "sf")

#crop extent to focal area
bbox <- st_bbox(scan_prepped %>% st_as_sf(coords = c("long", "lat"), crs = 4326))
bbox_expanded <- bbox + c(-0.1, -0.1, 0.1, 0.1)
coast_crop <- st_crop(coastline, bbox_expanded)

#buffer 400m offshore
coast_utm <- st_transform(coast_crop, crs = 32610)  # UTM zone 10N for Central CA
coast_buffer <- st_buffer(coast_utm, dist = 400)
coast_buffer_wgs84 <- st_transform(coast_buffer, crs = 4326)



################################################################################
#Build shiny with interpolation using IDW


ui <- fluidPage(
  titlePanel("Sea Otter Total IND Interpolation (by Quarter)"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "year", "Select Year:",
        choices = sort(unique(scan_prepped$year)),
        selected = 2024
      ),
      selectInput(
        "quarter", "Select Quarter:",
        choices = 1:4,
        selected = 3
      ),
      checkboxGroupInput(
        "behav", "Select Behavior(s):",
        choices = unique(scan_prepped$behav),
        selected = unique(scan_prepped$behav)
      )
    ),
    mainPanel(
      leafletOutput("map", height = 600)
    )
  )
)

server <- function(input, output, session) {
  
  # Filter data reactively
  filtered_data <- reactive({
    scan_prepped %>%
      filter(
        year == input$year,
        quarter == input$quarter,
        behav %in% input$behav
      ) %>%
      group_by(lat, long) %>%
      summarise(total_ind = sum(ind, na.rm = TRUE), .groups = "drop")
  })
  
  # Interpolation raster
  idw_raster <- reactive({
    dat <- filtered_data()
    if (nrow(dat) < 5) return(NULL)
    
    coordinates(dat) <- ~long + lat
    proj4string(dat) <- CRS("+proj=longlat +datum=WGS84")
    
    grid <- expand.grid(
      long = seq(min(dat$long) - 0.01, max(dat$long) + 0.01, length.out = 200),
      lat = seq(min(dat$lat) - 0.01, max(dat$lat) + 0.01, length.out = 200)
    )
    coordinates(grid) <- ~long + lat
    gridded(grid) <- TRUE
    proj4string(grid) <- CRS("+proj=longlat +datum=WGS84")
    
    idw_result <- idw(total_ind ~ 1, dat, newdata = grid)
    r <- raster(idw_result)
    
    # Mask raster to coastline buffer
    r_masked <- mask(r, as_Spatial(coast_buffer_wgs84))
    return(r_masked)
  })
  
  # Initialize map once
  output$map <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%
      setView(
        lng = mean(scan_prepped$long, na.rm = TRUE),
        lat = mean(scan_prepped$lat, na.rm = TRUE),
        zoom = 11
      )
  })
  
  # Update map when input changes
  observe({
    r <- idw_raster()
    if (is.null(r)) return()
    
    vals <- values(r)
    vals_valid <- vals[!is.na(vals) & vals >= 1]
    
    pal <- colorNumeric(
      palette = viridis::viridis(100, option = "D", direction = 1),
      domain = vals_valid,
      na.color = "#f0f9ff"
    )
    
    leafletProxy("map") %>%
      clearImages() %>%
      clearControls() %>%
      addRasterImage(
        r,
        colors = pal,
        opacity = 0.85,
        layerId = "raster"
      ) %>%
      addPolylines(data = coast_buffer_wgs84, color = "black", weight = 1.5) %>%
      addLegend(
        pal = pal,
        values = vals_valid,
        title = "Total IND",
        opacity = 1
      )
  })
  
  
}

shinyApp(ui, server)


################################################################################
#Build shiny with actual values

ui <- fluidPage(
  titlePanel("Sea Otter Total IND Interpolation (by Quarter)"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "year", "Select Year:",
        choices = sort(unique(scan_prepped$year)),
        selected = 2024
      ),
      selectInput(
        "quarter", "Select Quarter:",
        choices = 1:4,
        selected = 3
      ),
      checkboxGroupInput(
        "behav", "Select Behavior(s):",
        choices = unique(scan_prepped$behav),
        selected = unique(scan_prepped$behav)
      )
    ),
    mainPanel(
      leafletOutput("map", height = 600)
    )
  )
)

server <- function(input, output, session) {
  
  # Filter data reactively
  filtered_data <- reactive({
    scan_prepped %>%
      filter(
        year == input$year,
        quarter == input$quarter,
        behav %in% input$behav
      ) %>%
      group_by(lat, long) %>%
      summarise(total_ind = sum(ind, na.rm = TRUE), .groups = "drop")
  })
  
  # Interpolation raster
  idw_raster <- reactive({
    dat <- filtered_data()
    if (nrow(dat) < 5) return(NULL)
    
    coordinates(dat) <- ~long + lat
    proj4string(dat) <- CRS("+proj=longlat +datum=WGS84")
    
    grid <- expand.grid(
      long = seq(min(dat$long) - 0.01, max(dat$long) + 0.01, length.out = 200),
      lat = seq(min(dat$lat) - 0.01, max(dat$lat) + 0.01, length.out = 200)
    )
    coordinates(grid) <- ~long + lat
    gridded(grid) <- TRUE
    proj4string(grid) <- CRS("+proj=longlat +datum=WGS84")
    
    idw_result <- idw(total_ind ~ 1, dat, newdata = grid)
    r <- raster(idw_result)
    
    # Mask raster to coastline buffer
    r_masked <- mask(r, as_Spatial(coast_buffer_wgs84))
    return(r_masked)
  })
  
  # Initialize map once
  output$map <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%
      setView(
        lng = mean(scan_prepped$long, na.rm = TRUE),
        lat = mean(scan_prepped$lat, na.rm = TRUE),
        zoom = 11
      )
  })
  
  # Update map when input changes
  observe({
    dat <- scan_prepped %>%
      filter(
        year == input$year,
        quarter == input$quarter,
        behav %in% input$behav
      ) %>%
      group_by(lat, long) %>%
      summarise(total_ind = sum(ind, na.rm = TRUE), .groups = "drop")
    
    if (nrow(dat) == 0) return()
    
    pal <- colorNumeric(
      palette = "Blues",
      domain = dat$total_ind,
      na.color = "#f0f9ff"
    )
    
    leafletProxy("map") %>%
      clearMarkers() %>%
      clearControls() %>%
      addCircleMarkers(
        data = dat,
        lng = ~long,
        lat = ~lat,
        radius = ~sqrt(total_ind) * 2,
        fillColor = ~pal(total_ind),
        fillOpacity = 0.8,
        color = "#444",
        weight = 0.5,
        popup = ~paste("Total IND:", total_ind)
      ) %>%
      addPolylines(data = coast_buffer_wgs84, color = "black", weight = 1.5) %>%
      addLegend(
        pal = pal,
        values = dat$total_ind,
        title = "Total IND",
        opacity = 1
      )
  })
  
}

shinyApp(ui, server)



################################################################################
#Plot by recovery polygons

# Helper function to unzip KMZ and return the path to KML
unzip_kmz <- function(kmz_path, out_dir = tempdir()) {
  unzip(kmz_path, exdir = out_dir)
  list.files(out_dir, pattern = "\\.kml$", full.names = TRUE)
}

# Define the actual paths to your KMZ files
kmz_paths <- list(
  incipient = "/Users/jossmith/code_respositories/MBA_kelp_recovery_analyses/output/incipient_forests.kmz",
  persistent_barrens = "/Users/jossmith/code_respositories/MBA_kelp_recovery_analyses/output/persistent_barrens.kmz",
  persistent_forests = "/Users/jossmith/code_respositories/MBA_kelp_recovery_analyses/output/persistent_forests.kmz"
)

# Unzip and read each .kmz as sf
kml_layers <- lapply(kmz_paths, function(path) {
  kml_path <- unzip_kmz(path)
  st_read(kml_path[1], quiet = TRUE)
})

# Add polygon_type and combine into one sf
polygons_sf <- imap_dfr(kml_layers, ~ mutate(.x, polygon_type = .y))
polygons_sf <- st_make_valid(polygons_sf)

# Prepare scan_prepped (assumed to already exist)
# Convert scan_prepped to sf
scan_sf <- scan_prepped %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326, remove = FALSE)

# Spatial join: assign each point to a polygon (if any)
scan_joined <- st_join(scan_sf, polygons_sf, join = st_within)

# Filter only points that fall inside a polygon
scan_in_poly <- scan_joined %>%
  filter(!is.na(polygon_type))

# Summarize total IND by polygon_type, year, quarter
summary_tbl <- scan_in_poly %>%
  group_by(polygon_type, year, quarter) %>%
  summarise(total_ind = sum(ind, na.rm = TRUE), .groups = "drop")

# View result
print(summary_tbl)


