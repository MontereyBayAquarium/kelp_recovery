

rm(list=ls())


################################################################################
#Prep workspace

#required packages
require(librarian)
librarian::shelf(tidyverse, readxl, googledrive, readr, stringr, janitor, here)

drive_auth()

#set directories 
gdir <- "1vsT-_TrHs0A3xWBG7Gd3sNPXxP3a7wxX"
outdir <- here::here("output")

#load scan data from Google drive
scan_raw <- drive_ls(as_id(gdir), type = "csv")

#download scan data and compile into a single file
scan_dat <- purrr::map_dfr(scan_raw$name, function(file_name) {
  # Get the file info row
  file_row <- scan_raw %>% filter(name == file_name)
  
  # Download temporarily
  temp_path <- tempfile(fileext = ".csv")
  drive_download(file = file_row$id, path = temp_path, overwrite = TRUE)
  
  # Extract date from filename
  date_str <- str_extract(file_name, "\\d{8}")
  date_parsed <- as.Date(date_str, format = "%m%d%Y")
  
  # Read and add date, force all columns to character to avoid type mismatch
  read_csv(temp_path, col_types = cols(.default = "c")) %>%
    mutate(date = date_parsed)
})

################################################################################
#Tidy up dataframe
scan_build1 <- scan_dat %>%
  select(-1) %>%
  janitor::clean_names() %>%
  mutate(
    ind = as.numeric(ind),
    large = as.numeric(large),
    small = as.numeric(small),
    temperatur = as.numeric(temperatur),
    behav = as.factor(behav),
    canopy = as.factor(canopy),
    kelptype = as.factor(kelptype),
    seafix = as.factor(seafix),
    prey = as.factor(prey),
    prey2 = as.factor(prey2),
    prey3 = as.factor(prey3),
    prey4 = as.factor(prey4),
    observer1 = as.factor(observer1),
    observer2 = as.factor(observer2),
    visibility = as.factor(visibility),
    sky = as.factor(sky),
    wind = as.factor(wind),
    dir = as.factor(dir),
    seaop = as.factor(seaop),
    swell = as.factor(swell)
  ) %>%
  select(date, lat, long, ind, large, small, behav, canopy, kelptype, seafix,
         prey, prey2, prey3, prey4, visibility, sky, wind, dir, temp = temperatur,
         seaop, swell)

################################################################################
#Export
write.csv(scan_build1, file = file.path(outdir, "scans", "scans_data.csv"), row.names = F)


                       


















