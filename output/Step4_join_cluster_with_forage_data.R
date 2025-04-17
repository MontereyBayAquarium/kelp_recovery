

rm(list=ls())


################################################################################
#Prep workspace

#required packages
require(librarian)
librarian::shelf(tidyverse, sf, raster, shiny, tmap, scales, readxl, lubridate)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("figures")
output <- here::here("output")

#landsat data
clusters <- st_read(file.path(output, "/landsat/processed/named_clustersv2.geojson"))

#Forage data
forage_meta <- read_xlsx(file.path(basedir, "foraging_data/raw/ForageKey_lookup_table.xlsx"))
for_dat <- read_csv(file.path(basedir, "foraging_data/raw/Forage_data_2016to2023.csv"))
dives <- read_csv(file.path(basedir, "foraging_data/raw/Forage_dives_2016to2023.csv"))
index <- read_csv(file.path(basedir, "foraging_data/raw/Forage_index_2016to2023.csv"))


################################################################################
#Prep patch polygons

# Create a convex hull for each cluster while keeping the site type attribute
cluster_polygons <- clusters %>%
  filter(year == 2023) %>%
  group_by(site_num, site_type) %>%       
  summarize(geometry = st_combine(geometry), .groups = "drop") %>%  
  st_convex_hull()                                


################################################################################
#Prep forage data

#Step 1: prep metadata----------------------------------------------------------

meta_build1 <- forage_meta %>%
  janitor::clean_names() %>%
  #Sentence case
  mutate(category_code = ifelse(data_type == "Prey_vars",
                                str_to_sentence(category_code), category_code),
         category_long = str_to_sentence(category_long),
         description = str_to_sentence(description))

#Step 2: tidy data-----------==========-----------------------------------------

#clean up forage
for_dat_build1 <- for_dat %>% janitor::clean_names() %>%
  #correct inconsistencies
  mutate(qualifier = factor(tolower(qualifier)),
         pup_behav = factor(tolower(pup_behav)),
         mom_resp = factor(tolower(mom_resp)),
         outcome = factor(tolower(outcome)),
         prey = factor(prey),
         tooltype = factor(tooltype),
         pupsh = factor(pupsh),
         pup_behav = factor(pup_behav),
         stolenfr = factor(stolenfr),
         steal = factor(steal),
         sizecm = factor(sizecm)
  )

#clean up dives
dive_build1 <- dives %>% janitor::clean_names() %>%
  #correct inconsistencies
  mutate(where = tolower(where),
         #fix incorrect sign
         long_obs_deg = ifelse(long_obs_deg == 121,-121,long_obs_deg),
         #fix kelp type
         kelptype = ifelse(kelptype == "xx","x",kelptype),
         kelptype = factor(tolower(kelptype)),
         success = factor(tolower(success)),
         sb_area = factor(toupper(sb_area))
         
  ) 

#clean up index
index_build1 <- index %>% janitor::clean_names() %>%
  mutate(visib = factor(tolower(visib)),
         sex = factor(sex),
         ageclass = factor(ageclass),
         status = factor(status),
         pupsz = factor(pupsz),
         consort = factor(consort),
         obsbeg = factor(obsbeg),
         obsend = factor(obsend),
         daynight = factor(daynight),
         area = factor(area),
         sky = factor(sky),
         winddir = factor(winddir),
         seaopen = factor(seaopen),
         seafix = factor(seafix),
         swell = factor(swell)
  )


#Step 3: merge -----------------------------------------------------------------

forage_join <- for_dat_build1 %>%
  #filter(steal == "No")%>% 
  #select variables of interest
  dplyr::select(foragdiv_id, foragdata_id, preynum, prey, number,
                size, qualifier)%>%
  mutate_if(is.character, str_trim)

dive_join <- dive_build1 %>%
  dplyr::select(foragdiv_id, bout, subbout, lat, long, canopy,
                kelptype, divenum, dt, success) %>%
  mutate_if(is.character, str_trim) %>%
  #filter successful dives only
  filter(success == "y")

index_join <- index_build1 %>%
  dplyr::select(bout, date, otterno, sex, ageclass, status)%>%
  mutate_if(is.character, str_trim)


#merge forage data and bout info
data_build1 <- left_join(dive_join, forage_join, by="foragdiv_id")
anti_build <- anti_join(dive_join, forage_join, by="foragdiv_id")


#merge with index
data_build2 <- left_join(data_build1, index_join, by="bout")
anti_build2 <- anti_join(data_build1, index_join, by="bout")

#tidy
data_build3 <- data_build2 %>%
  mutate(
    date  = dmy(date),
    year  = year(date),
    month = month(date),
    day   = day(date)
  ) %>%
  select(
    year, month, day,
    foragdiv_id, foragdata_id,
    bout, subbout,
    lat, long, otterno,
    everything()
  ) %>%
  filter(year >= 2016)


################################################################################
#Merge patch ID and forage data

str(cluster_polygons)
str(data_build3)

#Step1: drop any forage obs with missing lat/long
points_clean <- data_build3 %>%
  filter(!is.na(long), !is.na(lat))

#Step2: convert forage obs to sf
points_sf <- st_as_sf(
  points_clean,
  coords = c("long", "lat"),
  crs    = st_crs(cluster_polygons)
)


#Step3: spatial join to grab site_type for any forage obs that fall inside
inside_join <- st_join(
  points_sf,
  cluster_polygons %>% select(site_type),
  left = TRUE,
  join = st_within
)

#Step4: for remaining points (NA site_type), find the nearest polygon:
na_idx       <- which(is.na(inside_join$site_type))
nearest_poly <- st_nearest_feature(inside_join[na_idx, ], cluster_polygons)

#Step5: fill in those NAs from the nearest polygon’s site_type
inside_join$site_type[na_idx] <- cluster_polygons$site_type[nearest_poly]

#check that it worked
nrow(inside_join)
nrow(points_clean) 

################################################################################
#Export

write_csv(inside_join, file = 
            file.path(basedir, 
                      "foraging_data/processed/recovery_patches/forage_cluster_2016_22.csv")) #last write 17 April 2023





