

rm(list=ls())


################################################################################
#Prep workspace

#required packages
require(librarian)
librarian::shelf(tidyverse, readxl)

#set directories 
basedir <- "Z:/kelp_recovery/data"
figdir <- here::here("figures")
output <- here::here("output")

#sofa data with cluster id
sofa_dat <- read_csv(file.path(basedir, "foraging_data/processed/recovery_patches/forage_cluster_2016_24.csv"))

#restructure sofa
sofa_dat$Site <- sofa_dat$site_type

#subset columns
sofa_dat <- sofa_dat %>% dplyr::select(-geometry, -site_type)

################################################################################
#Export

write_csv(sofa_dat, file = 
            file.path(basedir, 
                      "sofa_data/intermediate/recovery_patches/sofa_input_forage_cluster_2016_24.xlsx")) #last write 18 April 2025 LK


