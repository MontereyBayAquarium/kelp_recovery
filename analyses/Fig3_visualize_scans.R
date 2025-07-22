

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
