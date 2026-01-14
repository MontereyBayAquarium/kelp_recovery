
#jogsmith@ucsc.edu

rm(list = ls())

#prepare raw dissection data for plotting. Rows should be indiviudal dissected
#urchins for the model predicted site types. Variables from dissection data 
#include animal mass, gonad mass, gonad index. Variables from field data include
#size frequency

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


#load LDA-predicted patch types
#lda_patch <- load(here::here("output","lda_patch_transitionsv2.rda"))
lda_patch <- load(here::here("output","lda_patch_transitionsv5.rda"))

################################################################################
#Step 1: prep authoritative patch type table

site_tab <- transitions_tbl_constrained %>%
  pivot_longer(
    cols = starts_with("patch_"),
    names_to = "year",
    names_prefix = "patch_",
    values_to = "site_type_predicted"
  ) %>%
  mutate(
    year = as.integer(year),
    site_type_predicted = as.factor(site_type_predicted)
  ) 

################################################################################
#Step 2: prep dissection data
dissect_build1 <- dissection_orig %>%
  mutate(
    year = year(survey_date)
  ) %>%
  #drop field that we don't need
  select(-test_height_mm, test_diameter_mm)


dissect_with_pred <- dissect_build1 %>%
  left_join(
    site_tab,
    by = c("site_official" = "site", "site_type_official"="site_type", "zone" = "zone", "year" = "year")
  )


################################################################################
#Step 2: prepare benthic survey data
View(urchin_sizefq)

sizefq_with_pred <- urchin_sizefq %>%
  mutate(year = year(survey_date))%>%
  left_join(
    site_tab,
    by = c("site", "site_type", "zone" = "zone", "year" = "year")
  )


################################################################################
#Step 3: explore

#gonad index
ggplot(
  dissect_with_pred,
  aes(x = gonad_index, y = site_type_predicted, fill = site_type_predicted)
) +
  geom_density_ridges(alpha = 0.7, color = "white", linewidth = 0.3) +
  facet_wrap(~ species, scales = "free_x") +
  labs(x = "Gonad Index", y = NULL) +
  theme_bw()


#gonad mass
ggplot(
  dissect_with_pred,
  aes(x = gonad_mass_g, y = site_type_predicted, fill = site_type_predicted)
) +
  geom_density_ridges(alpha = 0.7, color = "white", linewidth = 0.3) +
  facet_wrap(~ species, scales = "free_x") +
  labs(x = "Gonad Mass (g)", y = NULL) +
  theme_bw()



#size


sizefq_with_pred %>%
  filter(!is.na(site_type_predicted), !is.na(size_cm), !is.na(count)) %>%
  group_by(species, site_type_predicted, size_cm) %>%
  summarise(n = sum(count), .groups = "drop") %>%
  ggplot(aes(x = size_cm, y = n, color = site_type_predicted)) +
  geom_step(linewidth = 1) +
  facet_wrap(~ species, scales = "free_y") +
  labs(x = "Size (cm)", y = "Count", title = "Size-frequency (1 cm bins)") +
  theme_bw()

ggplot(
  size_box,
  aes(x = size_vec, color = site_type_predicted)
) +
  stat_ecdf(linewidth = 1) +
  facet_wrap(~ species, scales = "free_x") +
  labs(
    x = "Size (cm)",
    y = "Cumulative proportion"
  ) +
  theme_bw()































