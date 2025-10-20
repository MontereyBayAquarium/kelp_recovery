
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


#load site patches
site_patches <- st_read(here::here("output","gis_data","processed","site_patch_polygons.shp"))

#load LDA-predicted patch types
lda_patch <- load(here::here("output","lda_patch_transitionsv2.rda"))

# read CA state
ca_state <- st_read("/Volumes/enhydra/data/kelp_recovery/gis_data/raw/CA_state/ca_boundary_wgs84.shp", quiet=TRUE) |> st_transform(4326)


################################################################################
#Step 1: prep dissection data
dissect_build1 <- dissection_orig %>%
  mutate(
    year = year(survey_date)
  ) %>%
  group_by(year, site_official, site_type_official, zone, species) %>%
  summarise(
    #mean_gonad_mass_g = mean(gonad_mass_g, na.rm = TRUE),
    mean_gonad_index  = mean(gonad_index,  na.rm = TRUE),
    sd_gonad_index = sd(gonad_index,  na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = species,
    values_from = c(sd_gonad_index, mean_gonad_index, n),
    names_sep   = "_"
  ) %>% select(-mean_gonad_index_red_urchin, -n_red_urchin, -sd_gonad_index_red_urchin)

################################################################################
#Step 1: extrapolate urchin densities to total surveyed area, then calculate
#sea urchin biomass

#Summarize density to zone level (and add year)
urchin_zone_density <- quad_data %>%
  mutate(year = year(survey_date)) %>%
  group_by(site, site_type, zone, year) %>%
  summarise(
    mean_density = mean(purple_urchin_densitym2, na.rm = TRUE),
    se_density   = sd(purple_urchin_densitym2, na.rm = TRUE) / sqrt(n()),
    n_quadrats   = n(),
    total_urchins_zone = mean_density * 80,
    .groups = "drop"
  )

#Size-frequency (purple only, keep year)
size_dist_clean <- urchin_sizefq %>%
  mutate(year = year(survey_date)) %>%
  filter(species == "Purple") %>%
  group_by(site, site_type, zone, year, size_cm) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  group_by(site, site_type, zone, year) %>%
  mutate(prob = count / sum(count)) %>%
  filter(!is.na(prob) & prob > 0)

#Fallback distribution (global purple urchin)
fallback_dist <- size_dist_clean %>%
  group_by(size_cm) %>%
  summarise(prob = sum(count) / sum(size_dist_clean$count), .groups = "drop")

#Join and simulate individuals (keeping year)
simulated_sizes <- urchin_zone_density %>%
  left_join(size_dist_clean, by = c("site", "site_type", "zone", "year")) %>%
  group_by(site, site_type, zone, year) %>%
  reframe(
    simulated_sizes_cm = {
      n_urchins <- round(first(total_urchins_zone))
      probs <- if (all(is.na(prob))) fallback_dist$prob else prob
      sizes <- if (all(is.na(size_cm))) fallback_dist$size_cm else size_cm
      sample(x = sizes, size = n_urchins, replace = TRUE, prob = probs)
    }
  )

################################################################################
#Step 2: build model to predict urchin biomass

urch_dat_orig <- dissection_orig %>%
  filter(species == "purple_urchin") %>%
  filter(!(is.na(test_diameter_mm) | is.na(animal_24hr_mass_g))) %>%
  #drop some outliers
  filter(animal_24hr_mass_g < 100) %>%
  filter(!(animal_24hr_mass_g >20 & test_diameter_mm <30)) %>%
  filter(!(animal_24hr_mass_g >38 & test_diameter_mm <40)) %>%
  filter(!(animal_24hr_mass_g <20 & test_diameter_mm > 45))
  

# derive parameters
#take a look
plot(urch_dat_orig$test_diameter_mm, urch_dat_orig$animal_24hr_mass_g)


# Initial estimates based on data
a_init <- -20
b_init <- 10
c_init <- 0.03

# Fit the biomass_fun model to the data with initial estimates
set.seed(1985)
fit <- nls(animal_24hr_mass_g ~ a + b * exp(c * test_diameter_mm), 
           data = urch_dat_orig,
           start = list(a = a_init, b = b_init, c = c_init))


# Extract the estimated parameters
a_est <- coef(fit)["a"]
b_est <- coef(fit)["b"]
c_est <- coef(fit)["c"]

# Print the estimated parameters
cat("Estimated Parameters:\n")
cat("a:", a_est, "\n")
cat("b:", b_est, "\n")
cat("c:", c_est, "\n")


#determine fit
# Calculate the predicted values from the model
predicted_values <- fitted(fit)

# Calculate the residuals
residuals <- urch_dat_orig$animal_24hr_mass_g - predicted_values

# Calculate the RSS (Residual Sum of Squares)
rss <- sum(residuals^2)

# Calculate the TSS (Total Sum of Squares)
mean_soft_mass <- mean(urch_dat_orig$animal_24hr_mass_g)
tss <- sum((urch_dat_orig$animal_24hr_mass_g - mean_soft_mass)^2)

# Calculate R-squared
r_squared <- 1 - (rss / tss)

# Print the R-squared value
cat("R-squared:", r_squared, "\n")


#plot
base_theme <-  theme(axis.text=element_text(size=12, color = "black"),
                     axis.title=element_text(size=12,color = "black"),
                     plot.tag=element_text(size=9,color = "black"),
                     plot.title=element_text(size=12,color = "black", face = "bold"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key.size = unit(0.3, "cm"), 
                     #legend.key = element_rect(fill = "white"), # Set it to transparent
                     legend.spacing.y = unit(0.1, "cm"),  
                     legend.text=element_text(size=8,color = "black"),
                     legend.title=element_blank(),
                     #legend.key.height = unit(0.1, "cm"),
                     #legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=10, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())

# Generate equation expression
equation_text <- substitute(italic(y) == a + b %*% e^(c * italic(x)) * "," ~ italic(R)^2 ~ "=" ~ r2,
                            list(a = round(a_est, 2), 
                                 b = round(b_est, 2), 
                                 c = round(c_est, 2),
                                 r2 = round(r_squared, 2)))

# Sample size
n <- nrow(urch_dat_orig)
sample_size_text <- paste("n =", n)

# Create the plot
g <- ggplot(urch_dat_orig, aes(x = test_diameter_mm, y = animal_24hr_mass_g)) +
  geom_point() +
  geom_line(aes(y = a_est + b_est * exp(c_est * test_diameter_mm)), color = "purple", size = 1) +
  labs(x = "Test Diameter (mm)", y = "Purple sea urchin biomass (g)") +
  theme_bw() + 
  base_theme +
  annotate("text", x = min(urch_dat_orig$test_diameter_mm), y = max(urch_dat_orig$animal_24hr_mass_g), 
           label = as.expression(equation_text), hjust = 0, vjust = 1, size = 5, color = "black") +
  annotate("text", x = min(urch_dat_orig$test_diameter_mm), y = max(urch_dat_orig$animal_24hr_mass_g) - 5, 
           label = sample_size_text, hjust = 0, vjust = 1, size = 5, color = "black")

g

#ggsave(g, file = file.path(figdir, "purple_urchin_td_biomass.png"), width = 6.5,
 #      height = 6.5, units = "in")

################################################################################
#Step 3: apply function to sampled urchin sizes to determine biomass

simulated_biomass <- simulated_sizes %>%
  unnest(simulated_sizes_cm) %>%
  mutate(
    test_diameter_mm = simulated_sizes_cm * 10,  # cm → mm
    biomass_g = -14.2 + 7.44 * exp(0.04 * test_diameter_mm),
    biomass_g = if_else(biomass_g < 0.5, 0.5, biomass_g)  # floor small/negative values
  )


################################################################################
#Step 4: convert to total biomass at each site and then calculate gonad mass

site_biomass <- simulated_biomass %>%
                group_by(year, site, site_type, zone) %>%
                summarize(total_biomass_g = sum(biomass_g),
                          n_biomass = n(),
                          sd_biomass = sd(biomass_g))


gonad_summary <- dissect_build1 %>%
  rename(
    site = site_official,
    site_type = site_type_official,
    mean_gi = mean_gonad_index_purple_urchin,
    sd_gi = sd_gonad_index_purple_urchin,
    n_gi = n_purple_urchin
  )

biomass_gonad <- site_biomass %>%
  left_join(gonad_summary, by = c("site", "site_type", "zone", "year"))

#calculate total goand mass
biomass_gonad <- site_biomass %>%
  left_join(gonad_summary, by = c("site", "site_type", "zone", "year")) %>%
  mutate(
    # compute mean biomass per individual
    mean_biomass_g = total_biomass_g / n_biomass,
    
    # compute gonad metrics
    mean_gonad_mass_g  = (mean_gi / 100) * mean_biomass_g,
    total_gonad_mass_g = (mean_gi / 100) * total_biomass_g,
    se_gonad_mass_g    = (sd_gi / 100) * mean_biomass_g / sqrt(n_gi)
  )


ggplot(biomass_gonad, aes(x = site, y = total_gonad_mass_g, fill = site_type)) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_wrap(~ year) +
  labs(
    x = "Site",
    y = expression("Total gonad mass (g·80 m"^-2*")"),
    fill = "Site type",
    title = "Total purple urchin gonad mass by site and year"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )


################################################################################
#Step 5: Average benthoc data to site, zone, site_type for each year 
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
            .groups = "drop") %>%
  mutate(year = year(survey_date))%>%
  #join gonad summary
  left_join(biomass_gonad, by = c("year","site","site_type","zone"))



#intermediate step: prepare site metadata table

site_table <- quad_zone %>% 
  group_by(site, site_type, zone) %>%
  distinct(latitude, longitude)

#write_csv(site_table, here::here("output","site_meta_data","site_table.csv"))


################################################################################
#Step 6: assign model-predicted patch types

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


#check
ggplot(quad_zone_with_pred %>% filter(total_gonad_mass_g < 4000), 
       aes(x = pred_patch, y = total_gonad_mass_g, fill = pred_patch)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  labs(
    x = "Predicted patch type",
    y = expression("Total gonad mass (g·80 m"^-2*")"),
    title = "Total purple urchin gonad mass across patch types"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 13, face = "bold")
  )


ggplot(quad_zone_with_pred, aes(x = site, y = total_gonad_mass_g, fill = pred_patch)) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_wrap(~ year) +
  labs(
    x = "Site",
    y = expression("Total gonad mass (g·80 m"^-2*")"),
    fill = "Site type",
    title = "Total purple urchin gonad mass by site and year"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )



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
  filter(!(is.na(pred_patch))) %>%
  select(patch_id, latitude, longitude, survey_date, year, site, site_type,
         pred_patch, patch_cat, zone, total_biomass_g, n_biomass, sd_biomass, sd_gi, 
         mean_gi, n_gi, mean_biomass_g, mean_gonad_mass_g, total_gonad_mass_g, 
         se_gonad_mass_g, everything())

ggplot(quad_build3) +
  geom_sf(aes(fill = pred_patch), color = "black") +
  geom_sf(data = st_centroid(quad_build3), color = "black", size = 1) +   # overlay points
  facet_wrap(~patch_cat, nrow=1)+
  theme_minimal() +
  labs(
    title = "Predicted Patch Type by Independent Polygon (2024)",
    fill = "Predicted Patch"
  )


#save(quad_build3, file = here::here("output","survey_data","processed","zone_level_data4.rda"))

################################################################################
#Step 6: prepare scan data for plot

# Convert to sf using lat/lon
scan_sf <- st_as_sf(
  scan_orig,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
) %>% filter(year == 2024 | year == 2025)