################################################################################
# Supplementary Figure S1: Sea otter occupancy, behavior, and focal bouts
# (aggregated across years, focal bouts grouped by prey type)
# Joshua G. Smith – UCSC Nearshore Ecology Research Group
################################################################################

# --- 1. Load packages ---------------------------------------------------------
rm(list = ls())
require(librarian)
librarian::shelf(
  tidyverse, janitor, lubridate, sf,
  gridExtra, RColorBrewer, here
)

# --- 2. Load data -------------------------------------------------------------
load(here::here("output", "survey_data", "processed", "zone_level_data4.rda"))
scan_orig   <- read_csv(file.path(here::here("output", "scans", "scans_data.csv")))
forage_orig <- read_csv("/Volumes/enhydra/data/foraging_data/processed/foraging_data_2024_2025_processed.csv")

years_keep <- c(2024, 2025)
patch_cols <- c("Barren"="#7570B3", "Incipient"="#D95F02", "Forest"="#1B9E77")

################################################################################
# 3. Process behavioral scan data
################################################################################

scan_clean <- scan_orig %>%
  clean_names() %>%
  mutate(
    lat  = as.numeric(lat),
    long = as.numeric(long),
    year = year(date),
    quarter = lubridate::quarter(date),
    behavior_group = ifelse(behav == "foraging", "Foraging", "Other")
  ) %>%
  drop_na(lat, long, ind, behavior_group)

scan_sf <- st_as_sf(scan_clean, coords = c("long","lat"), crs = st_crs(quad_build3))

pts_joined <- st_join(scan_sf, quad_build3, join = st_intersects, left = TRUE) %>%
  mutate(
    year = lubridate::year(date),
    pred_patch = case_when(
      str_detect(tolower(pred_patch),"bar")   ~ "BAR",
      str_detect(tolower(pred_patch),"incip") ~ "INCIP",
      str_detect(tolower(pred_patch),"for")   ~ "FOR",
      TRUE ~ toupper(pred_patch)
    )
  ) %>%
  filter(year %in% years_keep)

pts_behav <- pts_joined %>%
  st_drop_geometry() %>%
  mutate(
    pred_patch = case_when(
      str_detect(tolower(pred_patch),"bar")   ~ "Barren",
      str_detect(tolower(pred_patch),"incip") ~ "Incipient",
      str_detect(tolower(pred_patch),"for")   ~ "Forest",
      TRUE ~ NA_character_
    ),
    pred_patch = factor(pred_patch, levels=c("Barren","Incipient","Forest")),
    behav = str_to_title(str_trim(behav))
  ) %>%
  filter(!is.na(pred_patch), !is.na(behav), !is.na(ind), ind > 0)

################################################################################
# 4. Otter occupancy by patch type
################################################################################

patch_area_df <- quad_build3 %>%
  st_make_valid() %>%
  mutate(area_m2 = as.numeric(st_area(geometry))) %>%
  st_drop_geometry() %>%
  distinct(patch_id, area_m2)

occupancy_df <- pts_behav %>%
  group_by(year, patch_id, pred_patch) %>%
  summarise(mean_otters = mean(ind, na.rm = TRUE), .groups = "drop") %>%
  left_join(patch_area_df, by = "patch_id") %>%
  mutate(otters_per_km2 = mean_otters / (area_m2 / 1e6))

occupancy_summary <- occupancy_df %>%
  group_by(pred_patch) %>%
  summarise(
    mean_density = mean(otters_per_km2, na.rm = TRUE),
    se_density   = sd(otters_per_km2, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

################################################################################
# 5. Behavioral composition (pooled across years)
################################################################################

behavior_summary <- pts_behav %>%
  count(pred_patch, behav) %>%
  group_by(pred_patch) %>%
  mutate(prop = n / sum(n))

################################################################################
# 6. Focal bouts (foraging data, grouped prey)
################################################################################

# Identify focal bouts (>3 successful dives per prey per bout)
forage_bouts <- forage_orig %>%
  filter(success == "y") %>%
  group_by(year, bout, prey) %>%
  summarise(n_success = n_distinct(foragdiv_id), .groups = "drop") %>%
  mutate(focal_bout = n_success > 3) %>%
  filter(focal_bout)

# Get one coordinate per bout–prey–year
forage_locs <- forage_orig %>%
  semi_join(forage_bouts, by = c("year", "bout", "prey")) %>%
  group_by(year, bout, prey) %>%
  slice(1) %>%
  ungroup()

forage_sf <- st_as_sf(forage_locs, coords = c("long","lat"), crs = 4326, remove = FALSE)
quad3_same_crs <- st_transform(quad_build3, st_crs(forage_sf))

# Spatial join to assign patch type
focal_joined <- st_join(forage_sf, quad3_same_crs, join = st_intersects, left = TRUE) %>%
  mutate(
    pred_patch = case_when(
      str_detect(tolower(pred_patch),"bar")   ~ "Barren",
      str_detect(tolower(pred_patch),"incip") ~ "Incipient",
      str_detect(tolower(pred_patch),"for")   ~ "Forest",
      TRUE ~ NA_character_
    )
  ) %>%
  st_drop_geometry() %>%
  filter(!is.na(pred_patch))

# Group prey types into broader categories
focal_joined <- focal_joined %>%
  mutate(prey_group = case_when(
    prey %in% c("cam","clm","mus") ~ "Bivalve",
    prey %in% c("cra","kcr","can","rcr") ~ "Cancer crab",
    prey %in% c("gas","teg","whe") ~ "Gastropod",
    prey %in% c("pur","red","urc") ~ "Sea urchin",
    TRUE ~ "Other invertebrate"
  ))

# Summarize focal bouts by prey_group × patch × year
focal_summary <- focal_joined %>%
  rename(year = year.x) %>%
  group_by(year, pred_patch, prey_group) %>%
  summarise(n_bouts = n_distinct(bout), .groups = "drop") %>%
  group_by(pred_patch, prey_group) %>%
  summarise(
    mean_bouts = mean(n_bouts),
    se_bouts   = sd(n_bouts) / sqrt(n()),
    .groups = "drop"
  ) %>%
  group_by(pred_patch) %>%
  mutate(prop = mean_bouts / sum(mean_bouts))

prey_cols <- c(
  "Bivalve"           = "#1F78B4",  # teal blue
  "Cancer crab"       = "#FDBF6F",  # golden orange
  "Gastropod"         = "forestgreen",  # lavender
  "Sea urchin"        = "purple",  # deep rose
  "Other invertebrate"= "#B2B2B2"   # neutral gray
)


################################################################################
# 7. Plot panels
################################################################################

################################################################################
# Custom theme (matches Josh's usual style)
################################################################################

# Custom theme
my_theme <- theme(
  axis.text.x      = ggplot2::element_text(size = 8, color = "black"),
  axis.text.y      = ggplot2::element_text(size = 8, color = "black"),
  axis.title       = ggplot2::element_text(size = 10, color = "black"),
  legend.text      = ggplot2::element_text(size = 8, color = "black"),
  legend.title     = ggplot2::element_text(size = 8, color = "black"),
  plot.tag         = ggplot2::element_text(size = 10, color = "black"),
  panel.grid       = ggplot2::element_blank(),
  panel.background = ggplot2::element_blank(),
  axis.line        = ggplot2::element_line(colour = "black"),
  legend.key       = ggplot2::element_blank(),
  plot.title = element_text(face = "plain", size = 11, hjust = 0),
  legend.background = ggplot2::element_rect(
    fill = scales::alpha("blue", 0)
  ),
  strip.text       = ggplot2::element_text(
    size = 10, face = "bold", color = "black", hjust = 0
  ),
  strip.background = ggplot2::element_blank()
)

# --- Panel A: Otter occupancy -------------------------------------------------
gA <- ggplot(occupancy_summary, aes(x = pred_patch, y = mean_density, fill = pred_patch)) +
  geom_col(color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = mean_density - se_density, ymax = mean_density + se_density), width = 0.2) +
  scale_fill_manual(values = patch_cols) +
  labs(
    x = "Patch type",
    y = expression("Mean sea otter density (individuals km"^{-2}*")"),
    title = "A. Sea otter occupancy by patch type"
  ) +
  theme_bw()+
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  ) + my_theme

# --- Panel B: Behavioral composition -----------------------------------------
gB <- ggplot(behavior_summary, aes(x = pred_patch, y = prop, fill = behav)) +
  geom_col(color = "black", width = 0.7, position = "stack") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 8, "Set2"), name = "Behavior") +
  labs(
    x = "Patch type",
    y = "Proportion of observed behaviors",
    title = "B. Behavioral composition by patch type"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    plot.title = element_text(face = "plain", size = 13, hjust = 0)
  ) + my_theme

# --- Panel C: Focal bouts by prey group --------------------------------------
gC <- ggplot(focal_summary, aes(x = pred_patch, y = prop, fill = prey_group)) +
  geom_col(color = "black", width = 0.7, position = "stack") +
  scale_fill_manual(values = prey_cols, name = "Prey group") +
  labs(
    x = "Patch type",
    y = "Proportion of focal bouts (>3 successful dives)",
    title = "C. Bouts by prey and patch type"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "plain", size = 13, hjust = 0)
  ) + my_theme

################################################################################
# 8. Combine panels into Supplementary Figure S1
################################################################################

combined_fig <- gridExtra::grid.arrange(
  gA, gB, gC,
  ncol = 3,
  widths = c(1, 1.1, 1.1)
)



#gridExtra::grid.arrange(
#  gA, gB, gC,
#  layout_matrix = rbind(c(1, 2), c(3, 3)),
#  heights = c(2, 1.2)
#)

################################################################################
# 9. Save figure
################################################################################

ggsave(
  here("figures", "Fig5_patch_utilization.png"),
  combined_fig,
  width =14, height = 5, dpi = 600, bg = "white"
)
