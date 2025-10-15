
# RF model of sea otter foraging intensity across benthic patch types
# Joshua G. Smith – UCSC Nearshore Ecology Research Group


################################################################################
#load packages and data

rm(list = ls())
require(librarian)
librarian::shelf(
  tidyverse, janitor, lubridate, sf,
  randomForest, gridExtra, tidytext,
  viridis, here
)

datdir <- here::here("output")
load(here::here("output", "survey_data", "processed", "zone_level_data3.rda"))
scan_orig <- read_csv(file.path(here::here("output","scans","scans_data.csv")))
dissection_data <- read_csv("/Volumes/enhydra/data/kelp_recovery/MBA_kelp_forest_database/processed/dissection/dissection_data_cleanedv2.csv")

################################################################################
#Prep dissection data

dissect_build1 <- dissection_data %>%
  mutate(
    species = str_to_lower(species),        # make all lowercase
    species = str_trim(species),            # remove extra spaces
    species = case_when(
      species %in% c("red_urchin", "red urchin", "red_urchins") ~ "red_urchin",
      species %in% c("purple_urchin", "purple_urchins", "purple urchin") ~ "purple_urchin",
      TRUE ~ species))%>%
  mutate(year = year(date_collected)) %>%
  group_by(year, site_number, site_type, zone, species) %>%
  summarise(
    mean_gonad_mass_g = mean(gonad_mass_g, na.rm = TRUE),
    mean_gonad_index  = mean(gonad_index, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = species,
    values_from = c(mean_gonad_mass_g, mean_gonad_index, n),
    names_sep = "_"
  )

#Join with quad data
quad_build4 <- quad_build3 %>%
                mutate(year = year(survey_date))%>%
                left_join(dissect_build1, by = c("year","site"="site_number","site_type","zone")) 


################################################################################
#Prep RF model

predictors <- c(
  "purple_urchin_densitym2",
  "purple_urchin_conceiledm2",
  "mean_gonad_index_purple_urchin",
  "red_urchin_densitym2",
  "red_urchin_conceiledm2",
  #"mean_gonad_index_red_urchin",
 # "relief_cm",
#  "risk_index",
  "macro_stipe_density_20m2",
  "density20m2_nerlue"
)

center_effects <- TRUE
years_keep     <- c(2024, 2025)
patch_colors   <- c("BAR"="purple", "INCIP"="orange", "FOR"="forestgreen")


#clean up scan data
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

# Convert to sf
scan_sf <- scan_clean %>%
  st_as_sf(coords = c("long", "lat"), crs = st_crs(quad_build4))

# Join scans to patches
pts_joined <- st_join(scan_sf, quad_build4, join = st_intersects, left = TRUE) %>%
  mutate(
    year = lubridate::year(date),
    pred_patch = case_when(
      str_detect(tolower(pred_patch), "bar")   ~ "BAR",
      str_detect(tolower(pred_patch), "incip") ~ "INCIP",
      str_detect(tolower(pred_patch), "for")   ~ "FOR",
      TRUE ~ toupper(pred_patch)
    )
  ) %>%
  filter(year %in% years_keep)

#build patch-year level foraging data
# Average number of foraging individuals per scan within each patch-year
avg_foraging_by_patchyr <- pts_joined %>%
  st_drop_geometry() %>%
  mutate(is_foraging = tolower(behav) == "foraging",
         foraging_count_this_scan = ifelse(is_foraging, ind, 0)) %>%
  group_by(patch_id, year, pred_patch) %>%
  summarise(
    avg_foraging = mean(foraging_count_this_scan, na.rm = TRUE),
    n_scans      = n(),
    .groups = "drop"
  )

# aggregate benthic predictors
quad_year <- quad_build4 %>%
  mutate(
    year = lubridate::year(survey_date),
    pred_patch = case_when(
      str_detect(tolower(pred_patch), "bar")   ~ "BAR",
      str_detect(tolower(pred_patch), "incip") ~ "INCIP",
      str_detect(tolower(pred_patch), "for")   ~ "FOR",
      TRUE ~ toupper(pred_patch)
    )
  ) %>%
  filter(year %in% years_keep) %>%
  st_drop_geometry() %>%
  group_by(patch_id, year, pred_patch) %>%
  summarise(across(all_of(predictors), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# merge and prep rf data
rf_data <- avg_foraging_by_patchyr %>%
  inner_join(quad_year, by = c("patch_id", "year", "pred_patch")) %>%
  mutate(avg_foraging_log = log1p(avg_foraging)) %>%
  drop_na(avg_foraging_log, all_of(predictors))

# run RF model
patches <- sort(unique(rf_data$pred_patch))
data_imp <- data_r2 <- data_marg <- NULL

for (i in seq_along(patches)) {
  
  set.seed(1985)
  patch_do <- patches[i]
  sdata <- rf_data %>% filter(pred_patch == patch_do) %>% as.data.frame()
  
  if (nrow(sdata) < max(10, length(predictors) * 2)) {
    warning(sprintf("Skipping %s: not enough rows (%d)", patch_do, nrow(sdata)))
    next
  }
  
  message("Running RF for ", patch_do, " (n=", nrow(sdata), ")")
  
  fmla <- as.formula(paste("avg_foraging_log ~", paste(predictors, collapse = " + ")))
  rf_fit <- randomForest(fmla, data = sdata, ntree = 1501, importance = TRUE)
  
  preds <- predict(rf_fit, sdata)
  obs   <- sdata$avg_foraging_log
  r2    <- 1 - sum((obs - preds)^2) / sum((obs - mean(obs))^2)
  r2_df <- tibble(pred_patch = patch_do, r2 = r2)
  
  imp_df <- randomForest::importance(rf_fit, type = 2, scale = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("variable") %>%
    mutate(pred_patch = patch_do)
  
  # partial effects
  vars <- predictors
  for (j in seq_along(vars)) {
    parts <- randomForest::partialPlot(
      x = rf_fit, pred.data = sdata, x.var = vars[j], plot = FALSE
    )
    eff <- parts$y
    if (center_effects) eff <- eff - mean(eff, na.rm = TRUE)
    df <- tibble(pred_patch = patch_do, variable = vars[j],
                 value = parts$x, effect = eff)
    if (j == 1) marg_effects <- df else marg_effects <- bind_rows(marg_effects, df)
  }
  
  if (i == 1) {
    data_imp <- imp_df; data_marg <- marg_effects; data_r2 <- r2_df
  } else {
    data_imp <- bind_rows(data_imp, imp_df)
    data_marg <- bind_rows(data_marg, marg_effects)
    data_r2 <- bind_rows(data_r2, r2_df)
  }
}

# format for plotting
data_imp1 <- data_imp %>%
  rename(importance = IncNodePurity) %>%
  group_by(pred_patch) %>%
  mutate(importance_scaled = importance / max(importance, na.rm = TRUE)) %>%
  ungroup()

variable_importance_rank <- data_imp1 %>%
  group_by(variable) %>%
  summarise(average_importance = mean(importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(average_importance)) %>%
  mutate(rank = row_number())

data_marg2 <- data_marg %>% left_join(variable_importance_rank, by = "variable")

################################################################################
#Plot
my_theme <- theme(
  axis.text  = element_text(size = 8),
  axis.title = element_text(size = 9),
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 9),
  plot.tag = element_text(size = 10, face = "bold"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.text = element_text(size = 9, face = "bold")
)

g1 <- ggplot(
  data_imp1,
  aes(y = importance_scaled,
      x = tidytext::reorder_within(variable, -importance_scaled, pred_patch),
      fill = pred_patch)
) +
  facet_wrap(~ pred_patch, nrow = 1, scales = "free_x") +
  geom_col() +
  geom_text(
    data = data_r2,
    mapping = aes(x = Inf, y = 0.95, label = paste0("R² = ", round(r2, 2))),
    inherit.aes = FALSE, hjust = 1.05, size = 3
  ) +
  tidytext::scale_x_reordered() +
  scale_fill_manual(values = patch_colors) +
  labs(x = "", y = "Node impurity (scaled within patch)", tag = "A") +
  theme_bw() + my_theme +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))

g2 <- ggplot(
  data_marg2,
  aes(x = value, y = effect, color = pred_patch)
) +
  facet_wrap(~ variable, ncol = 2, scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(size = 1) +
  scale_color_manual(values = patch_colors, name = "Patch type") +
  labs(x = "Predictor value",
       y = "Centered partial effect (on log(1 + avg foraging per scan))",
       tag = "B") +
  theme_bw() + my_theme +
  theme(legend.position = "top")


gridExtra::grid.arrange(g1, g2, heights = c(0.42, 0.58))



################################################################################
#Plot gonad index by patch type

#prep dissection data
dissect_raw <- dissection_data %>%
  mutate(
    species = str_to_lower(species),
    species = str_trim(species),
    species = case_when(
      species %in% c("red_urchin", "red urchins", "red urchin") ~ "red_urchin",
      species %in% c("purple_urchin", "purple urchins", "purple urchin") ~ "purple_urchin",
      TRUE ~ species
    ),
    year = lubridate::year(date_collected)
  )

quad_clean <- quad_build3 %>%
  st_drop_geometry() %>%
  mutate(
    pred_patch = case_when(
      str_detect(tolower(pred_patch), "bar")   ~ "BAR",
      str_detect(tolower(pred_patch), "incip") ~ "INCIP",
      str_detect(tolower(pred_patch), "for")   ~ "FOR",
      TRUE ~ toupper(pred_patch)
    ),
    year = lubridate::year(survey_date)
  ) %>%
  select(year, site, site_type, zone, pred_patch) %>%
  distinct()

#join patches with dissection data
dissect_joined <- dissect_raw %>%
  left_join(quad_clean,
            by = c("year", "site_number" = "site",
                   "site_type", "zone"))


# choose species
sp_use <- "purple_urchin"

dissect_sp <- dissect_joined %>%
  filter(species == sp_use, !is.na(pred_patch), !is.na(gonad_index)) %>%
  mutate(pred_patch = factor(pred_patch, levels = c("BAR", "INCIP", "FOR")))


# Kruskal–Wallis 
kruskal_test(gonad_index ~ pred_patch, data = dissect_sp)

# Pairwise Wilcoxon with FDR correction
pairwise_res <- dissect_sp %>%
  pairwise_wilcox_test(gonad_index ~ pred_patch, p.adjust.method = "fdr")
pairwise_res

patch_cols <- c("BAR"="purple", "INCIP"="orange", "FOR"="forestgreen")


ggplot(dissect_sp %>% filter(gonad_index < 20),
       aes(x = pred_patch, y = gonad_index, fill = pred_patch)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.12, height = 0,
    alpha = 0.4, size = 1.8, shape = 21,
    color = "black"          # 👈 black outline
  ) +
  scale_fill_manual(values = patch_cols) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = 20.5            
  ) +
  stat_compare_means(
    comparisons = list(c("BAR","INCIP"), c("BAR","FOR"), c("INCIP","FOR")),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE,
    label.y = c(18, 19, 20.5)   
  ) +
  labs(
    x = "Patch Type",
    y = "Gonad Index",
    title = paste("Gonad index by patch type")
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )


