################################################################################
# Sea otter foraging and energetics — single-row figure
# Panels: A = Foraging effort (SOFA posterior)
#          B = Urchin consumption rate (SOFA posterior)
#          C = Patch-type probability vs. otter urchin take (Random Forest PDP)
#
# Author: Joshua G. Smith — UCSC Nearshore Ecology Research Group
################################################################################

rm(list = ls())

require(librarian)
librarian::shelf(
  tidyverse, janitor, lubridate, sf, here,
  randomForest, pdp, ggridges, ggpubr
)

################################################################################
# 1. Load data -----------------------------------------------------------------
################################################################################

load(here::here("output", "survey_data", "processed", "zone_level_data4.rda"))
forage_orig <- read_csv("/Volumes/enhydra/data/foraging_data/processed/foraging_data_2024_2025_processed.csv")

years_keep   <- c(2024, 2025)
focal_months <- c(6, 7, 8, 9)
patch_colors <- c("BAR"="purple", "INCIP"="orange", "FOR"="forestgreen")

################################################################################
# 2. Convert purple urchin prey to biomass -------------------------------------
################################################################################

size_key <- expand_grid(size = 1:4, qualifier = c("a","b","c")) %>%
  mutate(
    size_class = paste0(size, qualifier),
    size_cm = (size - 1) * 5 + case_when(
      qualifier == "a" ~ 1.66,
      qualifier == "b" ~ 3.32,
      qualifier == "c" ~ 4.98
    )
  )

pur_forage <- forage_orig %>%
  group_by(bout) %>%
  mutate(dives_for_pur = n_distinct(foragdiv_id[str_detect(prey, regex("^pur$", ignore_case = TRUE))])) %>%
  ungroup() %>%
  filter(
    month %in% focal_months,
    str_detect(prey, regex("^pur$", ignore_case = TRUE))
  ) %>%
  left_join(size_key, by = c("size", "qualifier")) %>%
  mutate(
    test_diameter_mm = size_cm * 10,
    biomass_g = -14.2 + 7.44 * exp(0.04 * test_diameter_mm),
    biomass_g = if_else(biomass_g < 0.5, 0.5, biomass_g),
    total_biomass_g = biomass_g * number
  ) %>%
  filter(size_cm < 8)

################################################################################
# 3. Match foraging data to mapped patches -------------------------------------
################################################################################

pur_sf <- st_as_sf(pur_forage, coords = c("long","lat"), crs = 4326, remove = FALSE)
quad4_same_crs <- st_transform(quad_build3, st_crs(pur_sf)) %>%
  mutate(survey_year = year(survey_date), survey_doy = yday(survey_date))

pur_spatial <- st_join(
  pur_sf,
  quad4_same_crs %>% select(patch_id, pred_patch, survey_date, survey_year, survey_doy),
  join = st_intersects,
  left = TRUE
)

pur_patch <- pur_spatial %>%
  mutate(forag_doy = yday(date),
         diff_days = abs(forag_doy - survey_doy)) %>%
  group_by(foragdata_id) %>%
  slice_min(diff_days, with_ties = FALSE) %>%
  ungroup() %>%
  st_drop_geometry() %>%
  filter(survey_year %in% years_keep)

################################################################################
# 4. Summarise biomass and predictors ------------------------------------------
################################################################################

biomass_patch <- pur_patch %>%
  group_by(survey_year, patch_id, pred_patch) %>%
  summarise(total_biomass_g = sum(total_biomass_g, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(foraging_biomass_kg = total_biomass_g / 1000,
         year = survey_year)

predictors_focus <- c("purple_urchin_densitym2", "mean_gonad_mass_g")

patch_predictors <- quad_build3 %>%
  mutate(year = year(survey_date)) %>%
  st_drop_geometry() %>%
  group_by(patch_id, year, pred_patch) %>%
  summarise(across(all_of(predictors_focus), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

model_data <- left_join(biomass_patch, patch_predictors,
                        by = c("patch_id", "year", "pred_patch")) %>%
  filter(!is.na(foraging_biomass_kg)) %>%
  mutate(pred_patch = factor(pred_patch, levels = c("BAR","INCIP","FOR")))

################################################################################
# 5. Random Forest model -------------------------------------------------------
################################################################################

rf_dat <- model_data %>%
  select(pred_patch, purple_urchin_densitym2, mean_gonad_mass_g, foraging_biomass_kg) %>%
  drop_na()

set.seed(42)
rf_patch <- randomForest(
  pred_patch ~ purple_urchin_densitym2 + mean_gonad_mass_g + foraging_biomass_kg,
  data = rf_dat, ntree = 1000, mtry = 2, importance = TRUE
)

################################################################################
# 6. SOFA posteriors (Panels A & B) --------------------------------------------
################################################################################

sofa_eta <- tribble(
  ~patch_type, ~mean, ~sd,
  "BAR", 0.1126, 0.0220,
  "FOR", 0.2139, 0.0311,
  "INCIP", 0.0241, 0.0098
)
sofa_CR <- tribble(
  ~patch_type, ~mean, ~sd,
  "BAR", 13.0056, 1.4666,
  "FOR", 13.5477, 1.2696,
  "INCIP", 16.8748, 3.0836
)

expand_posterior <- function(tbl, n = 4000) {
  tbl %>% rowwise() %>% mutate(draws = list(rnorm(n, mean, sd))) %>% unnest(draws)
}

posterior_summary <- function(df) {
  df %>%
    group_by(patch_type) %>%
    summarise(
      med = median(draws),
      l80 = quantile(draws, 0.10),
      u80 = quantile(draws, 0.90),
      l90 = quantile(draws, 0.05),
      u90 = quantile(draws, 0.95),
      l95 = quantile(draws, 0.025),
      u95 = quantile(draws, 0.975),
      .groups = "drop"
    )
}

eta_draws <- expand_posterior(sofa_eta)
CR_draws  <- expand_posterior(sofa_CR)
eta_sum <- posterior_summary(eta_draws)
CR_sum  <- posterior_summary(CR_draws)

theme_post <- theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold", size = 11),
        legend.position = "none",
        axis.title.y = element_blank())

plot_posterior <- function(draw_df, sum_df, xlab, tag_title) {
  ggplot() +
    geom_density_ridges(
      data = draw_df,
      aes(x = draws, y = patch_type, fill = patch_type),
      alpha = 0.85, scale = 1.1, color = "black"
    ) +
    # 90% credible interval
    geom_errorbar(
      data = sum_df,
      aes(y = patch_type, xmin = l90, xmax = u90),
      orientation = "y",
      height = 0.25, linewidth = 1.2, color = "black"
    ) +
    # median point
    geom_point(
      data = sum_df,
      aes(y = patch_type, x = med),
      color = "white", fill = "black", shape = 21, size = 3, stroke = 1
    ) +
    scale_fill_manual(values = patch_colors) +
    labs(x = xlab, y = NULL, tag = tag_title) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text  = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = 9),
      legend.position = "none",
      axis.text.y = element_blank()
    )
}



pA <- plot_posterior(eta_draws, eta_sum, "Foraging effort (proportion of total)", "A")
pB <- plot_posterior(CR_draws, CR_sum, "Urchin consumption rate (g min⁻¹)", "B")

################################################################################
# 7. Random Forest PDP (Panel C) -----------------------------------------------
################################################################################

make_pdp <- function(var) {
  pd_all <- lapply(levels(model_data$pred_patch), function(cl) {
    pdp::partial(
      rf_patch,
      pred.var = var,
      which.class = cl,
      prob = TRUE,
      grid.resolution = 150
    ) %>%
      mutate(patch_type = cl)
  }) %>% bind_rows()
  
  ggplot(pd_all, aes_string(x = var, y = "yhat", color = "patch_type")) +
    geom_smooth(se = FALSE, method = "loess", span = 0.6, linewidth = 1.2) +
    scale_color_manual(
      values = patch_colors,
      name = "Patch type",
      labels = c("BAR" = "Barren", "INCIP" = "Incipient", "FOR" = "Forest")
    ) +
    labs(
      x = "Consumed urchin biomass (kg patch⁻¹ yr⁻¹)",
      y = "Predicted P(patch type)",
      tag = "C"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text  = element_text(color = "black"),
      legend.position = "right"
    )
}


pC <- make_pdp("foraging_biomass_kg")

################################################################################
# 8. Combine (A–C) -------------------------------------------------------------
################################################################################



energetics_combined <- pA + pB + pC +
  plot_layout(ncol = 3, widths = c(1, 1, 1.1)) &
  theme(plot.margin = unit(c(1, 1, 1, 1), "pt"))

energetics_combined

ggsave(
  filename = here::here("figures", "Fig4_energetics_combined.png"),
  plot = energetics_combined,
  width = 8,        # in inches
  height = 3,        # adjust as needed
  dpi = 600,         # high-res for publication
  bg = "white"       # ensures white background if saving to PNG
)



