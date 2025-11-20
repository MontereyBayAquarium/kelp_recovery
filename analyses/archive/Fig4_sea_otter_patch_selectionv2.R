################################################################################
# Sea otter foraging vs. patch type
# Row 1 (A–C): Partial dependence of patch type probability
# Row 2 (D–F): SOFA foraging effort, consumption, and energy intake
#
# Author: Joshua G. Smith — UCSC Nearshore Ecology Research Group
################################################################################

rm(list = ls())

require(librarian)
librarian::shelf(
  tidyverse, janitor, lubridate, sf, here,
  randomForest, pdp, patchwork, ggpubr, ggridges
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
  mutate(
    dives_for_pur = n_distinct(foragdiv_id[str_detect(prey, regex("^pur$", ignore_case = TRUE))]),
    focal_bout    = if_else(dives_for_pur > 5, "yes", "no")
  ) %>%
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
  mutate(
    survey_year = year(survey_date),
    survey_doy  = yday(survey_date)
  )

pur_spatial <- st_join(
  pur_sf,
  quad4_same_crs %>% select(patch_id, pred_patch, survey_date, survey_year, survey_doy),
  join = st_intersects,
  left = TRUE
)

pur_patch <- pur_spatial %>%
  mutate(
    forag_doy = yday(date),
    diff_days = abs(forag_doy - survey_doy)
  ) %>%
  group_by(foragdata_id) %>%
  slice_min(diff_days, with_ties = FALSE) %>%
  ungroup() %>%
  st_drop_geometry() %>%
  filter(survey_year %in% years_keep)

################################################################################
# 4. Aggregate urchin biomass per patch-year -----------------------------------
################################################################################

biomass_patch <- pur_patch %>%
  group_by(survey_year, patch_id, pred_patch) %>%
  summarise(
    total_biomass_g   = sum(total_biomass_g, na.rm = TRUE),
    n_foraging_obs    = n(),
    .groups = "drop"
  ) %>%
  mutate(foraging_biomass_kg = total_biomass_g / 1000) %>%
  rename(year = survey_year)

################################################################################
# 5. Merge with benthic patch predictors ---------------------------------------
################################################################################

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
# 6. Random Forest classifier ---------------------------------------------------
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
# 7. Partial dependence plots (A–C) --------------------------------------------
################################################################################

my_theme <- theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color="black"),
        legend.position = "none")

make_pdp <- function(var) {
  pd_all <- lapply(levels(rf_dat$pred_patch), function(cl) {
    pdp::partial(rf_patch, pred.var = var, which.class = cl, prob = TRUE, grid.resolution = 150) %>%
      mutate(patch_type = cl)
  }) %>% bind_rows()
  
  labels <- c(
    "purple_urchin_densitym2" = "Purple urchin density (No. m⁻²)",
    "mean_gonad_mass_g"       = "Mean gonad mass (g)",
    "foraging_biomass_kg"     = "Sea otter urchin take (kg patch⁻¹ yr⁻¹)"
  )
  
  ggplot(pd_all, aes_string(x = var, y = "yhat", color = "patch_type")) +
    geom_smooth(se = FALSE, method = "loess", span = 0.6, linewidth = 1.2) +
    scale_color_manual(values = patch_colors, name = "Patch type") +
    labs(x = labels[[var]], y = NULL) +
    my_theme
}

pA <- make_pdp("purple_urchin_densitym2") + labs(tag = "A", y = "Predicted P(patch type)")
pB <- make_pdp("mean_gonad_mass_g") + labs(tag = "B")
pC <- make_pdp("foraging_biomass_kg") + labs(tag = "C")

legend_plot <- make_pdp("purple_urchin_densitym2") +
  theme(legend.position = "bottom")
legend_grob <- get_legend(legend_plot)

top_row_main <- ggarrange(pA, pB, pC, ncol = 3, align = "hv")
top_row <- ggarrange(top_row_main, legend_grob, ncol = 1, heights = c(1, 0.12))

################################################################################
# 8. SOFA posterior ridgeplots (D–F) -------------------------------------------
################################################################################

# Posterior summaries (from your SOFA output)
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
sofa_ER <- tribble(
  ~patch_type, ~mean, ~sd,
  "BAR", 10.0613, 1.1407,
  "FOR", 10.4841, 0.9878,
  "INCIP", 13.0557, 2.3916
)

set.seed(123)
expand_posterior <- function(tbl, n = 4000) {
  tbl %>% rowwise() %>% mutate(draws = list(rnorm(n, mean, sd))) %>% unnest(draws)
}
eta_draws <- expand_posterior(sofa_eta)
CR_draws  <- expand_posterior(sofa_CR)
ER_draws  <- expand_posterior(sofa_ER)

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

eta_sum <- posterior_summary(eta_draws)
CR_sum  <- posterior_summary(CR_draws)
ER_sum  <- posterior_summary(ER_draws)

theme_post <- theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color="black"),
        plot.title = element_text(face="bold", size=11),
        legend.position="none",
        axis.title.y = element_blank())

plot_posterior <- function(draw_df, sum_df, xlab, tag_title) {
  ggplot(draw_df, aes(x = draws, y = patch_type, fill = patch_type)) +
    geom_density_ridges(alpha = 0.85, scale = 1.1, color = "black") +
    geom_errorbar(data = sum_df, aes(y = patch_type, xmin = l95, xmax = u95),
                  orientation = "y", height = 0.25, linewidth = 1, color = "black",
                  inherit.aes = FALSE) +
    geom_errorbar(data = sum_df, aes(y = patch_type, xmin = l90, xmax = u90),
                  orientation = "y", height = 0.25, linewidth = 1.8, color = "gray30",
                  alpha = 0.5, inherit.aes = FALSE) +
    geom_errorbar(data = sum_df, aes(y = patch_type, xmin = l80, xmax = u80),
                  orientation = "y", height = 0.25, linewidth = 2.2, color = "gray50",
                  alpha = 0.4, inherit.aes = FALSE) +
    geom_point(data = sum_df, aes(y = patch_type, x = med),
               color = "white", fill = "black", shape = 21, size = 3, stroke = 1,
               inherit.aes = FALSE) +
    scale_fill_manual(values = patch_colors) +
    labs(x = xlab, y = NULL, tag = tag_title) +
    theme_post+
    theme(axis.text.y = element_blank())
}

pD <- plot_posterior(eta_draws, eta_sum, "Foraging effort (proportion of total)", "D")
pE <- plot_posterior(CR_draws, CR_sum, "Urchin consumption rate (g min⁻¹)", "E")
pF <- plot_posterior(ER_draws, ER_sum, "Energy intake rate (kJ min⁻¹)", "F")

bottom_row <- ggarrange(pD, pE, pF, ncol = 3, align = "hv")

################################################################################
# 9. Combine and export --------------------------------------------------------
################################################################################

full_fig <- ggarrange(top_row, bottom_row, ncol = 1, heights = c(1.2, 1))

final_fig <- annotate_figure(
  full_fig,
  top = text_grob("Sea otter patch selection across benthic state space",
                  face = "bold", size = 14),
  bottom = text_grob(
    "Top (A–C): Partial dependence of predicted patch type probability (Random Forest).\nBottom (D–F): Posterior distributions (median ± 80/90/95% CrI) for SOFA-estimated foraging effort, consumption rate, and energy intake.",
    size = 10
  )
)

print(final_fig)

# ggsave(here::here("figures","Fig_patch_selection_RF_SOFA_posteriors_CrIs.png"),
#        final_fig, width = 10, height = 8, dpi = 600)
