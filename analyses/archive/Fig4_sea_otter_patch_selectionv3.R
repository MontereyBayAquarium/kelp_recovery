################################################################################
# Sea otter energetics — Panels A & B using SOFA posterior summaries
# Author: Joshua G. Smith — UCSC Nearshore Ecology Research Group
################################################################################

rm(list = ls())

require(librarian)
librarian::shelf(
  tidyverse, janitor, readxl, here,
  ggridges, patchwork
)

################################################################################
# 1. Load SOFA data ------------------------------------------------------------
################################################################################

sofa_path <- here::here("output", "sofa", "SOFA_Results_2025-10-31.xlsx")

eta_df <- read_xlsx(sofa_path, sheet = "propforagingeffort")
cr_df  <- read_xlsx(sofa_path, sheet = "meanbiomassconsumption")

################################################################################
# 2. Filter to purple urchins and map to patch types ---------------------------
################################################################################

recode_patches <- function(df) {
  df %>%
    filter(prey_name %in% c("urchin_b", "urchin_f", "urchin_i")) %>%
    mutate(
      patch_type = case_when(
        prey_name == "urchin_b" ~ "BAR",
        prey_name == "urchin_f" ~ "FOR",
        prey_name == "urchin_i" ~ "INCIP"
      ),
      patch_type = factor(patch_type, levels = c("FOR", "BAR", "INCIP"))
    )
}

eta_focus <- recode_patches(eta_df)
cr_focus  <- recode_patches(cr_df)

################################################################################
# 3. Remove low-value duplicate summary rows (keep true consumption rates) -----
################################################################################

cr_focus <- cr_focus %>%
  group_by(patch_type) %>%
  slice_max(mean, n = 1, with_ties = FALSE) %>%
  ungroup()

################################################################################
# 4. Expand posterior draws from mean + SD -------------------------------------
################################################################################

expand_posterior <- function(tbl, n = 4000) {
  tbl %>%
    rowwise() %>%
    mutate(draws = list(rnorm(n, mean, sd))) %>%
    unnest(draws)
}

eta_draws <- expand_posterior(eta_focus)
cr_draws  <- expand_posterior(cr_focus)

################################################################################
# 5. Summarise 90% credible intervals ------------------------------------------
################################################################################

posterior_summary <- function(df) {
  df %>%
    group_by(patch_type) %>%
    summarise(
      med = median(draws),
      l90 = quantile(draws, 0.05),
      u90 = quantile(draws, 0.95),
      .groups = "drop"
    )
}

eta_sum <- posterior_summary(eta_draws)
cr_sum  <- posterior_summary(cr_draws)

################################################################################
# 6. Plot function -------------------------------------------------------------
################################################################################

patch_colors <- c("BAR"="purple", "INCIP"="orange", "FOR"="forestgreen")

plot_posterior <- function(draw_df, sum_df, xlab, tag_title) {
  ggplot() +
    geom_density_ridges(
      data = draw_df,
      aes(x = draws, y = patch_type, fill = patch_type),
      alpha = 0.85, scale = 1.1, color = "black"
    ) +
    # 90% credible interval bars
    geom_errorbar(
      data = sum_df,
      aes(y = patch_type, xmin = l90, xmax = u90),
      orientation = "y",
      height = 0.25,
      linewidth = 1.2,
      color = "black"
    ) +
    # Median points
    geom_point(
      data = sum_df,
      aes(y = patch_type, x = med),
      color = "white", fill = "black",
      shape = 21, size = 3, stroke = 1
    ) +
    scale_fill_manual(values = patch_colors) +
    labs(x = xlab, y = NULL, tag = tag_title) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.y = element_blank(),
      legend.position = "none"
    )
}

################################################################################
# 7. Build plots ---------------------------------------------------------------
################################################################################

pA <- plot_posterior(eta_draws, eta_sum, "Foraging effort (proportion of total)", "A")
pB <- plot_posterior(cr_draws,  cr_sum,  "Urchin consumption rate (g min⁻¹)",    "B")

################################################################################
# 8. Combine A + B -------------------------------------------------------------
################################################################################

energetics_combined <- pA + pB +
  plot_layout(ncol = 2, widths = c(1, 1)) &
  theme(plot.margin = unit(c(2, 2, 2, 2), "pt"))

energetics_combined
