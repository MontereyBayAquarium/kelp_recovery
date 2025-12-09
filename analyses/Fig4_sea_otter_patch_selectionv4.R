################################################################################
# Sea otter energetics and prey use — Panels A & B
#   Panel A: SOFA posterior – PUR foraging effort across patch types
#   Panel B: Focal bouts by prey group and patch type
# Author: Joshua G. Smith — UCSC Nearshore Ecology Research Group
################################################################################

rm(list = ls())

require(librarian)
librarian::shelf(
  tidyverse, janitor, readxl, here,
  ggridges, patchwork, sf, lubridate
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
      # order for ridge y-axis (FOR, BAR, INCIP)
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

expand_posterior <- function(tbl, n = 5000) {
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
# 6. Plot function for posterior ridgelines ------------------------------------
################################################################################

patch_colors <- c(
  "BAR"   = "#7570B3",
  "INCIP" = "#D95F02",
  "FOR"   = "#1B9E77"
)

plot_posterior <- function(draw_df, sum_df, xlab, tag_title, show_legend = FALSE) {
  
  p <- ggplot() +
    geom_density_ridges(
      data = draw_df,
      aes(x = draws, y = patch_type, fill = patch_type),
      alpha = 0.85, scale = 1.1, color = "black"
    ) +
    geom_errorbar(
      data = sum_df,
      aes(y = patch_type, xmin = l90, xmax = u90),
      orientation = "y", height = 0.25,
      linewidth = 1.2, color = "black"
    ) +
    geom_point(
      data = sum_df,
      aes(y = patch_type, x = med),
      color = "white", fill = "black",
      shape = 21, size = 3, stroke = 1
    ) +
    scale_fill_manual(values = patch_colors, name = "Patch type") +
    labs(x = xlab, y = NULL, tag = tag_title) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid   = element_blank(),
      axis.text    = element_text(color = "black"),
      axis.text.y  = element_blank(),
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8),
      # match tag styling to Panel B
      plot.tag     = element_text(size = 10, color = "black")
    )
  
  # Legend toggle
  if (!show_legend) {
    p <- p + theme(legend.position = "right")
  }
  
  return(p)
}

################################################################################
# 7. Panel A: PUR foraging effort across patch types ---------------------------
################################################################################

pA <- plot_posterior(
  eta_draws,
  eta_sum,
  xlab      = "Foraging effort (proportion of sea \nurchins in patch vs. all other prey)",
  tag_title = "A"
)

################################################################################
# 8. Panel B: Focal bouts by prey group and patch type -------------------------
################################################################################

# Load quad/foraging data used for patch-level focal bout summaries
load(here::here("output", "survey_data", "processed", "zone_level_data4.rda"))
forage_orig <- readr::read_csv(
  "/Volumes/enhydra/data/foraging_data/processed/foraging_data_2024_2025_processed.csv"
)

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

# Convert to sf and join to patch polygons
forage_sf <- st_as_sf(
  forage_locs,
  coords = c("long", "lat"),
  crs    = 4326,
  remove = FALSE
)

quad3_same_crs <- st_transform(quad_build3, st_crs(forage_sf))

# Spatial join to assign patch type
focal_joined <- st_join(
  forage_sf,
  quad3_same_crs,
  join = st_intersects,
  left = TRUE
) %>%
  mutate(
    pred_patch = case_when(
      str_detect(tolower(pred_patch), "bar")   ~ "Barren",
      str_detect(tolower(pred_patch), "incip") ~ "Incipient",
      str_detect(tolower(pred_patch), "for")   ~ "Forest",
      TRUE ~ NA_character_
    )
  ) %>%
  st_drop_geometry() %>%
  filter(!is.na(pred_patch))

# Group prey types into broader categories
focal_joined <- focal_joined %>%
  mutate(prey_group = case_when(
    prey %in% c("cam", "clm", "mus")               ~ "Bivalve",
    prey %in% c("cra", "kcr", "can", "rcr")        ~ "Cancer crab",
    prey %in% c("gas", "teg", "whe")               ~ "Gastropod",
    prey %in% c("pur", "red", "urc")               ~ "Sea urchin",
    TRUE                                           ~ "Other invertebrate"
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
    .groups    = "drop"
  ) %>%
  group_by(pred_patch) %>%
  mutate(prop = mean_bouts / sum(mean_bouts)) %>%
  ungroup() %>%
  # enforce order: Forest, Barren, Incipient
  mutate(
    pred_patch = factor(
      pred_patch,
      levels = c("Forest", "Barren", "Incipient")
    )
  )

prey_cols <- c(
  "Bivalve"            = "#1F78B4",   # teal blue
  "Cancer crab"        = "#FDBF6F",   # golden orange
  "Gastropod"          = "forestgreen",
  "Sea urchin"         = "purple",
  "Other invertebrate" = "#B2B2B2"    # neutral gray
)

# Custom theme (matches Josh's general style)
my_theme <- theme(
  axis.text.x      = element_text(size = 8, color = "black"),
  axis.text.y      = element_text(size = 8, color = "black"),
  axis.title       = element_text(size = 10, color = "black"),
  legend.text      = element_text(size = 8, color = "black"),
  legend.title     = element_text(size = 8, color = "black"),
  plot.tag         = element_text(size = 10, color = "black"),
  panel.grid       = element_blank(),
  panel.background = element_blank(),
  axis.line        = element_line(colour = "black"),
  legend.key       = element_blank()
)

# Panel B: focal bouts by prey group
pB <- ggplot(
  focal_summary,
  aes(x = pred_patch, y = prop, fill = prey_group)
) +
  geom_col(color = "black", width = 0.7, position = "stack") +
  scale_fill_manual(values = prey_cols, name = "Prey group") +
  labs(
    x   = "Patch type",
    y   = "Proportion of focal bouts (>3 successful dives)",
    tag = "B"
  ) +
  theme_bw() +
  theme(
    panel.grid     = element_blank(),
    legend.position = "right"
  ) +
  my_theme

################################################################################
# 9. Combine Panel A + Panel B -------------------------------------------------
################################################################################

energetics_combined <- pA + pB +
  plot_layout(ncol = 2, widths = c(1, 1)) &
  theme(plot.margin = unit(c(2, 2, 2, 2), "pt"))

energetics_combined

ggsave(
  here::here("figures", "Fig4_effort_allocation.png"),
  energetics_combined,
  width = 8, height = 4, dpi = 600, bg = "white"
)
