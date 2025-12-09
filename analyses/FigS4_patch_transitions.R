
#jogsmith@ucsc.edu

rm(list = ls())


################################################################################
#Step 0: set paths and load data
require(librarian)
librarian::shelf(tidyverse, lubridate, sf, stringr, purrr, terra, janitor,
                 rnaturalearth, rnaturalearthdata)

datadir <- "/Volumes/enhydra/data/kelp_recovery/"
localdir <- here::here("output")

#load benthic survey data
load(file.path(localdir, "lda_patch_transitionsv3.rda"))


str(transitions_tbl_constrained)

################################################################################
#Figure S4


# Factor levels + pretty labels
patch_levels  <- c("BAR",  "FOR",   "INCIP")
patch_labels  <- c("Barren","Forest","Incipient")

# Build 3x3 transition matrix
trans_mat <- transitions_tbl_constrained %>%
  mutate(
    patch_2024 = factor(patch_2024,
                        levels = patch_levels,
                        labels = patch_labels),
    patch_2025 = factor(as.character(patch_2025),
                        levels = patch_levels,
                        labels = patch_labels)
  ) %>%
  count(patch_2024, patch_2025, name = "n") %>%
  tidyr::complete(patch_2025, patch_2024, fill = list(n = 0))



p_transitions <- ggplot(trans_mat, aes(x = patch_2024, y = patch_2025, fill = n)) +
  geom_tile(color = "grey40") +
  geom_text(aes(label = n), size = 4) +
  coord_equal() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    x = "Patch type in 2024",
    y = "Patch type in 2025"
  ) +
  theme_bw() +
  theme(
    panel.grid       = element_blank(),
    legend.position  = "none",      # REMOVE LEGEND ✔
    axis.text        = element_text(size = 9, color = "black"),
    axis.title       = element_text(size = 10),
    strip.background = element_blank()
  )

p_transitions



ggsave(p_transitions, file = file.path(here::here("figures","S4_patch_transitions.png")), width = 6.5,
       height = 6.5, units = "in")




