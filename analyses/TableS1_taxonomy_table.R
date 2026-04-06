################################################################################
# Supplementary Table: Variables recorded during diver-based benthic surveys
#
# Output:
#   tables/Table_SX_recorded_variables.html
#   tables/Table_SX_recorded_variables.csv
################################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

require(librarian)
librarian::shelf(
  tidyverse, here, janitor, gt, sf
)

################################################################################
# 1. Load data -----------------------------------------------------------------
################################################################################

load(here::here("output", "survey_data", "processed", "zone_level_data4.rda"))

if (!exists("quad_build3")) {
  stop("quad_build3 was not found in zone_level_data4.rda")
}

dat <- quad_build3 %>%
  sf::st_drop_geometry() %>%
  janitor::clean_names()

################################################################################
# 2. Taxon label lookup --------------------------------------------------------
################################################################################
# Edit any labels here if you want different taxonomic wording in the table.

taxon_lookup <- c(
  # Density variables (20 m2)
  "ptecal"     = "Pterygophora californica",
  "eisarb"     = "Eisenia arborea",
  "nerlue"     = "Nereocystis luetkeana",
  "lamset"     = "Laminaria setchellii",
  "cancer_spp" = "Cancer spp.",
  "lamstump"   = "Laminarian holdfast / stump",
  "macstump"   = "Macrocystis holdfast / stump",
  "purps_on_kelp" = "Purple urchins on kelp",
  
  # Recruits
  "lamr"  = "Laminaria spp. recruits",
  "macr"  = "Macrocystis pyrifera recruits",
  "macj"  = "Macrocystis pyrifera juveniles",
  "nerj"  = "Nereocystis luetkeana juveniles",
  "ptej"  = "Pterygophora californica juveniles",
  "lsetj" = "Laminaria setchellii juveniles",
  "eisj"  = "Eisenia arborea juveniles",
  
  # Cover variables
  "crustose_coralline"      = "Crustose coralline algae",
  "desmarestia_spp"         = "Desmarestia spp.",
  "dodecaceria_spp"         = "Dodecaceria spp.",
  "lam_holdfast_live"       = "Live laminaria spp. holdfast",
  "articulated_coralline"   = "Articulated coralline algae",
  "shell_debris"            = "Shell debris",
  "fleshy_red"              = "Fleshy red algae",
  "encrusting_red"          = "Encrusting red algae",
  "bare_rock"               = "Bare rock",
  "stephanocystis"          = "Stephanocystis spp.",
  "green_algae"             = "Green algae",
  "barnacle"                = "Barnacles",
  "phragmatopoma"           = "Phragmatopoma",
  "cup_corals"              = "Cup corals",
  "bare_sand"               = "Bare sand",
  "diopatra_chaetopterus"   = "Diopatra / Chaetopterus tubes",
  "colonial_tunicate"       = "Colonial tunicates",
  "tubeworm_other_solitary" = "Other solitary tubeworms",
  "dictyota_dictyopteris"   = "Dictyota / Dictyopteris",
  "sponge"                  = "Sponges",
  "other_anemone"           = "Other anemones",
  "hydroids"                = "Hydroids",
  "red_turf_2_cm"           = "Red turf algae (< 2 cm)",
  "dead_kelp_holdfast_any"  = "Dead kelp holdfast",
  "mac_holdfast_live"       = "Live Macrocystis holdfast",
  "salmacina"               = "Salmacina",
  "corynactis_californica"  = "Corynactis californica",
  "dictyoneurum_spp"        = "Dictyoneurum spp.",
  "solitary_tunicate"       = "Solitary tunicates",
  "bryozoan"                = "Bryozoans",
  "embedded_cucumber"       = "Embedded cucumbers",
  "sediment_mud"            = "Mud / fine sediment",
  "thylacodes_orange_rim"   = "Thylacodes spp.",
  "mussel"                  = "Mussels",
  "egregia_menziesii"       = "Egregia menziesii",
  "petaloconchus"           = "Petaloconchus"
)

################################################################################
# 3. Helper functions ----------------------------------------------------------
################################################################################

nice_title <- function(x) {
  x %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_squish() %>%
    stringr::str_to_sentence()
}

lookup_taxon <- function(x) {
  dplyr::coalesce(unname(taxon_lookup[x]), nice_title(x))
}

################################################################################
# 4. Physical + structural variables -------------------------------------------
################################################################################

physical_table <- tibble::tribble(
  ~Category,            ~Variable,                     ~`Taxon / group`,       ~Unit,
  "Physical reef",      "Vertical relief",             "Reef structure",       "cm",
  "Physical reef",      "Seafloor rugosity",           "Reef structure",       "index",
  "Urchin condition",   "Total urchin biomass",        "Purple sea urchin",    "g",
  "Urchin condition",   "Mean gonad index",            "Purple sea urchin",    "ratio",
  "Urchin condition",   "Mean individual biomass",     "Purple sea urchin",    "g",
  "Urchin condition",   "Mean gonad mass",             "Purple sea urchin",    "g",
  "Urchin condition",   "Total gonad mass",            "Purple sea urchin",    "g",
)

################################################################################
# 5. Densities of mobile invertebrates -----------------------------------------
################################################################################

density_m2_table <- tibble::tribble(
  ~raw,                            ~Category,         ~Variable,                           ~`Taxon / group`,          ~Unit,
  "purple_urchin_densitym2",       "Mobile invertebrates", "Density",                     "Purple sea urchin",       "ind. m²",
  "purple_urchin_conceiledm2",     "Mobile invertebrates", "Concealed density",           "Purple sea urchin",       "ind. m²",
  "red_urchin_densitym2",          "Mobile invertebrates", "Density",                     "Red sea urchin",          "ind. m²",
  "red_urchin_conceiledm2",        "Mobile invertebrates", "Concealed density",           "Red sea urchin",          "ind. m²",
  "tegula_densitym2",              "Mobile invertebrates", "Density",                     "Tegula spp.",             "ind. m²",
  "pomaulax_densitym2",            "Mobile invertebrates", "Density",                     "Pomaulax spp.",           "ind. m²"
) %>%
  dplyr::select(-raw)

################################################################################
# 6. Macroalgal / invertebrate density variables per 20 m2 ---------------------
################################################################################

dens20_vars <- names(dat) %>%
  stringr::str_subset("^density20m2_")

dens20_table <- tibble(raw = dens20_vars) %>%
  dplyr::mutate(
    code = stringr::str_remove(raw, "^density20m2_"),
    Category = "Large benthic taxa",
    Variable = "Density",
    `Taxon / group` = lookup_taxon(code),
    Unit = "ind. 20 m²"
  ) %>%
  dplyr::select(Category, Variable, `Taxon / group`, Unit)

################################################################################
# 7. Recruit / juvenile variables ----------------------------------------------
################################################################################

recruit_vars <- c("lamr", "macr", "macj", "nerj", "ptej", "lsetj", "eisj")
recruit_vars <- recruit_vars[recruit_vars %in% names(dat)]

recruit_table <- tibble(raw = recruit_vars) %>%
  dplyr::mutate(
    Category = "Recruitment",
    Variable = "Recruit / juvenile abundance",
    `Taxon / group` = lookup_taxon(raw),
    Unit = "count"
  ) %>%
  dplyr::select(Category, Variable, `Taxon / group`, Unit)

################################################################################
# 8. Cover variables -----------------------------------------------------------
################################################################################

cover_vars <- names(dat) %>%
  stringr::str_subset("^cov_")

cover_table <- tibble(raw = cover_vars) %>%
  dplyr::mutate(
    code = stringr::str_remove(raw, "^cov_"),
    Category = "Benthic cover",
    Variable = "Percent cover",
    `Taxon / group` = lookup_taxon(code),
    Unit = "% cover"
  ) %>%
  dplyr::select(Category, Variable, `Taxon / group`, Unit)

################################################################################
# 9. Combine full table --------------------------------------------------------
################################################################################

full_table <- dplyr::bind_rows(
  physical_table,
  density_m2_table,
  dens20_table,
  recruit_table,
  cover_table
) %>%
  dplyr::mutate(
    Category = factor(
      Category,
      levels = c(
        "Physical reef",
        "Mobile invertebrates",
        "Urchin condition",
        "Large benthic taxa",
        "Recruitment",
        "Benthic cover"
      )
    )
  ) %>%
  dplyr::arrange(Category, Variable, `Taxon / group`)

################################################################################
# 10. Write CSV ----------------------------------------------------------------
################################################################################

dir.create(here::here("tables"), recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  full_table,
  here::here("tables", "Table_S1_taxonomy_table.csv")
)

################################################################################
# 11. Make formatted gt table --------------------------------------------------
################################################################################

gt_tbl <- full_table %>%
  gt::gt(groupname_col = "Category") %>%
  gt::cols_label(
    Variable = "Variable",
    `Taxon / group` = "Taxon / group",
    Unit = "Unit"
  ) %>%
  gt::cols_align(
    align = "left",
    columns = c(Variable, `Taxon / group`)
  ) %>%
  gt::cols_align(
    align = "center",
    columns = Unit
  ) %>%
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_row_groups()
  ) %>%
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_column_labels(everything())
  ) %>%
  gt::tab_options(
    table.font.size = 11,
    data_row.padding = gt::px(4)
  )

gt_tbl

################################################################################
# 12. Save formatted table -----------------------------------------------------
################################################################################

gt::gtsave(
  gt_tbl,
  filename = here::here("tables", "Table_S1_taxonomy_table.html")
)

# Optional PNG export:
# install.packages("webshot2")
# gt::gtsave(
#   gt_tbl,
#   filename = here::here("tables", "Table_SX_recorded_variables.png")
# )