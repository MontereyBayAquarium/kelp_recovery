################################################################################
# Patch-State Habitat Classification 
# Joshua G. Smith; jogsmith@ucsc.edu
################################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

require(librarian)
librarian::shelf(
  tidyverse, janitor, lubridate, sf, here,
  randomForest, patchwork, ggrepel, grid
)

set.seed(1985)

################################################################################
# 1. Load Data -----------------------------------------------------------------
################################################################################

load(here::here("output", "survey_data", "processed", "zone_level_data4.rda"))

patch_colors <- c("BAR"="#7570B3","INCIP"="#D95F02","FOR"="#1B9E77")

# Pretty variable labels
pretty_labs <- c(
  cov_fleshy_red = "Fleshy red algae cover",
  cov_mac_holdfast_live = "Macrocystis cover",
  cov_phragmatopoma = "Phragmatopoma spp. cover",
  density20m2_macstump = "Macrocystis dead holdfast cover",
  macj = "Juvenile Macrocystis density",
  ptej = "Juvenile Pterygophora density",
  cov_diopatra_chaetopterus = "Diopatra spp. cover",
  density20m2_lamset = "Laminaria setchelli density",
  cov_crustose_coralline = "Crustose coralline cover",
  tegula_densitym2 = "Tegula spp. density",
  n_macro_plants_20m2 = "Macrocystis density",
  density20m2_ptecal = "Pterygophora density",
  macro_stipe_density_20m2 = "Macrocystis stipe density",
  cov_demarestia_spp = "Desmarestia spp. cover",
  density20m2_nerlue = "Nereocystis luetkeana density",
  cov_lam_holdfast_live = "Laminariales holdfast cover",
  cov_articulated_coralline = "Articualted coralline cover",
  pomaulax_densitym2 = "Pomaulax spp. density",
  cov_encrusting_red = "Encrusting red algae cover",
  relief_cm = "Vertical relief",
  lamr = "Laminariales recruit density",
  cov_bare_sand = "Band sand cover",
  risk_index = "Reef rugosity",
  cov_desmarestia_spp = "Demarestia spp. cover"
)

################################################################################
# 2. Helper Functions ----------------------------------------------------------
################################################################################

clean_state_label <- function(x){
  x_up <- toupper(x)
  dplyr::case_when(
    x_up %in% c("BAR","BARREN","URCHIN BARREN","BARRENS") ~ "BAR",
    x_up %in% c("INCIP","INCIPIENT","TRANSITION","TRANSITIONAL") ~ "INCIP",
    x_up %in% c("FOR","FOREST","KELP FOREST") ~ "FOR",
    TRUE ~ NA_character_
  )
}

mode_char <- function(x){
  x <- x[!is.na(x)]
  if(!length(x)) return(NA_character_)
  names(sort(table(x),decreasing=TRUE))[1]
}

asin_sqrt <- function(x){ p <- pmin(pmax(x,0),100)/100; asin(sqrt(p)) }

med_impute <- function(x){
  if(!is.numeric(x)) return(x)
  x[!is.finite(x)] <- NA_real_
  if(all(is.na(x))) return(x)
  x[is.na(x)] <- median(x,na.rm=TRUE)
  x
}

drop_nonfinite_cols <- function(df){
  keep <- vapply(df,function(col) is.numeric(col)&&any(is.finite(col)),logical(1))
  df[,keep,drop=FALSE]
}

################################################################################
# 3. Habitat-only RF: Train on 2024, Predict 2025 ------------------------------
################################################################################

dat_raw <- quad_build3 %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(year = lubridate::year(survey_date))

truth_all <- dat_raw %>%
  dplyr::mutate(site_type_clean = clean_state_label(site_type)) %>%
  dplyr::filter(year %in% c(2024,2025)) %>%
  dplyr::group_by(patch_id,site,zone,year) %>%
  dplyr::summarise(state = mode_char(site_type_clean), .groups="drop") %>%
  dplyr::mutate(state = factor(state, levels=c("BAR","INCIP","FOR")))

ban_regex <- "(urchin|gonad|biomass|forag|gi|density20m2_otter)"
num_cols_all <- names(dat_raw)[vapply(dat_raw, is.numeric, logical(1))]
allowed_numeric <- num_cols_all[!grepl(ban_regex,num_cols_all,ignore.case=TRUE)]
allowed_numeric <- setdiff(allowed_numeric, c("year","latitude","longitude"))

agg_hab <- dat_raw %>%
  dplyr::group_by(patch_id,site,zone,year) %>%
  dplyr::summarise(dplyr::across(dplyr::any_of(allowed_numeric),~mean(.x,na.rm=TRUE)),.groups="drop")

meta_hab <- agg_hab %>% dplyr::select(patch_id,site,zone,year)
X_raw <- agg_hab %>% dplyr::select(-patch_id,-site,-zone,-year)

cover_cols <- names(X_raw)[grepl("^cov_",names(X_raw))]
other_cols <- setdiff(names(X_raw), cover_cols)

X_tr <- X_raw %>%
  dplyr::mutate(dplyr::across(dplyr::any_of(cover_cols),asin_sqrt),
                dplyr::across(dplyr::any_of(other_cols),~log1p(pmax(.x,0))))

too_na <- vapply(X_tr,function(x)mean(!is.finite(x)|is.na(x)),numeric(1))
X_tr2 <- X_tr %>% dplyr::select(dplyr::any_of(names(too_na)[too_na<=0.5]))
X_tr2 <- X_tr2 %>% dplyr::select(where(is.numeric))
X_imp <- X_tr2 %>% dplyr::mutate(dplyr::across(dplyr::everything(), med_impute))
X_imp <- drop_nonfinite_cols(X_imp)
X_scaled <- scale(X_imp) %>% tibble::as_tibble()
colnames(X_scaled) <- colnames(X_imp)

hab_scaled_df <- dplyr::bind_cols(meta_hab, X_scaled)

model_df <- hab_scaled_df %>%
  dplyr::left_join(truth_all, by=c("patch_id","site","zone","year")) %>%
  dplyr::mutate(state = droplevels(state))

train_df <- dplyr::filter(model_df, year==2024, !is.na(state))
test_df  <- dplyr::filter(model_df, year==2025)

predictor_cols <- setdiff(colnames(train_df),
                          c("patch_id","site","zone","year","state"))

rf_train2024 <- randomForest(
  x=train_df[,predictor_cols],
  y=train_df$state,
  ntree=1500,
  mtry=max(2,floor(sqrt(length(predictor_cols)))),
  importance=TRUE,
  na.action=na.omit
)

cat("\n[PART A] OOB accuracy = ",
    1 - rf_train2024$err.rate[nrow(rf_train2024$err.rate),"OOB"], "\n")

# Predict 2025 patch types
test_pred_class <- predict(rf_train2024, newdata=test_df[,predictor_cols])

# ---------------------------------------------------------------------------
# CONSTRAINT: If 2024 state == FOR, then 2025 prediction cannot be INCIP.
# Force any FOR->INCIP predictions to FOR (forest stays forest).
# ---------------------------------------------------------------------------
prev_states <- truth_all %>%
  dplyr::filter(year == 2024) %>%
  dplyr::select(patch_id, state_2024 = state)

pred_2025 <- test_df %>%
  dplyr::select(patch_id,site,zone,year) %>%
  dplyr::mutate(predicted_state_2025 = as.character(test_pred_class)) %>%
  dplyr::left_join(prev_states, by = "patch_id") %>%
  dplyr::mutate(
    predicted_state_2025 = ifelse(
      state_2024 == "FOR" & predicted_state_2025 == "INCIP", "FOR", predicted_state_2025
    ),
    predicted_state_2025 = factor(predicted_state_2025, levels = c("BAR","FOR","INCIP"))
  ) %>%
  dplyr::select(-state_2024)

################################################################################
#    Exports patch shapes with patch_2024 and patch_2025
################################################################################

# 5a. Diver-called per patch (patch_2024)
patch_calls_tbl <- quad_build3 %>%
  dplyr::mutate(
    patch_call_clean = clean_state_label(site_type)
  ) %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(patch_id, site, zone) %>%
  dplyr::summarise(
    patch_2024 = mode_char(patch_call_clean),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    patch_2024 = as.character(patch_2024)
  )

# 5b. Patch geometry dissolved per patch_id
patch_geom_tbl <- quad_build3 %>%
  dplyr::select(patch_id, site, zone, geometry) %>%
  dplyr::group_by(patch_id, site, zone) %>%
  dplyr::summarise(
    geometry = sf::st_union(geometry),
    .groups  = "drop"
  )

# 5c. Predicted 2025 state (constrained)
patch_pred2025_tbl <- pred_2025 %>%
  dplyr::transmute(
    patch_id,
    patch_2025 = predicted_state_2025
  ) %>%
  dplyr::mutate(
    patch_2025 = factor(patch_2025, levels = c("BAR","FOR","INCIP"))
  )

# 5d. Combine and export both spatial and tabular forms
final_patch_sf <- patch_geom_tbl %>%
  dplyr::left_join(
    patch_calls_tbl,
    by = c("patch_id","site","zone"),
    relationship = "many-to-many",
  ) %>%
  dplyr::left_join(
    patch_pred2025_tbl,
    by = "patch_id",
    relationship = "many-to-many",
  ) %>%
  dplyr::mutate(
    patch_2024 = as.character(patch_2024),
    patch_2025 = factor(patch_2025, levels = c("BAR","FOR","INCIP"))
  ) %>%
  dplyr::distinct(patch_id, .keep_all = TRUE) %>%
  sf::st_as_sf()

# Non-spatial tibble (mirrors transitions_tbl_constrained)
transitions_tbl_constrained <- final_patch_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(
    site_type = factor(patch_2024, levels = c("BAR","FOR","INCIP")),
    patch_2024 = as.character(patch_2024),
    patch_2025 = factor(patch_2025, levels = c("BAR","FOR","INCIP"))
  ) %>%
  dplyr::select(site, zone, site_type, patch_2024, patch_2025)

# Confirm structure
str(transitions_tbl_constrained)

# Save to disk (.rda)
# save(transitions_tbl_constrained,
 #    file = here::here("output","lda_patch_transitionsv4.rda"))

# 5e. Write shapefile for GIS visualization (optional)
# out_dir <- here::here("output", "gis_data", "processed", "patch_state_RFsummary_2024_2025_shp")
# if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
# out_path <- file.path(out_dir, "patch_state_RFsummary_2024_2025.shp")
# sf::st_write(final_patch_sf, out_path, delete_dsn = TRUE)

################################################################################
# 4. Prepare PCA Dataset -------------------------------------------------------
################################################################################

habitat_df <- hab_scaled_df %>%
  dplyr::left_join(truth_all %>% dplyr::filter(year==2024) %>% dplyr::rename(state_2024=state),
                   by=c("patch_id","site","zone","year")) %>%
  dplyr::left_join(pred_2025, by=c("patch_id","site","zone","year")) %>%
  dplyr::mutate(state_final = dplyr::case_when(
    year==2024 ~ as.character(state_2024),
    year==2025 ~ as.character(predicted_state_2025),
    TRUE ~ NA_character_
  ),
  state_final = factor(state_final, levels=c("BAR","FOR","INCIP"))) %>%
  dplyr::filter(!is.na(state_final))

predictor_cols_pca <- setdiff(names(habitat_df),
                              c("patch_id","site","zone","year",
                                "state_2024","predicted_state_2025","state_final"))

pca_obj <- prcomp(habitat_df %>% dplyr::select(dplyr::all_of(predictor_cols_pca)) %>% as.matrix(),
                  center=TRUE, scale.=TRUE)

scores_df <- tibble::as_tibble(pca_obj$x[,1:2]) %>%
  dplyr::rename(PC1=PC1,PC2=PC2) %>%
  dplyr::bind_cols(habitat_df %>% dplyr::select(patch_id,site,zone,year,state_final))

loadings_df <- tibble::as_tibble(pca_obj$rotation[,1:2], rownames="variable") %>%
  dplyr::mutate(vec_len = sqrt(PC1^2+PC2^2))

################################################################################
# 5. Panel A: PCA Biplot -------------------------------------------------------
################################################################################

range_x <- range(scores_df$PC1, na.rm=TRUE)
range_y <- range(scores_df$PC2, na.rm=TRUE)
mul <- min(diff(range_x)/(max(loadings_df$PC1)-min(loadings_df$PC1)),
           diff(range_y)/(max(loadings_df$PC2)-min(loadings_df$PC2))) * 0.8

arrow_df <- loadings_df %>%
  dplyr::filter(!grepl(ban_regex, variable, ignore.case=TRUE)) %>%
  dplyr::arrange(dplyr::desc(vec_len)) %>%
  dplyr::slice_head(n=8) %>%
  dplyr::mutate(PC1=PC1*mul, PC2=PC2*mul,
                variable_pretty = dplyr::recode(variable, !!!pretty_labs, .default = variable))

my_theme <- theme(
  axis.text.x   = element_text(size=8, color="black"),
  axis.text.y   = element_text(size=8, color="black"),
  axis.title    = element_text(size=10, color="black"),
  legend.text   = element_text(size=8, color="black"),
  legend.title  = element_text(size=8, color="black"),
  plot.tag      = element_text(size=10, color="black"),
  panel.grid    = element_blank(),
  panel.background = element_blank(),
  axis.line     = element_line(colour="black"),
  legend.key    = element_blank(),
  legend.background = element_rect(fill=alpha('blue', 0)),
  strip.text    = element_text(size=10, face="bold", color="black", hjust=0),
  strip.background = element_blank()
)

p_A <- ggplot2::ggplot(scores_df, ggplot2::aes(x=PC1, y=PC2, color=state_final)) +
  ggplot2::stat_ellipse(ggplot2::aes(fill=state_final), geom="polygon", alpha=0.15, color=NA) +
  ggplot2::stat_ellipse(linewidth=0.6) +
  ggplot2::geom_point(ggplot2::aes(shape=state_final), size=3, alpha=0.85) +
  ggplot2::geom_segment(data=arrow_df, ggplot2::aes(x=0,y=0,xend=PC1,yend=PC2),
                        arrow=ggplot2::arrow(length=grid::unit(0.2,"cm")), color="black", linewidth=0.7) +
  ggrepel::geom_label_repel(
    data = arrow_df,
    ggplot2::aes(x = PC1, y = PC2, label = variable_pretty),
    size = 2, fontface = "plain",
    segment.color = NA,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_color_manual(values=patch_colors, labels=c("Barren","Forest","Incipient")) +
  ggplot2::scale_fill_manual(values=patch_colors, guide="none") +
  ggplot2::scale_shape_manual(values=c("BAR"=16,"FOR"=15,"INCIP"=17),
                              labels=c("Barren","Forest","Incipient")) +
  ggplot2::coord_equal(expand=TRUE) +
  ggplot2::labs(x="PC1", y="PC2") +
  ggplot2::theme_bw(base_size=12) +
  ggplot2::theme(panel.grid=element_blank(), legend.position="none")

################################################################################
# 6. Panel B: RF Variable Importance ------------------------------------------
################################################################################
rf_imp <- rf_train2024$importance %>%
  as.data.frame() %>%
  tibble::rownames_to_column("variable")

if ("MeanDecreaseGini" %in% names(rf_imp)) {
  imp_col <- "MeanDecreaseGini"
} else if ("IncNodePurity" %in% names(rf_imp)) {
  imp_col <- "IncNodePurity"
} else stop("Variable importance column not found in RF object.")

rf_imp <- rf_imp %>%
  dplyr::arrange(dplyr::desc(.data[[imp_col]])) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::rename(Importance = !!imp_col) %>%
  dplyr::mutate(variable_pretty = dplyr::recode(variable, !!!pretty_labs, .default = variable))

p_B <- rf_imp %>%
  dplyr::mutate(variable_pretty = forcats::fct_reorder(variable_pretty, Importance)) %>%
  ggplot2::ggplot(ggplot2::aes(x = variable_pretty, y = Importance)) +
  ggplot2::geom_col(fill = "grey40") +
  ggplot2::coord_flip() +
  ggplot2::theme_bw() +
  ggplot2::labs(x = NULL, y = "Node impurity") +
  ggplot2::theme(panel.grid = element_blank())

################################################################################
# 7. Panel C: Raw Observed Values for Top 15 Vars ------------------------------
################################################################################

top15_tbl <- rf_imp %>%
  dplyr::arrange(dplyr::desc(Importance)) %>%
  dplyr::mutate(variable_pretty = forcats::fct_inorder(variable_pretty))
panelC_levels <- top15_tbl$variable_pretty

# matches "Barren", "Forest", "Incipient"
patch_colors_named <- c(
  "Barren"   = "#7570B3",
  "Forest"   = "#1B9E77",
  "Incipient"= "#D95F02"
)

box_df <- dat_raw %>%
  # optional: drop extreme outliers
  dplyr::filter(density20m2_lamset < 3) %>%
  dplyr::left_join(truth_all, by = c("patch_id","site","zone","year")) %>%
  dplyr::select(state, dplyr::all_of(top15_tbl$variable)) %>%
  tidyr::pivot_longer(cols=-state, names_to="variable", values_to="value") %>%
  dplyr::mutate(
    variable_pretty = dplyr::recode(variable, !!!pretty_labs, .default = variable),
    variable_pretty = factor(variable_pretty, levels = panelC_levels),
    state_label = factor(
      state,
      levels = c("BAR","FOR","INCIP"),
      labels = c("Barren","Forest","Incipient")
    )
  )

p_C <- ggplot2::ggplot(box_df, ggplot2::aes(x=state_label, y=value,
                                            fill=state_label, color=state_label)) +
  ggplot2::geom_boxplot(outlier.shape=NA, alpha=0.75, linewidth=0.4) +
  ggplot2::facet_wrap(~variable_pretty, scales="free_y", ncol=5) +
  ggplot2::scale_fill_manual(values=patch_colors_named, name="Patch state") +
  ggplot2::scale_color_manual(values=patch_colors_named, guide="none") +
  ggplot2::theme_bw() + my_theme +
  ggplot2::theme(
    legend.position="bottom",
    strip.text=ggplot2::element_text(size=8, face="bold"),
    axis.title.x=ggplot2::element_blank(),
    axis.text.x=ggplot2::element_text(size=7, angle=20, hjust=1),
    axis.text.y=ggplot2::element_text(size=7),
    panel.spacing=grid::unit(0.4,"lines"),
    panel.border=ggplot2::element_rect(color="black", linewidth=0.5)
  ) +
  ggplot2::labs(y="Observed value", x=NULL)

################################################################################
# 8. Combine Panels ------------------------------------------------------------
################################################################################

top_row <- p_A + p_B + patchwork::plot_spacer() +
  patchwork::plot_layout(widths=c(1.2,0.6,0.2))

final_fig_pca <- (top_row) / p_C +
  patchwork::plot_layout(heights=c(1,1.25)) +
  patchwork::plot_annotation(tag_levels="A") &
  ggplot2::theme(plot.tag=ggplot2::element_text(size=10),
                 plot.margin=ggplot2::margin(t=5,r=10,b=5,l=10))

final_fig_pca

ggplot2::ggsave(
  filename = here::here("figures", "Fig2_patch_drivers.png"),
  plot     = final_fig_pca,
  width    = 8.5,
  height   = 8,
  dpi      = 600,
  bg       = "white"
)
