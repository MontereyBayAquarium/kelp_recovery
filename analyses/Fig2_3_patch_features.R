################################################################################
# Predicting patch type and identifying features of incipient forests
# Author: Josh Smith
# Purpose: Build predictive models for patch type (BAR / FOR / INCIP)
#          and identify distinguishing features of incipient forests
################################################################################

rm(list = ls())

require(librarian)
librarian::shelf(
  tidyverse, MASS, lubridate, ranger, broom, pdp, rpart, rpart.plot, here,
  MLmetrics, glmnet, grid, vegan, patchwork
)

datadir <- "/Volumes/enhydra/data/kelp_recovery/MBA_kelp_forest_database/"
load(file.path(datadir, "processed/recovery/kelp_recovery_data.rda"))

################################################################################
# Step 1: Aggregate data to site-zone-site_type per year


kelp_avg <- kelp_data %>%
  dplyr::select(-macro_stipe_sd_20m2) %>%
  group_by(site, site_type, latitude, longitude, zone, survey_date) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-transect)

quad_avg <- quad_data %>%
  group_by(site, site_type, latitude, longitude, zone, survey_date) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-quadrat, -transect)

dat_agg <- inner_join(
  kelp_avg, quad_avg,
  by = c("site","site_type","latitude","longitude","zone","survey_date"),
  suffix = c("_kelp","_quad")
)

################################################################################
#Build predictor set

id_cols <- c("site","site_type","latitude","longitude","zone","survey_date",
             "transect","quadrat")

num_cols <- dat_agg %>% dplyr::select(where(is.numeric)) %>% names()
exclude_numeric <- c("latitude","longitude")
predictor_pool <- setdiff(num_cols, exclude_numeric)

base_2024 <- dat_agg %>%
  filter(year(survey_date) == 2024) %>%
  mutate(across(all_of(predictor_pool), ~ replace_na(., 0))) %>%
  dplyr::select(site, zone, site_type, survey_date, all_of(predictor_pool)) %>%
  mutate(site_type = factor(site_type))

base_2025 <- dat_agg %>%
  filter(year(survey_date) == 2025) %>%
  mutate(across(all_of(predictor_pool), ~ replace_na(., 0))) %>%
  dplyr::select(site, zone, site_type, survey_date, all_of(predictor_pool)) %>%
  mutate(site_type = factor(site_type, levels = levels(base_2024$site_type)))

cls <- levels(base_2024$site_type)

# Remove zero-variance predictors
nzv_cols <- caret::nearZeroVar(base_2024 %>% dplyr::select(all_of(predictor_pool)), saveMetrics = TRUE)
predictor_cols <- rownames(nzv_cols)[!nzv_cols$zeroVar & !nzv_cols$nzv]

################################################################################
#Fit LDA and Random Forest models 

#Train LDA on 2024 data
x_train <- base_2024 %>% dplyr::select(all_of(predictor_cols))
y_train <- droplevels(base_2024$site_type)

set.seed(1985)

# ---- LDA ----
lda_fit <- MASS::lda(y_train ~ ., data = data.frame(y_train, x_train))
#Warning is ok and expected
lda_pred_train <- predict(lda_fit)$class
lda_acc <- mean(lda_pred_train == y_train)
cat("\nLDA training accuracy:", round(lda_acc, 3), "\n")

# ---- Random Forest (ranger) ----
rf_fit <- ranger::ranger(
  formula = y_train ~ .,
  data = data.frame(y_train, x_train),
  num.trees = 500,
  mtry = floor(sqrt(ncol(x_train))),
  importance = "impurity",
  probability = TRUE,
  classification = TRUE
)
rf_acc <- 1 - rf_fit$prediction.error
cat("RF OOB accuracy:", round(rf_acc, 3), "\n")

# ---- Variable importance ----
rf_imp <- sort(rf_fit$variable.importance, decreasing = TRUE)
rf_imp_tbl <- tibble(feature = names(rf_imp), importance = rf_imp)
#write_csv(rf_imp_tbl, here::here("output","feature_importance_rf.csv"))

################################################################################
#Predict 2025 patch type probabilities

x_test <- base_2025 %>% dplyr::select(all_of(colnames(x_train)))

lda_post_2025 <- as_tibble(predict(lda_fit, newdata = x_test)$posterior)
rf_prob_2025  <- as_tibble(predict(rf_fit, data.frame(x_test))$predictions)

# Align factor levels
for (m in setdiff(cls, names(rf_prob_2025))) rf_prob_2025[[m]] <- 0
rf_prob_2025 <- rf_prob_2025[, cls]

# Ensemble average (LDA + RF)
ens_prob_2025 <- (lda_post_2025 + rf_prob_2025) / 2

################################################################################
#Build patch transition table 

prob_2025_lng <- bind_cols(
  base_2025 %>% dplyr::select(site, zone, site_type),
  ens_prob_2025
)

prob_2025_by_szst <- prob_2025_lng %>%
  group_by(site, zone, site_type) %>%
  summarise(across(all_of(cls), mean), .groups = "drop")

# LDA 2024 class predictions for constraint rule
lda_pred_2024 <- predict(lda_fit)$class
pred_2024_by_szst <- base_2024 %>%
  mutate(pred_2024 = lda_pred_2024) %>%
  group_by(site, zone, site_type) %>%
  summarise(patch_2024 = names(which.max(table(pred_2024))), .groups = "drop")

# Apply constraint: FOR cannot become INCIP
pred_2025_by_szst_constrained <- pred_2024_by_szst %>%
  inner_join(prob_2025_by_szst, by = c("site","zone","site_type")) %>%
  rowwise() %>%
  mutate(
    patch_2025 = {
      probs <- c_across(all_of(cls)); names(probs) <- cls
      if (patch_2024 == "FOR") {
        allowed <- c("FOR","BAR")
        allowed[which.max(probs[allowed])]
      } else {
        names(probs)[which.max(probs)]
      }
    }
  ) %>%
  ungroup() %>%
  mutate(patch_2025 = factor(patch_2025, levels = cls)) %>%
  dplyr::select(site, zone, site_type, patch_2024, patch_2025)

transitions_tbl_constrained <- pred_2025_by_szst_constrained
#save(transitions_tbl_constrained, file = here::here("output","lda_patch_transitionsv2.rda"))

################################################################################
#Feature interpretation – what defines incipient forests


incip_label <- if ("INCIP" %in% cls) "INCIP" else cls[grepl("INCIP", cls, ignore.case = TRUE)][1]

summ_inc <- base_2024 %>%
  mutate(is_incip = site_type == incip_label) %>%
  summarise(across(all_of(predictor_cols), list(
    med_incip = ~median(.x[is_incip], na.rm=TRUE),
    med_other = ~median(.x[!is_incip], na.rm=TRUE),
    s_md = ~{
      m1 <- mean(.x[is_incip], na.rm=TRUE)
      m0 <- mean(.x[!is_incip], na.rm=TRUE)
      s  <- sd(.x, na.rm=TRUE)
      (m1 - m0)/s
    }
  ), .names = "{.col}__{.fn}")) %>%
  pivot_longer(everything(),
               names_to = c("feature","stat"), names_sep="__",
               values_to = "value") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  arrange(desc(abs(s_md)))

#write_csv(summ_inc, here::here("output","incipient_signature_table.csv"))

# ---- Surrogate decision tree for simple rules ----
set.seed(42)
# Get predicted class labels from RF probability matrix
rf_pred_mat <- predict(rf_fit, data.frame(x_train))$predictions
rf_pred_labels <- colnames(rf_pred_mat)[max.col(rf_pred_mat, ties.method = "first")]

# Add predictions to a new data frame for the surrogate tree
surro_df <- base_2024 %>%
  mutate(pred = rf_pred_labels) %>%
  dplyr::select(all_of(predictor_cols), pred)


surro_tree <- rpart::rpart(as.factor(pred) ~ ., data = surro_df,
                           control = rpart::rpart.control(cp = 0.01, maxdepth = 4))
rules_txt <- rpart.plot::rpart.rules(surro_tree, style = "tallw")
capture.output(rules_txt, file = here::here("output","incipient_rules.txt"))


# ---- Partial dependence for top INCIP features ----
top_incip_feats <- head(summ_inc$feature, 6)
pdp_list <- lapply(top_incip_feats, function(v) {
  pdp::partial(rf_fit, pred.var = v,
               train = base_2024 %>% dplyr::select(all_of(predictor_cols)),
               which.class = incip_label, prob = TRUE)
})



################################################################################
# Build Figure 2

# Set labels
var_labels <- c(
  "purple_urchin_conceiledm2" = "Conceiled urchins",
  "cov_fleshy_red"            = "Fleshy red algae",
  "cov_crustose_coralline"    = "Crustose coralline algae",
  "purple_urchin_densitym2"   = "Purple urchin density",
  "n_macro_plants_20m2"       = "No. Macrocystis",
  "tegula_densitym2"          = "Tegula spp.",
  "cov_articulated_coralline" = "Articulated coralline algae",
  "cov_barnacle"              = "Barnacle cover",
  "cov_lam_holdfast_live"     = "Live kelp holdfast",
  "cov_phragmatopoma"         = "Phragmatopoma spp.",
  "cov_diopatra_chaetopterus" = "Diopatra chaetopterus",
  "cov_tubeworm_other_solitary" = "Tubeworm",
  "cov_bare_sand"             = "Bare sand",
  "nerj"                      = "Juvenile N. luetkeana",
  "lsetj"                     = "Juvenile L. setchellii",
  "lamr"                      = "Kelp recruits",
  "red_urchin_densitym2"      = "Red urchin density",
  "cov_encrusting_red"        = "Encrusting red algae",
  "cov_bare_rock"             = "Bare rock",
  "macro_stipe_density_20m2"  = "No. kelp stipes",
  "risk_index"                = "Rugosity",
  "relief_cm"                 = "Habitat relief (cm)",
  "density20m2_eisarb"        = "E. arborea",
  "density20m2_lamset"        = "L. setchellii",
  "density20m2_ptecal"        = "P. californica",
  "density20m2_nerlue"        = "N. luetkeana",
  "pomaulax_densitym2"        = "Pomaulax spp.",
  "macj"  = "Juvenile Macrocystis"
)

# set theme
my_theme <- theme(
  axis.text.x = element_text(size=8, color = "black"),
  axis.text.y = element_text(size=8, color = "black"),
  axis.title  = element_text(size=10,color = "black"),
  legend.text = element_text(size=8,color = "black"),
  legend.title= element_text(size=8,color = "black"),
  plot.tag = element_text(size = 10, color = "black"),
  # Gridlines
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_blank(), 
  axis.line = element_line(colour = "black"),
  # Legend
  legend.key = element_blank(),
  legend.background = element_rect(fill=alpha('blue', 0)),
  # Facets
  strip.text = element_text(size=10, face = "bold",color = "black", hjust=0),
  strip.background = element_blank()
)


# Get LDA predicted values
lda_pred <- predict(lda_fit)
lda_df <- data.frame(lda_pred$x, site_type = y_train)

# Build first panel
pA <- ggplot(lda_df, aes(LD1, LD2, color = site_type)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2) +
  scale_color_manual(values = c("BAR"="purple","INCIP"="orange","FOR"="forestgreen")) +
  theme_bw() + my_theme +
  theme(legend.position = "bottom") +
  labs(x = "LD1", y = "LD2", color = "Patch type")
# Removed coord_equal() to allow more height


# Overlay EnvFit vectors

#Get ord scores
ord_mat <- as.matrix(lda_pred$x[, 1:2])  # rows align with base_2024/x_train
drivers <- x_train %>% mutate(across(everything(), ~ as.numeric(scale(.x))))

#Fit envfit (permutation test for vector significance)
set.seed(123)
ef <- vegan::envfit(ord_mat, drivers, permutations = 999)

#Extract vector loadings, r2 and p-values
vec_raw <- as.data.frame(vegan::scores(ef, display = "vectors"))
vec_raw$variable <- rownames(vec_raw)
vec_raw$r2       <- ef$vectors$r      # squared correlation (explanatory power)
vec_raw$pval     <- ef$vectors$pvals  # permutation p-value

#mphasize variables that are both important and significant
#Here we keep the top 12 by r2 among those with p < 0.05
vec_keep <- vec_raw %>%
  filter(pval < 0.05) %>%
  arrange(desc(r2)) %>%
  slice_head(n = 12)

#Scale arrows to fit the LDA point cloud nicely
#Correct way to compute arrow scaling for envfit vectors
mul <- min(
  (max(ord_mat[,1]) - min(ord_mat[,1])) / (max(vec_keep$LD1) - min(vec_keep$LD1)),
  (max(ord_mat[,2]) - min(ord_mat[,2])) / (max(vec_keep$LD2) - min(vec_keep$LD2))
) * 0.8  # 0.8 = aesthetic shrink factor

vec_keep <- vec_keep %>%
  mutate(LD1 = LD1 * mul,
         LD2 = LD2 * mul)


vec_keep <- vec_keep %>%
  mutate(LD1 = LD1 * mul,
         LD2 = LD2 * mul)

#Apply labels when available
vec_keep <- vec_keep %>%
  mutate(label = ifelse(variable %in% names(var_labels),
                        unname(var_labels[variable]),
                        variable))

#Nudge labels a touch past the arrow tip (direction-aware)
lab_nudge <- 0.04 * max(dist(range(lda_df$LD1)), dist(range(lda_df$LD2)))

# --- Scale down the envfit vectors ---
shrink_factor <- 0.01  # try 0.3–0.5 if still too long

vec_keep <- vec_keep %>%
  mutate(
    LD1 = LD1 * mul * shrink_factor,
    LD2 = LD2 * mul * shrink_factor,
    xend = LD1,
    yend = LD2,
    xlab = LD1 + sign(LD1) * lab_nudge,
    ylab = LD2 + sign(LD2) * lab_nudge
  )



# Define color palette 
patch_colors <- c(
  "BAR"   = "#7570B3",
  "FOR"   = "#1B9E77",
  "INCIP" = "#D95F02"
)

# Rebuild p1 using your custom colors
p1 <- pA +
  scale_color_manual(
    values = patch_colors,
    name = "Patch type",
    labels = c("Barren", "Forest", "Incipient")  
  )+
  geom_segment(
    data = vec_keep,
    aes(x = 0, y = 0, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.18, "cm")),
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = vec_keep,
    aes(x = xend, y = yend, label = label),
    size = 3,
    fontface = "bold",
    min.segment.length = 0,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.25,
    segment.color = NA,
    inherit.aes = FALSE
  ) +
 # labs(
  #  caption = "Arrows: envfit vectors (length ~ explanatory power; labels shown for p < 0.05)"
  #) +
  theme(
    plot.caption = element_text(size = 8, color = "grey30"),
    legend.position = "none"
  )

# Feature importance barplot
p2 <- rf_imp_tbl %>%
  slice_max(importance, n = 15) %>%
  mutate(feature = ifelse(feature %in% names(var_labels), var_labels[feature], feature)) %>%
  ggplot(aes(x = reorder(feature, importance), y = importance)) +
  geom_col(fill = "grey40") +
  coord_flip() +
  theme_bw() + my_theme +
  labs(x = "", y = "Node impurity")

# Combine side-by-side with patchwork
p <- p1 + p2 +
  plot_layout(widths = c(1.2, 0.8)) +
  plot_annotation(tag_levels = "A")

p


#Build Panel C
#Identify and order the top 15 features (keep raw names for selection)
top15_tbl <- rf_imp_tbl %>%
  slice_max(importance, n = 15) %>%
  mutate(feature_label = ifelse(feature %in% names(var_labels),
                                var_labels[feature], feature)) %>%
  mutate(feature = forcats::fct_reorder(feature, importance, .desc = TRUE))

# 2. Create tidy data for boxplots, preserving raw feature names for select()
box_df <- base_2024 %>%
  dplyr::select(site_type, all_of(as.character(top15_tbl$feature))) %>%
  pivot_longer(
    cols = -site_type,
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    # Now apply readable labels after pivoting
    variable = factor(ifelse(variable %in% names(var_labels),
                             var_labels[variable], variable),
                      levels = top15_tbl$feature_label),
    site_type = factor(site_type,
                       levels = c("BAR", "FOR", "INCIP"),
                       labels = c("Barren", "Forest", "Incipient"))
  )

# 3. Define color palette
patch_colors_named <- c(
  "Barren"    = "#7570B3",
  "Forest"    = "#1B9E77",
  "Incipient" = "#D95F02"
)

# 4. Faceted boxplot ordered by node impurity
p_box <- ggplot(box_df, aes(x = site_type, y = value,
                            fill = site_type, color = site_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, linewidth = 0.4) +
  scale_fill_manual(values = patch_colors_named, name = "Patch type") +
  scale_color_manual(values = patch_colors_named, guide = "none") +
  facet_wrap(~ variable, scales = "free_y", ncol = 5) +
  theme_bw() + my_theme +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 8, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 7, angle = 20, hjust = 1),
    axis.text.y = element_text(size = 7),
    panel.spacing = unit(0.4, "lines"),
    panel.border = element_rect(color = "black", linewidth = 0.5)
  ) +
  labs(y = "Value")

# 5. Combine with panels A & B
final_fig <- (p / p_box) +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10))

final_fig



#ggsave(
#  here::here("figures", "Fig2_patch_features_with_boxplots.png"),
#  final_fig,
#  width = 9, height = 9, dpi = 600
#)



################################################################################
#Export figure

#ggsave(
#  here::here("figures","Fig2_patch_features.png"),
#  p,
#  width = 9, height = 5, dpi = 600
#)


################################################################################
#Build Figure 3


g1 <- summ_inc %>%
  slice_max(abs(s_md), n = 15) %>%
  mutate(feature = ifelse(feature %in% names(var_labels), var_labels[feature], feature)) %>%
  mutate(feature = fct_reorder(feature, s_md)) %>%
  ggplot(aes(x = feature, y = s_md, fill = s_md > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("TRUE"="#D95F02","FALSE"="gray40")) +
  theme_bw() + my_theme +
  labs(x = "Features that define \nincipient forests", y = "Node impurity") +
  guides(fill = "none")
g1



# ---- 2. Partial dependence plots (right column, 6 panels) ----
pdp_multi_list <- lapply(top_incip_feats, function(v) {
  pd_all <- lapply(cls, function(cl) {
    pdp::partial(
      rf_fit,
      pred.var = v,
      train = base_2024 %>% select(all_of(predictor_cols)),
      which.class = cl,
      prob = TRUE
    ) %>%
      mutate(site_type = cl)
  }) %>%
    bind_rows()
  
  x_label <- ifelse(v %in% names(var_labels), var_labels[[v]], v)
  
  ggplot(pd_all, aes_string(x = v, y = "yhat", color = "site_type")) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = patch_colors) +
    theme_bw() + my_theme +
    labs(x = x_label, color = "Patch type") +
    theme(
      legend.position = "none",
      axis.title.y = element_blank()
    )
})

# ---- Combine PDPs into one grid with shared legend ----
# ---- Combine PDPs into one grid with shared legend ----
pdp_grid <- patchwork::wrap_plots(
  lapply(top_incip_feats, function(v) {
    pd_all <- lapply(cls, function(cl) {
      pdp::partial(
        rf_fit,
        pred.var = v,
        train = base_2024 %>% select(all_of(predictor_cols)),
        which.class = cl,
        prob = TRUE
      ) %>% mutate(site_type = cl)
    }) %>% bind_rows()
    
    x_label <- ifelse(v %in% names(var_labels), var_labels[[v]], v)
    
    ggplot(pd_all, aes_string(x = v, y = "yhat", color = "site_type")) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(
        values = patch_colors,
        name = "Patch type",
        labels = c("Barren", "Forest", "Incipient")   # ← added custom labels
      ) +
      theme_bw() + my_theme +
      labs(x = x_label) +
      theme(axis.title.y = element_blank())
  }),
  ncol = 3, nrow = 2,
  guides = "collect"
) &
  theme(
    plot.tag = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8)
  )

# --- Make PDP grid core object (where tags will appear) ---
bottom_row_core <- pdp_grid

# --- Add shared y-axis label as an inset element (same method you used before) ---
bottom_row <- patchwork::inset_element(
  bottom_row_core,
  left = 0.06, right = 1, bottom = 0, top = 1,
  align_to = "full"
) + 
  patchwork::wrap_elements(
    full = grid::textGrob(
      "Predicted probability", rot = 90,
      gp = grid::gpar(fontsize = 10, col = "black", fontface = "plain")
    )
  )

# ---- Combine left column (g1) with right column (PDPs) ----
combined_final <- (g1 | bottom_row_core) +
  plot_layout(widths = c(0.7, 1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10))

combined_final


################################################################################
#Export
ggsave(
  here::here("figures","Fig3_incipient_features.png"),
  combined_final,
  width = 9, height = 5, dpi = 600
)
