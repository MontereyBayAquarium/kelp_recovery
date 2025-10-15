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
  MLmetrics, glmnet
)

datadir <- "/Volumes/enhydra/data/kelp_recovery/MBA_kelp_forest_database/"
load(file.path(datadir, "processed/recovery/kelp_recovery_data.rda"))

################################################################################
# Step 1: Aggregate data to site-zone-site_type per year


kelp_avg <- kelp_data %>%
  select(-macro_stipe_sd_20m2) %>%
  group_by(site, site_type, latitude, longitude, zone, survey_date) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  select(-transect)

quad_avg <- quad_data %>%
  group_by(site, site_type, latitude, longitude, zone, survey_date) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  select(-quadrat, -transect)

dat_agg <- inner_join(
  kelp_avg, quad_avg,
  by = c("site","site_type","latitude","longitude","zone","survey_date"),
  suffix = c("_kelp","_quad")
)

################################################################################
#Build predictor set

id_cols <- c("site","site_type","latitude","longitude","zone","survey_date",
             "transect","quadrat")

num_cols <- dat_agg %>% select(where(is.numeric)) %>% names()
exclude_numeric <- c("latitude","longitude")
predictor_pool <- setdiff(num_cols, exclude_numeric)

base_2024 <- dat_agg %>%
  filter(year(survey_date) == 2024) %>%
  mutate(across(all_of(predictor_pool), ~ replace_na(., 0))) %>%
  select(site, zone, site_type, survey_date, all_of(predictor_pool)) %>%
  mutate(site_type = factor(site_type))

base_2025 <- dat_agg %>%
  filter(year(survey_date) == 2025) %>%
  mutate(across(all_of(predictor_pool), ~ replace_na(., 0))) %>%
  select(site, zone, site_type, survey_date, all_of(predictor_pool)) %>%
  mutate(site_type = factor(site_type, levels = levels(base_2024$site_type)))

cls <- levels(base_2024$site_type)

# Remove zero-variance predictors
nzv_cols <- caret::nearZeroVar(base_2024 %>% select(all_of(predictor_pool)), saveMetrics = TRUE)
predictor_cols <- rownames(nzv_cols)[!nzv_cols$zeroVar & !nzv_cols$nzv]

################################################################################
#Fit LDA and Random Forest models 

#Train LDA on 2024 data
x_train <- base_2024 %>% select(all_of(predictor_cols))
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

x_test <- base_2025 %>% select(all_of(colnames(x_train)))

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
  base_2025 %>% select(site, zone, site_type),
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
  select(site, zone, site_type, patch_2024, patch_2025)

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
  select(all_of(predictor_cols), pred)


surro_tree <- rpart::rpart(as.factor(pred) ~ ., data = surro_df,
                           control = rpart::rpart.control(cp = 0.01, maxdepth = 4))
rules_txt <- rpart.plot::rpart.rules(surro_tree, style = "tallw")
capture.output(rules_txt, file = here::here("output","incipient_rules.txt"))


# ---- Partial dependence for top INCIP features ----
top_incip_feats <- head(summ_inc$feature, 9)
pdp_list <- lapply(top_incip_feats, function(v) {
  pdp::partial(rf_fit, pred.var = v,
               train = base_2024 %>% select(all_of(predictor_cols)),
               which.class = incip_label, prob = TRUE)
})



################################################################################
# visualize



rf_imp_tbl %>%
  slice_max(importance, n = 15) %>%
  ggplot(aes(x = reorder(feature, importance), y = importance)) +
  geom_col(fill = "darkolivegreen4") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top 15 predictors distinguishing patch types",
    x = "Feature",
    y = "Variable importance (RF impurity)"
  )


summ_inc %>%
  slice_max(abs(s_md), n = 20) %>%
  ggplot(aes(x = reorder(feature, s_md), y = s_md)) +
  geom_col(aes(fill = s_md > 0)) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "forestgreen", "FALSE" = "firebrick")) +
  theme_minimal() +
  labs(
    title = "Features elevated (+) or reduced (–) in Incipient forests",
    x = "Feature",
    y = "Standardized mean difference (INCIP vs others)"
  )


library(patchwork)
pdp_plots <- lapply(seq_along(pdp_list), function(i) {
  autoplot(pdp_list[[i]]) +
    labs(
      title = paste0("INCIP probability vs ", top_incip_feats[i]),
      y = "Predicted INCIP probability",
      x = top_incip_feats[i]
    ) +
    theme_minimal(base_size = 11)
})
patchwork::wrap_plots(pdp_plots)



lda_pred <- predict(lda_fit)
lda_df <- data.frame(lda_pred$x, site_type = y_train)

ggplot(lda_df, aes(LD1, LD2, color = site_type)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Linear Discriminant space of patch types",
    x = "LD1 (main separation axis)",
    y = "LD2"
  )





