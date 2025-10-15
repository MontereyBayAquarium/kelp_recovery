
#jogsmith@ucsc.edu

rm(list = ls())


################################################################################
#Step 0: set paths and load data
require(librarian)
librarian::shelf(tidyverse, MASS, lubridate, caret, ranger)

datadir <- "/Volumes/enhydra/data/kelp_recovery/MBA_kelp_forest_database/"
load(file.path(datadir, "processed/recovery/kelp_recovery_data.rda"))

################################################################################
#Step 1: Average to site, zone, site_type for each year 
kelp_avg <- kelp_data %>%
  dplyr::select(-macro_stipe_sd_20m2) %>%
  dplyr::group_by(site, site_type, latitude, longitude, zone, survey_date) %>%
  dplyr::summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  dplyr::select(-transect)

quad_avg <- quad_data %>%
  dplyr::group_by(site, site_type, latitude, longitude, zone, survey_date) %>%
  dplyr::summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  dplyr::select(-quadrat, -transect)

dat_agg <- kelp_avg %>%
  dplyr::inner_join(
    quad_avg,
    by = c("site","site_type","latitude","longitude","zone","survey_date"),
    suffix = c("_kelp","_quad")
  )

###############################################################################
#Step 2: seledct predictors and prepare training data

#select a stipitate-only model for predicting incipient forests
predictor_cols <- c(
  # "relief_cm", "risk_index",
  # "purple_urchin_densitym2", "red_urchin_densitym2",
  "macro_stipe_density_20m2",
  "density20m2_eisarb",
  "density20m2_ptecal", 
  "density20m2_lamset",
  "density20m2_nerlue"
  #"cov_desmarestia_spp", "cov_stephanocystis"
   #"cov_fleshy_red","cov_encrusting_red","cov_crustose_coralline","cov_bare_rock"
)

base_2024 <- dat_agg %>%
  dplyr::filter(lubridate::year(survey_date) == 2024) %>%
  dplyr::select(site, zone, site_type, survey_date, dplyr::all_of(predictor_cols)) %>%
  dplyr::mutate(across(where(is.numeric), ~tidyr::replace_na(.x, 0))) %>%
  dplyr::mutate(site_type = factor(site_type))

base_2025 <- dat_agg %>%
  dplyr::filter(lubridate::year(survey_date) == 2025) %>%
  dplyr::select(site, zone, site_type, survey_date, dplyr::all_of(predictor_cols)) %>%
  dplyr::mutate(across(where(is.numeric), ~tidyr::replace_na(.x, 0))) %>%
  dplyr::mutate(site_type = factor(site_type, levels = levels(base_2024$site_type)))

cls <- levels(base_2024$site_type) 


################################################################################
#Step 3: train LDA on 2024, then apply to 2025

set.seed(1985)

X_train <- scale(dplyr::select(base_2024, dplyr::all_of(predictor_cols)))
train_df <- data.frame(site_type = base_2024$site_type, X_train)

mu <- attr(X_train, "scaled:center")
sd <- attr(X_train, "scaled:scale")

X_test  <- scale(dplyr::select(base_2025, dplyr::all_of(predictor_cols)), center = mu, scale = sd)
test_df <- data.frame(site_type = base_2025$site_type, X_test)

lda_fit <- MASS::lda(site_type ~ ., data = train_df)

# Honest 2024 preds cross validation on train data
lda_cv_2024   <- MASS::lda(site_type ~ ., data = train_df, CV = TRUE)
lda_pred_2024 <- lda_cv_2024$class

# Get 2025 posterior probs for ensemble
lda_post_2025 <- as_tibble(predict(lda_fit, newdata = test_df)$posterior) %>%
  dplyr::select(dplyr::all_of(cls))


###############################################################################
#Step 4: train a RF model on 2024, then apply to 2025


ctrl <- caret::trainControl(
  method = "repeatedcv", number = 10, repeats = 3,
  sampling = "up"  
)

p <- length(predictor_cols)
grid <- expand.grid(
  mtry = unique(pmax(1, round(c(sqrt(p), p/3, p/2)))),
  splitrule = "gini",
  min.node.size = c(1, 3, 5)
)

rf_fit <- caret::train(
  site_type ~ .,
  data      = base_2024 %>% dplyr::select(site_type, dplyr::all_of(predictor_cols)),
  method    = "ranger",
  trControl = ctrl,
  tuneGrid  = grid,
  num.trees = 1000,
  importance= "impurity"
)

# RF class probabilities for 2025 via underlying ranger model
align_probs <- function(probs_df, cls_levels) {
  probs_df <- as.data.frame(probs_df)
  miss <- setdiff(cls_levels, colnames(probs_df))
  for (m in miss) probs_df[[m]] <- 0
  probs_df <- probs_df[, cls_levels, drop = FALSE]
  tibble::as_tibble(probs_df)
}

rf_prob_2025_raw <- predict(
  rf_fit$finalModel,
  data        = base_2025 %>% dplyr::select(dplyr::all_of(predictor_cols)),
  type        = "response",
  probability = TRUE
)$predictions %>% as.data.frame()

rf_prob_2025 <- align_probs(rf_prob_2025_raw, cls)

################################################################################
#Step 5: build a probability ensemble for identifying best predictive model


ens_prob_2025 <- (lda_post_2025 + rf_prob_2025) / 2

# attach probs to 2025 rows
prob_2025_lng <- dplyr::bind_cols(
  base_2025 %>% dplyr::select(site, zone, site_type),
  ens_prob_2025
)

# aggregate to site / zone / site_type level
prob_2025_by_szst <- prob_2025_lng %>%
  dplyr::group_by(site, zone, site_type) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(cls), mean), .groups = "drop")

# Collapse 2024 to one label per site / zone / site_type
mode_factor <- function(x, levels_keep) {
  xt <- table(x)
  # tie-break using factor level order
  lv <- levels_keep[which.max(sapply(levels_keep, function(L) ifelse(L %in% names(xt), xt[[L]], 0)))]
  factor(lv, levels = levels_keep)
}

pred_2024_by_szst <- base_2024 %>%
  dplyr::mutate(pred_2024 = lda_pred_2024) %>%
  dplyr::group_by(site, zone, site_type) %>%
  dplyr::summarise(patch_2024 = mode_factor(pred_2024, cls), .groups = "drop")

#Apply constrain -- FOR cannot become incipient
pred_2025_by_szst_constrained <- pred_2024_by_szst %>%
  dplyr::inner_join(prob_2025_by_szst, by = c("site","zone","site_type")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    patch_2025 = {
      probs <- c_across(dplyr::all_of(cls)); names(probs) <- cls
      if (patch_2024 == "FOR") {
        allowed <- c("FOR","BAR")
        allowed[which.max(probs[allowed])]
      } else {
        names(probs)[which.max(probs)]
      }
    }
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(patch_2025 = factor(patch_2025, levels = cls)) %>%
  dplyr::select(site, zone, site_type, patch_2024, patch_2025)

# extract final patch transition table
transitions_tbl_constrained <- pred_2025_by_szst_constrained
transitions_tbl_constrained

save(transitions_tbl_constrained, file = here::here("output","lda_patch_transitions.rda"))

################################################################################
# Step 6: Explore accuracy summaries

# LDA labels on 2025 (argmax posterior)
lda_pred_labels_2025 <- factor(
  colnames(lda_post_2025)[max.col(lda_post_2025, ties.method = "first")],
  levels = cls
)
# RF labels on 2025
rf_pred_labels_2025 <- predict(
  rf_fit, newdata = base_2025 %>% dplyr::select(dplyr::all_of(predictor_cols))
)

# Ensemble labels on 2025 (argmax)
ens_pred_labels_2025 <- factor(
  colnames(ens_prob_2025)[max.col(ens_prob_2025, ties.method = "first")],
  levels = cls
)

lda_cm <- caret::confusionMatrix(lda_pred_labels_2025, base_2025$site_type)
rf_cm  <- caret::confusionMatrix(rf_pred_labels_2025,  base_2025$site_type)
ens_cm <- caret::confusionMatrix(ens_pred_labels_2025, base_2025$site_type)

cat("\n--- LDA ---\n"); print(lda_cm)
cat("\n--- RF (ranger) ---\n"); print(rf_cm)
cat("\n--- Ensemble ---\n"); print(ens_cm)

lda_cm$overall["Accuracy"]; rf_cm$overall["Accuracy"]; ens_cm$overall["Accuracy"]



