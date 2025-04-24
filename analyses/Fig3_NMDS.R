

rm(list=ls())

require(librarian)

librarian::shelf(tidyverse, ggplot2, RColorBrewer, vegan, grid, tidytext, forcats)

datdir <- "/Volumes/seaotterdb$/kelp_recovery/data/MBA_kelp_forest_database"
figdir <- here::here("figures")

quad_dat <- read_csv(file.path(datdir,"processed/recovery/recovery_quad.csv"))

urch_size_dat <- read_csv(file.path(datdir,"processed/recovery/recovery_urch_sizefq.csv")) 

kelp_dat <- read_csv(file.path(datdir,"processed/recovery/recovery_kelpswath.csv")) 

gonad_dat <- read_csv(file.path(datdir, "processed/dissection_data_recovery.csv")) 

################################################################################
#Prep derived datasets

#Step 1: aggregate quad data to transect level
quad_transect <- quad_dat %>%
                  #drop substrate for now, we'll handle this later
                  dplyr::select(-substrate, -quadrat)%>%
                  group_by(site, site_type, survey_date, latitude, longitude,
                           zone, transect) %>%
                  dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

#check nrow
nrow(quad_transect)
nrow(kelp_dat)
#wooo they match!

#Step 2: aggregate dissection data to transect level
gonad_zone <- gonad_dat %>%
                  #drop substrate for now, we'll handle this later
                  dplyr::select(site, site_type, zone, species, test_diameter_mm,
                         gonad_index)%>%
                  group_by(site, site_type, zone, species) %>%
                  dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
                pivot_wider(
                  names_from = species,
                  values_from = c(test_diameter_mm, gonad_index)
                )

#step 3: join data 
merge_dat <- quad_transect %>%
              left_join(., gonad_zone, by = c("site", "site_type", "zone")) %>%
              left_join(., kelp_dat, by = c("survey_date","site","site_type", 
                                            "zone","transect","latitude", "longitude"
                                            )) %>%
  #reorder
  dplyr::select(site, site_type, survey_date, latitude, longitude, zone, transect,
         depth_m, relief_cm, risk_cm, risk_index, purple_urchin_densitym2, 
         purple_urchin_densitym2, red_urchin_densitym2, red_urchin_conceiledm2, 
         tegula_densitym2, pomaulax_densitym2, test_diameter_mm_purple_urchin,
         test_diameter_mm_red_urchin, gonad_index_purple_urchin, 
         gonad_index_red_urchin, n_macro_plants_20m2, macro_stipe_density_20m2, 
         macro_stipe_sd_20m2, density20m2_purps_on_kelp, density20m2_ptecal, 
         density20m2_eisarb, density20m2_nerlue, density20m2_lamset,  
         density20m2_cancer_spp, density20m2_lamstump,density20m2_macstump,
         everything())


################################################################################
#Analyses

# Function to clean up variable names:
nice_name <- function(x) {
  # Remove the "cov_" prefix, if present
  x <- gsub("^cov_", "", x)
  # Remove the "density20m2_" prefix, if present (for density vectors)
  x <- gsub("^density20m2_", "", x)
  x <- gsub("m2", "", x)
  # Replace underscores with spaces
  x <- gsub("_", " ", x)
  # Convert to lowercase and then capitalize the first letter
  x <- tolower(x)
  x <- paste0(toupper(substr(x,1,1)), substr(x,2, nchar(x)))
  return(x)
}


#Helper Function for Pairwise PERMANOVA-----------------------------------------

pairwise_adonis2 <- function(dist_matrix, groups, perm = 999, p.adjust.method = "bonferroni") {
  groups <- as.factor(groups)
  pair_list <- combn(levels(groups), 2, simplify = FALSE)
  results <- data.frame()
  dmat <- as.matrix(dist_matrix)  # Convert dist object to full matrix
  for (pair in pair_list) {
    idx <- groups %in% pair
    dpair <- as.dist(dmat[idx, idx])
    ad_res <- adonis2(dpair ~ groups[idx], permutations = perm)
    res <- data.frame(Group1 = pair[1],
                      Group2 = pair[2],
                      F.Model = ad_res$F[1],
                      R2 = ad_res$R2[1],
                      p.value = ad_res$`Pr(>F)`[1])
    results <- rbind(results, res)
  }
  results$p.adjusted <- p.adjust(results$p.value, method = p.adjust.method)
  return(results)
}

#Subset and Clean Cover Data
cover_data <- merge_dat[, grep("^cov_", colnames(merge_dat))]
complete_idx_cover <- complete.cases(cover_data)
cover_data_complete <- cover_data[complete_idx_cover, ]
meta_cover <- merge_dat[complete_idx_cover, ]
meta_cover$site_type <- as.factor(meta_cover$site_type)
non_constant_cover <- sapply(cover_data_complete, function(x) sd(x, na.rm = TRUE)) > 0
cover_data_complete <- cover_data_complete[, non_constant_cover]

#Run NMDS on raw cover data-----------------------------------------------------

#NMDS on Raw Cover Data Using Bray–Curtis Distance
set.seed(123)
nmds_cover <- metaMDS(cover_data_complete, distance = "bray", trymax = 100)
cat("NMDS stress (Cover):", nmds_cover$stress, "\n")
scores_cover <- as.data.frame(scores(nmds_cover, display = "sites"))
scores_cover$site_type <- meta_cover$site_type

##Plot NMDS for Cover Data: Option to show only centroids
centroids_cover <- scores_cover %>%
  group_by(site_type) %>%
  summarize(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))

#Global and Pairwise PERMANOVA on Cover Data
dist_cover <- vegdist(cover_data_complete, method = "bray")
perm_cover <- adonis2(dist_cover ~ meta_cover$site_type, permutations = 999)
print(perm_cover)
pairwise_cover <- pairwise_adonis2(dist_cover, meta_cover$site_type, perm = 999)
print(pairwise_cover)

#SIMPER Analysis for Cover Data
simper_cover <- simper(cover_data_complete, group = meta_cover$site_type, permutations = 999)
simper_summary_cover <- summary(simper_cover)
print(simper_summary_cover)

#Overlay Vectors on Cover NMDS Plot
# Physical Vectors (Blue)
physical_vars <- meta_cover[, c("relief_cm", "risk_index")]
envfit_phys <- envfit(nmds_cover, physical_vars, permutations = 999)
phys_scores <- as.data.frame(scores(envfit_phys, display = "vectors"))
phys_scores$variable <- rownames(phys_scores)
phys_sig <- phys_scores[envfit_phys$vectors$pvals < 0.05, ]
cat("Significant physical vectors (blue):\n")
print(phys_sig)

#Biological Vectors (Black) from cover data
envfit_bio <- envfit(nmds_cover, cover_data_complete, permutations = 999)
bio_scores <- as.data.frame(scores(envfit_bio, display = "vectors"))
bio_scores$variable <- rownames(bio_scores)
bio_sig <- bio_scores[envfit_bio$vectors$pvals < 0.05, ]
cat("Significant biological vectors (black):\n")
print(bio_sig)

# Update vector names using nice_name() for better display
phys_sig$nice_variable <- nice_name(phys_sig$variable)
bio_sig$nice_variable  <- nice_name(bio_sig$variable)


#Run NMDS on density data-------------------------------------------------------

density_vars <- c(
  "red_urchin_densitym2",
  "red_urchin_conceiledm2",
  "tegula_densitym2",
  "pomaulax_densitym2",
  "test_diameter_mm_purple_urchin",
  "test_diameter_mm_red_urchin",
  "gonad_index_purple_urchin",
  "gonad_index_red_urchin",
  "n_macro_plants_20m2",
  "macro_stipe_density_20m2",
  #"density20m2_purps_on_kelp",
  "density20m2_ptecal",
  "density20m2_eisarb",
  "density20m2_nerlue",
  "density20m2_lamset",
  #"density20m2_cancer_spp",
  "density20m2_lamstump",
  "density20m2_macstump",
  "purple_urchin_densitym2",
  "purple_urchin_conceiledm2"
)

density_data <- merge_dat[, density_vars]
complete_idx_density <- complete.cases(density_data)
density_data_complete <- density_data[complete_idx_density, ]
meta_density <- merge_dat[complete_idx_density, ]
meta_density$site_type <- as.factor(meta_density$site_type)
non_constant_density <- sapply(density_data_complete, function(x) sd(x, na.rm = TRUE)) > 0
density_data_complete <- density_data_complete[, non_constant_density]

## NMDS on Raw Density Data Using Euclidean Distance
set.seed(123)
nmds_density <- metaMDS(density_data_complete, distance = "euclidean", trymax = 100)
cat("NMDS stress (Density):", nmds_density$stress, "\n")
scores_density <- as.data.frame(scores(nmds_density, display = "sites"))
scores_density$site_type <- meta_density$site_type


## Global and Pairwise PERMANOVA on Density Data
dist_density <- dist(density_data_complete, method = "euclidean")
perm_density <- adonis2(dist_density ~ meta_density$site_type, permutations = 999)
print(perm_density)
pairwise_density <- pairwise_adonis2(dist_density, meta_density$site_type, perm = 999)
print(pairwise_density)

## (Optional) SIMPER Analysis for Density Data
simper_density <- simper(density_data_complete, group = meta_density$site_type, permutations = 999)
simper_summary_density <- summary(simper_density)
print(simper_summary_density)

## Overlay Vectors on Density NMDS Plot
# Physical Vectors (Blue)
physical_vars_density <- meta_density[, c("relief_cm", "risk_index")]
envfit_phys_density <- envfit(nmds_density, physical_vars_density, permutations = 999)
phys_scores_density <- as.data.frame(scores(envfit_phys_density, display = "vectors"))
phys_scores_density$variable <- rownames(phys_scores_density)
phys_sig_density <- phys_scores_density[envfit_phys_density$vectors$pvals < 0.05, ]
cat("Significant physical vectors for density (blue):\n")
print(phys_sig_density)

# Biological Vectors (Black) using density data
envfit_bio_density <- envfit(nmds_density, density_data_complete, permutations = 999)
bio_scores_density <- as.data.frame(scores(envfit_bio_density, display = "vectors"))
bio_scores_density$variable <- rownames(bio_scores_density)
bio_sig_density <- bio_scores_density[envfit_bio_density$vectors$pvals < 0.05, ]
cat("Significant biological vectors for density (black):\n")
print(bio_sig_density)

# Use nice_name() to clean up names; this will now remove "density20m2_" if present
phys_sig_density$nice_variable <- nice_name(phys_sig_density$variable)
bio_sig_density$nice_variable  <- nice_name(bio_sig_density$variable)



################################################################################
#Plot

my_theme <-  theme(axis.text.x=element_text(size=18, color = "black"),
                   axis.text.y=element_text(size=18, color = "black"),
                   axis.title=element_text(size=20,color = "black"),
                   legend.text=element_text(size=18,color = "black"),
                   legend.title=element_text(size=18,color = "black"),
                   plot.tag=element_text(size=18,color = "black"),
                   # Gridlines
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   legend.key = element_rect(fill=alpha('blue', 0)),
                   legend.background = element_rect(fill=alpha('blue', 0)),
                   #facets
                   strip.text = element_text(size=20, face = "bold",color = "black", hjust=0),
                   strip.background = element_blank())

p1 <- ggplot() +
  # individual points
  geom_point(
    data = scores_cover %>% filter(site_type %in% c("FOR", "BAR")),
    aes(x = NMDS1, y = NMDS2, color = site_type, shape = site_type),
    size = 2, alpha = 0.2
  ) +
  # ellipses
  stat_ellipse(
    data = scores_cover %>% filter(site_type %in% c("FOR", "BAR")),
    aes(x = NMDS1, y = NMDS2, color = site_type),
    type = "norm", linetype = 1, size = 1
  ) +
  # centroids
  geom_point(
    data = centroids_cover %>% filter(site_type %in% c("FOR", "BAR")),
    aes(x = NMDS1, y = NMDS2, color = site_type, shape = site_type),
    size = 4, alpha = 1
  ) +
  scale_color_manual(values = c(
    "FOR"   = "#1B9E77",
    "INCIP" = "#D95F02",
    "BAR"   = "#7570B3"
  )) +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    color = "Site Type",
    shape = "Site Type"
  ) +
  theme_bw() +
  my_theme

p1

#ggsave(p1,  filename=file.path(figdir, "Fig6_NMDS_BAR_FOR.png"), width = 10, height = 7.5, units = "in",
#       bg = "white", dpi = 600) 



p2 <- ggplot(scores_cover, aes(x = NMDS1, y = NMDS2)) +
  # individual points
  geom_point(
    data = scores_cover %>% filter(site_type %in% c("FOR", "BAR", "INCIP")),
    aes(x = NMDS1, y = NMDS2, color = site_type, shape = site_type),
    size = 2, alpha = 0.2
  ) +
  stat_ellipse(aes(color = site_type), type = "norm", linetype = 1, size = 1) +
  # geom_point(aes(color = site_type, shape = site_type), size = 3, alpha = 0.5) +
  geom_point(data = centroids_cover, aes(x = NMDS1, y = NMDS2, color = site_type, shape = site_type),
             size = 4, alpha = 1) +
  scale_color_manual(values = c(
    "FOR" = "#1B9E77",
    "INCIP"  = "#D95F02",
    "BAR" = "#7570B3"
  )) +
  labs(title = "",
       x = "NMDS1", y = "NMDS2", color = "Site Type", shape = "Site Type") +
  theme_bw() + my_theme

p2


#ggsave(p2,  filename=file.path(figdir, "Fig7_NMDS_BAR_FOR_INCIP.png"), width = 10, height = 7.5, units = "in",
#       bg = "white", dpi = 600) 





# Define multipliers for envfit vectors:
cover_multiplier <- 1.2   # For cover data vectors
density_multiplier <- 0.26 # For density data vectors

# Increase vector lengths using cover_multiplier
phys_sig$NMDS1 <- phys_sig$NMDS1 * cover_multiplier
phys_sig$NMDS2 <- phys_sig$NMDS2 * cover_multiplier
bio_sig$NMDS1  <- bio_sig$NMDS1 * cover_multiplier
bio_sig$NMDS2  <- bio_sig$NMDS2 * cover_multiplier

p3 <- p2 +
  geom_segment(data = phys_sig, inherit.aes = FALSE,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), color = "blue", size = 1) +
  geom_text(data = phys_sig, inherit.aes = FALSE,
            aes(x = NMDS1, y = NMDS2, label = nice_variable),
            color = "blue", vjust = 1.5, size = 3) +
  geom_segment(data = bio_sig, inherit.aes = FALSE,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), color = "gray70", size = 1) +
  geom_text(data = bio_sig, inherit.aes = FALSE,
            aes(x = NMDS1, y = NMDS2, label = nice_variable),
            color = "black", vjust = -0.5, size = 3)
p3



#Plot density-------------------------------------------------------------------


## For Density, show only centroids
centroids_density <- scores_density %>%
  group_by(site_type) %>%
  summarize(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))

p4 <- ggplot(scores_density %>% filter(site_type %in% c("FOR", "BAR")), aes(x = NMDS1, y = NMDS2)) +
  stat_ellipse(aes(color = site_type), type = "norm", linetype = 1, size = 1) +
  # Uncomment the next line to show individual points:
  # geom_point(aes(color = site_type, shape = site_type), size = 3, alpha = 0.5) +
  geom_point(data = centroids_density %>% filter(site_type %in% c("FOR", "BAR")), aes(x = NMDS1, y = NMDS2, color = site_type, shape = site_type),
             size = 4, alpha = 1) +
  scale_color_manual(values = c(
    "FOR" = "#1B9E77",
    "INCIP"  = "#D95F02",
    "BAR" = "#7570B3"
  )) +
  labs(title = "",
       x = "NMDS1", y = "NMDS2", color = "Site Type", shape = "Site Type") +
  theme_bw() + my_theme
p4


#ggsave(p4,  filename=file.path(figdir, "Fig8_NMDS_density_FOR_BAR.png"), width = 10, height = 7.5, units = "in",
#       bg = "white", dpi = 600) 



p5 <- ggplot(scores_density, aes(x = NMDS1, y = NMDS2)) +
  stat_ellipse(aes(color = site_type), type = "norm", linetype = 1, size = 1) +
  # Uncomment the next line to show individual points:
  # geom_point(aes(color = site_type, shape = site_type), size = 3, alpha = 0.5) +
  geom_point(data = centroids_density, aes(x = NMDS1, y = NMDS2, color = site_type, shape = site_type),
             size = 4, alpha = 1) +
  scale_color_manual(values = c(
    "FOR" = "#1B9E77",
    "INCIP"  = "#D95F02",
    "BAR" = "#7570B3"
  )) +
  labs(title = "",
       x = "NMDS1", y = "NMDS2", color = "Site Type", shape = "Site Type") +
  theme_bw() + my_theme
p5


#ggsave(p5,  filename=file.path(figdir, "Fig9_density_BAR_FOR_INCIP.png"), width = 10, height = 7.5, units = "in",
#       bg = "white", dpi = 600) 



# Multiply density vectors by density_multiplier
phys_sig_density$NMDS1 <- phys_sig_density$NMDS1 * density_multiplier
phys_sig_density$NMDS2 <- phys_sig_density$NMDS2 * density_multiplier
bio_sig_density$NMDS1 <- bio_sig_density$NMDS1 * density_multiplier
bio_sig_density$NMDS2 <- bio_sig_density$NMDS2 * density_multiplier

p6 <- p5 +
  geom_segment(data = phys_sig_density, inherit.aes = FALSE,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), color = "indianred", size = 1) +
  geom_text(data = phys_sig_density, inherit.aes = FALSE,
            aes(x = NMDS1, y = NMDS2, label = nice_variable),
            color = "indianred", vjust = 1.5, size = 3) +
  geom_segment(data = bio_sig_density, inherit.aes = FALSE,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), color = "gray70", size = 1,
               alpha = 0.5) +
  #geom_text(data = bio_sig_density, inherit.aes = FALSE,
   #         aes(x = NMDS1, y = NMDS2, label = nice_variable),
    #        color = "black", vjust = -0.5, 
     #       nudge_x = 0.001, nudge_y = 0.001, size = 5)+
  #coord_cartesian(clip = "off") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
p6

#ggsave(p6,  filename=file.path(figdir, "Fig10_density_vectors_unlabeled.png"), width = 10, height = 7.5, units = "in",
 #      bg = "white", dpi = 600) 


################################################################################
#Testing predictor models


library(randomForest)


# -----------------------------------------------------------------------------
# 1. Prepare df_rf with ONLY your six predictors + response
# -----------------------------------------------------------------------------
df_rf <- meta_cover %>%
  data.frame()%>%
  filter(site_type %in% c("BAR", "INCIP", "FOR")) %>%
  select(
    site_type,
    relief_cm,
    risk_index,
    purple_urchin_densitym2,
    purple_urchin_conceiledm2,
    red_urchin_densitym2,
    red_urchin_conceiledm2
  ) %>%
  na.omit() %>%
  mutate(site_type = factor(site_type, levels = c("BAR","INCIP","FOR"))) 


# Double‐check that df_rf has exactly 7 columns:
print(colnames(df_rf))
#> [1] "site_type"                   "relief_cm"                  
#> [3] "risk_index"                  "purple_urchin_densitym2"    
#> [5] "purple_urchin_conceiledm2"   "red_urchin_densitym2"       
#> [7] "red_urchin_conceiledm2"

# -----------------------------------------------------------------------------
# 2. Fit the Random Forest using an EXPLICIT formula
# -----------------------------------------------------------------------------
set.seed(1985)
rf_fit <- randomForest(
  site_type ~ relief_cm + risk_index +
    purple_urchin_densitym2 + purple_urchin_conceiledm2 +
    red_urchin_densitym2    + red_urchin_conceiledm2,
  data       = df_rf,
  ntree      = 1501,
  importance = TRUE
)


# -----------------------------------------------------------------------------
# 3a. Compute OOB accuracy per class
# -----------------------------------------------------------------------------
# rf_fit$confusion has a column "class.error" with OOB error rates
class_error <- rf_fit$confusion[, "class.error"]
class_acc   <- 1 - class_error
df_acc      <- tibble(
  site_type = names(class_acc),
  accuracy  = class_acc
) %>%
  mutate(site_type = factor(site_type, levels = levels(df_rf$site_type)))

# -----------------------------------------------------------------------------
# 3b. Prepare importance data (ranked within each facet)
# -----------------------------------------------------------------------------
imp_full <- importance(rf_fit) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  select(variable, all_of(levels(df_rf$site_type))) %>%
  pivot_longer(
    cols      = -variable,
    names_to  = "site_type",
    values_to = "importance"
  ) %>%
  mutate(
    site_type = factor(site_type, levels = levels(df_rf$site_type)),
    variable  = reorder_within(variable, -importance, site_type)
  )

# -----------------------------------------------------------------------------
# 3c. Plot with per‐panel OOB accuracy
# -----------------------------------------------------------------------------

# 1) pivot into long form
imp_long <- importance(rf_fit) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  select(variable, BAR, INCIP, FOR) %>%
  pivot_longer(-variable, names_to="site_type", values_to="importance") %>%
  mutate(site_type = factor(site_type, levels=c("BAR","INCIP","FOR")))

# 2) recode to friendly labels
imp_long <- imp_long %>%
  mutate(variable = fct_recode(variable,
                               "Relief"            = "relief_cm",
                               "Risk index"        = "risk_index",
                               "Purple density"    = "purple_urchin_densitym2",
                               "Purple concealed"  = "purple_urchin_conceiledm2",
                               "Red density"       = "red_urchin_densitym2",
                               "Red concealed"     = "red_urchin_conceiledm2"
  ))

# 3) reorder **descending** within each site_type
imp_long <- imp_long %>%
  mutate(variable = reorder_within(variable, -importance, site_type))

# 4) plot
imp_long %>%
  ggplot(aes(x = variable, y = importance, fill = site_type)) +
  geom_col() +
  facet_wrap(~ site_type, scales = "free_x", nrow = 1) +
  tidytext::scale_x_reordered() +
  scale_fill_manual(values = site_type_colors) +
  geom_text(
    data = df_acc,
    aes(x = Inf, y = Inf, label = paste0("OOB acc = ", round(accuracy,2))),
    inherit.aes = FALSE, hjust = 1.1, vjust = 1.1,
    fontface = "bold", size = 3
  ) +
  labs(
    x     = "Variable",
    y     = "Class-specific Importance",
    title = "Variable Importance by Site Type\n(with OOB accuracy)"
  ) +
  theme_bw(base_size = 14) +
  my_theme +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )















# -----------------------------------------------------------------------------
# 4. Multi‐class partial‐dependence curves for all three site types
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 4. Multi‐class partial‐dependence curves with pretty variable names
# -----------------------------------------------------------------------------

# define your pretty labels
var_labels <- c(
  "relief_cm"                 = "Relief",
  "risk_index"                = "Risk index",
  "purple_urchin_densitym2"   = "Purple density",
  "purple_urchin_conceiledm2" = "Purple concealed",
  "red_urchin_densitym2"      = "Red density",
  "red_urchin_conceiledm2"    = "Red concealed"
)

# a) define predictors and compute baseline at median of each
predictors <- setdiff(names(df_rf), "site_type")
baseline   <- df_rf %>%
  summarise(across(all_of(predictors), median, na.rm = TRUE))

# b) class labels
classes <- levels(df_rf$site_type)

# c) build tibble of predicted probabilities
pd_all <- map_dfr(predictors, function(var) {
  vals    <- seq(min(df_rf[[var]], na.rm=TRUE),
                 max(df_rf[[var]], na.rm=TRUE),
                 length.out = 50)
  newdata <- map_df(vals, ~ mutate(baseline, !!var := .x))
  probmat <- predict(rf_fit, newdata=newdata, type="prob")
  as_tibble(probmat) %>%
    mutate(value = vals, variable = var) %>%
    pivot_longer(cols     = all_of(classes),
                 names_to = "site_type",
                 values_to= "prob")
})

# recode raw variable names to pretty ones
pd_all <- pd_all %>%
  mutate(variable = recode(variable, !!!var_labels))

# site‐type colors
site_type_colors <- c(
  "FOR"   = "#1B9E77",
  "INCIP" = "#D95F02",
  "BAR"   = "#7570B3"
)

# d) plot
pd_all %>%
  ggplot(aes(x=value, y=prob, color=site_type)) +
  facet_wrap(~variable, scales="free", ncol=2) +
  geom_line(size=1) +
  labs(
    x     = "Predictor value",
    y     = "Predicted probability",
    color = "Site type",
    title = "Partial Dependence of P(site_type) on Each Predictor"
  ) +
  scale_color_manual(values=site_type_colors) +
  theme_bw(base_size=12) +
  my_theme
