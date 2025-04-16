

rm(list=ls())

require(librarian)

librarian::shelf(tidyverse, ggplot2, RColorBrewer, vegan, grid)

datdir <- "/Volumes/seaotterdb$/kelp_recovery/data/MBA_kelp_forest_database"

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

my_theme <- theme(axis.text = element_text(size = 7, color = "black"),
                  axis.title = element_text(size = 8, color = "black"),
                  legend.text = element_text(size = 7, color = "black"),
                  legend.title = element_text(size = 8, color = "black"),
                  strip.text = element_text(size = 8, hjust = 0, face = "bold", color = "black"),
                  strip.background = element_blank(),
                  plot.title = element_text(size = 9, color = "black"),
                  plot.tag = element_text(size = 9, color = "black", face = 'bold'),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  axis.line = element_line(colour = "black"),
                  legend.key.size = unit(0.5, "cm"),
                  legend.background = element_rect(fill = alpha('blue', 0)))

p1 <- ggplot(scores_cover %>%
               filter(site_type %in% c("FOR", "BAR"))
             , aes(x = NMDS1, y = NMDS2)) +
  stat_ellipse(aes(color = site_type), type = "norm", linetype = 1, size = 1) +
  # geom_point(aes(color = site_type, shape = site_type), size = 3, alpha = 0.5) +
  geom_point(data = centroids_cover %>%
               filter(site_type %in% c("FOR", "BAR"))
             , aes(x = NMDS1, y = NMDS2, color = site_type, shape = site_type),
             size = 4, alpha = 1) +
  scale_color_manual(values = c(
    "FOR" = "#1B9E77",
    "INCIP"  = "#D95F02",
    "BAR" = "#7570B3"
  )) +
  labs(title = "",
       x = "NMDS1", y = "NMDS2", color = "Site Type", shape = "Site Type") +
  theme_bw() + my_theme

p1



p2 <- ggplot(scores_cover, aes(x = NMDS1, y = NMDS2)) +
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

p4 <- ggplot(scores_density, aes(x = NMDS1, y = NMDS2)) +
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
p4




# Multiply density vectors by density_multiplier
phys_sig_density$NMDS1 <- phys_sig_density$NMDS1 * density_multiplier
phys_sig_density$NMDS2 <- phys_sig_density$NMDS2 * density_multiplier
bio_sig_density$NMDS1 <- bio_sig_density$NMDS1 * density_multiplier
bio_sig_density$NMDS2 <- bio_sig_density$NMDS2 * density_multiplier

p5 <- p4 +
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
  geom_text(data = bio_sig_density, inherit.aes = FALSE,
            aes(x = NMDS1, y = NMDS2, label = nice_variable),
            color = "black", vjust = -0.5, 
            nudge_x = 0.001, nudge_y = 0.001, size = 3)+
  #coord_cartesian(clip = "off") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
p5

################################################################################
#Toy with 3D nMDS

# Full 3D NMDS with convex‐hull bubbles and all (or more) vectors

# Assumes you’ve defined:
#   density_data_complete, meta_density,
#   physical_vars_density, density_multiplier,
#   and your nice_name() function

library(vegan)
library(plotly)
library(geometry)
library(dplyr)

# 1) 3‑dimensional NMDS
set.seed(123)
nmds3 <- metaMDS(
  density_data_complete,
  distance = "euclidean",
  k        = 3,
  trymax   = 100
)

# 2) Extract site scores and add site_type
scr3 <- as.data.frame(scores(nmds3, display = "sites"))
scr3$site_type <- meta_density$site_type

# 3) Compute convex‐hull centroids (optional)
centroids3 <- scr3 %>%
  group_by(site_type) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    NMDS3 = mean(NMDS3)
  )

# 4) Build matrices for manual 3D env‑fit
scores_mat <- as.matrix(scr3[, c("NMDS1","NMDS2","NMDS3")])
phys_mat   <- as.matrix(physical_vars_density)
bio_mat    <- as.matrix(density_data_complete)

# 5) Compute Pearson correlations of each predictor with each axis
cor_phys <- cor(phys_mat, scores_mat, use = "pairwise.complete.obs")
cor_bio  <- cor(bio_mat,  scores_mat, use = "pairwise.complete.obs")

# 6) Scale them by your multiplier
ar_phys3 <- cor_phys * density_multiplier
ar_bio3  <- cor_bio  * density_multiplier

# 7) Inspect vector lengths (overall 3D correlation)
len_phys <- sqrt(rowSums(ar_phys3^2))
len_bio  <- sqrt(rowSums(ar_bio3^2))
print(len_phys)
print(len_bio)

# 8) Decide which vectors to show
# Option A: keep only those above a lowered threshold (e.g. 0.1)
sig_phys <- ar_phys3[len_phys > 0.1, , drop = FALSE]
sig_bio  <- ar_bio3 [len_bio  > 0.1, , drop = FALSE]

# Option B: to show all vectors, uncomment these lines:
# sig_phys <- ar_phys3
# sig_bio  <- ar_bio3

# 9) Colors for site types
type_colors <- c("FOR" = "#1B9E77", "INCIP" = "#D95F02", "BAR" = "#7570B3")

# 10) Start empty Plotly canvas
p3d <- plot_ly()

# 11) Add convex‐hull mesh “bubbles” per site_type
for (st in unique(scr3$site_type)) {
  pts  <- filter(scr3, site_type == st)[, c("NMDS1","NMDS2","NMDS3")]
  hull <- convhulln(as.matrix(pts), options = "Qt")
  p3d <- p3d %>% add_trace(
    type       = "mesh3d",
    x          = pts$NMDS1,
    y          = pts$NMDS2,
    z          = pts$NMDS3,
    i          = hull[,1] - 1,
    j          = hull[,2] - 1,
    k          = hull[,3] - 1,
    opacity    = 0.2,
    color      = type_colors[st],
    name       = st,
    showlegend = TRUE
  )
}

# 12) (Optional) add centroids as markers
p3d <- p3d %>% add_trace(
  data       = centroids3,
  x          = ~NMDS1, y = ~NMDS2, z = ~NMDS3,
  type       = "scatter3d",
  mode       = "markers",
  marker     = list(size = 8, symbol = "diamond", opacity = 1),
  inherit    = FALSE,
  showlegend = FALSE
)

# 13) Add physical vectors as red lines
for (i in seq_len(nrow(sig_phys))) {
  df_vec <- data.frame(
    x = c(0, sig_phys[i,1]),
    y = c(0, sig_phys[i,2]),
    z = c(0, sig_phys[i,3])
  )
  p3d <- p3d %>% add_trace(
    data       = df_vec,
    x          = ~x, y = ~y, z = ~z,
    type       = "scatter3d",
    mode       = "lines",
    line       = list(color = "indianred", width = 4),
    inherit    = FALSE,
    showlegend = FALSE
  )
}

# 14) Add labels for physical vectors
df_phys_lbl <- data.frame(
  x   = sig_phys[,1],
  y   = sig_phys[,2],
  z   = sig_phys[,3],
  txt = nice_name(rownames(sig_phys))
)
p3d <- p3d %>% add_trace(
  data         = df_phys_lbl,
  x            = ~x, y = ~y, z = ~z, text = ~txt,
  type         = "scatter3d",
  mode         = "text",
  textposition = "top right",
  inherit      = FALSE,
  showlegend   = FALSE
)

# 15) Add biological vectors as gray lines
for (i in seq_len(nrow(sig_bio))) {
  df_vec <- data.frame(
    x = c(0, sig_bio[i,1]),
    y = c(0, sig_bio[i,2]),
    z = c(0, sig_bio[i,3])
  )
  p3d <- p3d %>% add_trace(
    data       = df_vec,
    x          = ~x, y = ~y, z = ~z,
    type       = "scatter3d",
    mode       = "lines",
    line       = list(color = "gray70", width = 2),
    inherit    = FALSE,
    showlegend = FALSE
  )
}

# 16) Add labels for biological vectors
df_bio_lbl <- data.frame(
  x   = sig_bio[,1],
  y   = sig_bio[,2],
  z   = sig_bio[,3],
  txt = nice_name(rownames(sig_bio))
)
p3d <- p3d %>% add_trace(
  data         = df_bio_lbl,
  x            = ~x, y = ~y, z = ~z, text = ~txt,
  type         = "scatter3d",
  mode         = "text",
  textposition = "bottom left",
  inherit      = FALSE,
  showlegend   = FALSE
)

# 17) Final layout
p3d %>% layout(
  scene = list(
    xaxis = list(title = "NMDS1"),
    yaxis = list(title = "NMDS2"),
    zaxis = list(title = "NMDS3")
  )
)
