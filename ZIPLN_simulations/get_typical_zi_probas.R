############################### Useful libraries ###############################
library(PLNmodels)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(blockmodels)
set.seed(1)
############################### Data loading ###################################
load("data/NORTHERN_RANGE_DATASET_COUNTS.Rdata")
NR_COUNTS <- NORTHERN_RANGE_DATASET_COUNTS %>% as_tibble()
NR_COUNTS <- NR_COUNTS[, -21]
NORTHERN_RANGE_DATASET_COUNTS <- NORTHERN_RANGE_DATASET_COUNTS[, -21]
NR_ZEROS <- NORTHERN_RANGE_DATASET_COUNTS == 0

load("data/NORTHERN_RANGE_DATASET_COVARIATES.Rdata")
NR_COVARIATES <-
  NORTHERN_RANGE_DATASET_COVARIATES %>% as_tibble() %>%
  rename(latitude = LATITUDE,
         longitude = LONGITUDE,
         long_lat = LONGLAT,
         coarse_gravel = coarse.gravel, # We do not include geographical coordinates here
         fine_gravel = fine.gravel,
         leaf_litter = leaf.litter,
         time_step = TimeStep, year = YEAR, month = MONTH,
         stream = STREAM,
         disturbance = DISTURBANCE, # human activity / binary categorical variable
         altitude = ALTITUDE,
         season = SEASON) %>%
  mutate(
    # season = case_match(time_step,
    #   c(1, 5,  9, 13, 17) ~ "dry_start",
    #   c(2, 6, 10, 14, 18) ~ "dry_end",
    #   c(3, 7, 11, 15, 19) ~ "wet_start",
    #   c(4, 8, 12, 16    ) ~ "wet_end"
    # ),
    season = as.character(season),
    year = as.character(year),
    disturbance = ifelse(disturbance == 1, "yes", "no")
  ) %>%
  ## flow is not reliable
  ## time_step + month are redundant with season + year
  ## site is redudant with disturbance + stream
  ## long_lat is the 1st PCA axis of latitude + longitude
  dplyr::select(-latitude, -longitude, - long_lat, -site, -time_step,
                -month, -flow) %>% # we remove geographical coordinates
  relocate(altitude, width, depth, volume, garbage,
           conductivity, O2, pH, temperature, turbidity, # site features
           coarse_gravel, fine_gravel, leaf_litter, cobble, sand, silt, boulders, canopy, # soil
           season, year, # sampling time
           stream, disturbance # categorical: place + human activity
  )
# turbidity: ordinal variable indicating the level of turbidity
# season:-> dry-start/dry-end, wet-start/wet-end, (january-> may: dry)
quali_ind <- seq(ncol(NR_COVARIATES) - 3, ncol(NR_COVARIATES))

############################### Preparing the covar for ZIPLN ##################

nb_pc <- ncol(NR_COVARIATES) - 4
NR_PCA <-
  NR_COVARIATES %>%
  PCA(quali.sup = quali_ind, ncp = nb_pc, scale.unit = TRUE, graph = FALSE)

NR_PCA_SCORES <- setNames(data.frame(NR_PCA$ind$coord), paste0("PC",1:nb_pc))

NR_DATA <- prepare_data(
  counts     = NR_COUNTS,
  covariates = cbind(NR_COVARIATES, NR_PCA_SCORES)
)

############################### Running ZIPLN ##################################

model_ZIPLN_covar  <- ZIPLNnetwork(Abundance ~ PC1 + PC2 + disturbance | stream,  data = NR_DATA,
                                   control = ZIPLNnetwork_param(min_ratio = 0.05,
                                                                penalize_diagonal = FALSE))
model_ZIPLN_covar_BIC <- model_ZIPLN_covar$getBestModel("BIC")
zi_proba_crossed <- model_ZIPLN_covar_BIC$model_par$Pi

model_ZIPLN_site  <- ZIPLNnetwork(Abundance ~ PC1 + PC2 + disturbance,  data = NR_DATA,
                                  zi = "row",
                                  control = ZIPLNnetwork_param(min_ratio = 0.05,
                                                               penalize_diagonal = FALSE))
model_ZIPLN_site_BIC <- model_ZIPLN_site$getBestModel("BIC")
zi_proba_sites <- model_ZIPLN_site_BIC$model_par$Pi[,1]

model_ZIPLN_species  <- ZIPLNnetwork(Abundance ~ PC1 + PC2 + disturbance,  data = NR_DATA,
                                     zi = "col",
                                     control = ZIPLNnetwork_param(min_ratio = 0.05,
                                                                  penalize_diagonal = FALSE))
model_ZIPLN_species_BIC <- model_ZIPLN_species$getBestModel("BIC")
zi_proba_species <- model_ZIPLN_species_BIC$model_par$Pi[1,]
############################### ZI PROBA WITH ZERO COUNTS ######################
############################### Sites ZI proba #################################
# Grouping by stream
streams <- sub("_.*", "", rownames(NORTHERN_RANGE_DATASET_COUNTS))
row_zero_prop <- rowMeans(NORTHERN_RANGE_DATASET_COUNTS == 0)
mean_prop_per_stream <- tapply(row_zero_prop, streams, mean)

# Leaving sites separated
mean_prop_per_site <- rowMeans(NORTHERN_RANGE_DATASET_COUNTS == 0)
cl_sites <- kmeans(mean_prop_per_site, centers = 3)

############################### Species ZI proba ###############################
mean_prop_per_species <- colMeans(NORTHERN_RANGE_DATASET_COUNTS == 0)
cl_species <- kmeans(mean_prop_per_species, centers = 4)

############################### Crossed sites-species ZI proba #################
cl_sites_species <- BM_bernoulli("LBM", NR_ZEROS)
cl_sites_species$estimate()
cl_sites_species_best_ICL <- cl_sites_species$memberships[[12]]

sites_clusters <- apply(cl_sites_species_best_ICL$Z1, c(1), which.max)
species_clusters <- apply(cl_sites_species_best_ICL$Z2, c(1), which.max)
list_sites_labels <- sort(unique(sites_clusters))
list_species_labels <- sort(unique(species_clusters))
sites_ind <- sapply(sites_clusters, function(cl) list_sites_labels == cl)
species_ind <- sapply(species_clusters, function(cl) list_species_labels == cl)


zero_per_cluster <- sites_ind %*% NR_ZEROS %*% t(species_ind)

elements_per_site_clusters <- unlist(lapply(1:length(list_sites_labels),
                                            f <- function(i) length(which(sites_clusters == i))))
elements_per_species_clusters <- unlist(lapply(1:length(list_species_labels),
                                               f <- function(i) length(which(species_clusters == i))))
n_element_per_double_cluster <- outer(elements_per_site_clusters, elements_per_species_clusters)

zero_prop_per_cluster <- zero_per_cluster / n_element_per_double_cluster

zero_prop_per_cluster_logit <- apply(zero_prop_per_cluster, c(1, 2),
                                     f <- function(p) log(p/(1 - p)))


############################### ZI PROBA WITH PLN ZI PROBA #####################
############################### Sites ZI proba [PLN] ###########################
cl_sites_pln <- kmeans(zi_proba_sites, centers = 3)

############################### Species ZI proba [PLN]##########################
cl_species_pln <- kmeans(zi_proba_species, centers = 4)

############################### Crossed sites-species ZI proba [PLN] ###########
cl_sites_species_pln <- BM_gaussian("LBM", zi_proba_crossed)
cl_sites_species_pln$estimate()
cl_sites_species_pln_best_ICL <- cl_sites_species_pln$memberships[[19]]

sites_clusters_pln <- apply(cl_sites_species_pln_best_ICL$Z1, c(1), which.max)
species_clusters_pln <- apply(cl_sites_species_pln_best_ICL$Z2, c(1), which.max)
list_sites_labels_pln <- sort(unique(sites_clusters_pln))
list_species_labels_pln <- sort(unique(species_clusters_pln))
sites_ind_pln <- sapply(sites_clusters, function(cl) list_sites_labels_pln == cl)
species_ind_pln <- sapply(species_clusters_pln, function(cl) list_species_labels_pln == cl)


zero_per_cluster_pln <- sites_ind_pln %*% zi_proba_crossed %*% t(species_ind_pln)

elements_per_site_clusters_pln <- unlist(lapply(1:length(list_sites_labels_pln),
                                                f <- function(i) length(which(sites_clusters_pln == i))))
elements_per_species_clusters_pln <- unlist(lapply(1:length(list_species_labels_pln),
                                                   f <- function(i) length(which(species_clusters_pln == i))))
n_element_per_double_cluster_pln <- outer(elements_per_site_clusters_pln, elements_per_species_clusters_pln)

zero_prop_per_cluster_pln <- zero_per_cluster_pln / n_element_per_double_cluster_pln

zero_prop_per_cluster_logit_pln <- apply(zero_prop_per_cluster_pln, c(1, 2),
                                         f <- function(p){
                                           p <- min(0.999, max(p, 0.001))
                                           log(p/(1 - p))
                                         })




