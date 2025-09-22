############################### Useful libraries ###############################
library(blockmodels)
set.seed(1)
################################# Data loading #################################
load("../data/NORTHERN_RANGE_DATASET_COUNTS.Rdata")
NORTHERN_RANGE_DATASET_COUNTS <- NORTHERN_RANGE_DATASET_COUNTS[, -21]
NR_ZEROS <- NORTHERN_RANGE_DATASET_COUNTS == 0


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
proba_per_site_clusters <- elements_per_site_clusters / sum(elements_per_site_clusters)

elements_per_species_clusters <- unlist(lapply(1:length(list_species_labels),
                                               f <- function(i) length(which(species_clusters == i))))
proba_per_species_clusters <- elements_per_species_clusters / sum(elements_per_species_clusters)

n_element_per_double_cluster <- outer(elements_per_site_clusters, elements_per_species_clusters)

zero_prop_per_cluster <- zero_per_cluster / n_element_per_double_cluster

zero_prop_per_cluster_logit <- apply(zero_prop_per_cluster, c(1, 2),
                                     f <- function(p){
                                       p <- min(0.999, max(p, 0.001))
                                       log(p/(1 - p))
                                     } )

verif <- apply(zero_prop_per_cluster_logit, c(1, 2),
               f <- function(x) exp(x) / (1 + exp(x)))

