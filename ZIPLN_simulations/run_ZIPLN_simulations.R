############################ Loading useful libraries ###########################
source("ZIPLN_simulations.R")
set.seed(2)
############################ Reference ZI values from real data - ZERO COUNTS ##
##### SITES #####
# Stream-aggregated
# n_mode_zi_proba_stream_ref = 3
# zi_mode_values_stream_ref = c(0.5, 0.7, 0.8)
# proba_mode_zi_values_stream_ref = c(0.25, 0.44, 0.31)
#
# # Site-specific
# n_mode_zi_proba_site_ref = 3
# zi_mode_values_site_ref = c(0.5, 0.7, 0.8)
# proba_mode_zi_values_site_ref = c(0.35, 0.45, 0.2)
#
# ##### SPECIES #####
# n_mode_zi_proba_species_ref = 4
# zi_mode_values_species_ref = c(0.13, 0.5, 0.8, 0.95)
# proba_mode_zi_values_species_ref = c(0.15, 0.3, 0.15, 0.4)
#
# ##### SITES X SPECIES #####
#
#
# block_values_ref <- matrix(c(0.7, 7, 6.1, 2.9,
#                              -1.4, 7, -1.5, 0.4,
#                              3.7, 2.9, 4.7, 7,
#                              -3.1, -1.2, -3.2, 0.1,
#                              -2.3, 4, 2.7, -3.2,
#                              -0.2, -3.7, 0.6, 0.8,
#                              0.9, 1.8, 1.7, 7,
#                              -0.7, 0.8, 0.05, 3.9), nrow = 4)
# row_clusters_proba_ref <- c(0.27, 0.28, 0.37, 0.08)
# col_clusters_proba_ref <- c(0.20, 0.10, 0.25, 0.15, 0.10, 0.05, 0.05, 0.10)

############################ Reference ZI values from real data - ZIPLN ##
##### SITES #####
n_mode_zi_proba_sites_ref = 3
zi_mode_values_sites_ref = c(0.2, 0.4, 0.6)
proba_mode_zi_values_sites_ref = c(0.3, 0.3, 0.4)

##### SPECIES #####
n_mode_zi_proba_species_ref = 4
zi_mode_values_species_ref = c(0, 0.3, 0.6, 0.9)
proba_mode_zi_values_species_ref = c(0.4, 0.15, 0.15, 0.3)

##### SITES X SPECIES #####


block_values_ref <- matrix(c(-0.8, 1.5, 3.2, -0.8, -7,
                             -2.2, 7, -1, -2.7, 7,
                             -0.7, 7, -1, -0.6, -7,
                             -0.7, 7, 7, -3.6, -7,
                              7, 7, 7, -1.7, -7,
                              7, 7, 0.8, -0.4, -7,
                              -2, -4.8, -4.4, -1.9, -7,
                              7, -2, 7, -1, -7,
                              7, 7, 7, 0.6, -7,
                              7, 7, 7, -0.1, -7,
                              -0.1, 7, 7, 0, -7,
                              7, 7, 7, 0.2, -7 ), nrow = 5)
row_clusters_proba_ref <- c(0.125, 0.125, 0.125, 0.125, 0.5)
col_clusters_proba_ref <- rep(0.0833, 12)

############################ Simulations parameters ############################
n_simu = 50
n_list = c(300)
p_list = c(20) #, 100)
omega_structure_list = c("erdos_renyi", "community", "preferential_attachment")
zi_type_list =  c("sites", "species")#, "covar")

n_mode_zi_proba_sites_list = c(3, 3, 3)
zi_mode_values_sites_list = list(0.1 * zi_mode_values_sites_ref,
                                 0.5 * zi_mode_values_sites_ref,
                                 zi_mode_values_sites_ref)
proba_mode_zi_sites_list = list(proba_mode_zi_values_sites_ref,
                                proba_mode_zi_values_sites_ref,
                                proba_mode_zi_values_sites_ref)


n_mode_zi_proba_species_list = c(4, 4, 4)
zi_mode_values_species_list = list(0.1 * zi_mode_values_species_ref,
                                   0.5 * zi_mode_values_species_ref,
                                   zi_mode_values_species_ref)
proba_mode_zi_species_list = list(proba_mode_zi_values_species_ref)

block_values_list = list(0.1 * block_values_ref,
                         0.5 * block_values_ref)#,
                         # block_values_ref) #, matrix(c(-10,-10, -10, 2,-10,0, -10,2,0), nrow = 3))
row_clusters_proba_list = list(row_clusters_proba_ref, row_clusters_proba_ref)#, row_clusters_proba_ref)#, c(0.4, 0.4, 0.2))
col_clusters_proba_list = list(col_clusters_proba_ref, col_clusters_proba_ref)#, col_clusters_proba_ref)#, c(0.75, 0.125, 0.125))


# n_mode_zi_proba_species_list = NULL ; zi_mode_values_species_list = NULL ; proba_mode_zi_species_list = NULL
# block_values_list = NULL ; row_clusters_proba_list = NULL ; col_clusters_proba_list = NULL

############################ Running simulations ###############################
res <- grid_ZIPLN_simulation(n_simu, n_list, p_list, omega_structure_list,
                             zi_type_list, n_mode_zi_proba_sites_list = n_mode_zi_proba_sites_list,
                             zi_mode_values_sites_list = zi_mode_values_sites_list,
                             proba_mode_zi_sites_list = zi_mode_values_sites_list,
                             n_mode_zi_proba_species_list,
                             zi_mode_values_species_list,
                             proba_mode_zi_species_list, block_values_list,
                             row_clusters_proba_list, col_clusters_proba_list,
                             min_X = 0,  max_X = 10, SNR = 0.75, min_X0 = 0,
                             max_X0 = 10, max_X0B0 = -0.2,
                             mc.cores = max(1, parallel::detectCores() - 2))

############################ Saving the results and parameters #################
write.csv(res, "ZIPLN_simulations_res/ZIPLN_simus_zi_from_real_zi_proba_PLNonly_BIConly_1.csv")


all_params <- list(n_simu = n_simu, n_list = n_list, p_list = p_list,
                   omega_structure_list = omega_structure_list,
                   zi_type_list = zi_type_list,
                   n_mode_zi_proba_sites_list = n_mode_zi_proba_sites_list,
                   zi_mode_values_sites_list = zi_mode_values_sites_list,
                   proba_mode_zi_sites_list = proba_mode_zi_sites_list,
                   n_mode_zi_proba_species_list = n_mode_zi_proba_species_list,
                   zi_mode_values_species_list = zi_mode_values_species_list,
                   proba_mode_zi_species_list = proba_mode_zi_species_list)
writeLines(capture.output(str(all_params)), "ZIPLN_simulations_res/ZIPLN_simus_zi_from_real_zi_proba_PLNonly_BIConly_1_parameters.txt")

############################ Debugging bits ####################################

PLN_formula <- "Abundance ~ 0 + V1"
ZIPLN_formula <- "Abundance ~ 0 + V1" # | 0 + VZI1"
# PLN_formula_ZIvar <- "Abundance ~ 0 " #+ V1 + VZI1"



# simu_params = list(n = setting$n,
#                    p = setting$p,
#                    d = 1,
#                    omega_structure = setting$omega_structure,
#                    zi_type = setting$zi_type,
#                    zi_covar_cluster = TRUE,
#                    min_X = min_X, max_X = max_X, SNR = SNR,
#                    n_mode_zi_proba = setting$n_mode_zi_proba,
#                    zi_mode_values = setting$zi_mode_values,
#                    proba_mode_zi = setting$proba_mode_zi,
#                    block_values = setting$block_values,
#                    row_clusters = NULL,
#                    col_clusters = NULL,
#                    row_clusters_proba = setting$row_clusters_proba,
#                    col_clusters_proba = setting$col_clusters_proba,
#                    X0 = NULL, B0 = NULL,
#                    min_X0 = min_X0, max_X0 = max_X0,
#                    max_X0B0 = max_X0B0)



# n = 300 ; p = 20 ; zi_type = "sites"; block_values = block_values_ref
# row_clusters_proba = row_clusters_proba_ref ; col_clusters_proba = col_clusters_proba_ref
#
# simu_params = list(n = n,
#                    p = p,
#                    d = 1,
#                    omega_structure = "erdos_renyi",
#                    zi_type = "sites",
#                    zi_covar_cluster = TRUE,
#                    min_X = 0, max_X = 10, SNR = 0.75,
#                    n_mode_zi_proba = 3,
#                    zi_mode_values = 0.5 * zi_mode_values_sites_ref,
#                    proba_mode_zi = zi_mode_values_sites_ref,
#                    block_values = NULL,
#                    row_clusters = NULL,
#                    col_clusters = NULL,
#                    row_clusters_proba = NA,
#                    col_clusters_proba = NA,
#                    X0 = NULL, B0 = NULL,
#                    min_X0 = 0, max_X0 = 10,
#                    max_X0B0 = 0.2)
#
# res <- one_ZIPLN_simulation(1, "sites_1", simu_params,
#                             PLN_formula, ZIPLN_formula,
#                             PLN_formula_ZIvar = NA)
#


