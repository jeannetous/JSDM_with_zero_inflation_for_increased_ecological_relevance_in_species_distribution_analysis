############################ Loading useful libraries ###########################
source("PLN_simulations.R")
set.seed(2)


############################ Reference ZI values from real data - ZIPLN ########
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
n_simu = 3
n_list = c(300)
p_list = c(20) #, 100)
omega_structure_list = c("erdos_renyi")# , "community", "preferential_attachment")
zi_type_list =  c("covar")#c("sites", "species")

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
                         0.5 * block_values_ref,
                         block_values_ref)
row_clusters_proba_list = list(row_clusters_proba_ref, row_clusters_proba_ref, row_clusters_proba_ref)
col_clusters_proba_list = list(col_clusters_proba_ref, col_clusters_proba_ref, col_clusters_proba_ref)


# n_mode_zi_proba_species_list = NULL ; zi_mode_values_species_list = NULL ; proba_mode_zi_species_list = NULL
# block_values_list = NULL ; row_clusters_proba_list = NULL ; col_clusters_proba_list = NULL

############################ Running simulations ###############################
res <- grid_PLN_simulation(n_simu, n_list, p_list, add_intercept = TRUE,
                           omega_structure_list, zi_type_list,
                           n_mode_zi_proba_sites_list = n_mode_zi_proba_sites_list,
                           zi_mode_values_sites_list = zi_mode_values_sites_list,
                           proba_mode_zi_sites_list = zi_mode_values_sites_list,
                           n_mode_zi_proba_species_list,
                           zi_mode_values_species_list,
                           proba_mode_zi_species_list, block_values_list,
                           row_clusters_proba_list, col_clusters_proba_list,
                           mean_X = 0,  sd_X = 1, SNR = 0.75, mean_X0 = 0,
                           sd_X0 = 10, max_X0B0 = -0.2,
                           mc.cores = max(1, parallel::detectCores() - 2))

############################ Saving the results and parameters #################
# write.csv(res, "ZIPLN_simulations_res/ZIPLN_simus_zi_from_real_zi_proba_PLNonly_BIConly_2.csv")


all_params <- list(n_simu = n_simu, n_list = n_list, p_list = p_list,
                   add_intercept = "TRUE",
                   omega_structure_list = omega_structure_list,
                   zi_type_list = zi_type_list,
                   n_mode_zi_proba_sites_list = n_mode_zi_proba_sites_list,
                   zi_mode_values_sites_list = zi_mode_values_sites_list,
                   proba_mode_zi_sites_list = proba_mode_zi_sites_list,
                   n_mode_zi_proba_species_list = n_mode_zi_proba_species_list,
                   zi_mode_values_species_list = zi_mode_values_species_list,
                   proba_mode_zi_species_list = proba_mode_zi_species_list)
# writeLines(capture.output(str(all_params)), "ZIPLN_simulations_res/ZIPLN_simus_zi_from_real_zi_proba_PLNonly_BIConly_2_parameters.txt")

############################ Debugging bits ####################################

PLN_formula <- "Abundance ~ 1 + V1"
PLN_formula_ZIvar <- NA
# ZIPLN_formula <- "Abundance ~ 0 + V1" # | 0 + VZI1"
# PLN_formula_ZIvar <- "Abundance ~ 0 " #+ V1 + VZI1"

n = 300 ; p = 20
mean_X = 0;   sd_X = 10; SNR = 0.75
omega_structure = "erdos_renyi"

zi_type = "sites"
zi_config = "sites_1"
n_mode_zi_proba = 3
zi_mode_values =  c(0.2, 0.4, 0.6)
proba_mode_zi = c(0.3, 0.3, 0.4)

simu_params = list(n = n,
                   p = p,
                   d = 1,
                   add_intercept = TRUE,
                   omega_structure = omega_structure,
                   zi_type = zi_type,
                   zi_covar_cluster = TRUE,
                   mean_X = 0, sd_X = 10, SNR = 0.75,
                   n_mode_zi_proba = n_mode_zi_proba,
                   zi_mode_values = zi_mode_values,
                   proba_mode_zi = proba_mode_zi,
                   block_values = NULL,
                   row_clusters = NULL,
                   col_clusters = NULL,
                   row_clusters_proba = NULL,
                   col_clusters_proba = NULL,
                   X0 = NULL, B0 = NULL,
                   mean_X0 = 0, sd_X0 = 10,
                   max_X0B0 = 0.2)

# test <- one_PLN_simulation(simu = 1, zi_config, simu_params,
#                           PLN_formula, PLN_formula_ZIvar = NA)
