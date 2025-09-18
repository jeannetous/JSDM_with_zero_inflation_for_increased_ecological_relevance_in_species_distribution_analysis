source("ZIPLN_simulations.R")

# n = 30 ; p = 10;  d = 1
# omega_structure = "erdos_renyi"
# zi_type = "covar"
# zi_covar_cluster = T
# min_X = 0 ;  max_X = 10 ; SNR = 0.75
# n_mode_zi_proba = c(3) ;  zi_mode_values = c(0, 0.3, 0.8)
# block_values = matrix(c(-10,-10, 2,-10,-10,2), nrow = 2)
# row_clusters = NULL
# col_clusters = NULL

# params <- generate_all_ZIPLN_parameters(n, p, d, omega_structure = omega_structure,
#                                         zi_type = zi_type,
#                                         zi_covar_cluster = zi_covar_cluster,
#                                         min_X = min_X, max_X = max_X, SNR = SNR,
#                                         n_mode_zi_proba = n_mode_zi_proba,
#                                         zi_mode_values = zi_mode_values,
#                                         block_values = block_values,
#                                         row_clusters = row_clusters,
#                                         col_clusters = col_clusters,
#                                         X0 = NULL, B0 = NULL,
#                                         min_X0 = 0, max_X0 = 10,
#                                         max_X0B0 = -0.2)
#
# Y <- simulate_ZIPLN_data(params)


# simu = 1
# simu_params =list(n = n, p = p, d = d, omega_structure = omega_structure,
#                   zi_type = zi_type,
#                   zi_covar_cluster = zi_covar_cluster,
#                   min_X = min_X, max_X = max_X, SNR = SNR,
#                   n_mode_zi_proba = n_mode_zi_proba,
#                   zi_mode_values = zi_mode_values,
#                   block_values = block_values,
#                   row_clusters = row_clusters,
#                   col_clusters = col_clusters,
#                   X0 = NULL, B0 = NULL,
#                   min_X0 = 0, max_X0 = 10,
#                   max_X0B0 = -0.2)
# PLN_formula = "Abundance ~ V1"
# ZIPLN_formula = "Abundance ~ V1 | VZI1"
# PLN_formula_ZIvar = "Abundance ~ V1 + VZI1"

# res <- one_ZIPLN_simulation(1, simu_params,
#                             PLN_formula, ZIPLN_formula,
#                             PLN_formula_ZIvar = PLN_formula_ZIvar)

# res <- multiple_ZIPLN_simulation(3, simu_params,
#                                  PLN_formula, ZIPLN_formula,
#                                  PLN_formula_ZIvar = NULL)

n_simu = 2
n_list = c(100, 300)
p_list = c(20)#, 100)
omega_structure_list = c("erdos_renyi")#, "community")
zi_type_list =  c("covar") #c("species", "sites")#, "covar")
n_mode_zi_proba_sites_list = c(2, 3)
zi_mode_values_sites_list = list(c(0.1, 0.3), c(0, 0.2, 0.4))
n_mode_zi_proba_species_list = c(2, 3)
zi_mode_values_species_list = list(c(0.1, 0.3), c(0, 0.2, 0.4))
block_values_list = list(matrix(c(-10,-10,2,-10,-10,2), nrow = 2))#, matrix(c(-10,-10, -10, 2,-10,0, -10,2,0), nrow = 3))


res <- grid_ZIPLN_simulation(n_simu, n_list, p_list, omega_structure_list,
                             zi_type_list, n_mode_zi_proba_sites_list,
                             zi_mode_values_sites_list, n_mode_zi_proba_species_list,
                             zi_mode_values_species_list, block_values_list,
                             min_X = 0,  max_X = 10, SNR = 0.75, min_X0 = 0,
                             max_X0 = 10, max_X0B0 = -0.2,
                             mc.cores = max(1, parallel::detectCores() - 2))
# write.csv(res, "ZIPLN_simus_test2.csv")



# simu_params = list(n = n,
#                    p = p,
#                    d = 1,
#                    omega_structure = omega_structure,
#                    zi_type = zi_type,
#                    zi_covar_cluster = TRUE,
#                    min_X = 0, max_X = 10, SNR = 0.75,
#                    n_mode_zi_proba = NULL,
#                    zi_mode_values = NULL,
#                    block_values = block_values,
#                    row_clusters = NULL,
#                    col_clusters = NULL,
#                    X0 = NULL, B0 = NULL,
#                    min_X0 = 0, max_X0 = 10,
#                    max_X0B0 = -0.2)
#
PLN_formula <- "Abundance ~ 0 + V1"
ZIPLN_formula <- "Abundance ~ 0 + V1 | 0 + VZI1"
PLN_formula_ZIvar <- "Abundance ~ 0 + V1 + VZI1"



# setting <- settings[1,]
# simu_params = list(n =  setting$n,
#                    p =  setting$p,
#                    d = 1,
#                    omega_structure =  setting$omega_structure,
#                    zi_type =  setting$zi_type,
#                    zi_covar_cluster = TRUE,
#                    min_X = min_X, max_X = max_X, SNR = SNR,
#                    n_mode_zi_proba =  setting$n_mode_zi_proba,
#                    zi_mode_values =  setting$zi_mode_values,
#                    block_values =  setting$block_values[[1]],
#                    row_clusters = NULL,
#                    col_clusters = NULL,
#                    X0 = NULL, B0 = NULL,
#                    min_X0 = min_X0, max_X0 = max_X0,
#                    max_X0B0 = max_X0B0)
#
# res <-  multiple_ZIPLN_simulations(
#   n_simu = 2, simu_params,
#   PLN_formula = PLN_formula,
#   ZIPLN_formula = ZIPLN_formula,
#   PLN_formula_ZIvar = PLN_formula_ZIvar,
# )
