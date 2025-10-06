# source("PLN_simulations.R")
source("ZIPLN_simulations.R")


n = 300 ; p = 20
mean_X = 0;   sd_X = 10; SNR = 2
omega_structure = "erdos_renyi"

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

zi_type = "covar"
n_mode_zi_proba = 3
zi_mode_values =  c(0.2, 0.4, 0.6)
proba_mode_zi = c(0.3, 0.3, 0.4)
block_values_list = list(block_values_ref)
row_clusters_proba_list = list(row_clusters_proba_ref)
col_clusters_proba_list = list(col_clusters_proba_ref)


simu_params = list(n = n,
                   p = p,
                   d = 1,
                   add_intercept = TRUE,
                   omega_structure = omega_structure,
                   zi_type = zi_type,
                   zi_covar_cluster = TRUE,
                   mean_X = mean_X, sd_X = sd_X, SNR = SNR,
                   n_mode_zi_proba = n_mode_zi_proba,
                   zi_mode_values = zi_mode_values,
                   proba_mode_zi = proba_mode_zi,
                   block_values = block_values_ref,
                   row_clusters = NULL,
                   col_clusters = NULL,
                   row_clusters_proba = row_clusters_proba_ref,
                   col_clusters_proba = col_clusters_proba_ref,
                   X0 = NULL, B0 = NULL,
                   mean_X0 = 0, sd_X0 = 10,
                   sd_X0B0 = 0.2)

for(i in 3:200){
  print(i)
  set.seed(i)
  res <- one_ZIPLN_simulation(simu = 1, zi_config = "sites_1", simu_params,
                              PLN_formula = "Abundance ~ 1 + V1",
                              ZIPLN_formula = "Abundance ~ 1 + V1 | 0 + VZI1",
                              PLN_formula_ZIvar = "Abundance ~ 1 + V1 + VZI1")
  print(res)
}

# res <- one_PLN_simulation(simu = 1, zi_config = "sites_1", simu_params,
#                           PLN_formula = "Abundance ~ 1 + V1",
#                           PLN_formula_ZIvar = NA)



