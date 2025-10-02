source("PLN_simulations.R")
source("ZIPLN_simulations.R")


n = 300 ; p = 20
min_X = 0;   max_X = 10; SNR = 2
omega_structure = "erdos_renyi"

zi_type = "sites"
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
                   min_X = min_X, max_X = max_X, SNR = SNR,
                   n_mode_zi_proba = n_mode_zi_proba,
                   zi_mode_values = zi_mode_values,
                   proba_mode_zi = proba_mode_zi,
                   block_values = NULL,
                   row_clusters = NULL,
                   col_clusters = NULL,
                   row_clusters_proba = NULL,
                   col_clusters_proba = NULL,
                   X0 = NULL, B0 = NULL,
                   min_X0 = 0, max_X0 = 10,
                   max_X0B0 = 0.2)


set.seed(139)
# res <- one_PLN_simulation(simu = 1, zi_config = "sites_1", simu_params,
#                           PLN_formula = "Abundance ~ 1 + V1",
#                           PLN_formula_ZIvar = NA)
res <- one_ZIPLN_simulation(simu = 1, zi_config = "sites_1", simu_params,
                          PLN_formula = "Abundance ~ 1 + V1",
                          ZIPLN_formula = "Abundance ~ 1 +V1",
                          PLN_formula_ZIvar = NA)


print(res)
