source("ZIPLN_simulate_data.R")

n = 30 ; p = 10;  d = 1
omega_structure = "erdos_renyi"
zi_type = "sites"
zi_covar_cluster = T
min_X = 0 ;  max_X = 10 ; SNR = 0.75
n_mode_proba = c(3) ;  zi_mode_values = c(0, 0.3, 0.8)
block_values = matrix(c(-10,-10, 2,-10,-10,2), nrow = 2)
row_clusters = NULL
col_clusters = NULL

params <- generate_all_ZIPLN_parameters(n, p, d, omega_structure = omega_structure,
                                        zi_type = zi_type,
                                        zi_covar_cluster = zi_covar_cluster,
                                        min_X = min_X, max_X = max_X, SNR = SNR,
                                        n_mode_proba = n_mode_proba,
                                        zi_mode_values = zi_mode_values,
                                        block_values = block_values,
                                        row_clusters = row_clusters,
                                        col_clusters = col_clusters,
                                        X0 = NULL, B0 = NULL,
                                        min_X0 = 0, max_X0 = 10,
                                        max_X0B0 = -0.2)

Y <- simulate_ZIPLN_data(params)


simu = 1
simu_params =list(n = n, p = p, d = d, omega_structure = omega_structure,
                  zi_type = zi_type,
                  zi_covar_cluster = zi_covar_cluster,
                  min_X = min_X, max_X = max_X, SNR = SNR,
                  n_mode_proba = n_mode_proba,
                  zi_mode_values = zi_mode_values,
                  block_values = block_values,
                  row_clusters = row_clusters,
                  col_clusters = col_clusters,
                  X0 = NULL, B0 = NULL,
                  min_X0 = 0, max_X0 = 10,
                  max_X0B0 = -0.2)
PLN_formula = "Abundance ~ V1"
ZIPLN_formula = "Abundance ~ V1"
# PLN_formula_ZIvar = NULL

res <- one_ZIPLN_simulation(1, simu_params,
                            PLN_formula, ZIPLN_formula,
                            PLN_formula_ZIvar = NULL)
