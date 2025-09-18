library(PLNmodels)
library(tidyr)
library(dplyr)
source("ZIPLN_measures.R")
source("ZIPLN_simulate_data.R")

one_ZIPLN_simulation <- function(simu, simu_params,
                                 PLN_formula, ZIPLN_formula,
                                 PLN_formula_ZIvar = NA){

  params <- do.call(generate_all_ZIPLN_parameters, simu_params)
  Y <- simulate_ZIPLN_data(params)

  if(!is.null(params$zi_params)){
    X <- data.frame(params$X, params$zi_params$X0)
  }else{X <- params$X}
  simu_data <- prepare_data(Y, X)
  params$Y <- simu_data$Abundance



  ########################## Running PLN model #################################
  myPLN <- PLNnetwork(as.formula(PLN_formula), simu_data,
                      control = PLNnetwork_param(min_ratio = 0.05))
  PLN_StARS_measures <- get_measures(myPLN, params, model_selection = "StARS",
                                     stability = 0.8)
  PLN_BIC_measures <- get_measures(myPLN, params, model_selection = "BIC",
                                   AUC = PLN_StARS_measures[["AUC"]])


  ############### Running PLN model with ZI covar, if applicable ###############
  zi <- ifelse(simu_params$zi_type == "sites", "row",
               ifelse(simu_params$zi_type == "species", "col", "single") )
  if(!is.na(PLN_formula_ZIvar)){
    myPLN_ZIvar <- PLNnetwork(as.formula(PLN_formula_ZIvar), simu_data,
                        control = PLNnetwork_param(min_ratio = 0.05))
    PLN_ZIvar_StARS_measures <- get_measures(myPLN_ZIvar, params, model_selection = "StARS",
                                       stability = 0.8)
    PLN_ZIvar_BIC_measures <- get_measures(myPLN_ZIvar, params, model_selection = "BIC",
                                     AUC = PLN_ZIvar_StARS_measures[["AUC"]])
  }else{
    PLN_ZIvar_StARS_measures <- NULL ; PLN_ZIvar_BIC_measures <- NULL
  }


  ######################### Running ZIPLN model ################################
  myZIPLN <- ZIPLNnetwork(as.formula(ZIPLN_formula), simu_data, zi = "row",
                          control = ZIPLNnetwork_param(min_ratio = 0.05))
  ZIPLN_StARS_measures <- get_measures(myZIPLN, params, model_selection = "StARS",
                                       stability = 0.8)
  ZIPLN_BIC_measures <- get_measures(myZIPLN, params, model_selection = "BIC",
                                     AUC = ZIPLN_StARS_measures[["AUC"]])

  ################# Merging all the measures in one data frame #################
  measure_rows <- list(c(method = "PLN", PLN_StARS_measures),
                        c(method = "PLN", PLN_BIC_measures),
                        c(method = "ZIPLN", ZIPLN_StARS_measures),
                        c(method = "ZIPLN", ZIPLN_BIC_measures)
                        )

  if(!is.na(PLN_formula_ZIvar)){
    measure_rows <- c(measure_rows,
                      list(c(method = "PLN_ZIvar", PLN_ZIvar_StARS_measures),
                           c(method = "PLN_ZIvar", PLN_ZIvar_BIC_measures)))
  }
  res <- as.data.frame(cbind(simu = simu, n = simu_params$n, p = simu_params$p,
                             zi_type = simu_params$zi_type, do.call(rbind, measure_rows)))
  return(res)
}


multiple_ZIPLN_simulations <- function(n_simu, simu_params,
                                       PLN_formula, ZIPLN_formula,
                                       PLN_formula_ZIvar = NULL,
                                       min_X = 0,  max_X = 10, SNR = 0.75,
                                       mc.cores = max(1, parallel::detectCores() - 2)){
  cat("Settings: (n, p, omega structure, zi type) = (",simu_params$n, simu_params$p,
      simu_params$omega_structure, simu_params$zi_type, ")\n")

  multiple_res <- parallel::mclapply(1:n_simu,
                                     one_ZIPLN_simulation,
                                     simu_params = simu_params,
                                     PLN_formula = PLN_formula,
                                     ZIPLN_formula = ZIPLN_formula,
                                     PLN_formula_ZIvar = PLN_formula_ZIvar,
                                     mc.cores = mc.cores)
  res <- do.call(rbind, multiple_res)
  res
}

grid_ZIPLN_simulation <- function(n_simu, n_list, p_list, omega_structure_list,
                                  zi_type_list, n_mode_zi_proba_sites_list,
                                  zi_mode_values_sites_list, n_mode_zi_proba_species_list,
                                  zi_mode_values_species_list, block_values_list,
                                  min_X = 0,  max_X = 10, SNR = 0.75,
                                  min_X0 = 0, max_X0 = 10, max_X0B0 = -0.2,
                                  mc.cores = max(1, parallel::detectCores() - 2)) {
  settings <- expand.grid(n = n_list,
                          p = p_list,
                          omega_structure = omega_structure_list,
                          zi_type = zi_type_list,
                          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
                          )

  sites_zi_pairs <- tibble(n_mode_zi_proba = n_mode_zi_proba_sites_list,
                           zi_mode_values  = zi_mode_values_sites_list)
  species_zi_pairs <- tibble(n_mode_zi_proba = n_mode_zi_proba_species_list,
                             zi_mode_values  = zi_mode_values_species_list)


  settings <- settings %>%
    rowwise() %>%
    do({
      row <- .
      if (row$zi_type == "covar") {
        # Keep row with NA for E and F
        tibble(n = row$n, p = row$p, omega_structure = row$omega_structure,
               zi_type = row$zi_type, n_mode_zi_proba = NA, zi_mode_values = NA)
      } else if (row$zi_type == "sites") {
        # Expand with lookup1
        cbind(row[1:4], sites_zi_pairs)
      } else if (row$zi_type == "species") {
        cbind(row[1:4], species_zi_pairs)
      }
    }) %>%
    ungroup()

  block_values_rows <- settings %>%
                       filter(zi_type == "covar") %>%
                       mutate(block_values = list(block_values_list)) %>%
                       unnest(cols = c(block_values))

  settings <- settings %>% filter(zi_type != "covar")
  settings <- settings %>% mutate(block_values = NA)
  settings <- rbind(settings, block_values_rows)

  settings$PLN_formula <- "Abundance ~ 0 + V1"
  settings$ZIPLN_formula <- "Abundance ~ 0 + V1"
  settings[settings$zi_type == "covar",]$ZIPLN_formula <- "Abundance ~ 0 + V1 | 0 + VZI1"
  settings$PLN_formula_ZIvar <- NA
  settings[settings$zi_type == "covar",]$PLN_formula_ZIvar <- "Abundance ~ 0 + V1 + VZI1"

  settings$n_simu <- n_simu

  final_res <- purrr::pmap(settings, f <- function(n, p, omega_structure, zi_type,
                                                   n_mode_zi_proba, zi_mode_values,
                                                   block_values, PLN_formula,
                                                   ZIPLN_formula, PLN_formula_ZIvar,
                                                   n_simu){

                                          simu_params = list(n = n,
                                                             p = p,
                                                             d = 1,
                                                             omega_structure = omega_structure,
                                                             zi_type = zi_type,
                                                             zi_covar_cluster = TRUE,
                                                             min_X = min_X, max_X = max_X, SNR = SNR,
                                                             n_mode_zi_proba = n_mode_zi_proba,
                                                             zi_mode_values = zi_mode_values,
                                                             block_values = block_values,
                                                             row_clusters = NULL,
                                                             col_clusters = NULL,
                                                             X0 = NULL, B0 = NULL,
                                                             min_X0 = min_X0, max_X0 = max_X0,
                                                             max_X0B0 = max_X0B0)
                                          multiple_ZIPLN_simulations(
                                            n_simu = n_simu, simu_params,
                                            PLN_formula = PLN_formula,
                                            ZIPLN_formula = ZIPLN_formula,
                                            PLN_formula_ZIvar = PLN_formula_ZIvar,
                                          )})
  final_res <- do.call(rbind, final_res) %>% as_tibble()
  final_res
}
