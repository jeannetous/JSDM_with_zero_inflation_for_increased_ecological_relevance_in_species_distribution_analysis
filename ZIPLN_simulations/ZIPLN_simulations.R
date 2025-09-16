library(PLNmodels)
source("ZIPLN_measures.R")
source("ZIPLN_simulate_data.R")

one_ZIPLN_simulation <- function(simu, simu_params,
                                 PLN_formula, ZIPLN_formula,
                                 PLN_formula_ZIvar = NULL){
  params <- do.call(generate_all_ZIPLN_parameters, simu_params)
  Y <- simulate_ZIPLN_data(params)
  if(!is.null(params$zi_params)){
    X <- cbind(params$X, params$zi_params$X0_num)
  }
  simu_data <- prepare_data(Y, params$X)

  ########################## Running PLN model #################################
  myPLN <- PLNnetwork(as.formula(PLN_formula), simu_data,
                      control = PLNnetwork_param(min_ratio = 0.05))
  PLN_StARS_measures <- get_measures(myPLN, params, model_selection = "StARS",
                                     stability = 0.8)
  PLN_BIC_measures <- get_measures(myPLN, params, model_selection = "BIC",
                                   AUC = PLN_StARS_measures[["AUC"]])

  ############### Running PLN model with ZI covar, if applicable ###############
  if(!is.null(PLN_formula_ZIvar)){
    myPLN_ZIvar <- PLNnetwork(as.formula(PLN_formula_ZIvar), simu_data,
                        control = PLNnetwork_param(min_ratio = 0.05))
    PLN_ZIvar_StARS_measures <- get_measures(myPLN_ZIvar, params, model_selection = "StARS",
                                       stability = 0.8)
    PLN_ZIvar_BIC_measures <- get_measures(myPLN_ZIvar, params, model_selection = "BIC",
                                     AUC = PLN_ZIvar_BIC_measures[["AUC"]])
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
  res <- as.data.frame(
    cbind(simu = simu, n = simu_params$n, p = simu_params$p,
          zi_type = simu_params$zi_type,
          omega_structure = simu_params$omega_structure,
          rbind(c(method = "PLN", PLN_StARS_measures),
                c(method = "PLN", PLN_BIC_measures),
                c(method = "PLN_ZIvar", PLN_ZIvar_StARS_measures),
                c(method = "PLN_ZIvar", PLN_ZIvar_BIC_measures),
                c(method = "ZIPLN", ZIPLN_StARS_measures),
                c(method = "ZIPLN", ZIPLN_BIC_measures)
          )
    )
  )
  return(res)
}
