library(Metrics)

#' @description computes recall, fallout, precision and f1-score for
#' an inferred network (precision matrices here) given the true one
#' @param omega_true true precision matrix
#' @param omega_estimate estimated precision matrix
roc_metrics <- function(omega_true, omega_estimate){

  diag(omega_true) <- 0 ; diag(omega_estimate) <- 0

  true.nzero <- which(omega_true != 0)
  true.zero  <- which(omega_true == 0)

  nzero <- which(omega_estimate != 0)
  zero  <- which(omega_estimate == 0)

  TP <- 0.5 * sum(nzero %in% true.nzero)
  TN <- 0.5 * (sum(zero %in%  true.zero) - nrow(omega_true)) # removing diagonal values that do not count
  FP <- 0.5 * sum(nzero %in% true.zero)
  FN <- 0.5 * sum(zero %in%  true.nzero)

  recall    <- TP/(TP + FN)
  fallout   <- FP/(FP + TN)
  precision <- TP/(TP + FP)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  recall[TP + FN == 0] <- NA
  fallout[TN + FP == 0] <- NA
  precision[TP + FP == 0] <- NA

  res <-  round(c(fallout,recall,precision, f1_score), 3)
  res[is.nan(res)] <- 0
  names(res) <- c("fallout","recall", "precision", "f1_score")

  return(res)
}

#' @description computes the AUC given the list of recall and fallout values
#' @param recall list of recall values
#' @param fallout list of corresponding fallout values
auc <- function(recall, fallout){
  return(sum(diff(fallout) * (recall[-1] + recall[-length(recall)]) / 2))
}

#' @description given a precision matrix and a PLN model (collection of PLN fit
#' for different penalties), computes the AUC associated to the model
#' @param omega_true true precision matrix
#' @param PLN_model fitted PLN model (with multiple penalties)
get_auc <- function(omega_true, PLN_model){
  fallout <- c() ; recall <- c()
  for(pen in PLN_model$penalties){
    omega_estimate <- PLN_model$getModel(pen)$model_par$Omega
    res <- roc_metrics(omega_true, omega_estimate)
    if(!is.na(res[["fallout"]]) && !is.na(res[["recall"]])){
      fallout <- c(fallout, res[["fallout"]]) ; recall <- c(recall, res[["recall"]])
    }
  }
  if(pen == max(PLN_model$penalties)){recall <- rev(recall) ; fallout <- rev(fallout)}
  # One value of fallout may correspond to different recall values depending on the penalty
  fallout_unique <- unique(fallout) ; recall_unique <- c()
  for(x in fallout_unique){
    recall_unique <- c(recall_unique, max(recall[which(fallout == x)]))
  }
  return(auc(recall_unique, fallout_unique))
}

#' @description plots the ROC curve plot (True Positive Rate as a function
#' of the False Positive Rate) of a model with different penalties
#' @param omega_true true precision matrix
#' @param PLN_model fitted PLN model (with multiple penalties)
plot_roc_curve <- function(omega_true, PLN_model){
  fallout <- c() ; recall <- c()
  for(pen in PLN_model$penalties){
    omega_estimate <- PLN_model$getModel(pen)$model_par$Omega
    res <- roc_metrics(omega_true, omega_estimate)
    if(!is.na(res[["fallout"]]) && !is.na(res[["recall"]])){
      fallout <- c(fallout, res[["fallout"]]) ; recall <- c(recall, res[["recall"]])
    }
  }
  if(pen == max(PLN_model$penalties)){recall <- rev(recall) ; fallout <- rev(fallout)}
  # One value of fallout may correspond to different recall values depending on the penalty
  plot(recall, fallout)
}

#' @description computes a collection of measures associated to a given PLN model
#' @param PLN_model fitted PLN model (with multiple penalties)
#' @param params true parameters under which the data was simulated
#' @param model_selection model selection criterion to compute the measures for (BIC, ICL, StARS)
#' @param stability when model_selection = StARS, level of stability to use for stability selection
#' @param AUC AUC if already known, to avoid recomputing it if it's already been done with another model_selection criterion
get_measures <- function(PLN_model, params, model_selection = NULL,
                         stability = 0.8, AUC = NULL) {
  # Select best sparsity level according to the chosen criterion

  if(is.numeric(model_selection)){
    model <- PLN_model$getModel(model_selection)
  }else{
    if(model_selection == "StARS"){
      model <- PLN_model$getBestModel(model_selection, stability)
    }else{model <- PLN_model$getBestModel(model_selection)}
  }

  omega_hat <- model$model_par$Omega
  omega_rmse <- Metrics::rmse(omega_hat, params$Omega)
  ## get metrics
  if(!is.null(AUC)) AUC = get_auc(params$Omega, PLN_model)
  res <- c(
    criterion = model_selection,
    omega_rmse = round(omega_rmse, 2),
    AUC = get_auc(params$Omega, PLN_model),
    rmse_fit = rmse(model$fitted, params$Y),
    roc_metrics(params$Omega, omega_hat)
  )
  res
}
