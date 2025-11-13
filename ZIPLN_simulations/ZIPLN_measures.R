library(Metrics)

#' @description computes recall, fallout, precision and f1-score for
#' an inferred network (precision matrices here) given the true one
#' @param omega_true true precision matrix
#' @param omega_estimate estimated precision matrix
roc_metrics <- function(omega_true, omega_estimate) {

  diag(omega_true) <- 0 ; p <- nrow(omega_true)
  roc <- function(theta) {
    diag(theta) <- 0

    nzero <- which(theta != 0)
    zero  <- which(theta == 0)

    true.nzero <- which(omega_true != 0)
    true.zero  <- which(omega_true == 0)

    TP <- sum(nzero %in% true.nzero)
    TN <- sum(zero %in%  true.zero) - p
    FP <- sum(nzero %in% true.zero)
    FN <- sum(zero %in%  true.nzero)
    recall    <- TP/(TP + FN) ## also recall and sensitivity
    fallout   <- FP/(FP + TN) ## also 1 - specificit
    precision <- TP/(TP + FP) ## also PPR
    f1_score  <- 2 * (precision * recall) / (precision + recall)
    recall[TP + FN == 0] <- NA
    fallout[TN + FP == 0] <- NA
    precision[TP + FP == 0] <- NA

    res <-  round(c(fallout,recall,precision, f1_score),3)
    res[is.nan(res)] <- 0
    names(res) <- c("fallout","recall", "precision", "f1_score")
    res
  }

  if (is.list(omega_estimate)) {
    return(as.data.frame(do.call(rbind, lapply(omega_estimate, roc))))
  } else {
    return(roc(omega_estimate))
  }
}

#' @description computes AUC from the roc measures
#' @param roc measures as computed by function roc_metrics, named list that contains
#' a list of fallout and a list of recall values
#' @param threshold from which the list of recall / fallout values should be cut
#' before only adding a (1, 1) point to the ROC curve
perf_auc <- function(roc, threshold = 1) {
  cut <- (roc$fallout < threshold) & (roc$recall < threshold)
  fallout <- c(0, roc$fallout[cut], threshold)
  recall  <- c(0, roc$recall[cut] , threshold)
  dx <- diff(fallout)
  res <- sum(c(recall[-1]*dx, recall[-length(recall)]*dx))/2
  res <- ifelse(is.character(res), NA, res)
  res
}


#' @description given a precision matrix and a PLN model (collection of PLN fit
#' for different penalties), computes the AUC associated to the model
#' @param omega_true true precision matrix
#' @param PLN_model fitted PLN model (with multiple penalties)
get_auc <- function(omega_true, PLN_model){
  roc <- roc_metrics(omega_true,
                     lapply(PLN_model$models, function(model) model$model_par$Omega))

  return(perf_auc(roc))
}

#' @description plots the ROC curve plot (True Positive Rate as a function
#' of the False Positive Rate) of a model with different penalties
#' @param omega_true true precision matrix
#' @param PLN_model fitted PLN model (with multiple penalties)
plot_roc_curve <- function(omega_true, PLN_model){
  roc <- roc_metrics(omega_true,
                     lapply(PLN_model$models, function(model) model$model_par$Omega))
  plot(roc$recall, roc$fallout)
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
