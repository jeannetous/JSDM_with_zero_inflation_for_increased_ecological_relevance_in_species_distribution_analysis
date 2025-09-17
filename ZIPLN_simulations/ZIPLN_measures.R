library(Metrics)

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

auc <- function(recall, fallout){
  return(sum(diff(fallout) * (recall[-1] + recall[-length(recall)]) / 2))
}

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
  ## get metrics
  if(!is.null(AUC)) AUC = get_auc(params$Omega, PLN_model)
  res <- c(
    criterion = model_selection,
    AUC = get_auc(params$Omega, PLN_model),
    rmse_fit = rmse(model$fitted, params$Y),
    roc_metrics(params$Omega, omega_hat)
  )
  res
}
