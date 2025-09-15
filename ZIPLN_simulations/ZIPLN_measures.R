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

  recall    <- TP/(TP + FN) ## also recall and sensitivity
  fallout   <- FP/(FP + TN) ## also 1 - specificit
  precision <- TP/(TP + FP) ## also PPR
  recall[TP + FN == 0] <- NA
  fallout[TN + FP == 0] <- NA
  precision[TP + FP == 0] <- NA

  res <-  round(c(fallout,recall,precision), 3)
  res[is.nan(res)] <- 0
  names(res) <- c("fallout","recall", "precision")

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

get_measures <- function(PLN_model, param, model_selection = NULL,
                         stability = 0.8, fixed_blocks = FALSE) {
  # Select best sparsity level according to the chosen criterion

  if(is.numeric(model_selection)){
    model <- PLN_model$getModel(model_selection)
  }else{
    if(model_selection == "StARS"){
      model <- PLN_model$getBestModel(model_selection, stability)
    }else{model <- PLN_model$getBestModel(model_selection)}
  }

  # Get best permutation of Omega according to rmse when possible
  omega_hat <- model$model_par$Omega

  ## get metrics
  res <- c(
    criterion = model_selection,
    fixed_blocks = fixed_blocks,
    AUC = get_auc(param$Omega, PLN_model),
    rmse_fit = rmse(model$fitted, param$Y),
    roc_metrics(param$Omega, omega_hat, best_perm)
  )
  res
}
