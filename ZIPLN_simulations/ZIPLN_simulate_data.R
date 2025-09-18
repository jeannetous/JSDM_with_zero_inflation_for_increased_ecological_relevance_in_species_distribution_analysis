source("ZIPLN_generate_simulation_parameters.R")
library(MASS)

# to document !!
generate_all_ZIPLN_parameters <- function(n, p, d, omega_structure = "erdos_renyi",
                                          zi_type = c("covar", "sites", "species"),
                                          zi_covar_cluster = FALSE,
                                          min_X = 0, max_X = 10, SNR = 0.75,
                                          v = 0.3, u = 0.1, n_mode_zi_proba = 2,
                                          zi_mode_values = NULL,
                                          proba_mode_zi = NULL,
                                          block_values = NULL,
                                          row_clusters = NULL,
                                          col_clusters = NULL,
                                          X0 = NULL, B0 = NULL,
                                          min_X0 = 0, max_X0 = 10,
                                          max_X0B0 = -0.2){
  Omega <- generate_omega(p, omega_structure, v, u)
  Sigma <- chol2inv(chol(Omega))
  X <- generate_X(n, d, min_X, max_X)
  B <- generate_B(p, X, Sigma, SNR)

  if(zi_type == "covar"){
    if(is.null(X0)){
      if(zi_covar_cluster){
        zi_params <- generate_X0_B0_cluster(n, p, block_values,
                                            row_clusters, col_clusters)
        B0 <- zi_params$B0 ; X0 <- zi_params$X0_num
        zi_params <- list(X0 = zi_params$X0, B0 = B0, X0_num = zi_params$X0_num)
      }else{
        X0 <- generate_X(n, d, min_X0, max_X0)
        colnames(X0) <- unlist(lapply(1:ncol(X0), f <- function(x) paste0("VZI", as.character(x))))
        B0 <- generate_B0(X0, max_X0B0)
        zi_params <- list(X0 = X0, B0 = B0)
      }
    }else{zi_params <- list(X0 = X0, B0 = B0)}
    }else{zi_params <- NULL}
  zi_proba <- generate_zi_proba(n, p, zi_type,
                                n_mode_zi_proba, zi_mode_values, proba_mode_zi,
                                X0, B0)
  return(list(Omega = Omega, Sigma = Sigma, X = X, B = B,
              zi_params = zi_params, zi_proba = zi_proba))
}

simulate_ZIPLN_data_fixed_parameters <- function(Sigma, X, B, zi_proba){
  n <- nrow(X) ; p <- nrow(Sigma)
  Y = matrix(rep(1, n*p), nrow=n)
  W = mvrnorm(n, mu=matrix(rep(0, p), p, 1), Sigma = Sigma)
  for(j in 1:p){
    for(i in 1:n){ Y[i, j] = rpois(1, exp(t(X[i,]) %*% B[,j] + W[i,j])[1,])}
  }
  Y <- add_zero_inflation(Y, zi_proba)
  return(Y)
}

simulate_ZIPLN_data <- function(params){
  Y <- simulate_ZIPLN_data_fixed_parameters(params$Sigma, params$X,
                                            params$B, params$zi_proba)
  return(Y)
}
