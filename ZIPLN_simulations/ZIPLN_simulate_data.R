source("ZIPLN_generate_simulation_parameters.R")
library(MASS)

#' @description generates a named list with all the required parameters to simulate
#' data under the ZIPLN model
#' @param n number of rows in the Abundance matrix
#' @param p number of columns in the Abundance matrix
#' @param d number of covariates in the covariates matrix
#' @param omega_structure network structure for the precision matrix (erdos_renyi,
#' community or preferential_attachment)
#' @param zi_type type of zero-inflation
#' @param zi_covar_cluster boolean, if zi_type = covar, whether there should be
#' a division of the ZI values into rows and column clusters
#' @param min_X minimum value for X, either one single value for X, or a list of length d for each dimension
#' @param max_X maximum value for X, either one single value for X, or a list of length d for each dimension
#' @param SNR signal to noise ratio, ratio between Sigma's variance and that of XB
#' @param v calibration parameter to get Omega from a graph
#' @param u calibration parameter to get Omega from a graph
#' @param n_mode_zi_proba if zi_type = sites or species, number of different zi probabilities
#' @param zi_mode_values if zi_type = sites or species, list of zi probabilities of length n_mode_zi_proba
#' @param proba_mode_zi list of probabilities of having each ZI contained in zi_mode_values
#' @param block_values if zi_type = "covar" and  zi_covar_cluster = TRUE, values
#' of X0 %*% B0 expected for each pair (row_cluster, col_cluster)
#' @param X0 optional, if zi_type = covar, list of ZI covariates
#' @param B0 optional, regression matrix for the ZI covariates, required if X0 is not null
#' @param min_X0 minimum value for X0, either one single value for X0, or a list of length d for each dimension,
#' applied only if zi_covar_cluster = FALSE (otherwise X0 is discrete)
#' @param max_X0 maximum value for X0, either one single value for X0, or a list of length d for each dimension
#' applied only if zi_covar_cluster = FALSE (otherwise X0 is discrete)
#' @param max_X0B0 maximum value for the mean of each column of X0 %*% B0
#' applied only if zi_covar_cluster = FALSE (otherwise X0 is discrete)
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
                                          row_clusters_proba = NULL,
                                          col_clusters_proba = NULL,
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
                                            row_clusters, col_clusters,
                                            row_clusters_proba, col_clusters_proba)
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

#' @description simulates data under the ZIPLN model for fixed parameters
#' @param Sigma variance-covariance matrix of the model
#' @param X covariates matrix
#' @param B regression coefficient matrix
#' @param zi_proba matrix of zero-inflation probabilities
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

#' @description simulates data under the ZIPLN model for a list of fixed parameters
#' @param params named list of parameters for the ZIPLN model
simulate_ZIPLN_data <- function(params){
  Y <- simulate_ZIPLN_data_fixed_parameters(params$Sigma, params$X,
                                            params$B, params$zi_proba)
  return(Y)
}
