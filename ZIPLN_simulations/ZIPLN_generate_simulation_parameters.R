library(igraph)

####################### Functions to generate Omega ############################

# Erdos-Reyni
erdos_reyni_graph <- function(p, prob = 0.5){
  as_adjacency_matrix(sample_gnp(p, prob))
}

# Preferential attachment
preferential_attachment_graph <- function(p){
  as_adjacency_matrix(sample_pa(p, directed = FALSE))
}

# Community structure
community_graph <- function(p, prob = c(1/2,1/4,1/4), prob_in = 0.5, prob_out = 0.1) {
  pref_mat <- matrix(prob_out, length(prob), length(prob))
  diag(pref_mat) <- prob_in
  graph_mat <- as_adjacency_matrix(sample_sbm(p,
                                              pref.matrix = pref_mat,
                                              block.sizes = c(rmultinom(1, p, prob)) ))
  graph_mat
}


generate_omega <- function(p, omega_structure, v = 0.3, u = 0.1){
  cond <- FALSE
  while(!cond){
    if(omega_structure == "erdos_reyni") G <- erdos_reyni_graph(p)
    if(omega_structure == "preferential_attachment") G <- preferential_attachment_graph(p)
    if(omega_structure == "community") G <- community_graph(p)

    # Ensuring that the network is not empty for AUC to make sense
    if(max(G) == 0){
      off_diag_indices <- which(row(matrix(1:p, p, p)) != col(matrix(1:p, p, p)), arr.ind = TRUE)
      selected_index <- off_diag_indices[sample(nrow(off_diag_indices), 1), ]
      G[selected_index[["row"]], selected_index[["col"]]] <- 1
      G[selected_index[["col"]], selected_index[["row"]]] <- 1
    }
    omega_tilde <- G * v
    omega <- omega_tilde + diag(abs(min(eigen(omega_tilde)$values)) + u, p, p)
    # Ensuring that the network is not full for AUC to make sense
    if(min(omega) > 0){ # Ensuring that the network has 0s for AUC to make sense
      off_diag_indices <- which(row(matrix(1:p, p, p)) != col(matrix(1:p, p, p)), arr.ind = TRUE)
      selected_index <- off_diag_indices[sample(nrow(off_diag_indices), 1), ]
      omega[selected_index[["row"]], selected_index[["col"]]] <- 0
      omega[selected_index[["col"]], selected_index[["row"]]] <- 0
    }
    cond <- all(eigen(omega)$values > 0)
  }
  as.matrix(omega)
}


################## Functions to add zero-inflation in the data #################
# %% max_X0B0 controls for the column mean
generate_B0_simple <- function(X0, max_X0B0 = -0.2){
  d <- ncol(X0)
  B0 <- matrix(rep(1, d*p), nrow=d)
  for(dim in 1:d){B0[dim,] = runif(p, min=-1, max = 1)}
  correcting_factors <- unlist(lapply(colMeans(X0 %*% B0),
                                      f <- function(x){ifelse(x <= max_X0B0, 1,
                                                              max_X0B0 / x)}))
  B0 <- sweep(B0, 2, correcting_factors, `*`)
}

# %% C gives the expected block values of X0 %*% B0
# %% X0 and B0 are divided in the clusters given by row_clusters and col_clusters
generate_X0_B0_cluster <- function(n, p, block_values,
                                row_clusters = NULL, col_clusters = NULL) {

  a <- nrow(block_values) ; b <- ncol(block_values)
  X0 <- generate_discrete_X(n, 1, a)
  X0_num <- model.matrix(~ . - 1, data = as.data.frame(X0))
  if(is.null(row_clusters)) row_clusters <- sort(rep(1:a, length.out = n))
  if(is.null(col_clusters)) col_clusters <- sort(rep(1:b, length.out = p))

  M <- block_values[row_clusters, col_clusters]
  M <- apply(M, c(1,2), f <- function(x){rnorm(1, x, 0.05)})

  # Solve X0 %*% B0 = M for B0 using least squares
  B0 <- solve(t(X0_num) %*% X0_num, t(X0_num) %*% M)
  return(list("X0" = X0, "B0" = B0))
}

generate_zi_proba <- function(n, p, n_mode_proba = c(1), zi_mode_values = NULL,
                              zi_type, X0 = NULL, B0 = NULL){
  if(zi-type == "covar"){
    X0B0 <- X0 %*% B0
    zi_proba <- exp(X0B0) / (1 + exp(X0B0))
    zi_proba <- apply(zi_proba, c(1, 2), f <- function(x) min(1, max(0, x)))
  }else{
    breaks <- cumsum(n_mode_proba)
    if(zi-type == "sites"){groups <- cut(1:n, c(0, round(breaks * n)), labels = FALSE)}
    if(zi-type == "species"){groups <- cut(1:p, c(0, round(breaks * p)), labels = FALSE)}
    zi_proba_list <- unlist(lapply(1:length(n_mode_proba),
                                   f <- function(i){unlist(lapply(rnorm(table(groups)[[i]],
                                                                        mean = zi_mode_values, sd = 0.05),
                                                                  f <- function(x){return(min(1, max(x, 0)))}))}))
    if(zi-type == "sites"){zi_proba <- matrix(rep(zi_proba_list, p), nrow = n, byrow = F)}
    if(zi-type == "species"){zi_proba <- matrix(rep(zi_proba_list, n), nrow = n, byrow = T)}
  }
  return(zi_proba)
}

add_zero_inflation <- function(Y, zi_proba){
  Z <- apply(zi_proba, c(1, 2), f <- function(x) rbinom(1,1,x))
  Y[Z == 1] <- 0
  return(Y)
}

##################### Functions to generate other parameters ###################
generate_D <- function(p, min_D = 0.5, max_D = 1.5){
  D <- matrix(rep(0, p*p), nrow = p)
  diag(D) <- runif(p, min_D, max_D)
  return(D)
}

generate_X <- function(n, d, min_X = 0, max_X = 10){
  if(length(min_X == 1)) min_X <- rep(min_X, d)
  if(length(max_X == 1)) max_X <- rep(max_X, d)
  X = matrix(rep(1, n * d), nrow=n)
  for(dim in 1:d){X[,dim] = runif(n, min=min_X[[dim]], max = max_X[[dim]])}
  return(X)
}

generate_discrete_X <- function(n, d, n_cat_values){
  X = matrix(rep(1, n * d), nrow=n)
  for(dim in 1:d){X[,dim] = sort(rep(1:n_cat_values[[dim]], length.out = n))}
  X <- apply(X, c(1,2), as.character) # for X values to be treated as categorical variables
  return(X)
}

generate_B <- function(p, X, Sigma, SNR = 0.75){
  d <- ncol(X)
  B <- matrix(rep(1, d*p), nrow=d)
  for(dim in 1:d){B[dim,] = runif(p, min=0, max = 1)}
  correcting_factor <- SNR * var(as.vector(Sigma)) / (var(as.vector(X %*% B)))
  B <- sqrt(correcting_factor) * B
  return(B)
}

generate_parameters <- function(X, p, omega_structure) {
  Omega <- generate_omega(p, omega_structure)
  Sigma <- chol2inv(chol(Omega))
  list(
    B = generate_B(p, X, Sigma),
    D = generate_D(p),
    Omega = Omega,
    Sigma = Sigma
  )
}
