library(MASS)
# library(PLNmodels)
library(pheatmap)
library(paletteer)

# Function to get partial correlations from precision matrix
partial_correlations <- function(Omega){
  p <- nrow(Omega)
  partial_cor <- matrix(rep(0, p*p), nrow = p)
  diag(partial_cor) <- 1
  for(i in 1:(p-1)){
    for(j in (i + 1):p){
      partial_cor[i, j] <- -Omega[i, j] / sqrt(Omega[i, i] * Omega[j, j])
    }
  }
  partial_cor[lower.tri(partial_cor)] <- t(partial_cor)[lower.tri(partial_cor)]
  return(partial_cor)
}

# Generating simulated data
set.seed(1)
n = 500
p = 3

X <- matrix(0.01 * as.vector(-200:299), nrow = n)
B <- matrix(as.vector(c(2, 2, 1.7)), nrow = 1)

R <- matrix(c(1, -0.05, -0.1,
              -0.05, 1, 0.9,
              -0.1, 0.9, 1), nrow = 3, byrow = TRUE)
variances <- c(0.5, 0.5, 0.5)
D <- diag(sqrt(variances))
Sigma <- D %*% R %*% D
Omega <- solve(Sigma)

Y = matrix(rep(1, n*p), nrow=n)
Z = mvrnorm(n, mu=matrix(rep(0, p), p, 1), Sigma=Sigma)
for(j in 1:p){
  for(i in 1:n){ Y[i, j] = rpois(1, exp(t(X[i,]) %*% B[,j] + Z[i,j])[1,])}
}

Yres = matrix(rep(1, n*p), nrow=n)
for(j in 1:p){
  for(i in 1:n){ Yres[i, j] =  rpois(1, exp( Z[i,j]))}
}

colnames(Y) <- c("A", "B", "C")
simu_data <- prepare_data(Y, X)

# Visualising the abundances
log_Y <- log(1 + Y)
color_palette <- paletteer_c("ggthemes::Orange-Gold", n = 30)
palette_colors <- paletteer::paletteer_c("ggthemes::Orange-Gold", n = 29)  # 29 instead of 30
color_palette <- c("white", palette_colors)
max_val <- max(log_Y, na.rm = TRUE)

breaks <- seq(0, max_val, length.out = length(color_palette) + 1)
pheatmap(
  log_Y,
  color = color_palette,
  fontsize_col = 8,
  angle_col = 45,
  show_rownames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  legend = TRUE,
  width = 7,
  height = 6
)


# Running regular PLN model
myPLN <- PLNnetwork(Abundance ~ 0 + V1, simu_data, penalties = c(1e-9))
myPLN$getModel(0)$plot_network()

# Adding zero-inflation

## Species-dependent
pi = matrix(c(rep(0.7, n), rep(0.3, n), rep(0.7, n)), ncol = 3)
Y_zi <- Y
Y_zi[matrix(data = rbinom(n * p, size = 1, prob = pi),
         nrow = n,ncol = p) == 1] <- 0


log_Y_zi <- log(1 + Y_zi)
color_palette <- paletteer_c("ggthemes::Orange-Gold", n = 30)
palette_colors <- paletteer::paletteer_c("ggthemes::Orange-Gold", n = 29)  # 29 instead of 30
color_palette <- c("white", palette_colors)
max_val <- max(log_Y_zi, na.rm = TRUE)
pheatmap(
  log_Y_zi,
  color = color_palette,
  fontsize_col = 8,
  angle_col = 45,
  show_rownames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  legend = TRUE,
  width = 7,
  height = 6
)

simu_data_zi <- prepare_data(Y_zi, X)
myPLN_zi <- PLNnetwork(Abundance ~ 0 + V1, simu_data_zi, penalties = c(1e-9))
myPLN_zi$getModel(0)$plot_network()
myZIPLN_zi <- ZIPLNnetwork(Abundance ~ 0 + V1 | 1, simu_data_zi, penalties = c(1e-9))
myZIPLN_zi$getModel(0)$plot_network()

## Site-dependent
pi = pi <- matrix(rep(runif(n, 0, 0.8), 3), ncol = 3)
Y_zi <- Y
Y_zi[matrix(data = rbinom(n * p, size = 1, prob = pi),
            nrow = n,ncol = p) == 1] <- 0

log_Y_zi <- log(1 + Y_zi)
color_palette <- paletteer_c("ggthemes::Orange-Gold", n = 30)
palette_colors <- paletteer::paletteer_c("ggthemes::Orange-Gold", n = 29)  # 29 instead of 30
color_palette <- c("white", palette_colors)
max_val <- max(log_Y_zi, na.rm = TRUE)
pheatmap(
  log_Y_zi,
  color = color_palette,
  fontsize_col = 8,
  angle_col = 45,
  show_rownames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  legend = TRUE,
  width = 7,
  height = 6
)

simu_data_zi <- prepare_data(Y_zi, X)
myPLN_zi <- PLNnetwork(Abundance ~ 0 + V1, simu_data_zi, penalties = c(1e-9))
myPLN_zi$getModel(0)$plot_network()
myZIPLN_zi <- ZIPLNnetwork(Abundance ~ 0 + V1, simu_data_zi, zi = "row", penalties = c(1e-9))
myZIPLN_zi$getModel(0)$plot_network()
# Covariate dependent

generate_pi <- function(p, X0 = NA, max_X0B0 = -0.2){
  d <- ncol(X0)
  B0 <- matrix(rep(1, d*p), nrow=d)
  for(dim in 1:d){B0[dim,] = runif(p, min=-1, max = 1)}
  correcting_factors <- unlist(lapply(colMeans(X0 %*% B0),
                                      f <- function(x){ifelse(x <= max_X0B0, 1,
                                                              max_X0B0 / x)}))
  B0 <- sweep(B0, 2, correcting_factors, `*`)
  X0B0 <- X0 %*% B0
  pi <- exp(X0B0) / (1 + exp(X0B0))
  return(list("pi"= pi, "B0" = B0))
}

pi_param <- generate_pi(p, X0 = X, 0.2)
pi <- pi_param$pi
Y_zi <- Y
Y_zi[matrix(data = rbinom(n * p, size = 1, prob = pi),
            nrow = n,ncol = p) == 1] <- 0

max_val <- max(log_Y_zi, na.rm = TRUE)
pheatmap(
  log(1 + Y_zi),
  color = color_palette,
  fontsize_col = 8,
  angle_col = 45,
  show_rownames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  legend = TRUE,
  width = 7,
  height = 6
)

simu_data_zi <- prepare_data(Y_zi, X)
myPLN_covar_zi <- PLNnetwork(Abundance ~ 0 + V1, simu_data_zi, penalties = c(1e-9))
myPLN_covar_zi$getModel(0)$plot_network()
myZIPLN_covar_zi <- ZIPLNnetwork(Abundance ~ 0 + V1 | V1, simu_data_zi, penalties = c(1e-9))
myZIPLN_covar_zi$getModel(0)$plot_network()
