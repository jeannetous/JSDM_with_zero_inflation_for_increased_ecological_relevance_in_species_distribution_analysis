n = 300 ; p = 20
min_X = 0;   max_X = 10; SNR = 0.75
omega_structure = "community"

# zi_type = "sites"
# n_mode_zi_proba = 3
# zi_mode_values =  c(0.2, 0.4, 0.6)
# proba_mode_zi = c(0.3, 0.3, 0.4)

zi_type = "species"
n_mode_zi_proba = 4
zi_mode_values = c(0, 0.3, 0.6, 0.9)
proba_mode_zi = c(0.4, 0.15, 0.15, 0.3)

simu_params = list(n = n,
                   p = p,
                   d = 1,
                   omega_structure = omega_structure,
                   zi_type = zi_type,
                   zi_covar_cluster = TRUE,
                   min_X = 0, max_X = 10, SNR = 0.75,
                   n_mode_zi_proba = n_mode_zi_proba,
                   zi_mode_values = zi_mode_values,
                   proba_mode_zi = proba_mode_zi,
                   block_values = NULL,
                   row_clusters = NULL,
                   col_clusters = NULL,
                   row_clusters_proba = NULL,
                   col_clusters_proba = NULL,
                   X0 = NULL, B0 = NULL,
                   min_X0 = 0, max_X0 = 10,
                   max_X0B0 = 0.2)

set.seed(1)

PLN_formula <- "Abundance ~ 0 + V1"
ZIPLN_formula <- "Abundance ~ 0 + V1" # | 0 + VZI1"
# PLN_formula_ZIvar <- "Abundance ~ 0 " #+ V1 + VZI1


plot_network = function(Omega,
                        type  = "partial_corr",
                        output = "igraph",
                        edge.color      = c("#F8766D", "#00BFC4"),
                        remove.isolated = FALSE,
                        node.labels     = NULL,
                        layout          = layout_in_circle,
                        plot = TRUE) {

  net <- - Omega / tcrossprod(sqrt(diag(Omega))); diag(net) <- 1
  colnames(net) <- unlist(lapply(1:ncol(net), f <- function(i) paste0("species_", i)))
  if (output == "igraph") {

    G <-  graph_from_adjacency_matrix(net, mode = "undirected", weighted = TRUE, diag = FALSE)

    if (!is.null(node.labels)) {
      igraph::V(G)$label <- node.labels
    } else {
      igraph::V(G)$label <- colnames(net)
    }
    ## Nice nodes
    V.deg <- degree(G)/sum(degree(G))
    igraph::V(G)$label.cex <- V.deg / max(V.deg) + .5
    igraph::V(G)$size <- V.deg * 100
    igraph::V(G)$label.color <- rgb(0, 0, .2, .8)
    igraph::V(G)$frame.color <- NA
    ## Nice edges
    igraph::E(G)$color <- ifelse(igraph::E(G)$weight > 0, edge.color[1], edge.color[2])
    if (type == "support"){igraph::E(G)$width <- abs(igraph::E(G)$weight)
    }else{igraph::E(G)$width <- 15*abs(igraph::E(G)$weight)}

    if (remove.isolated) {
      G <- delete.vertices(G, which(degree(G) == 0))
    }
    if (plot) plot(G, layout = layout)
  }
  if (output == "corrplot") {
    if (plot) {
      if (ncol(net) > 100)
        colnames(net) <- rownames(net) <- rep(" ", ncol(net))
      G <- net
      diag(net) <- 0
      corrplot(as.matrix(net), method = "color", is.corr = FALSE, tl.pos = "td", cl.pos = "n", tl.cex = 0.5, type = "upper")
    } else  {
      G <- net
    }
  }
  invisible(G)
}

get_partial_corr <- function(Omega){
  p <- nrow(Omega)
  net <- matrix(0, p, p)
  for(i in 1:(p -1)){
    for(j in ((i + 1):p)){
      net[i, j] <- - Omega[i, j] / sqrt(Omega[i, i] * Omega[j, j])
      net[j, i] <- net[i, j]
    }
  }
  diag(net) <- 1
  return(net)
}
