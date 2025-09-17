library(dplyr)
library(tidyr)

# Example input data
df <- data.frame(
  A = c("a1", "a2", "a3"),
  B = c("b1", "b2", "b3"),
  C = c(0, 1, 2)
)
df <- settings

# Example lists
E_list  <- c("e1", "e2")
F_list  <- c("f1", "f2")

E2_list <- c("e3", "e4", "e5")
F2_list <- c("f3", "f4", "f5")

# Build lookup tables
lookup1 <- data.frame(E = E_list, F = F_list)
lookup2 <- data.frame(E = E2_list, F = F2_list)

# Expand df
df_out <- df %>%
  rowwise() %>%
  do({
    row <- .
    if (row$zi_type == "covar") {
      # Keep row with NA for E and F
      tibble(n = row$n, p = row$p, omega_structure = row$omega_structure,zi_type = row$zi_type, E = NA, F = NA)
    } else if (row$zi_type == "species") {
      # Expand with lookup1
      cbind(row[1:4], lookup1)
    } else if (row$zi_type == "sites") {
      # Expand with lookup2
      cbind(row[1:4], lookup2)
    }
  }) %>%
  ungroup()
