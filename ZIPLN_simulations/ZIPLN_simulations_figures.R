library(ggplot2)
library(viridis)
library(tidyverse)

res1 <- read.csv("ZIPLN_simulations_res/ZIPLN_simu_ref_BIC_2.csv")

res <- res1[!grepl("^Error", res1$simu), ]
# res <- read.csv("ZIPLN_simus_test_3.csv")
res$AUC <- as.numeric(res$AUC) #; res <- res[!is.na(res$AUC),]
res$f1_score <- as.numeric(res$f1_score)
res$omega_rmse <- as.numeric(res$omega_rmse)
res$precision <- as.numeric(res$precision)
# res_auc <- res[res$criterion == "BIC",]
res$zi_strength <- as.numeric(sub(".*_", "", res$zi_config))
res <- res[res$n %in% c(100, 300) & res$method != "neighborhood_selection_network",]

################################ Figures zi_type ~ strength ####################
auc_zi_type_strength <- ggplot(res, aes(x = zi_strength, y = AUC, fill = method)) +
  geom_violin() +
  facet_grid(zi_strength ~ zi_type) +
  labs(title = "AUC (zi type ~ zi strength)", x = "", y = "AUC") +
  theme_minimal() +  # Apply a clean theme
  scale_fill_viridis_d() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  )

rmse_omega_zi_type_strength <- ggplot(res, aes(x = zi_strength, y = omega_rmse, fill = method)) +
  geom_violin() +
  facet_grid(zi_strength ~ zi_type) +
  labs(title = "Omega RMSE (zi type ~ zi strength)", x = "", y = "Omega RMSE") +
  theme_minimal() +  # Apply a clean theme
  scale_fill_viridis_d() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  ) +
  ylim(0, 1) 

f1score_zi_type_strength <- ggplot(res, aes(x = zi_strength, y = f1_score, fill = method)) +
  geom_violin() +
  facet_grid(zi_strength ~ zi_type) +
  labs(title = "F1 score (zi type ~ zi strength)", x = "", y = "F1 score") +
  theme_minimal() +  # Apply a clean theme
  scale_fill_viridis_d() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  )

precision_zi_type_strength <- ggplot(res, aes(x = zi_strength, y = precision, fill = method)) +
  geom_violin() +
  facet_grid(zi_strength ~ zi_type) +
  labs(title = "Precision (zi type ~ zi strength)", x = "", y = "Precision") +
  theme_minimal() +  # Apply a clean theme
  scale_fill_viridis_d() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  )

fit_zi_type_strength <- ggplot(res, aes(x = zi_strength, y = rmse_fit, fill = method)) +
  geom_violin() +
  facet_grid(zi_strength ~ zi_type) +
  labs(title = "RMSE fit (zi type ~ zi strength)", x = "", y = "RMSE Fit") +
  theme_minimal() +  # Apply a clean theme
  scale_fill_viridis_d() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  )


################################ Figures n ~ p #################################
auc_n_p <- ggplot(res, aes(x = n, y = AUC, fill = method)) +
  geom_violin() +
  facet_grid(n ~ p) +
  labs(title = "AUC (n ~ p)", x = "", y = "AUC") +
  theme_minimal() +  # Apply a clean theme
  scale_fill_viridis_d() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  )

rmse_omega_n_p <- ggplot(res, aes(x = n, y = omega_rmse, fill = method)) +
  geom_violin() +
  facet_grid(n ~ p) +
  labs(title = "Omega RMSE (n ~ p)", x = "", y = "Omega RMSE") +
  theme_minimal() +  # Apply a clean theme
  scale_fill_viridis_d() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  ) +
  ylim(0, 1) 
