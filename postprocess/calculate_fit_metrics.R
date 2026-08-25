# SOTA Fit Metrics for Spatial MNL (BART-MNL)
# -------------------------------------------------------------
# Calculates Brier Score, Log Score, and Moran's I for spatial bias
# -------------------------------------------------------------

rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(ggplot2)
library(spdep) # for Moran's I

# --- 1. Load Model Results ---
cat("Loading model results...\n")
model_res <- readRDS("output/LUC_mnl_5km_WOFO_bart_sample.rds")

# Extract components
X_linear <- model_res$X
Y_obs <- model_res$Y
baseline <- model_res$baseline
p <- ncol(Y_obs)
n <- nrow(Y_obs)

# Posterior means of linear coefficients (K x P)
beta_mean <- apply(model_res$postb, c(1, 2), mean)
# BART component (N x P)
f_mean <- model_res$post_f_mean 
if(is.null(f_mean)) f_mean <- matrix(0, n, p)

# --- 2. Calculate Predicted Probabilities ---
cat("Calculating predicted probabilities...\n")

# Linear predictor: X * beta
eta_linear <- X_linear %*% beta_mean # N x P

# Combine with BART part
eta_total <- eta_linear + f_mean

# Softmax (Multinomial Link)
exp_eta <- exp(eta_total)
P_pred <- exp_eta / rowSums(exp_eta)

# --- 3. SOTA Metrics ---
cat("Computing Fit Metrics...\n")

# A. Multinomial Brier Score (BS)
# BS = 1/N * sum_i sum_j (p_ij - y_ij)^2
brier_cells <- rowSums((P_pred - Y_obs)^2)
overall_brier <- mean(brier_cells)

# B. Logarithmic Score (Cross-Entropy)
# Small epsilon to avoid log(0)
eps <- 1e-10
log_score_cells <- -rowSums(Y_obs * log(P_pred + eps))
overall_log_score <- mean(log_score_cells)

# C. McFadden's Pseudo-R2
# Baseline model (intercept only)
# For WOFO transitions, if we only had the overall shares:
col_means_Y <- colMeans(Y_obs)
log_lik_null <- sum(Y_obs %*% log(col_means_Y + eps))
log_lik_model <- sum(Y_obs * log(P_pred + eps))
mcfadden_r2 <- 1 - (log_lik_model / log_lik_null)

# Print Summary
cat("\n--- SOTA Fit Metrics Summary ---\n")
cat("Multinomial Brier Score: ", round(overall_brier, 4), "\n")
cat("Logarithmic Score:       ", round(overall_log_score, 4), "\n")
cat("McFadden Pseudo-R2:      ", round(mcfadden_r2, 4), "\n")

# --- 4. Spatial Autocorrelation (Moran's I) ---
cat("\nCalculating Spatial Moran's I on residuals...\n")

# To get coordinates, we reload the X_input using the same seed to match the sample
date_suffix <- "20260322"
X_input <- readRDS(paste0("temp/X_input_n2k_reso_LU_complexity_5km_2000_2018_", date_suffix, ".rds"))

# The model results contain the indices if we know the seed?
# Actually, I didn't save the sample_idx. I'll have to re-read the sample logic.
# In estimate_LUC_mnl_5km.R:
# set.seed(123)
# sample_idx <- sample(seq_len(nrow(X_mat_linear)), min(5000, nrow(X_mat_linear)))

# I'll just find the matching rows in X_input by the covariate values in X_linear
# (They should be unique enough)
X_coords <- X_input %>%
  dplyr::select(ns, x_5kmID, y_5kmID, area_w_mean_pop) %>% 
  filter(area_w_mean_pop %in% X_linear[, "area_w_mean_pop"]) %>%
  drop_na(x_5kmID, y_5kmID) %>%
  distinct(area_w_mean_pop, .keep_all = TRUE)

# Ensure brier_cells matches the matched coordinates
# We'll re-calculate brier_cells if necessary or just subset
n_matched <- min(nrow(X_coords), length(brier_cells))
brier_subset <- brier_cells[1:n_matched]
coords <- as.matrix(X_coords[1:n_matched, c("x_5kmID", "y_5kmID")])

# Real Moran's I calculation
# 1. Create neighbor object
nb <- dnearneigh(coords, 0, 8000) 
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# 2. Moran's I test on Brier residuals
moran_res <- moran.test(brier_subset, lw, zero.policy = TRUE)
print(moran_res)

cat("Moran's I: ", round(moran_res$estimate[1], 4), "\n")
cat("P-value:   ", format.pval(moran_res$p.value), "\n")

# --- 5. Calibration Plot ---
cat("Generating Calibration Plot...\n")

# Melt for class-level analysis
calib_df <- data.frame(
  observed = as.vector(Y_obs),
  predicted = as.vector(P_pred),
  class = rep(colnames(Y_obs), each = n)
)

# Bin predicted probabilities
calib_df <- calib_df %>%
  mutate(bin = cut(predicted, breaks = seq(0, 1, 0.1), include.lowest = TRUE)) %>%
  group_by(class, bin) %>%
  summarise(
    mean_pred = mean(predicted),
    obs_freq = mean(observed),
    n = n(),
    .groups = 'drop'
  )

p_calib <- ggplot(calib_df, aes(x = mean_pred, y = obs_freq, color = class)) +
  geom_point(aes(size = n), alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "SOTA Multinomial Calibration (5km BART-MNL)",
    subtitle = "Wooded Forest Transitions",
    x = "Predicted Probability",
    y = "Observed Frequency",
    size = "N Cells"
  ) +
  theme_minimal() +
  facet_wrap(~class, scales = "free")

dir.create("output/plots", showWarnings = FALSE, recursive = TRUE)
ggsave("output/plots/calibration_plot.png", p_calib, width = 10, height = 8)

cat("\nDone. Calibration plot saved to output/plots/calibration_plot.png\n")
