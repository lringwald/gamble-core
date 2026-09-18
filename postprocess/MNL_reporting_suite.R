# =============================================================================
# MNL_reporting_suite.R
# 
# A comprehensive suite for visualizing and reporting results from Hybrid
# Hierarchical + BART Multinomial Logit models.
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(patchwork)
library(dbarts)
library(terra)
library(abind)
library(scales)

# --- 1. CORE PREDICTION HELPER ---

#' Predict Shares from Hybrid MNL Model
#' @param res Combined result list from mnlogit_rcpp
#' @param X_linear Matrix of linear covariates
#' @param X_bart Matrix of BART covariates
#' @param n_samples Number of posterior samples to use for prediction (speed/memory trade-off)
predict_mnl_hybrid <- function(res, X_linear, X_bart, n_samples = 100) {
  # Check if model has an intercept and X_linear is missing it
  k_expected <- dim(res$postb_pooled)[1]
  if (!is.null(X_linear)) {
    X_linear <- as.matrix(X_linear)
    if (ncol(X_linear) == k_expected - 1) {
      X_linear <- cbind(X_linear, intercept = 1)
    }
  }
  
  n_obs <- if (!is.null(X_linear)) nrow(X_linear) else nrow(X_bart)
  k <- k_expected
  p <- dim(res$postb_pooled)[2]
  baseline <- res$baseline
  X_bart <- as.matrix(X_bart)
  pp <- (1:(p+1))[-baseline]

  
  # Thin samples if requested
  p_all <- p + 1
  
  # Average over n_samples
  U_sum <- matrix(0, n_obs, p_all)
  n_available <- if (!is.null(res$tree_store)) length(res$tree_store) else dim(res$postb_pooled)[3]
  n_to_sample <- min(n_samples, n_available)
  
  if (n_to_sample == 0) stop("No posterior samples found in model object.")
  
  s_indices <- sample(seq_len(n_available), n_to_sample)
    
  for (iter in 1:n_to_sample) {
    s <- s_indices[iter]
    U <- matrix(0, n_obs, p_all)
    
    # 1. Linear Component (EU-mean used for global prediction suite)
    U[, pp] <- X_linear %*% res$postb_pooled[, , s]
    
    # 2. BART Component
    if (!is.null(res$tree_store)) {
      state_m <- res$tree_store[[s]]
      if (is.character(state_m)) {
        # DISK STORAGE (Batched)
        batch_size <- if (!is.null(res$bart_batch_size)) res$bart_batch_size else 50
        batch_idx <- ceiling(s / batch_size)
        sample_in_batch <- (s - 1) %% batch_size + 1
        curr_batch <- qs2::qs_read(res$tree_store[batch_idx])
        state_m <- curr_batch[[sample_in_batch]]
      }
      
      if (!is.null(state_m)) {
        for (ip in seq_along(pp)) {
          if (is.data.frame(state_m[[ip]])) {
            # SLIM TREE STORAGE (Optimized for disk space)
            # We use the C++ predictor for speed
            bart_preds <- predict_slim_bart_cpp(as.matrix(X_bart), state_m[[ip]])
            if (iter == 1 && ip == 1) {
            }
            U[, pp[ip]] <- U[, pp[ip]] + bart_preds
          } else {
            # FULL SERIALIZED STORAGE (Includes full sampler state)
            sampler <- unserialize(state_m[[ip]])
            U[, pp[ip]] <- U[, pp[ip]] + sampler$predict(as.matrix(X_bart))
          }
        }
      }
    } else if (!is.null(res$post_f_mean)) {
      # Fallback to mean (only if dimensions match, e.g. spatial maps of training data)
      if (nrow(res$post_f_mean) == n_obs) {
        U[, pp] <- U[, pp] + res$post_f_mean
      }
    }
    
    # Softmax
    exp_U <- exp(U)
    shares <- exp_U / rowSums(exp_U)
    U_sum <- U_sum + shares
  }
    
  return(U_sum / n_to_sample)
}

# --- 2. CONVERGENCE DIAGNOSTICS ---

plot_mnl_convergence_summary <- function(res_list, name, output_dir) {
  # This expects the list of chains before combination
  # We reuse the logic from assess_convergence in the main script but make it prettier
  
  # [Implementation node: assess_convergence already does a good job, 
  # this function would consolidate them into a report page]
  cat("Consolidating convergence diagnostics...\n")
}

# --- 3. PARAMETRIC EFFECTS ---

plot_mnl_parametric_report <- function(res_combined, cov_names, cat_names, output_dir) {
  source(file.path(if (dir.exists("postprocess")) "postprocess" else ".", "MNL_parameter_heatplot.R"))
  source(file.path(if (dir.exists("postprocess")) "postprocess" else ".", "MNL_re_viz_utils.R"))
  
  # 1. Heatplot
  p_heat <- MNL_parameter_heatplot(
    draws = res_combined$postb_pooled,
    cov_names = cov_names,
    cat_names = cat_names,
    exclude_intercept = TRUE,
    title = "Hybrid Model: Parametric (EU-Mean) Effects"
  )
  ggsave(file.path(output_dir, "parametric_heatplot.png"), p_heat, width = 12, height = 10)
  
  # 2. RE Profiles (if hierarchical)
  if (length(dim(res_combined$postb)) == 4) {
    p_re <- plot_mnl_hierarchical_profile(
      res_combined = res_combined,
      re_idx = 1:dim(res_combined$postb)[1],
      cov_names = cov_names,
      cat_names = cat_names,
      group_names = res_combined$nuts0_names %||% as.character(1:dim(res_combined$postb)[3]),
      title = "Hierarchical Response Profiles (NUTS0)"
    )
    ggsave(file.path(output_dir, "hierarchical_profiles.png"), p_re, width = 14, height = 10)
  }
}

# --- 4. BART INTERPRETABILITY ---

plot_mnl_bart_pdp <- function(res, X_linear_ref, X_bart, var_name, 
                              cat_names, output_dir, n_points = 50) {
  cat("Generating PDP for:", var_name, "...\n")
  
  target_idx <- which(colnames(X_bart) == var_name)
  if (length(target_idx) == 0) return(NULL)
  
  v_range <- seq(min(X_bart[, target_idx]), max(X_bart[, target_idx]), length.out = n_points)
  
  # Reference point: mean of all other variables
  X_bart_ref <- matrix(colMeans(X_bart), n_points, ncol(X_bart), byrow = TRUE)
  colnames(X_bart_ref) <- colnames(X_bart)
  X_bart_ref[, target_idx] <- v_range
  
  X_linear_ref_mat <- matrix(colMeans(X_linear_ref), n_points, ncol(X_linear_ref), byrow = TRUE)
  
  # Predict shares across range
  pred_shares <- predict_mnl_hybrid(res, X_linear_ref_mat, X_bart_ref, n_samples = 50)
  
  df_pdp <- as.data.frame(pred_shares)
  colnames(df_pdp) <- cat_names
  df_pdp$Value <- v_range
  
  df_pdp_long <- df_pdp %>%
    pivot_longer(cols = -Value, names_to = "Category", values_to = "Share")
  
  p <- ggplot(df_pdp_long, aes(x = Value, y = Share, color = Category)) +
    geom_line(linewidth = 1.2) +
    theme_minimal() +
    labs(title = paste("Partial Dependence Plot:", var_name),
         subtitle = "Effect of driver on predicted land-use shares (BART Mixed Model)",
         x = var_name, y = "Predicted Share") +
    scale_y_continuous(labels = percent)
  
  ggsave(file.path(output_dir, paste0("pdp_", var_name, ".png")), p, width = 8, height = 6)
  return(p)
}

# --- 5. SPATIAL VISUALIZATION ---

plot_mnl_spatial_maps <- function(pred_shares, coords, cat_names, output_dir, res_km = 10) {
  cat("Generating spatial maps...\n")
  
  df_maps <- cbind(coords, as.data.frame(pred_shares))
  colnames(df_maps)[(ncol(coords)+1):ncol(df_maps)] <- cat_names
  
  # For each category, create a map
  map_list <- list()
  for (cat in cat_names) {
    p <- ggplot(df_maps, aes(x = X, y = Y, fill = .data[[cat]])) +
      geom_tile() +
      scale_fill_viridis_c(option = "magma", labels = percent, name = "Share") +
      coord_fixed() +
      theme_void() +
      labs(title = cat) +
      theme(legend.position = "right")
    
    ggsave(file.path(output_dir, paste0("map_", tolower(cat), ".png")), p, width = 8, height = 7)
    map_list[[cat]] <- p
  }
  
  # Combine into mosaic
  p_mosaic <- wrap_plots(map_list, ncol = 3) + 
    plot_annotation(title = "Spatial Distribution of Predicted Land-Use Shares",
                    subtitle = paste0("Hybrid BART-MNL Model (", res_km, "km Resolution)"))
  
  ggsave(file.path(output_dir, "spatial_mosaic.png"), p_mosaic, width = 18, height = 12)
}

# --- 6. FIT ANALYSIS ---

plot_mnl_fit_analysis <- function(Y_obs, Y_pred, cat_names, output_dir) {
  cat("Performing fit analysis...\n")
  
  # Ensure dimensions match
  n_obs_pred <- nrow(Y_pred)
  n_cats_pred <- ncol(Y_pred)
  
  if (length(cat_names) != n_cats_pred) {
    cat_names <- paste0("Cat_", 1:n_cats_pred)
  }
  
  # Flatten and align
  df_fit <- data.frame(
    Observed  = as.vector(Y_obs[, 1:n_cats_pred, drop = FALSE]),
    Predicted = as.vector(Y_pred),
    Category  = rep(cat_names, each = n_obs_pred)
  )
  
  p_calib <- ggplot(df_fit, aes(x = Observed, y = Predicted, color = Category)) +
    geom_point(alpha = 0.1, size = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    facet_wrap(~Category, scales = "free") +
    theme_minimal() +
    labs(title = "Model Calibration: Observed vs Predicted Shares",
         subtitle = "Closer to diagonal = Better fit")
  
  ggsave(file.path(output_dir, "fit_calibration.png"), p_calib, width = 12, height = 10)
}

# --- MAIN REPORTING FUNCTION ---

generate_full_mnl_report <- function(res, X_linear, X_bart, Y_obs, coords, 
                                     cov_names, cat_names, model_label,
                                     output_root = "output/plots/reporting_suite") {
  
  out_dir <- file.path(output_root, model_label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("\n=== GENERATING FULL BART-MNL REPORT [", model_label, "] ===\n")
  
  # 1. Parametric Report
  plot_mnl_parametric_report(res, cov_names, cat_names, out_dir)
  
  # 2. Predicted Shares
  # Using a decent sample size for the final report
  Y_pred <- predict_mnl_hybrid(res, X_linear, X_bart, n_samples = 100)
  
  # 3. Fit Analysis
  plot_mnl_fit_analysis(Y_obs, Y_pred, cat_names, out_dir)
  
  # 4. Spatial Maps
  plot_mnl_spatial_maps(Y_pred, coords, cat_names, out_dir)
  
  # 5. BART PDPs
  # Only for the continuous variables actually passed to BART
  bart_vars <- colnames(X_bart)
  for (v in bart_vars) {
    try(plot_mnl_bart_pdp(res, X_linear, X_bart, v, cat_names, out_dir))
  }
  
  cat("\n=== Report Complete. Files saved to:", out_dir, "===\n")
}
