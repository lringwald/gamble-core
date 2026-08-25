# =============================================================================
# mnlogit_rcpp.R  —  R wrapper that delegates the hot path to C++ via Rcpp
#
# Usage:
#   Rcpp::sourceCpp("mnlogit_gibbs_core.cpp")   # one-time compile
#   source("mnlogit_rcpp.R")
#   fit <- mnlogit_rcpp(X, Y, ...)               # same API as before
# =============================================================================

library(progressr)
library(MASS)
library(TruncatedNormal)
library(dbarts)
library(loo)
library(matrixStats)
library(Matrix)
library(pg)
library(BayesLogit)
library(Rcpp)
library(qs2)
library(RcppArmadillo)

aux_path <- "codes/mnl_aux_func.R"
if (!file.exists(aux_path)) aux_path <- "mnl_aux_func.R"
source(aux_path)

# --- Compile C++ (idempotent — skips if already loaded) ---
# The exists() guard alone can find stale function objects from a previous
# R session whose compiled pointers are now NULL.  We also do a test call
# to detect that case and force recompilation.
.needs_cpp_compile <- TRUE
if (exists("chol_sample_precision_cpp", mode = "function")) {
  .needs_cpp_compile <- tryCatch(
    {
      # tiny test: 2x2 identity precision, zero rhs → should return a numeric(2)
      chol_sample_precision_cpp(diag(2), rep(0, 2))
      FALSE # call succeeded → no recompile needed
    },
    error = function(e) TRUE
  ) # stale pointer → recompile
}
if (.needs_cpp_compile) {
  cpp_path <- "mnlogit_gibbs_core.cpp"
  if (!file.exists(cpp_path)) cpp_path <- "codes/mnlogit_gibbs_core.cpp"

  # Parallel-safe compilation: Use a lock file to ensure only one process compiles at a time
  lock_file <- paste0(cpp_path, ".lock")

  # Simple retry loop with backoff
  max_retries <- 120 # Wait up to 2 minutes
  retry_count <- 0

  while (file.exists(lock_file) && retry_count < max_retries) {
    Sys.sleep(1 + runif(1, 0, 1)) # Random sleep 1-2s
    retry_count <- retry_count + 1
  }

  if (retry_count >= max_retries) {
    warning("Compilation lock timeout for ", cpp_path, ". Proceeding anyway...")
  }

  # Create lock
  cat(Sys.getpid(), file = lock_file)

  # Execute compilation
  tryCatch({
    Rcpp::sourceCpp(cpp_path)
  }, finally = {
    # Always remove lock
    if (file.exists(lock_file)) {
      pid_in_lock <- tryCatch(as.integer(readLines(lock_file, warn = FALSE, n = 1)), error = function(e) 0)
      if (pid_in_lock == Sys.getpid()) {
        unlink(lock_file)
      }
    }
  })
}
rm(.needs_cpp_compile)




# =============================================================================
# Main sampler — Rcpp-accelerated
# =============================================================================
mnlogit_rcpp <- function(X, Y, intercept = FALSE, baseline = ncol(Y),
                         niter = 1000, nburn = 500, A0 = 2.0,
                         empirical_intercept_prior = FALSE,
                         y_weight = NULL, use_re = FALSE, group_idx = NULL,
                         use_bart = FALSE, bart_idx = 1:ncol(X), n_trees_bart = 50, n_threads_bart = 1,
                         store_bart_trees = TRUE, do_slim_trees = TRUE,
                         save_bart_to_disk = FALSE,
                         save_posterior_to_disk = FALSE,
                         disk_path = tempdir(),
                         bart_batch_size = 50,
                         posterior_batch_size = 50,
                         bart_warmup = 0,
                         bart_tempering = FALSE,
                         bart_base = 0.95,
                         bart_power = 2.0,
                         re_cor_threshold = 0.85,
                         re_cor_min_groups = 0.5,
                         store_f = FALSE,
                         symmetric = FALSE,
                         positive_constraints = NULL, negative_constraints = NULL,
                         use_ncp = TRUE, calc_loo = FALSE,
                         support_prior_strength = 0, re_asis = FALSE,
                         use_tempering = FALSE, tempering_T0 = 0.05,
                         loo_scaling_factor = 1, chain_id = NULL, progress_cb = NULL,
                         re_idx = 1:ncol(X), linear_idx = 1:ncol(X), prior_mu = NULL, prior_V_inv = NULL,
                         use_horseshoe = FALSE, horseshoe_idx = NULL,
                         equation_specific_hs = FALSE,
                         estimate_c2 = FALSE, slab_df = 4, slab_s2 = 4, slab_s2_re = 4,
                         p0_mu = NULL, tau0_mu = NULL,
                         standardize = TRUE, method = c("standardize", "center", "scale", "QR", "none"),
                         gamma_matched_pg = FALSE,
                         init_state = NULL, prior_a_re = 0.01, prior_b_re = 0.01,
                         use_half_cauchy_re = TRUE, re_scale_A = 1.0) {
  # --- 1. SETUP (identical to original) ---
  if (niter <= nburn) stop("niter must be > nburn.")
  if (use_bart && nrow(X) > 50000 && nburn < 2000) {
    warning("For massive N with BART, consider nburn > 5000.")
  }
  if (any(!is.finite(Y))) stop("Y contains NA/Inf.")
  if (any(!is.finite(X))) stop("X contains NA/Inf.")
  if (any(rowSums(Y) <= 0)) stop("Zero/negative row sums in Y.")

  bart_disk_path <- disk_path
  posterior_disk_path <- disk_path

  if (is.null(y_weight)) {
    y_weight <- 1 / min(rowSums(Y))
    message("Auto y_weight: ", round(y_weight, 6))
  }
  if (missing(loo_scaling_factor) || loo_scaling_factor == 1) {
    loo_scaling_factor <- 1 / y_weight
  }

  nn_weighted <- as.vector(rowSums(Y) * y_weight)
  kappa_weighted <- (Y * y_weight) - (nn_weighted / 2)

  X <- as.matrix(X)
  if (intercept && !any(colSums(X) == nrow(X))) {
    X <- cbind(X, intercept = 1)
  }
  X_orig <- X

  # Convert character indices to integers early for robustness
  if (is.character(linear_idx)) linear_idx <- match(linear_idx, colnames(X))
  if (is.character(bart_idx)) bart_idx <- match(bart_idx, colnames(X))
  if (is.character(re_idx)) re_idx <- match(re_idx, colnames(X))
  if (is.character(horseshoe_idx)) horseshoe_idx <- match(horseshoe_idx, colnames(X))
  if (is.character(positive_constraints)) positive_constraints <- match(positive_constraints, colnames(X))
  if (is.character(negative_constraints)) negative_constraints <- match(negative_constraints, colnames(X))

  # Remove NAs from indexing (variables not found)
  linear_idx <- linear_idx[!is.na(linear_idx)]
  bart_idx <- bart_idx[!is.na(bart_idx)]

  n <- nrow(X)
  k <- ncol(X)
  p_all <- ncol(Y)
  p <- p_all - 1
  pp <- (1:p_all)[-baseline]

  # --- 2. PREPROCESSING ---
  if (missing(method)) {
    methods_chosen <- if (standardize) "standardize" else "none"
  } else {
    methods_chosen <- match.arg(method, several.ok = TRUE)
  }

  do_std <- "standardize" %in% methods_chosen
  do_cen <- "center" %in% methods_chosen || do_std
  do_scl <- "scale" %in% methods_chosen || do_std
  do_qr <- "QR" %in% methods_chosen

  # Identify continuous variables (needed for rescaling AND BART scaling)
  # Continuous variables are only candidates for centering/scaling if they are not naturally bounded in [-1, 1] (e.g. shares, contrasts, bounded indices)
  cont_idx <- apply(X, 2, function(x) {
    is.numeric(x) &&
      length(unique(x)) > 2 &&
      max(abs(x), na.rm = TRUE) > 1.0
  })
  X_rescaling <- matrix(0, 2, sum(cont_idx))
  rownames(X_rescaling) <- c("mean", "sd")
  X_rescaling[2, ] <- 1 # Default SD=1

  if (do_cen || do_scl) {
    if (sum(cont_idx) > 0) {
      if (do_cen) {
        X_rescaling[1, ] <- colMeans(X[, cont_idx, drop = FALSE])
        X[, cont_idx] <- sweep(X[, cont_idx, drop = FALSE], 2, X_rescaling[1, ], "-")
      }
      if (do_scl) {
        X_rescaling[2, ] <- apply(X[, cont_idx, drop = FALSE], 2, sd)
        X_rescaling[2, X_rescaling[2, ] == 0] <- 1
        X[, cont_idx] <- sweep(X[, cont_idx, drop = FALSE], 2, X_rescaling[2, ], "/")
      }
    }
  }

  # --- 2.05 Linear Design Matrix Selection & Index Remapping ---
  X_full_processed <- X
  X_full_for_bart <- X

  if (intercept) {
    # Ensure intercept is in linear_idx if it was added/exists
    int_pos <- which(colnames(X) == "intercept")
    if (length(int_pos) > 0 && !(int_pos %in% linear_idx)) {
      linear_idx <- c(linear_idx, int_pos)
    }
  }

  # Map indices from full X to linear subset (to maintain support for re_idx, horseshoe_idx, etc.)
  map_idx <- function(idx, full_idx) {
    if (is.null(idx)) {
      return(NULL)
    }
    m <- match(idx, full_idx)
    m[!is.na(m)]
  }

  re_idx <- map_idx(re_idx, linear_idx)
  horseshoe_idx <- map_idx(horseshoe_idx, linear_idx)
  positive_constraints <- map_idx(positive_constraints, linear_idx)
  negative_constraints <- map_idx(negative_constraints, linear_idx)

  # Redefine X and k for the linear part
  X <- X_full_processed[, linear_idx, drop = FALSE]
  k <- ncol(X)
  X_orig <- X_orig[, linear_idx, drop = FALSE]

  bart_scaling <- NULL
  if (use_bart) {
    # Robust bart_scaling: map bart_idx to X_rescaling columns
    # Note: X_rescaling only contains continuous variables
    bart_mu <- rep(0, length(bart_idx))
    bart_sd <- rep(1, length(bart_idx))

    # full_to_rescale maps full column index to X_rescaling column index
    full_to_rescale <- rep(NA, ncol(X_full_processed))
    full_to_rescale[cont_idx] <- 1:sum(cont_idx)

    for (i in seq_along(bart_idx)) {
      full_idx <- bart_idx[i]
      rescale_idx <- full_to_rescale[full_idx]
      if (!is.na(rescale_idx)) {
        bart_mu[i] <- X_rescaling[1, rescale_idx]
        bart_sd[i] <- X_rescaling[2, rescale_idx]
      }
    }
    bart_scaling <- list(mu = bart_mu, sd = bart_sd)
  }

  if (do_qr) {
    if (use_horseshoe || !is.null(positive_constraints) || !is.null(negative_constraints)) {
      warning(
        "QR decomposition ('method=QR') is active alongside Horseshoe shrinkage or sign constraints.\n",
        "Back-transformation (beta = R^-1 * theta) will destroy coefficient sparsity and sign consistency.\n",
        "Recommendation: Use 'method=c(\"center\", \"scale\")' instead of 'QR' for these configurations."
      )
    }
    # Apply QR only to the subsetted linear matrix
    qr_X <- qr(X)
    X_R_mat <- qr.R(qr_X)
    X <- qr.Q(qr_X)
  }

  # Update continuous variable tracking for back-transformation
  if (do_cen || do_scl) {
    # 1. Map each column of the FULL X to its position in the ORIGINAL X_rescaling
    full_to_rescale <- rep(NA, length(cont_idx))
    names(full_to_rescale) <- names(cont_idx)
    full_to_rescale[cont_idx] <- seq_len(sum(cont_idx))

    # 2. For the subsetted X (linear_idx), find which rescale columns we need
    needed_rescale_cols <- full_to_rescale[linear_idx]

    # Validation Check
    if (any(cont_idx[linear_idx] & is.na(needed_rescale_cols))) {
      offenders <- linear_idx[which(cont_idx[linear_idx] & is.na(needed_rescale_cols))]
      offender_names <- if (!is.null(colnames(X_full_processed))) colnames(X_full_processed)[offenders] else "unnamed"
      stop(sprintf(
        "Alignment Error: Continuous variables at indices [%s] (%s) have no corresponding rescaling factors in X_rescaling.",
        paste(offenders, collapse = ","), paste(offender_names, collapse = ",")
      ))
    }

    # Keep only those that are actually continuous
    needed_rescale_cols <- needed_rescale_cols[!is.na(needed_rescale_cols)]

    # 3. Subset X_rescaling and cont_idx
    X_rescaling <- X_rescaling[, needed_rescale_cols, drop = FALSE]
    cont_idx <- cont_idx[linear_idx]

    # Final Sanity Check
    if (sum(cont_idx) != ncol(X_rescaling)) {
      stop(sprintf(
        "Dimension Mismatch: %d continuous variables found but X_rescaling has %d columns.",
        sum(cont_idx), ncol(X_rescaling)
      ))
    }
  }

  # --- 2.1 BART Design Matrix Selection ---
  if (use_bart) {
    # Select subset of processed covariates for BART from the full matrix
    X_bart_final <- X_full_for_bart[, bart_idx, drop = FALSE]
    # Robust calibration references: use empirical range to ensure different terminal nodes
    bart_ref1_mat <- matrix(apply(X_bart_final, 2, min), 1, length(bart_idx))
    bart_ref2_mat <- matrix(apply(X_bart_final, 2, max), 1, length(bart_idx))
  }

  # --- 3. INITIALIZATION ---
  # Hoist matrix conversions and unnaming to avoid overhead in Gibbs loop
  cov_names_save <- colnames(X)
  int_idx <- integer(0)
  if (!is.null(cov_names_save)) {
    int_idx <- which(tolower(cov_names_save) == "intercept")
  }
  if (length(int_idx) == 0) {
    int_idx <- which(apply(X_orig, 2, function(col) all(col == 1)))
  }
  X <- as.matrix(unname(X))
  Xt <- as.matrix(unname(t(X)))
  kappa_weighted <- as.matrix(unname(kappa_weighted))

  curr_beta <- matrix(0, k, p)
  bart_shifts <- matrix(0, p, 1)
  prior_V_inv_originally_null <- is.null(prior_V_inv)
  if (is.null(prior_mu)) prior_mu <- matrix(0, k, p)
  if (!is.matrix(prior_mu)) prior_mu <- matrix(as.vector(prior_mu), k, p)
  bart_alpha <- 1.0 # NEW: Global tempering factor
  if (is.null(prior_V_inv)) {
    # Default: A0 for predictors, but allow wide prior (variance = 100.0, sd = 10) for the intercept (similar to brms/Stan)
    diag_vals <- rep(1 / A0, k)
    if (length(int_idx) > 0) {
      diag_vals[int_idx] <- 1 / 100.0
    } else if (intercept) {
      diag_vals[k] <- 1 / 100.0
    }
    prior_V_inv <- diag(diag_vals, k)
  }

  if (use_re && use_half_cauchy_re) {
    a_re <- matrix(re_scale_A^2 / 2, nrow = k, ncol = p)
  } else {
    a_re <- NULL
  }

  if (empirical_intercept_prior) {
    if (length(int_idx) > 0) {
      Y_smoothed <- colSums(Y) + 1.0
      empirical_log_odds <- log(Y_smoothed / Y_smoothed[baseline])
      empirical_mu <- empirical_log_odds[-baseline]

      for (idx in int_idx) {
        prior_mu[idx, ] <- empirical_mu
        if (prior_V_inv_originally_null) {
          prior_V_inv[idx, idx] <- 0.25
        }
      }
    }
  }

  prior_P <- as.matrix(unname(prior_V_inv))
  prior_Pb <- as.matrix(unname(prior_P %*% prior_mu)) # k x p, static

  # Baseline RE precision: set to effectively infinite (1e12) to lock fixed effects
  # Initialize at the half-Cauchy prior mode for re_scale_A for early iteration exploration
  prec_init <- 1 / (re_scale_A^2)
  prec_beta_pooled <- matrix(prec_init, k, p)
  sigma_beta_pooled <- matrix(re_scale_A, k, p)

  if (length(re_idx) < k) {
    fixed_idx <- setdiff(seq_len(k), re_idx)
    prec_beta_pooled[fixed_idx, ] <- 1e12
    sigma_beta_pooled[fixed_idx, ] <- 1e-6
  }

  # Precompute complement of re_idx to avoid recomputing in the Gibbs loops
  re_complement <- setdiff(seq_len(k), re_idx)

  mu_pooled <- curr_beta # Initial reference

  if (use_re) {
    groups <- unique(group_idx)
    n_groups <- length(groups)
    # Pre-convert to unnamed integer lists and matrices for Rcpp bridge speed
    idx_list <- lapply(groups, function(m) as.integer(unname(which(group_idx == m))))
    Xm_list <- lapply(idx_list, function(idx) as.matrix(unname(X[idx, , drop = FALSE])))
    Xmt_list <- lapply(Xm_list, function(m) as.matrix(unname(t(m))))

    # Precompute structural masking matrix for group-specific REs
    re_mask <- matrix(1.0, nrow = k, ncol = n_groups)
    y_mask <- matrix(1.0, nrow = p, ncol = n_groups)

    # Support-aware RE prior: per (covariate x group) SD multiplier in (0,1] = 1/sqrt((n_g/PR)^strength),
    # PR = within-group participation ratio (effective # obs informing that group's slope). Shrinks
    # information-sparse RE slopes (the graduated generalization of the A2 no-variation pin); column m
    # aligns with Xm_list[[m]] / re_mask column m. All-ones when off. Globally-constant cols (random
    # intercept) exempted. Applied as the SD scaler `sig_m = sig %*% re_mask %*% re_support` in the C++ core.
    re_support_mat <- matrix(1.0, nrow = k, ncol = n_groups)
    if (support_prior_strength > 0 && length(re_idx) >= 1) {
      re_glob_const <- apply(X[, re_idx, drop = FALSE], 2, function(cc) { s <- sd(cc); !is.finite(s) || s < 1e-8 })
      for (m in seq_len(n_groups)) {
        Xg <- Xm_list[[m]][, re_idx, drop = FALSE]; ng <- nrow(Xg)
        Xc <- sweep(Xg, 2, colMeans(Xg), "-")             # within-group centering -> slope information (PR is not location-invariant)
        s2 <- colSums(Xc^2); s4 <- colSums(Xc^4); pr <- s2^2 / pmax(s4, 1e-12)
        sd_fac <- 1 / sqrt(pmax((ng / pmax(pr, 1))^support_prior_strength, 1e-12))
        sd_fac[re_glob_const] <- 1.0
        re_support_mat[re_idx, m] <- sd_fac
      }
      message(sprintf("RE support prior (participation ratio, strength %.2f): min SD factor %.3f.", support_prior_strength, min(re_support_mat[re_idx, ])))
    }
    is_global_intercept <- sapply(seq_len(k), function(v) {
      all(X[, v] == 1.0)
    })
    group_sizes <- sapply(idx_list, length)

    for (v in seq_len(k)) {
      if (is_global_intercept[v]) next
      var_global <- var(X[, v])
      if (is.na(var_global) || var_global < 1e-12) {
        re_mask[v, ] <- 0.0
        next
      }
      for (m in 1:n_groups) {
        var_local <- var(Xm_list[[m]][, v])
        std_var <- var_local / var_global
        n_nonzero <- sum(Xm_list[[m]][, v] != 0.0)

        # Scale-dependent non-zero threshold (3 to 10 observations depending on group size)
        x_thresh <- min(10L, max(3L, as.integer(ceiling(0.02 * group_sizes[m]))))

        # Mask if local variance is less than 0.5% of global OR has sparse non-zero entries
        if (is.na(std_var) || std_var < 0.005 || n_nonzero < x_thresh) {
          re_mask[v, m] <- 0.0
        }
      }
    }

    # Y-side structural masking [p x n_groups] based on cell counts and standardized variance
    mean_pixel_total <- mean(rowSums(Y))

    for (ip in seq_len(p)) {
      j <- pp[ip]
      var_global <- var(Y[, j])
      if (is.na(var_global) || var_global < 1e-8) {
        y_mask[ip, ] <- 0.0
        next
      }

      for (m in seq_len(n_groups)) {
        group_obs <- Y[idx_list[[m]], j]

        # Calculate number of effective positive observations (scale/unit invariant)
        n_positive_eff <- sum(group_obs) / mean_pixel_total

        var_local <- var(group_obs)
        std_var <- var_local / var_global

        # Scale-dependent category count threshold (2 to 5 observations depending on group size)
        y_thresh <- min(5L, max(2L, as.integer(ceiling(0.01 * group_sizes[m]))))

        # Mask if less than threshold positive outcomes in group or local variance is zero/near-zero
        if (n_positive_eff < y_thresh || is.na(std_var) || std_var < 1e-4) {
          y_mask[ip, m] <- 0.0
        }
      }
    }

    # Baseline support check: if baseline has too few observations in group m, mask all REs in that group
    baseline_eff_counts <- sapply(idx_list, function(idx) sum(Y[idx, baseline]) / mean_pixel_total)

    for (m in seq_len(n_groups)) {
      # Scale-dependent baseline threshold (2 to 5 observations depending on group size)
      b_thresh <- min(5L, max(2L, as.integer(ceiling(0.01 * group_sizes[m]))))
      if (baseline_eff_counts[m] < b_thresh) {
        re_mask[, m] <- 0.0
        y_mask[, m] <- 0.0
      }
    }

    # Diagnostic summary
    n_y_masked <- sum(y_mask == 0)
    if (n_y_masked > 0) {
      message(sprintf(
        "Y-mask: %d of %d category-group combinations fully masked (constant share or sparse cell).",
        n_y_masked, p * n_groups
      ))
    }

    # =================================================================
    # AUTOMATIC RE COLLINEARITY SCREEN
    # Computes within-group pairwise correlations for all RE predictors.
    # For any pair exceeding cor_threshold in a majority of groups,
    # the variable with lower global variance is removed from re_mask
    # for ALL groups (structural collinearity → pooled-only effect).
    # Variables removed only in specific groups (local collinearity)
    # are masked just for those groups.
    # =================================================================

    if (length(re_idx) >= 2) {
      re_col_names <- colnames(X)[re_idx] # names for messaging
      n_re_vars <- length(re_idx)

      # --- Per-group correlation matrices ---
      # pairwise_hits[i, j] = number of groups where |cor(i,j)| > threshold
      pairwise_hits <- matrix(0L, n_re_vars, n_re_vars)
      pairwise_maxcor <- matrix(0.0, n_re_vars, n_re_vars)
      group_cor_list <- vector("list", n_groups)

      for (m in seq_len(n_groups)) {
        idx_m <- idx_list[[m]]
        if (length(idx_m) < 10L) {
          group_cor_list[[m]] <- matrix(NA_real_, n_re_vars, n_re_vars)
          next
        }
        Xm_re <- Xm_list[[m]][, re_idx, drop = FALSE]

        # Suppress warnings for zero-variance columns (already masked)
        cm <- suppressWarnings(cor(Xm_re))
        cm[is.na(cm)] <- 0.0
        diag(cm) <- 0.0 # ignore self-correlation
        group_cor_list[[m]] <- cm

        hit <- abs(cm) > re_cor_threshold
        pairwise_hits <- pairwise_hits + hit
        pairwise_maxcor <- pmax(pairwise_maxcor, abs(cm))
      }

      # Valid groups (enough obs) for denominator
      n_valid_groups <- sum(sapply(group_cor_list, function(cm) !all(is.na(cm))))
      hit_fraction <- pairwise_hits / max(n_valid_groups, 1L)

      # --- Identify collinear pairs ---
      # Work upper-triangle only to avoid double-processing
      globally_removed <- integer(0) # re_idx positions removed everywhere
      locally_masked <- list() # re_idx positions masked per group

      for (i in seq_len(n_re_vars - 1L)) {
        for (j in (i + 1L):n_re_vars) {
          if (hit_fraction[i, j] < re_cor_min_groups) next # not widespread enough

          # Decide which variable to drop: lower global variance loses
          var_i <- var(X[, re_idx[i]])
          var_j <- var(X[, re_idx[j]])
          drop_pos <- if (var_i <= var_j) i else j # position in re_idx
          keep_pos <- if (drop_pos == i) j else i

          is_global <- hit_fraction[i, j] >= 0.75 # >75% of groups → global drop

          if (is_global) {
            if (!(drop_pos %in% globally_removed)) {
              globally_removed <- c(globally_removed, drop_pos)
              message(sprintf(
                paste0(
                  "RE collinearity [GLOBAL]: removing '%s' from re_idx ",
                  "(r > %.2f in %.0f%% of groups; keeping '%s')."
                ),
                re_col_names[drop_pos],
                re_cor_threshold,
                hit_fraction[i, j] * 100,
                re_col_names[keep_pos]
              ))
            }
          } else {
            # Local: mask only in the offending groups
            bad_groups <- which(sapply(
              group_cor_list,
              function(cm) !is.na(cm[i, j]) && abs(cm[i, j]) > re_cor_threshold
            ))
            for (m in bad_groups) {
              if (re_mask[re_idx[drop_pos], m] != 0.0) {
                re_mask[re_idx[drop_pos], m] <- 0.0
                locally_masked[[length(locally_masked) + 1L]] <-
                  list(
                    var = re_col_names[drop_pos], group = m,
                    cor = group_cor_list[[m]][i, j]
                  )
              }
            }
            if (length(locally_masked) > 0L) {
              message(sprintf(
                paste0(
                  "RE collinearity [LOCAL]: masking '%s' in %d group(s) ",
                  "(r > %.2f; keeping '%s')."
                ),
                re_col_names[drop_pos],
                length(bad_groups),
                re_cor_threshold,
                re_col_names[keep_pos]
              ))
            }
          }
        }
      }

      # --- A2: within-group no-X-variation (the correlation screen MISSES constants:
      #     a constant column -> NA cor -> 0 -> never flagged). Mask such RE slopes locally.
      #     Exempt globally-constant columns: the intercept's random intercept is identified. ---
      re_global_const <- apply(X[, re_idx, drop = FALSE], 2, function(cc) { s <- sd(cc); !is.finite(s) || s < 1e-8 })
      n_a2 <- 0L
      for (m in seq_len(n_groups)) {
        sds_m <- apply(Xm_list[[m]][, re_idx, drop = FALSE], 2, sd)
        for (r in which((!is.finite(sds_m) | sds_m < 1e-8) & !re_global_const)) {
          if (re_mask[re_idx[r], m] != 0.0) { re_mask[re_idx[r], m] <- 0.0; n_a2 <- n_a2 + 1L }
        }
      }
      if (n_a2 > 0L) message(sprintf("RE no-variation screen (A2): masked %d (group x covariate) slopes with no within-group X variation.", n_a2))

      # --- Apply global removals to re_mask ---
      if (length(globally_removed) > 0L) {
        re_mask[re_idx[globally_removed], ] <- 0.0
        message(sprintf(
          "RE collinearity screen complete: %d variable(s) globally demoted to pooled-only effects.",
          length(unique(globally_removed))
        ))
      }

      # --- Summary table ---
      if (length(globally_removed) > 0L || length(locally_masked) > 0L) {
        cat("\n--- RE Collinearity Screen Summary ---\n")
        if (length(globally_removed) > 0L) {
          cat("Globally removed from RE:\n")
          cat(paste0("  ", re_col_names[unique(globally_removed)], "\n"))
        }
        if (length(locally_masked) > 0L) {
          cat("Locally masked (group-specific):\n")
          for (lm_entry in locally_masked) {
            cat(sprintf(
              "  %s in group %d (r = %.3f)\n",
              lm_entry$var, lm_entry$group, lm_entry$cor
            ))
          }
        }
        cat("--------------------------------------\n\n")
      } else {
        message("RE collinearity screen: no problematic pairs detected.")
      }

      rm(
        pairwise_hits, pairwise_maxcor, group_cor_list,
        hit_fraction, n_valid_groups
      )
    }

    # Precompute masked coordinates by group and category to avoid repeating which() in inner loops
    masked_by_group_cat <- lapply(seq_len(n_groups), function(m) {
      idx_m <- idx_list[[m]]
      lapply(seq_len(p), function(ip) {
        j <- pp[ip]
        base_masked <- which((re_mask[, m] == 0.0 | (y_mask[ip, m] == 0.0 & !is_global_intercept)) & seq_len(k) %in% re_idx)

        if (length(base_masked) == length(re_idx) || length(idx_m) < 10L) {
          return(base_masked)
        }

        # Dynamic outcome presence threshold — use absolute threshold
        # based on mean_pixel_total to avoid degeneracy when Y is stabilized
        # (where min(Y[,j]) ≈ stabilization epsilon, making relative threshold trivially pass)
        presence_thresh <- mean_pixel_total * 0.01
        group_obs <- Y[idx_m, j]
        is_j <- (group_obs > presence_thresh)

        # Find baseline outcome presence in this group
        is_base <- (Y[idx_m, baseline] > presence_thresh)

        additional_masks <- integer(0)

        if (any(is_j) && any(is_base)) {
          for (v in setdiff(re_idx, base_masked)) {
            if (is_global_intercept[v]) next
            vals <- Xm_list[[m]][, v]
            vals_1 <- vals[is_j]
            vals_base <- vals[is_base]

            # 1. Quantile Boundaries (Ignore extreme 5% tails)
            q1_low  <- quantile(vals_1, 0.05, names = FALSE)
            q1_high <- quantile(vals_1, 0.95, names = FALSE)
            qb_low  <- quantile(vals_base, 0.05, names = FALSE)
            qb_high <- quantile(vals_base, 0.95, names = FALSE)

            # 2. Local Variance Check (Loosened threshold)
            var_1 <- var(vals_1)
            if (is.na(var_1)) var_1 <- 0.0

            # 3. Zero-Bounded Functional Separation Check
            # If 95% of one group is practically zero, but the 95th percentile of the other is substantive
            base_is_zero <- (qb_high < 1e-3)
            out_is_zero  <- (q1_high < 1e-3)
            functional_separation <- (base_is_zero && q1_high > 0.05) || 
                                     (out_is_zero && qb_high > 0.05)

            # Trigger mask if distributions don't meaningfully overlap, 
            # if variance is critically low, or if zero-bounded separation exists.
            if (q1_low >= qb_high || qb_low >= q1_high || var_1 < 1e-4 || functional_separation) {
              additional_masks <- c(additional_masks, v)
            }
          }
        } else if (any(is_j) && !any(is_base)) {
          # If there is no baseline at all in this group, every RE for this category is unidentifiable
          additional_masks <- c(additional_masks, setdiff(re_idx, base_masked))
        } else if (!any(is_j)) {
          # If the target category is effectively absent (failed presence_thresh),
          # mask all covariates to prevent MCMC drift. We leave the intercept unmasked 
          # so it can cleanly drift to -Inf and correctly predict probability ~ 0.
          cov_idx <- setdiff(re_idx, base_masked)
          cov_idx <- cov_idx[!is_global_intercept[cov_idx]]
          additional_masks <- c(additional_masks, cov_idx)
        }

        sort(unique(c(base_masked, additional_masks)))
      })
    })

    curr_beta_c <- array(0, c(k, p, n_groups))
    mu_pooled <- matrix(0, k, p)
    # group_sizes already computed at line 399 — no recomputation needed

    # ---------------------------------------------------------------
    # PRECOMPUTE separation-aware masks (Issues 1-3 from Doc 10)
    # These are static (masked_by_group_cat never changes), so compute
    # once here instead of rebuilding every iteration in the Gibbs loop.
    # ---------------------------------------------------------------

    # 3D category-specific mask for precision updates
    re_mask_cube <- array(1.0, dim = c(k, p, n_groups))
    for (m_d in seq_len(n_groups)) {
      for (ip_d in seq_len(p)) {
        re_mask_cube[, ip_d, m_d] <- re_mask[, m_d]
        masked_vars <- masked_by_group_cat[[m_d]][[ip_d]]
        if (length(masked_vars) > 0L) {
          re_mask_cube[masked_vars, ip_d, m_d] <- 0.0
        }
      }
    }

    # Issue 3: sep_mask_by_ip — per-category separation mask matrices
    # for the ASIS interweaving step. Each element is a k × n_groups
    # matrix that can be applied via element-wise multiplication.
    sep_mask_by_ip <- lapply(seq_len(p), function(ip) {
      m_mat <- matrix(1.0, k, n_groups)
      for (m in seq_len(n_groups)) {
        sv <- masked_by_group_cat[[m]][[ip]]
        if (length(sv) > 0L) m_mat[sv, m] <- 0.0
      }
      m_mat
    })
  }

  # Pre-compute group-size weights for RE BART intercept absorption
  if (use_re && use_bart) {
    weights <- group_sizes / sum(group_sizes)
  }

  if (use_re && use_ncp) {
    z_c <- array(0, c(k, p, n_groups))
  }

  # --- Horseshoe Setup ---
  if (use_horseshoe) {
    hs_idx <- if (is.null(horseshoe_idx)) 1:k else horseshoe_idx
    k_hs <- length(hs_idx)

    # NEW: Tau0 Calibration (Piironen-Vehtari recommended)
    # sigma_util ≈ pi/sqrt(3) for Logit link
    sigma_util <- pi / sqrt(3)
    p0_m <- p0_mu %||% max(1, round(k_hs * 0.3))

    # Ensure p0 < k_hs to avoid Inf or negative global shrinkage
    p0_m <- min(p0_m, k_hs - 0.5)

    # Effective N: pooled uses total observations (Piironen-Vehtari recommendation).
    # Using n_groups systematically under-shrinks when n_groups << n.
    N_eff_pooled <- n

    tau0_pooled <- tau0_mu %||% ((p0_m / (k_hs - p0_m)) * (sigma_util / sqrt(N_eff_pooled)))

    cat(sprintf(
      "Horseshoe Calibration: tau0_pooled=%.4f (p0=%.1f)\n",
      tau0_pooled, p0_m
    ))

    if (equation_specific_hs) {
      hs_lambda2 <- matrix(1, k, p)
      hs_nu <- matrix(1, k, p)
    } else {
      hs_lambda2 <- rep(1, k)
      hs_nu <- rep(1, k)
    }
    hs_tau2 <- if (equation_specific_hs) rep(tau0_pooled^2, p) else tau0_pooled^2
    hs_xi <- if (equation_specific_hs) rep(1, p) else 1
    hs_c2 <- slab_s2
    hs_zeta <- 1
    # Pre-allocate precision corrections for C++
    hs_prec_mat <- matrix(0, k, p)
  } else {
    hs_prec_mat <- matrix(0, k, p)
  }
  # Support-aware FIXED-effect prior: global participation-ratio shrinkage of sparse covariates'
  # fixed effects (mirrors CLR's FE support_factor); multiplies the horseshoe precision rows below.
  fe_support <- fe_support_factor(X, support_prior_strength)

  if (use_bart) {
    curr_f <- matrix(0, n, p)
    bart_samplers <- list()
    for (ip in 1:p) {
      j <- pp[ip]
      ctrl <- dbarts::dbartsControl(
        n.samples = 1, n.burn = 0,
        n.trees = n_trees_bart,
        n.threads = n_threads_bart,
        keepTrees = TRUE,
        n.chains = 1, updateState = TRUE
      )
      bart_samplers[[ip]] <- dbarts::dbarts(
        X_bart_final, as.numeric(kappa_weighted[, j] / 0.25),
        control = ctrl, sigma = 1.0,
        tree.prior = dbarts:::cgm(base = bart_base, power = bart_power),
        resid.prior = dbarts:::chisq(df = 1e10, quant = 0.5)
      )
    }
  }

  has_constraints <- !is.null(positive_constraints) || !is.null(negative_constraints)
  if (has_constraints) {
    pv <- diag(MASS::ginv(prior_P))
    psd <- sqrt(pmax(pv, 1e-12))
    start_lo <- rep(-Inf, k)
    start_hi <- rep(Inf, k)
    if (!is.null(positive_constraints)) start_lo[positive_constraints] <- -2 * psd[positive_constraints]
    if (!is.null(negative_constraints)) start_hi[negative_constraints] <- 2 * psd[negative_constraints]
    decay_rate <- 5
  }

  # =========================================================================
  # HOT-START OVERRIDES
  # =========================================================================
  if (!is.null(init_state)) {
    cat(">>> Overriding initial values with hot-start state...\n")

    # --- Dimension Safety Check ---
    if (use_re && !is.null(init_state$beta)) {
      expected_dim <- as.integer(c(k, p, n_groups))
      actual_dim <- dim(init_state$beta)
      if (!identical(expected_dim, actual_dim)) {
        warning(sprintf(
          paste0(
            "Hot-start SKIPPED: beta dimension mismatch.\n",
            "  Expected: [%s]\n  Got:      [%s]\n",
            "  Likely cause: data filtering changed group count or covariate set."
          ),
          paste(expected_dim, collapse = ","),
          paste(actual_dim, collapse = ",")
        ))
        init_state <- NULL # abort hot-start, proceed cold
      }
    }
    if (!is.null(init_state) && !is.null(init_state$mu)) {
      if (!identical(dim(init_state$mu), as.integer(c(k, p)))) {
        warning("Hot-start SKIPPED: mu dimension mismatch.")
        init_state <- NULL
      }
    }

    if (!is.null(init_state$beta)) {
      if (use_re) curr_beta_c <- init_state$beta else curr_beta <- init_state$beta
    }
    if (!is.null(init_state$mu)) mu_pooled <- init_state$mu
    if (!is.null(init_state$prec_beta)) prec_beta_pooled <- init_state$prec_beta
    if (!is.null(init_state$z_c)) z_c <- init_state$z_c
    if (!is.null(init_state$sigma_re)) sigma_beta_pooled <- init_state$sigma_re
    if (!is.null(init_state$a_re)) a_re <- init_state$a_re

    if (use_horseshoe && !is.null(init_state$horseshoe)) {
      hs <- init_state$horseshoe
      if (!is.null(hs$lambda2)) hs_lambda2 <- hs$lambda2
      if (!is.null(hs$tau2)) hs_tau2 <- hs$tau2
      if (!is.null(hs$nu)) hs_nu <- hs$nu
      if (!is.null(hs$xi)) hs_xi <- hs$xi
      if (!is.null(hs$zeta)) hs_zeta <- hs$zeta
      if (!is.null(hs$c2)) hs_c2 <- hs$c2
    }

    if (use_bart && !is.null(init_state$bart_states)) {
      for (ip in 1:p) {
        if (!is.null(init_state$bart_states[[ip]])) {
          bart_samplers[[ip]]$setState(init_state$bart_states[[ip]])
        }
      }
    }
  }

  # --- 4. STORAGE & METADATA ---
  nretain <- niter - nburn

  if (save_bart_to_disk || save_posterior_to_disk) {
    dir.create(disk_path, recursive = TRUE, showWarnings = FALSE)
    meta_file <- file.path(disk_path, "model_metadata.qs")
    # Always save metadata to ensure it matches the current run's subsetting
    # We save it once per directory, but since linear_idx might change,
    # we should ideally overwrite if it exists?
    # No, usually one directory = one model config.
    if (!file.exists(meta_file)) {
      meta <- list(
        cov_names = if (!is.null(cov_names_save)) cov_names_save else paste0("V", 1:k),
        cat_names = if (!is.null(colnames(Y))) colnames(Y)[-baseline] else paste0("C", 1:p),
        baseline_name = if (!is.null(colnames(Y))) colnames(Y)[baseline] else "Base",
        k = k,
        p = p,
        use_re = use_re,
        X_rescaling = if (do_cen || do_scl) X_rescaling else NULL,
        do_cen = do_cen,
        do_scl = do_scl,
        cont_idx = cont_idx
      )
      qs2::qs_save(meta, meta_file)
    }
  }

  if (save_posterior_to_disk) {
    postb_total <- NULL
    postb_pooled <- NULL
  } else {
    postb_total <- if (use_re) array(0, c(k, p, n_groups, nretain)) else array(0, c(k, p, nretain))
    postb_pooled <- array(0, c(k, p, nretain))
  }
  if (use_bart) {
    post_f_sum <- matrix(0, n, p)
    post_f_sum_sq <- matrix(0, n, p)
    if (store_f) post_f <- array(0, c(n, p, nretain))
  }
  post_log_lik <- if (save_posterior_to_disk) numeric(0) else numeric(nretain)
  post_ll_pw <- if (calc_loo && !save_posterior_to_disk) matrix(0, nretain, n) else NULL
  post_kappa_pooled <- if (use_horseshoe && !save_posterior_to_disk) matrix(0, k_hs * p, nretain) else NULL
  post_c2 <- if (use_horseshoe && !save_posterior_to_disk) numeric(nretain) else NULL
  post_sigma_re <- if (use_re && !save_posterior_to_disk) matrix(0, k * p, nretain) else NULL
  tree_store <- if ((store_bart_trees || save_bart_to_disk) && use_bart && !save_posterior_to_disk) vector("list", nretain) else NULL

  # --- Batched Disk Buffers ---
  if (save_bart_to_disk && use_bart && !save_posterior_to_disk) {
    bart_batch_buffer <- list()
    bart_batch_files <- character()
    c_label <- if (is.null(chain_id)) "0" else as.character(chain_id)
    old_files <- list.files(bart_disk_path, pattern = sprintf("^bart_batch_[0-9]+_chain_%s\\.qs$", c_label), full.names = TRUE)
    if (length(old_files) > 0) unlink(old_files)
  }

  if (save_posterior_to_disk) {
    posterior_batch_buffer <- list()
    posterior_batch_files <- character()
    c_label <- if (is.null(chain_id)) "0" else as.character(chain_id)
    old_files <- list.files(posterior_disk_path, pattern = sprintf("^posterior_batch_[0-9]+_chain_%s\\.qs$", c_label), full.names = TRUE)
    if (length(old_files) > 0) unlink(old_files)
  }

  # Fix #3: Explicit initialization avoids fragile exists() scoping that could
  # pick up stale variables from a prior call in the parent environment.
  postb_total_std <- NULL
  postb_pooled_std <- NULL

  # Fix #2: Precompute back-transformation indices for disk path.
  # These use dedicated names (bt_*) to avoid shadowing the BART intercept
  # int_idx set on line 312, which is still needed inside Section E.
  if (save_posterior_to_disk && (do_cen || do_scl) && sum(cont_idx) > 0) {
    bt_is_cont <- which(cont_idx)
    bt_int_idx <- which(apply(X_orig, 2, var) == 0)
  }

  # Convert pp to integer vector for C++ (already 1-based, which C++ code expects)
  pp_int <- as.integer(pp)
  group_idx_0 <- if (use_re) as.integer(match(group_idx, groups) - 1L) else integer(0)

  if (is.null(chain_id)) {
    cat("Sampling (Rcpp-accelerated)...\n")
    pb <- utils::txtProgressBar(min = 0, max = niter, style = 3)
  } else {
    cat(sprintf("Chain %d: Sampling (Rcpp-accelerated)...\n", chain_id))
  }

  # =====================================================================
  # GIBBS LOOP
  # =====================================================================
  burn_buffer <- max(50, floor(0.05 * nburn))
  # TEMPERING SETUP
  if (use_tempering) {
    cat(sprintf("Likelihood Tempering Enabled: linear from T=%.3f to 1.0 over first half of burnin (%d iters)\n", tempering_T0, floor(nburn/4)))
  }
  nburn_half <- max(1L, floor(nburn / 4))
  for (iter in 1:niter) {
    if (use_tempering && iter <= nburn_half && nburn > 0) {
      temp_iter <- tempering_T0 + (1.0 - tempering_T0) * (iter / nburn_half)
    } else {
      temp_iter <- 1.0
    }
    nn_weighted_iter <- nn_weighted * temp_iter
    kappa_weighted_iter <- kappa_weighted * temp_iter
    # --- BART Tempering ---
    if (use_bart) {
      if (iter <= bart_warmup) {
        bart_alpha <- 0
      } else if (bart_tempering && iter <= nburn) {
        # Linear ramp from 0.1 to 1.0 during remaining burn-in
        # Buffer: ensures at least 10% of burn-in is spent at full alpha
        tempering_buffer <- floor(nburn / 10)
        ramp_end <- nburn - tempering_buffer

        if (iter <= ramp_end && ramp_end > bart_warmup) {
          progress <- (iter - bart_warmup) / (ramp_end - bart_warmup)
          bart_alpha <- 0.1 + 0.9 * max(0, min(1, progress))
        } else {
          bart_alpha <- 1.0
        }
      } else {
        bart_alpha <- 1.0
      }
    }

    # --- Annealing walls ---
    if (has_constraints) {
      if (iter <= nburn) {
        df <- exp(-decay_rate * (iter / nburn))
        cur_lo <- start_lo * df
        cur_hi <- start_hi * df
        if (iter > nburn * 0.95) {
          if (!is.null(positive_constraints)) cur_lo[positive_constraints] <- 0
          if (!is.null(negative_constraints)) cur_hi[negative_constraints] <- 0
        }
      } else {
        cur_lo <- start_lo
        cur_hi <- start_hi
        if (!is.null(positive_constraints)) cur_lo[positive_constraints] <- 0
        if (!is.null(negative_constraints)) cur_hi[negative_constraints] <- 0
      }
    }

    # =================================================================
    # A. UPDATE UTILITIES + C_J  —  delegated to C++
    # =================================================================
    f_bart_mat <- if (use_bart) bart_alpha * curr_f else matrix(0, n, p)

    if (use_re) {
      uc <- update_utilities_and_cj_re(
        X, curr_beta_c, f_bart_mat, pp_int, group_idx_0,
        baseline, p_all, use_bart
      )
      U <- uc$U
      c_j_mat <- uc$c_j_mat
    } else {
      # POOLED: full C++ path
      uc <- update_utilities_and_cj(
        X, curr_beta, f_bart_mat, pp_int,
        baseline, p_all, use_bart
      )
      U <- uc$U
      c_j_mat <- uc$c_j_mat
    }

    # =================================================================
    # B. POLYA-GAMMA DRAWS  —  batched in R (already compiled C under the hood)
    # =================================================================
    all_z <- as.vector(U[, pp, drop = FALSE] - c_j_mat)
    all_omega <- fast_rpg(n * p, rep(nn_weighted_iter, p), all_z, gamma_matched_pg = gamma_matched_pg)
    omega <- matrix(pmax(all_omega, 1e-6), n, p)

    # =================================================================
    # C. COEFFICIENT SAMPLING  —  delegated to C++
    # =================================================================
    if (!use_re && !has_constraints) {
      # >>> FAST PATH: entire category loop in one C++ call <<<
      curr_beta <- gibbs_step_pooled(
        X, Xt, kappa_weighted_iter, omega, c_j_mat,
        prior_P, hs_prec_mat, prior_Pb, pp_int,
        f_bart_mat, use_bart
      )
    } else if (use_re && !has_constraints) {
      # >>> RE path in C++ <<<
      sigma_mat <- matrix(
        1 / sqrt(pmax(prec_beta_pooled, 1e-8)),
        nrow = k, ncol = p
      )

      if (use_ncp) {
        # NOTE (Issue 6, Doc 10): gibbs_step_re_ncp uses re_mask/y_mask only.
        # Dynamic separation masks (masked_by_group_cat) are applied one
        # iteration later in ASIS step C.5. This one-step lag is acceptable —
        # the ASIS correction is deterministic and takes effect before the
        # next iteration's utility update.
        re_res <- gibbs_step_re_ncp(
          X, Xt, kappa_weighted_iter,
          omega, c_j_mat,
          prior_P, hs_prec_mat, prior_Pb,
          mu_pooled, sigma_mat,
          z_c, pp_int, idx_list,
          Xm_list, Xmt_list, n_groups, f_bart_mat,
          use_bart, as.integer(re_idx - 1L),
          re_mask, y_mask, as.integer(is_global_intercept),
          re_support_mat,
          if (use_half_cauchy_re && !is.null(a_re)) a_re else matrix(1, k, p),   # half-Cauchy aux for RE-scale ASIS
          isTRUE(re_asis) && use_half_cauchy_re && !is.null(a_re)                 # RE-scale ASIS interweave
        )

        curr_beta_c <- re_res$beta_c
        mu_pooled <- re_res$mu
        z_c <- re_res$z_c
      } else {
        re_res <- gibbs_step_re(
          X, Xt, kappa_weighted_iter, omega, c_j_mat,
          prior_P, hs_prec_mat, prior_Pb, mu_pooled, prec_beta_pooled,
          pp_int, idx_list, Xm_list, Xmt_list,
          n_groups, f_bart_mat, use_bart,
          as.integer(re_idx - 1L), # 0-based for C++
          re_mask, y_mask, as.integer(is_global_intercept)
        )
        curr_beta_c <- re_res$beta_c
        mu_pooled <- re_res$mu
      }
    } else {
      # =================================================================
      # C. CONSTRAINED PATH
      # Branches are evaluated ONCE per iteration, not per category.
      # sigma_mat and re_complement are hoisted outside the ip loop.
      # =================================================================

      if (!use_re) {
        # -----------------------------------------------------------------
        # C1. POOLED CONSTRAINED (no RE)
        # -----------------------------------------------------------------
        for (ip in 1:p) {
          j <- pp[ip]
          om_p <- omega[, ip]
          c_j <- c_j_mat[, ip]
          target <- kappa_weighted_iter[, j] + om_p * c_j
          if (use_bart) target <- target - om_p * (bart_alpha * curr_f[, ip])

          P <- prior_P + weighted_crossprod(Xt, X, om_p)
          if (use_horseshoe) diag(P) <- diag(P) + hs_prec_mat[, ip]
          Pb <- prior_Pb[, ip] + Xt %*% target

          lo <- rep(-Inf, k)
          hi <- rep(Inf, k)
          if (!is.null(positive_constraints)) lo[positive_constraints] <- cur_lo[positive_constraints]
          if (!is.null(negative_constraints)) hi[negative_constraints] <- cur_hi[negative_constraints]

          curr_beta[, ip] <- sample_tmvn_precision_gibbs_cpp(
            P       = P,
            Pb      = Pb,
            lo      = lo,
            hi      = hi,
            init    = curr_beta[, ip],
            n_steps = if (iter <= nburn) 2L else 1L
          )
        }
      } else if (use_ncp) {
        # -----------------------------------------------------------------
        # C2. CONSTRAINED RE — NCP PATH
        # beta_c[, m] = mu + sigma * z_c[, m],  z_c ~ N(0, I)
        # Constraint beta >= L  =>  z >= (L - mu) / sigma
        # mu bound accumulated as max over groups of (L - sigma * z_g)
        # -----------------------------------------------------------------
        sigma_mat <- matrix(1 / sqrt(pmax(prec_beta_pooled, 1e-8)), nrow = k, ncol = p)


        for (ip in 1:p) {
          j <- pp[ip]
          om_p <- omega[, ip]
          c_j <- c_j_mat[, ip]
          sig <- sigma_mat[, ip]

          P_mu_acc <- matrix(0, k, k)
          Pb_mu_acc <- rep(0, k)

          # Initialise mu walls at the hard constraint bounds
          mu_lo <- rep(-Inf, k)
          mu_hi <- rep(Inf, k)
          if (!is.null(positive_constraints)) mu_lo[positive_constraints] <- cur_lo[positive_constraints]
          if (!is.null(negative_constraints)) mu_hi[negative_constraints] <- cur_hi[negative_constraints]

          for (m in 1:n_groups) {
            idx <- idx_list[[m]]
            if (length(idx) == 0) next
            Xm <- Xm_list[[m]]
            Xmt <- Xmt_list[[m]]
            om_m <- om_p[idx]

            # Localize RE scales with group-specific structural mask
            sig_m <- sig * re_mask[, m]

            # Base precision term for this group (used for both z and mu)
            Xmt_om_Xm <- Xmt %*% (Xm * om_m)

            # Working residual after removing mu contribution
            r_z <- kappa_weighted_iter[idx, j] + om_m * (c_j[idx] - Xm %*% mu_pooled[, ip])
            if (use_bart) r_z <- r_z - om_m * (bart_alpha * curr_f[idx, ip])

            # CORRECT (Issue 5, Doc 10): element-wise Hadamard, not matrix product.
            # NCP precision = I + diag(sig) %*% X'OmX %*% diag(sig)
            #               = I + outer(sig, sig) * X'OmX
            P_z <- diag(1, k) + outer(sig_m, sig_m) * Xmt_om_Xm
            Pb_z <- sig_m * as.vector(Xmt %*% r_z)

            # Correct bound transformation: z >= (L - mu) / sigma_m
            lo_z <- rep(-Inf, k)
            hi_z <- rep(Inf, k)
            if (!is.null(positive_constraints)) {
              raw_lo <- (cur_lo[positive_constraints] - mu_pooled[positive_constraints, ip]) /
                pmax(sig_m[positive_constraints], 1e-8)
              lo_z[positive_constraints] <- pmax(raw_lo, -1e6)
            }
            if (!is.null(negative_constraints)) {
              raw_hi <- (cur_hi[negative_constraints] - mu_pooled[negative_constraints, ip]) /
                pmax(sig_m[negative_constraints], 1e-8)
              hi_z[negative_constraints] <- pmin(raw_hi, 1e6)
            }

            z_c[, ip, m] <- sample_tmvn_precision_gibbs_cpp(
              P       = P_z,
              Pb      = Pb_z,
              lo      = lo_z,
              hi      = hi_z,
              init    = z_c[, ip, m],
              n_steps = if (iter <= nburn) 2L else 1L
            )

            # Recover beta_c and hard clip for one-iteration lag safety
            curr_beta_c[, ip, m] <- mu_pooled[, ip] + sig_m * z_c[, ip, m]
            if (!is.null(positive_constraints)) {
              curr_beta_c[positive_constraints, ip, m] <- pmax(
                curr_beta_c[positive_constraints, ip, m], cur_lo[positive_constraints]
              )
            }
            if (!is.null(negative_constraints)) {
              curr_beta_c[negative_constraints, ip, m] <- pmin(
                curr_beta_c[negative_constraints, ip, m], cur_hi[negative_constraints]
              )
            }

            # Explicitly pin fixed effects and masked random effects
            if (length(re_complement) > 0) {
              curr_beta_c[re_complement, ip, m] <- mu_pooled[re_complement, ip]
              z_c[re_complement, ip, m] <- 0
            }
            masked_re <- masked_by_group_cat[[m]][[ip]]
            if (length(masked_re) > 0) {
              curr_beta_c[masked_re, ip, m] <- mu_pooled[masked_re, ip]
              z_c[masked_re, ip, m] <- 0
            }

            # Accumulate tighter mu walls from this group's z draw
            if (!is.null(positive_constraints)) {
              mu_lo[positive_constraints] <- pmax(
                mu_lo[positive_constraints],
                cur_lo[positive_constraints] - sig_m[positive_constraints] * z_c[positive_constraints, ip, m]
              )
            }
            if (!is.null(negative_constraints)) {
              mu_hi[negative_constraints] <- pmin(
                mu_hi[negative_constraints],
                cur_hi[negative_constraints] - sig_m[negative_constraints] * z_c[negative_constraints, ip, m]
              )
            }

            # z-residualised sufficient stats for mu update
            r_mu <- kappa_weighted_iter[idx, j] + om_m * (c_j[idx] - Xm %*% (sig_m * z_c[, ip, m]))
            if (use_bart) r_mu <- r_mu - om_m * (bart_alpha * curr_f[idx, ip])

            P_mu_acc <- P_mu_acc + Xmt_om_Xm
            Pb_mu_acc <- Pb_mu_acc + as.vector(Xmt %*% r_mu)
          }

          # mu draw — uses walls accumulated from ALL groups
          P_mu <- prior_P + P_mu_acc
          diag(P_mu) <- diag(P_mu) + hs_prec_mat[, ip]

          mu_pooled[, ip] <- sample_tmvn_precision_gibbs_cpp(
            P       = P_mu + diag(1e-10, k),
            Pb      = prior_Pb[, ip] + Pb_mu_acc,
            lo      = mu_lo, # accumulated wall, not raw cur_lo
            hi      = mu_hi,
            init    = mu_pooled[, ip],
            n_steps = if (iter <= nburn) 3L else 2L
          )
        }
      } else {
        # -----------------------------------------------------------------
        # C3. CONSTRAINED RE — CP PATH (fallback when use_ncp = FALSE)
        # -----------------------------------------------------------------


        for (ip in 1:p) {
          j <- pp[ip]
          om_p <- omega[, ip]
          c_j <- c_j_mat[, ip]

          resid_Pb_sum <- rep(0, k)
          P_fixed_sum <- matrix(0, k, k)

          for (m in 1:n_groups) {
            idx <- idx_list[[m]]
            if (length(idx) == 0) next
            Xm <- Xm_list[[m]]
            Xmt <- Xmt_list[[m]]
            om_m <- om_p[idx]

            target_m <- kappa_weighted_iter[idx, j] + om_m * c_j[idx]
            if (use_bart) target_m <- target_m - om_m * (bart_alpha * curr_f[idx, ip])

            P_m <- weighted_crossprod(Xmt, Xm, om_m)
            P <- diag(prec_beta_pooled[, ip], k) + P_m
            Pb <- as.vector(prec_beta_pooled[, ip] * mu_pooled[, ip]) + Xmt %*% target_m

            lo <- rep(-Inf, k)
            hi <- rep(Inf, k)
            if (!is.null(positive_constraints)) lo[positive_constraints] <- cur_lo[positive_constraints]
            if (!is.null(negative_constraints)) hi[negative_constraints] <- cur_hi[negative_constraints]

            curr_beta_c[, ip, m] <- sample_tmvn_precision_gibbs_cpp(
              P       = P,
              Pb      = Pb,
              lo      = lo,
              hi      = hi,
              init    = curr_beta_c[, ip, m],
              n_steps = if (iter <= nburn) 2L else 1L
            )
            if (length(re_complement) > 0) curr_beta_c[re_complement, ip, m] <- mu_pooled[re_complement, ip]
            masked_re <- masked_by_group_cat[[m]][[ip]]
            if (length(masked_re) > 0) curr_beta_c[masked_re, ip, m] <- mu_pooled[masked_re, ip]

            # Residual-based accumulation for mu update
            P_fixed_sum <- P_fixed_sum + P_m
            delta_g <- curr_beta_c[, ip, m] - mu_pooled[, ip]
            resid_target <- target_m - (Xm %*% delta_g) * om_m
            resid_Pb_sum <- resid_Pb_sum + as.vector(Xmt %*% resid_target)
          }

          cur_prior_P <- prior_P
          if (use_horseshoe) {
            cur_prior_P[hs_idx, hs_idx] <- cur_prior_P[hs_idx, hs_idx] +
              diag(hs_prec_mat[hs_idx, ip], length(hs_idx))
          }

          P_mu <- cur_prior_P + P_fixed_sum + diag(1e-10, k)
          Pb_mu <- as.vector(cur_prior_P %*% prior_mu[, ip]) + resid_Pb_sum

          lo <- rep(-Inf, k)
          hi <- rep(Inf, k)
          if (!is.null(positive_constraints)) lo[positive_constraints] <- cur_lo[positive_constraints]
          if (!is.null(negative_constraints)) hi[negative_constraints] <- cur_hi[negative_constraints]

          mu_pooled[, ip] <- sample_tmvn_precision_gibbs_cpp(
            P       = P_mu,
            Pb      = Pb_mu,
            lo      = lo,
            hi      = hi,
            init    = mu_pooled[, ip],
            n_steps = if (iter <= nburn) 3L else 2L
          )
        }
      }
    }

    # =================================================================
    # C.5 ASIS INTERWEAVING (NCP -> CP -> NCP)
    # Improves mixing for mu by sampling it in the Centered Parameterization
    # =================================================================
    if (use_re && use_ncp) {
      for (ip in 1:p) {
        # 1. Gather sufficient stats from total effects
        prec_re <- pmax(prec_beta_pooled[, ip], 1e-12)

        ym <- y_mask[ip, ]
        re_ym <- re_mask

        # Scale only the non-intercept rows by the group-specific y_mask for this category
        if (any(!is_global_intercept)) {
          re_ym[!is_global_intercept, ] <- t(t(re_ym[!is_global_intercept, , drop = FALSE]) * ym)
        }

        # --- CRITICAL PATCH: inject per-group separation masks -----------
        # masked_by_group_cat contains variables pinned to mu by perfect
        # separation or collinearity. Counting them as active groups would
        # feed mu back into itself, inflating precision and locking chains.
        # Uses precomputed sep_mask_by_ip (Issue 3, Doc 10) — single
        # vectorized matrix multiply replaces O(n_groups) loop per iteration.
        re_ym <- re_ym * sep_mask_by_ip[[ip]]
        # ----------------------------------------------------------------

        n_active_vec_ip <- rowSums(re_ym)

        # Matrix conversion ensures 2D layout even for single group/category cases
        beta_ip_mat <- matrix(curr_beta_c[, ip, ], nrow = k, ncol = n_groups)
        beta_sum_active <- rowSums(beta_ip_mat * re_ym)

        # 2. Build full CP Precision and Pb using active groups
        P_cp <- prior_P
        diag(P_cp)[re_idx] <- diag(P_cp)[re_idx] + n_active_vec_ip[re_idx] * prec_re[re_idx]
        if (use_horseshoe) diag(P_cp) <- diag(P_cp) + hs_prec_mat[, ip]

        Pb_cp <- prior_Pb[, ip]
        Pb_cp[re_idx] <- Pb_cp[re_idx] + prec_re[re_idx] * beta_sum_active[re_idx]

        # 3. Condition on non-RE (Fixed) elements
        if (length(re_complement) > 0) {
          P_cond <- P_cp[re_idx, re_idx, drop = FALSE]
          Pb_cond <- Pb_cp[re_idx] - P_cp[re_idx, re_complement, drop = FALSE] %*% mu_pooled[re_complement, ip]
        } else {
          P_cond <- P_cp
          Pb_cond <- Pb_cp
        }

        # 4. Map Constraints to the subsetted RE indices
        lo_k <- rep(-Inf, k)
        hi_k <- rep(Inf, k)
        if (has_constraints) {
          if (!is.null(positive_constraints)) lo_k[positive_constraints] <- cur_lo[positive_constraints]
          if (!is.null(negative_constraints)) hi_k[negative_constraints] <- cur_hi[negative_constraints]
        }
        lo_cond <- lo_k[re_idx]
        hi_cond <- hi_k[re_idx]

        # 5. Resample mu for RE indices in CP space
        asis_n_steps <- if (nburn == 0L || iter > nburn) 1L else
                        if (iter <= max(1L, floor(nburn * 0.1))) 5L else
                        if (iter <= max(2L, floor(nburn * 0.5))) 3L else 2L

        mu_pooled[re_idx, ip] <- sample_tmvn_precision_gibbs_cpp(
          P       = P_cond,
          Pb      = Pb_cond,
          lo      = lo_cond,
          hi      = hi_cond,
          init    = mu_pooled[re_idx, ip],
          n_steps = asis_n_steps
        )

        # 6. Map deterministically back to NCP (z_c)
        sig <- 1 / sqrt(prec_re)
        for (m in 1:n_groups) {
          z_c[re_idx, ip, m] <- (curr_beta_c[re_idx, ip, m] - mu_pooled[re_idx, ip]) / sig[re_idx]

          # Apply structural masking to z_c and beta_c
          active_masked_vars <- masked_by_group_cat[[m]][[ip]]
          if (length(active_masked_vars) > 0) {
            z_c[active_masked_vars, ip, m] <- 0.0
            curr_beta_c[active_masked_vars, ip, m] <- mu_pooled[active_masked_vars, ip]
          }

          # Maintain fixed effect consistency
          if (length(re_complement) > 0) {
            z_c[re_complement, ip, m] <- 0
            curr_beta_c[re_complement, ip, m] <- mu_pooled[re_complement, ip]
          }
        }
      }
    }

    # =================================================================
    # E. BART UPDATE (stays in R — dbarts has its own C++)
    # NOTE: Placed BEFORE Sections D/D2 so that BART intercept absorption
    #       is reflected in the subsequent sigma/horseshoe updates.
    # =================================================================
    if (use_bart && iter > bart_warmup) {
      for (ip in 1:p) {
        j <- pp[ip]
        lp <- if (use_re) rowSums(X * t(curr_beta_c[, ip, group_idx])) else X %*% curr_beta[, ip]

        y_star <- (kappa_weighted_iter[, j] / omega[, ip]) + c_j_mat[, ip] - lp
        # Robustness: Cap y_star to prevent exploding intercepts during burn-in
        y_star <- pmax(pmin(y_star, 100), -100)

        # NOTE: BART sampler sees 'raw' residuals, its own scale is internal.
        # We temper BARTs CONTRIBUTION to the utility, not the sampler itself.

        bart_samplers[[ip]]$setResponse(as.numeric(y_star))
        bart_samplers[[ip]]$setWeights(as.numeric(omega[, ip]))
        curr_f[, ip] <- bart_samplers[[ip]]$run(0L, 1L)$train

        # Robustness: Only perform this after some burn-in to avoid exploding intercepts from initial y_star spikes
        if (length(int_idx) > 0 && iter > burn_buffer) {
          shift_val <- mean(curr_f[, ip]) # Capture the mean BEFORE centering
          curr_f[, ip] <- curr_f[, ip] - shift_val

          # Bug 1 Fix: Only absorb the portion that actually affected the utility
          actual_utility_shift <- bart_alpha * shift_val

          if (use_re) {
            curr_beta_c[int_idx, ip, ] <- curr_beta_c[int_idx, ip, ] + actual_utility_shift
            # Robust weighted average for intercept(s)
            if (length(int_idx) == 1) {
              mu_pooled[int_idx, ip] <- sum(weights * curr_beta_c[int_idx, ip, ])
            } else {
              # Handle multiple constant columns if they exist
              for (ii in seq_along(int_idx)) {
                mu_pooled[int_idx[ii], ip] <- sum(weights * curr_beta_c[int_idx[ii], ip, ])
              }
            }
          } else {
            curr_beta[int_idx, ip] <- curr_beta[int_idx, ip] + actual_utility_shift
          }
          # Bug B Fix: Trees are calibrated to shift_val (what left curr_f),
          # not actual_utility_shift (what entered beta). Do NOT use +=.
          bart_shifts[ip, 1] <- shift_val
        } else {
          shift_val <- 0
          bart_shifts[ip, 1] <- 0
        }
      }
    }

    # =================================================================
    # D. SIGMA UPDATE (RE only)  —  C++ when unconstrained
    # NOTE: Placed AFTER Section E so that BART intercept absorption is
    #       reflected in curr_beta_c / mu_pooled before computing RE precision.
    # =================================================================
    if (use_re) {
      # Always use reconstructed total effects (curr_beta_c) to update precision/shrinkage,
      # as z_c in NCP mode is standardized and doesn't reflect the actual scale of signal.
      target_for_prec <- curr_beta_c

      # Use 3D category-specific mask

      if (use_half_cauchy_re) {
        re_prec <- update_re_precision_hc(
          beta_c       = target_for_prec,
          mu_pooled    = mu_pooled,
          re_idx       = as.integer(re_idx - 1L),
          n_groups     = n_groups,
          re_mask      = re_mask_cube,
          y_mask       = y_mask,
          is_intercept = as.integer(is_global_intercept),
          prec_prev    = prec_beta_pooled,
          a_aux_prev   = a_re,
          re_scale_A   = re_scale_A
        )
        a_re <- re_prec$a_aux
      } else {
        re_prec <- update_re_precision(
          target_for_prec,
          mu_pooled,
          as.integer(re_idx - 1L),
          n_groups,
          re_mask_cube,
          y_mask,
          as.integer(is_global_intercept),
          prior_a_re,
          prior_b_re
        )
      }
      prec_beta_pooled <- re_prec$prec
      sigma_beta_pooled <- re_prec$sigma
    }

    # =================================================================
    # D2. HORSESHOE UPDATE
    # NOTE: Placed AFTER Section E for the same reason as D.
    # =================================================================
    if (use_horseshoe) {
      hs <- update_horseshoe(
        if (use_re) mu_pooled else curr_beta, hs_lambda2, hs_nu, hs_tau2, hs_xi,
        hs_c2, hs_zeta, k, p, hs_idx,
        estimate_c2 = estimate_c2, slab_df = slab_df, slab_s2 = slab_s2,
        equation_specific = equation_specific_hs
      )
      hs_lambda2 <- hs$lambda2
      hs_nu <- hs$nu
      hs_tau2 <- hs$tau2
      hs_xi <- hs$xi
      hs_c2 <- hs$c2
      hs_zeta <- hs$zeta

      hs_tau2_mat <- if (equation_specific_hs) matrix(hs_tau2, k_hs, p, byrow = TRUE) else hs_tau2
      hs_lam2_sub <- if (equation_specific_hs) hs_lambda2[hs_idx, , drop = FALSE] else hs_lambda2[hs_idx]
      var_eff <- (hs_c2 * hs_lam2_sub * hs_tau2_mat) /
        (hs_c2 + hs_lam2_sub * hs_tau2_mat)
      hs_prec_mat[hs_idx, ] <- 1 / pmax(var_eff, 1e-12)
      if (support_prior_strength > 0) hs_prec_mat[hs_idx, ] <- hs_prec_mat[hs_idx, ] * fe_support[hs_idx]   # FE support-aware shrinkage
    }

    # =================================================================
    # F. STORE POSTERIOR  —  log-likelihood via C++
    # =================================================================
    if (iter > nburn) {
      s <- iter - nburn

      # Compute log-likelihood (needed for both RAM and Disk)
      ll <- compute_loglik(U, Y, as.numeric(y_weight))

      if (!save_posterior_to_disk) {
        if (use_re) {
          postb_total[, , , s] <- curr_beta_c
          postb_pooled[, , s] <- mu_pooled
          post_sigma_re[, s] <- as.vector(sigma_beta_pooled)
        } else {
          postb_total[, , s] <- curr_beta
          postb_pooled[, , s] <- curr_beta
        }

        if (use_horseshoe) {
          # Level 1: Pooled Mean Kappa
          hs_tau2_mat_store <- if (equation_specific_hs) matrix(hs_tau2, k_hs, p, byrow = TRUE) else hs_tau2
          hs_lam2_store <- if (equation_specific_hs) hs_lambda2[hs_idx, , drop = FALSE] else hs_lambda2[hs_idx]
          kappa_pooled <- 1 / (1 + hs_lam2_store * hs_tau2_mat_store)
          post_kappa_pooled[, s] <- as.vector(kappa_pooled)
          post_c2[s] <- hs_c2
        }
        post_log_lik[s] <- ll$total
        if (calc_loo) post_ll_pw[s, ] <- ll$pointwise
      }
      if (use_bart) {
        post_f_sum <- post_f_sum + bart_alpha * curr_f
        post_f_sum_sq <- post_f_sum_sq + (bart_alpha * curr_f)^2
        if (store_f) post_f[, , s] <- bart_alpha * curr_f
        if (store_bart_trees || save_bart_to_disk || save_posterior_to_disk) {
          state_m <- vector("list", p)
          for (ip in 1:p) {
            if (do_slim_trees) {
              # SLIM TREES: Extract structures and calibrate
              df <- bart_samplers[[ip]]$getTrees()

              # 1. CALIBRATION: Recover hidden dbarts Y-scaling and offset
              # (Uses standardized trees and empirical reference points)
              s1 <- as.numeric(bart_samplers[[ip]]$predict(bart_ref1_mat))
              s2 <- as.numeric(bart_samplers[[ip]]$predict(bart_ref2_mat))
              t1 <- as.numeric(predict_slim_bart_cpp(bart_ref1_mat, df))
              t2 <- as.numeric(predict_slim_bart_cpp(bart_ref2_mat, df))

              t_diff <- t1 - t2
              y_scale <- if (abs(t_diff) > 1e-10) (s1 - s2) / t_diff else 1.0

              y_offset <- s1 - (t1 * y_scale)

              # 2. TRANSFORM SPLIT POINTS: Move to raw scale (if standardized)
              if (!is.null(bart_scaling)) {
                split_rows <- which(df$var > 0)
                if (length(split_rows) > 0) {

                  v_indices <- df$var[split_rows]
                  df$value[split_rows] <- df$value[split_rows] * bart_scaling$sd[v_indices] + bart_scaling$mu[v_indices]

                }
              }

              # 3. APPLY CALIBRATION: Scale leaves and bake in offset
              df$value[df$var == -1] <- df$value[df$var == -1] * y_scale
              total_offset <- y_offset - as.numeric(bart_shifts[ip, 1])
              df$value[df$tree == 1 & df$var == -1] <- df$value[df$tree == 1 & df$var == -1] + total_offset

              # LIVE CHECK: test_pred should match target_pred at ref1
              test_pred <- as.numeric(predict_slim_bart_cpp(bart_ref1_mat, df))
              target_pred <- s1 - as.numeric(bart_shifts[ip, 1])

              if (abs(test_pred - target_pred) > 1e-5) {
                warning(sprintf("Calibration mismatch in sample %d cat %d: Slim=%.4f, Target=%.4f", s, ip, test_pred, target_pred))
              }

              state_m[[ip]] <- df
            } else {
              # FULL STORAGE: Serialized sampler state
              invisible(bart_samplers[[ip]]$state)
              state_m[[ip]] <- serialize(bart_samplers[[ip]], NULL)
            }
          }

          if (save_bart_to_disk && !save_posterior_to_disk) {
            # Add to buffer and write in batches
            bart_batch_buffer[[length(bart_batch_buffer) + 1]] <- state_m
            if (length(bart_batch_buffer) >= bart_batch_size || s == nretain) {
              b_idx <- length(bart_batch_files) + 1
              c_label <- if (is.null(chain_id)) "0" else as.character(chain_id)
              f_name <- file.path(bart_disk_path, sprintf("bart_batch_%d_chain_%s.qs", b_idx, c_label))
              qs2::qs_save(bart_batch_buffer, f_name, compress_level = 1)
              bart_batch_files[b_idx] <- f_name
              # Store the filename in tree_store for the first sample of this batch
              tree_store[[s - length(bart_batch_buffer) + 1]] <- f_name
              bart_batch_buffer <- list() # Clear RAM
            }
          } else if (!save_posterior_to_disk) {
            # RAM Storage
            tree_store[[s]] <- state_m
          }
        }
      }

      if (save_posterior_to_disk) {
        # Construct the whole state sample
        state_sample <- list(
          iter = iter,
          beta = if (use_re) curr_beta_c else curr_beta,
          mu = if (use_re) mu_pooled else curr_beta,
          sigma_re = if (use_re) sigma_beta_pooled else NULL,
          log_lik = ll$total,
          log_lik_pw = if (calc_loo) ll$pointwise else NULL,
          horseshoe = if (use_horseshoe) {
            list(
              lambda2 = hs_lambda2, tau2 = hs_tau2,
              c2 = hs_c2,
              nu = hs_nu, xi = hs_xi, zeta = hs_zeta
            )
          } else {
            NULL
          },
          bart_trees = if (use_bart) state_m else NULL
        )

        # BACK-TRANSFORMATION (Single Sample)
        if (!("none" %in% methods_chosen)) {
          # 1. Reverse QR
          if (do_qr) {
            if (use_re) {
              for (m in 1:n_groups) {
                state_sample$beta[, , m] <- backsolve(X_R_mat, state_sample$beta[, , m, drop = FALSE])
              }
              state_sample$mu <- backsolve(X_R_mat, state_sample$mu, drop = FALSE)
            } else {
              state_sample$beta <- backsolve(X_R_mat, state_sample$beta)
              state_sample$mu <- state_sample$beta
            }
          }
          # 2. Reverse Standardization (uses precomputed bt_is_cont / bt_int_idx
          #    from line ~579 to avoid shadowing the BART int_idx on line 312)
          if ((do_cen || do_scl) && sum(cont_idx) > 0) {
            if (use_re) {
              if (do_scl) {
                for (ci in seq_along(bt_is_cont)) {
                  state_sample$beta[bt_is_cont[ci], , ] <- state_sample$beta[bt_is_cont[ci], , ] / X_rescaling[2, ci]
                  state_sample$mu[bt_is_cont[ci], ] <- state_sample$mu[bt_is_cont[ci], ] / X_rescaling[2, ci]
                }
              }
              if (do_cen && length(bt_int_idx) > 0) {
                for (m in 1:n_groups) {
                  shift <- colSums(state_sample$beta[bt_is_cont, , m, drop = FALSE] * X_rescaling[1, ])
                  for (ii in bt_int_idx) state_sample$beta[ii, , m] <- state_sample$beta[ii, , m] - shift
                }
                shift_p <- colSums(state_sample$mu[bt_is_cont, , drop = FALSE] * X_rescaling[1, ])
                for (ii in bt_int_idx) state_sample$mu[ii, ] <- state_sample$mu[ii, ] - shift_p
              }
            } else {
              if (do_scl) {
                for (ci in seq_along(bt_is_cont)) state_sample$beta[bt_is_cont[ci], ] <- state_sample$beta[bt_is_cont[ci], ] / X_rescaling[2, ci]
              }
              if (do_cen && length(bt_int_idx) > 0) {
                shift <- colSums(state_sample$beta[bt_is_cont, , drop = FALSE] * X_rescaling[1, ])
                for (ii in bt_int_idx) state_sample$beta[ii, ] <- state_sample$beta[ii, ] - shift
              }
              state_sample$mu <- state_sample$beta
            }
          }
        }

        # Add to buffer
        posterior_batch_buffer[[length(posterior_batch_buffer) + 1]] <- state_sample
        if (length(posterior_batch_buffer) >= posterior_batch_size || s == nretain) {
          b_idx <- length(posterior_batch_files) + 1
          c_label <- if (is.null(chain_id)) "0" else as.character(chain_id)
          f_name <- file.path(posterior_disk_path, sprintf("posterior_batch_%d_chain_%s.qs", b_idx, c_label))
          qs2::qs_save(posterior_batch_buffer, f_name, compress_level = 1)
          posterior_batch_files[b_idx] <- f_name
          posterior_batch_buffer <- list() # Clear RAM
        }
      }

      if (!save_posterior_to_disk) {
        # ll was already computed and stored at the top of the loop (F. STORE POSTERIOR)
      }
    }

    if (!is.null(progress_cb)) {
      if (iter %% 10 == 0 || iter == niter) {
        phase <- if (iter <= nburn) "[Burn-in]" else "[Sampling]"
        p_msg <- sprintf("Chain %s: Iteration %d / %d %s", ifelse(is.null(chain_id), "?", chain_id), iter, niter, phase)
        progress_cb(message = p_msg)
      } else {
        progress_cb() # advance without changing message
      }
    } else if (is.null(chain_id)) {
      utils::setTxtProgressBar(pb, iter)
    } else if (iter %% 100 == 0 || iter == niter) {
      phase <- if (iter <= nburn) "[Burn-in]" else "[Sampling]"
      cat(sprintf("Chain %d: Iteration %d / %d %s\n", chain_id, iter, niter, phase))
    }
  }

  # --- Flush remaining BART ensembles ---
  if (save_bart_to_disk && use_bart && !save_posterior_to_disk && length(bart_batch_buffer) > 0) {
    b_idx <- length(bart_batch_files) + 1
    c_label <- if (is.null(chain_id)) "0" else as.character(chain_id)
    f_name <- file.path(bart_disk_path, sprintf("bart_batch_%d_chain_%s.qs", b_idx, c_label))
    qs2::qs_save(bart_batch_buffer, f_name, compress_level = 1)
    bart_batch_files[b_idx] <- f_name
    bart_batch_buffer <- list()
  }

  # --- Flush remaining Posterior samples ---
  if (save_posterior_to_disk && length(posterior_batch_buffer) > 0) {
    b_idx <- length(posterior_batch_files) + 1
    c_label <- if (is.null(chain_id)) "0" else as.character(chain_id)
    f_name <- file.path(posterior_disk_path, sprintf("posterior_batch_%d_chain_%s.qs", b_idx, c_label))
    qs2::qs_save(posterior_batch_buffer, f_name, compress_level = 1)
    posterior_batch_files[b_idx] <- f_name
    posterior_batch_buffer <- list()
  }

  if (is.null(chain_id)) {
    close(pb)
    cat("\nSampling finished.\n")
  } else {
    cat(sprintf("Chain %d: Sampling finished.\n", chain_id))
  }

  # --- Capture and Save Final State for Hot-Start ---
  horseshoe_state <- NULL
  if (use_horseshoe) {
    horseshoe_state <- list(
      lambda2 = hs_lambda2, tau2 = hs_tau2, nu = hs_nu, xi = hs_xi, zeta = hs_zeta,
      c2 = hs_c2
    )
  }

  final_state_raw <- list(
    beta = if (use_re) curr_beta_c else curr_beta,
    mu = if (use_re) mu_pooled else NULL,
    prec_beta = if (use_re) prec_beta_pooled else NULL,
    z_c = if (use_re && use_ncp) z_c else NULL,
    sigma_re = if (use_re) sigma_beta_pooled else NULL,
    a_re = if (use_re && use_half_cauchy_re) a_re else NULL,
    horseshoe = horseshoe_state
  )
  if (use_bart) {
    final_state_raw$bart_states <- lapply(bart_samplers, function(s) s$state)
  }

  if (save_posterior_to_disk) {
    c_label <- if (is.null(chain_id)) "1" else as.character(chain_id)
    f_name_fs <- file.path(posterior_disk_path, sprintf("final_state_chain_%s.qs", c_label))
    qs2::qs_save(final_state_raw, f_name_fs, compress_level = 1)
  }

  # --- 5. BACK-TRANSFORMATION (Only if not already done via Disk path) ---
  if (!save_posterior_to_disk && !("none" %in% methods_chosen)) {
    if (is.null(chain_id)) cat("Back-transforming...\n") else cat(sprintf("Chain %d: Back-transforming...\n", chain_id))

    # 1. Reverse QR FIRST
    if (do_qr) {
      if (use_re) {
        # Vectorized QR backsolve across all groups and samples by flattening extra dimensions
        dim_total <- dim(postb_total)
        dim(postb_total) <- c(dim_total[1], prod(dim_total[-1]))
        postb_total <- backsolve(X_R_mat, postb_total)
        dim(postb_total) <- dim_total

        dim_pooled <- dim(postb_pooled)
        dim(postb_pooled) <- c(dim_pooled[1], prod(dim_pooled[-1]))
        postb_pooled <- backsolve(X_R_mat, postb_pooled)
        dim(postb_pooled) <- dim_pooled
      } else {
        dim_total <- dim(postb_total)
        dim(postb_total) <- c(dim_total[1], prod(dim_total[-1]))
        postb_total <- backsolve(X_R_mat, postb_total)
        dim(postb_total) <- dim_total
        postb_pooled <- postb_total
      }
    }

    # Save standardized parameters before unscaling
    postb_total_std <- postb_total
    postb_pooled_std <- postb_pooled

    # 2. Reverse Standardization LATER
    if ((do_cen || do_scl) && sum(cont_idx) > 0) {
      is_cont <- which(cont_idx)
      int_idx <- which(apply(X_orig, 2, var) == 0) # Find constant columns (usually intercept)

      if (use_re) {
        if (do_scl) {
          for (ci in seq_along(is_cont)) {
            postb_total[is_cont[ci], , , ] <- postb_total[is_cont[ci], , , ] / X_rescaling[2, ci]
            postb_pooled[is_cont[ci], , ] <- postb_pooled[is_cont[ci], , ] / X_rescaling[2, ci]
          }
        }
        if (do_cen && length(int_idx) > 0) {
          # Fully vectorized shift subtraction
          shift <- colSums(postb_total[is_cont, , , , drop = FALSE] * X_rescaling[1, ], dims = 1)
          for (ii in int_idx) {
            target_dim <- dim(postb_total[ii, , , ])
            if (is.null(target_dim)) target_dim <- length(postb_total[ii, , , ])
            postb_total[ii, , , ] <- postb_total[ii, , , ] - array(shift, dim = target_dim)
          }
          shift_p <- colSums(postb_pooled[is_cont, , , drop = FALSE] * X_rescaling[1, ], dims = 1)
          for (ii in int_idx) {
            target_dim <- dim(postb_pooled[ii, , ])
            if (is.null(target_dim)) target_dim <- length(postb_pooled[ii, , ])
            postb_pooled[ii, , ] <- postb_pooled[ii, , ] - array(shift_p, dim = target_dim)
          }
        }
      } else {
        if (do_scl) {
          for (ci in seq_along(is_cont)) {
            postb_total[is_cont[ci], , ] <- postb_total[is_cont[ci], , ] / X_rescaling[2, ci]
          }
        }
        if (do_cen && length(int_idx) > 0) {
          shift <- colSums(postb_total[is_cont, , , drop = FALSE] * X_rescaling[1, ], dims = 1)
          for (ii in int_idx) {
            target_dim <- dim(postb_total[ii, , ])
            if (is.null(target_dim)) target_dim <- length(postb_total[ii, , ])
            postb_total[ii, , ] <- postb_total[ii, , ] - array(shift, dim = target_dim)
          }
        }
        postb_pooled <- postb_total
      }
    }
  }

  # --- 6. ATTACH NAMES ---
  cov_names <- cov_names_save
  cat_names <- colnames(Y)[pp]

  if (!save_posterior_to_disk) {
    if (use_re) {
      g_names <- as.character(seq_len(n_groups))
      dimnames(postb_total) <- list(cov_names, cat_names, g_names, NULL)
      dimnames(postb_pooled) <- list(cov_names, cat_names, NULL)
      if (!is.null(postb_total_std)) {
        dimnames(postb_total_std) <- list(cov_names, cat_names, g_names, NULL)
        dimnames(postb_pooled_std) <- list(cov_names, cat_names, NULL)
      }
    } else {
      dimnames(postb_total) <- list(cov_names, cat_names, NULL)
      dimnames(postb_pooled) <- list(cov_names, cat_names, NULL)
      if (!is.null(postb_total_std)) {
        dimnames(postb_total_std) <- list(cov_names, cat_names, NULL)
        dimnames(postb_pooled_std) <- list(cov_names, cat_names, NULL)
      }
    }
  }

  # --- 7. SUMMARY & DIAGNOSTICS ---
  diagnostics <- list(
    loglik = mean(post_log_lik)
  )
  res <- list(
    postb = postb_total, postb_total = postb_total, postb_pooled = postb_pooled,
    postb_total_std = postb_total_std,
    postb_pooled_std = postb_pooled_std,
    post_sigma_re = post_sigma_re,
    sigma_beta_pooled = sigma_beta_pooled,
    post_log_lik = post_log_lik, post_log_lik_pointwise = post_ll_pw,
    diagnostics = diagnostics, baseline = baseline,
    tree_store = if (save_posterior_to_disk) "Integrated in posterior_store" else (if (save_bart_to_disk) bart_batch_files else tree_store),
    posterior_store = if (save_posterior_to_disk) posterior_batch_files else NULL,
    bart_batch_size = if (save_bart_to_disk || save_posterior_to_disk) bart_batch_size else NULL,
    post_f_mean = if (use_bart) post_f_sum / nretain else NULL,
    post_f_sd = if (use_bart) sqrt(pmax(post_f_sum_sq / nretain - (post_f_sum / nretain)^2, 0)) else NULL,
    post_f = if (use_bart && store_f) post_f else NULL,
    bart_scaling = if (use_bart & !do_slim_trees) bart_scaling else NULL,
    bart_idx = if (use_bart & !do_slim_trees) bart_idx else NULL
  )

  # --- Capture Final State for Hot-Start ---
  res$final_state <- final_state_raw

  if (calc_loo && !save_posterior_to_disk) {
    res$waic <- loo::waic(post_ll_pw)
    res$loo <- loo::loo(post_ll_pw)
  }
  if (use_horseshoe) {
    res$horseshoe <- list(
      lambda2 = hs_lambda2, tau2 = hs_tau2,
      c2 = hs_c2,
      post_kappa_pooled = post_kappa_pooled,
      post_c2 = post_c2,
      nu = hs_nu, xi = hs_xi, zeta = hs_zeta
    )
  }
  return(res)
}

# =============================================================================
# Recovery Utility: Reconstructs the full RAM-based posterior list from disk batches
# =============================================================================
recover_mnlogit_posterior <- function(path_or_files, chain_id = NULL) {
  # 1. Handle Input: File List or Directory Path
  if (length(path_or_files) == 1 && dir.exists(path_or_files)) {
    cat(sprintf("Scanning directory for posterior batches: %s\n", path_or_files))

    # Define pattern: if chain_id is provided, filter for it
    pattern <- if (is.null(chain_id)) {
      "posterior_batch_.*\\.qs$"
    } else {
      sprintf("posterior_batch_.*_chain_%s\\.qs$", as.character(chain_id))
    }

    posterior_files <- list.files(path_or_files, pattern = pattern, full.names = TRUE)

    # Sort files numerically to ensure iteration order
    # (Extracts the batch number from 'posterior_batch_X_...')
    batch_nums <- as.integer(gsub(".*batch_([0-9]+)_.*", "\\1", basename(posterior_files)))
    posterior_files <- posterior_files[order(batch_nums)]
  } else {
    posterior_files <- path_or_files
  }

  if (is.null(posterior_files) || length(posterior_files) == 0) {
    if (!is.null(chain_id)) {
      message(sprintf("No posterior files found for chain %s.", chain_id))
      return(NULL)
    }
    stop("No posterior files provided or found.")
  }

  cat(sprintf("Recovering posterior from %d batch files...\n", length(posterior_files)))

  # 2. Load all batches into a flat list of samples
  all_samples <- list()
  for (f in posterior_files) {
    if (!file.exists(f)) {
      warning("File not found, skipping: ", f)
      next
    }
    batch_data <- qs2::qs_read(f)
    all_samples <- c(all_samples, batch_data)
  }

  nretain <- length(all_samples)
  if (nretain == 0) {
    return(NULL)
  }

  # 2. Inspect first sample to determine dimensions and features
  s1 <- all_samples[[1]]
  beta_dim <- dim(s1$beta)
  k <- beta_dim[1]
  p <- beta_dim[2]
  use_re <- length(beta_dim) == 3
  n_groups <- if (use_re) beta_dim[3] else NULL

  has_sigma_re <- !is.null(s1$sigma_re)
  has_loo <- !is.null(s1$log_lik_pw)
  has_hs <- !is.null(s1$horseshoe)
  has_bart <- !is.null(s1$bart_trees)

  cat(sprintf(
    "Dimensions detected: %d covariates, %d categories, %s\n",
    k, p, ifelse(use_re, paste0(n_groups, " groups"), "pooled model")
  ))

  # 3. Pre-allocate RAM storage
  res_list <- list()
  res_list$postb_total <- if (use_re) array(0, c(k, p, n_groups, nretain)) else array(0, c(k, p, nretain))
  res_list$postb_pooled <- array(0, c(k, p, nretain))
  res_list$post_log_lik <- numeric(nretain)

  if (has_sigma_re) res_list$post_sigma_re <- matrix(0, length(s1$sigma_re), nretain)
  if (has_loo) res_list$post_log_lik_pointwise <- matrix(0, nretain, length(s1$log_lik_pw))
  if (has_bart) res_list$tree_store <- vector("list", nretain)

  if (has_hs) {
    k_hs <- dim(s1$horseshoe$lambda2)[1] %||% length(s1$horseshoe$lambda2)
    res_list$horseshoe <- list(
      post_kappa_pooled = matrix(0, k_hs * p, nretain),
      post_kappa_re     = if (use_re && !is.null(s1$horseshoe$phi2)) matrix(0, k_hs * p, nretain) else NULL,
      post_c2           = numeric(nretain),
      post_c2_re        = if (use_re && !is.null(s1$horseshoe$c2_re)) numeric(nretain) else NULL
    )
  }

  # 4. Fill storage arrays
  for (s in 1:nretain) {
    samp <- all_samples[[s]]

    if (use_re) {
      res_list$postb_total[, , , s] <- samp$beta
      res_list$postb_pooled[, , s] <- samp$mu
    } else {
      res_list$postb_total[, , s] <- samp$beta
      res_list$postb_pooled[, , s] <- samp$mu
    }

    res_list$post_log_lik[s] <- samp$log_lik

    if (has_sigma_re) res_list$post_sigma_re[, s] <- as.vector(samp$sigma_re)
    if (has_loo) res_list$post_log_lik_pointwise[s, ] <- samp$log_lik_pw
    if (has_bart) res_list$tree_store[[s]] <- samp$bart_trees

    if (has_hs) {
      k_hs <- dim(samp$mu)[1]
      p_hs <- dim(samp$mu)[2]

      # Level 1: Finnish-regularized pooled kappa (matches RAM-path formula)
      hs_c2_rec <- samp$horseshoe$c2 %||% 15
      hs_tau2_mat_rec <- if (length(samp$horseshoe$tau2) > 1) matrix(samp$horseshoe$tau2, k_hs, p_hs, byrow = TRUE) else samp$horseshoe$tau2
      lam2_tau2 <- samp$horseshoe$lambda2 * hs_tau2_mat_rec
      var_eff_rec <- (hs_c2_rec * lam2_tau2) / (hs_c2_rec + lam2_tau2)
      res_list$horseshoe$post_kappa_pooled[, s] <- as.vector(1 - var_eff_rec / lam2_tau2)
      res_list$horseshoe$post_c2[s] <- hs_c2_rec

      if (!is.null(res_list$horseshoe$post_kappa_re)) {
        # Level 2: Finnish-regularized RE deviation kappa
        hs_c2_re_rec <- samp$horseshoe$c2_re %||% 3
        hs_psi2_mat_rec <- if (length(samp$horseshoe$psi2) > 1) matrix(samp$horseshoe$psi2, k_hs, p_hs, byrow = TRUE) else samp$horseshoe$psi2
        phi2_psi2 <- samp$horseshoe$phi2 * hs_psi2_mat_rec
        var_re_rec <- (hs_c2_re_rec * phi2_psi2) / (hs_c2_re_rec + phi2_psi2)
        res_list$horseshoe$post_kappa_re[, s] <- as.vector(1 - var_re_rec / phi2_psi2)
        res_list$horseshoe$post_c2_re[s] <- hs_c2_re_rec
      }
    }
  }

  # 5. Restore metadata and naming
  base_dir <- if (dir.exists(path_or_files[1])) path_or_files[1] else dirname(path_or_files[1])
  meta_file <- file.path(base_dir, "model_metadata.qs")
  meta <- if (length(meta_file) == 1 && file.exists(meta_file)) qs2::qs_read(meta_file) else NULL

  if (!is.null(meta)) {
    # Preferred: Use explicitly saved metadata
    cov_names <- meta$cov_names
    cat_names <- meta$cat_names
  } else {
    # Fallback: Try to get from first sample
    cov_names <- dimnames(s1$beta)[[1]] %||% paste0("V", k)
    cat_names <- dimnames(s1$beta)[[2]] %||% paste0("C", p)
  }

  if (use_re) {
    dimnames(res_list$postb_total) <- list(cov_names, cat_names, NULL, NULL)
    dimnames(res_list$postb_pooled) <- list(cov_names, cat_names, NULL)
  } else {
    dimnames(res_list$postb_total) <- list(cov_names, cat_names, NULL)
    dimnames(res_list$postb_pooled) <- list(cov_names, cat_names, NULL)
  }

  res_list$postb <- res_list$postb_total
  res_list$var_names <- cov_names
  res_list$cat_names <- cat_names

  # 6. Apply back-transformation if metadata contains rescaling info
  if (!is.null(meta) && !is.null(meta$X_rescaling)) {
    cat("Back-transforming posterior to physical scale based on metadata...\n")
    X_rescaling <- meta$X_rescaling
    cont_idx <- meta$cont_idx
    do_cen <- meta$do_cen
    do_scl <- meta$do_scl

    if ((do_cen || do_scl) && sum(cont_idx) > 0) {
      is_cont <- which(cont_idx)
      int_idx <- which(cov_names == "intercept")

      if (use_re) {
        if (do_scl) {
          for (ci in seq_along(is_cont)) {
            res_list$postb_total[is_cont[ci], , , ] <- res_list$postb_total[is_cont[ci], , , ] / X_rescaling[2, ci]
            res_list$postb_pooled[is_cont[ci], , ] <- res_list$postb_pooled[is_cont[ci], , ] / X_rescaling[2, ci]
          }
        }
        if (do_cen && length(int_idx) > 0) {
          # Vectorized shift: sum across all continuous variables in one shot
          # (matches main sampler approach — clearer and faster than sequential per-variable loop)
          shift <- colSums(res_list$postb_total[is_cont, , , , drop = FALSE] * X_rescaling[1, ], dims = 1)
          for (ii in int_idx) {
            target_dim <- dim(res_list$postb_total[ii, , , ])
            if (is.null(target_dim)) target_dim <- length(res_list$postb_total[ii, , , ])
            res_list$postb_total[ii, , , ] <- res_list$postb_total[ii, , , ] - array(shift, dim = target_dim)
          }
          shift_p <- colSums(res_list$postb_pooled[is_cont, , , drop = FALSE] * X_rescaling[1, ], dims = 1)
          for (ii in int_idx) {
            target_dim <- dim(res_list$postb_pooled[ii, , ])
            if (is.null(target_dim)) target_dim <- length(res_list$postb_pooled[ii, , ])
            res_list$postb_pooled[ii, , ] <- res_list$postb_pooled[ii, , ] - array(shift_p, dim = target_dim)
          }
        }
      } else {
        if (do_scl) {
          for (ci in seq_along(is_cont)) res_list$postb_total[is_cont[ci], , ] <- res_list$postb_total[is_cont[ci], , ] / X_rescaling[2, ci]
        }
        if (do_cen && length(int_idx) > 0) {
          shift <- colSums(res_list$postb_total[is_cont, , , drop = FALSE] * X_rescaling[1, ], dims = 1)
          for (ii in int_idx) {
            target_dim <- dim(res_list$postb_total[ii, , ])
            if (is.null(target_dim)) target_dim <- length(res_list$postb_total[ii, , ])
            res_list$postb_total[ii, , ] <- res_list$postb_total[ii, , ] - array(shift, dim = target_dim)
          }
        }
        res_list$postb_pooled <- res_list$postb_total
      }
    }
  }

  # NEW: Recover final state if it exists (for hot-starting)
  fs_file <- file.path(base_dir, sprintf("final_state_chain_%s.qs", if (is.null(chain_id)) "1" else as.character(chain_id)))
  if (file.exists(fs_file)) {
    res_list$final_state <- qs2::qs_read(fs_file)
  }

  cat("Recovery complete.\n")
  return(res_list)
}

# --- Multi-Chain Recovery Helper ---
# Reconstructs a list of chain results from a directory
recover_mnlogit_chains <- function(path) {
  if (!dir.exists(path)) stop("Directory does not exist: ", path)

  # Identify all unique chain IDs in the directory
  files <- list.files(path, pattern = "posterior_batch_.*_chain_.*\\.qs$")
  if (length(files) == 0) stop("No chain batches found in ", path)

  chain_ids <- unique(gsub(".*_chain_([0-9]+)\\.qs$", "\\1", files))
  cat(sprintf("Found %d chains in directory: %s\n", length(chain_ids), paste(chain_ids, collapse = ", ")))

  res_chains <- lapply(chain_ids, function(id) {
    recover_mnlogit_posterior(path, chain_id = id)
  })

  names(res_chains) <- paste0("chain_", chain_ids)
  return(res_chains)
}
