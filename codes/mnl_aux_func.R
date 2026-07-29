# =============================================================================
# mnl_aux_func.R — Shared auxiliary functions for Bayesian MNL Samplers
#
# Canonical source for:
#   - fast_rpg()           : Bias-corrected Polya-Gamma sampler
#   - rig()                : Safe inverse-gamma draw
#   - update_horseshoe()   : Finnish Horseshoe (pooled)
#
# Source this file at the top of mnlogit.R, mnlogit_fast.R, mnlogit_rcpp.R.
# =============================================================================
library(pg) # For high-performance Polya-Gamma sampling (hybrid exact/approx)
library(Matrix) # For nearPD in ensure_pd

# --- Null-coalescing operator ---
`%||%` <- function(a, b) if (!is.null(a)) a else b
# --- Bias-Corrected PG Sampler ---
# Normal approx for PG(h, z): mean = h/(2z)*tanh(z/2),
#                              var  = h*(sinh(z)-z)/(4z^3*cosh^2(z/2))
#
# pg_approx_thresh controls the boundary between exact and approximate draws:
#   - Default 170: strict mathematical safety (KL ~ 0)
#   - Set to 1.5 for speed-optimized samplers (KL < 0.005 nats for h >= 1.5)
# Optimized Polya-Gamma Sampler using the 'pg' package
# Handles fractional weights (h) and provides hybrid exact/approx logic for performance.
fast_rpg <- function(n, h, z, gamma_matched_pg = FALSE) {
  # Clean inputs
  z[!is.finite(z)] <- 0
  h <- as.vector(h)

  if (gamma_matched_pg) {
    # Moment-matched Gamma approximation: Captures mean, variance, and skewness
    # Handle z=0 safely (pg_mean/pg_var return NaN at 0)
    mu <- ifelse(abs(z) < 1e-10, h / 4, pg_mean(h, z))
    v <- ifelse(abs(z) < 1e-10, h / 24, pg_var(h, z))

    # Moment matching for Gamma(shape, rate):
    # rate = mu / v, shape = mu^2 / v
    rate <- mu / v
    shape <- (mu^2) / v

    res <- rgamma(length(mu), shape = shape, rate = rate)
    return(pmax(as.vector(res), 1e-12))
  }

  # pg::rpg_hybrid is vectorized and handles fractional h natively.
  # It is ~200x faster for large h (aggregated data).
  return(pmax(as.vector(rpg_hybrid(h, z)), 1e-12))
}

# --- Cholesky-based sampler from precision parameterization ---
# Replaces: V <- safe_invert(V_inv); beta <- mvrnorm(1, V %*% rhs, V)
# With:     L <- chol(V_inv); beta <- solve(L, solve(t(L), rhs) + rnorm(k))
# This avoids the explicit matrix inverse entirely.
chol_sample_precision <- function(P, Pb, k) {
  L <- tryCatch(chol(P), error = function(e) NULL)
  if (is.null(L)) {
    L <- tryCatch(chol(P + diag(1e-6, k)), error = function(e) NULL)
    if (is.null(L)) {
      V <- MASS::ginv(P)
      return(MASS::mvrnorm(1, V %*% Pb, (V + t(V)) / 2))
    }
  }
  mu <- backsolve(L, forwardsolve(t(L), Pb))
  mu + backsolve(L, rnorm(k))
}

# --- Truncated version for sign-constrained parameters ---
chol_mean_and_var <- function(P, Pb, k) {
  L <- tryCatch(chol(P), error = function(e) NULL)
  if (is.null(L)) {
    L <- tryCatch(chol(P + diag(1e-6, k)), error = function(e) NULL)
    if (is.null(L)) {
      V <- MASS::ginv(P)
      V <- ensure_pd(V)
      return(list(mu = as.vector(V %*% Pb), V = V))
    }
  }
  mu <- backsolve(L, forwardsolve(t(L), Pb))
  L_inv <- backsolve(L, diag(k))
  V <- ensure_pd(tcrossprod(L_inv))
  list(mu = as.vector(mu), V = V)
}

# --- Safe Inverse-Gamma draw (for Horseshoe updates) ---
rig <- function(shape, rate) {
  x <- rgamma(1, shape = shape, rate = rate)
  if (!is.finite(x) || x <= 0) {
    return(1)
  }
  1 / x
}

# --- Finnish Horseshoe Update: Pooled with Estimation Toggle ---
update_horseshoe <- function(curr_beta, hs_lambda2, hs_nu, hs_tau2, hs_xi,
                             hs_c2, hs_zeta, k, p, hs_idx,
                             estimate_c2 = TRUE, slab_df = 4, slab_s2 = 4,
                             equation_specific = FALSE) {
  # 1. Update Local Shrinkage
  if (equation_specific) {
    # Vectorized Local Shrinkage (Equation-Specific)
    n_hs <- length(hs_idx)
    # Update nu
    nu_rates <- 1 + 1 / hs_lambda2[hs_idx, , drop = FALSE]
    hs_nu[hs_idx, ] <- 1 / rgamma(n_hs * p, shape = 1, rate = nu_rates)
    # Update lambda2
    sq_beta <- curr_beta[hs_idx, , drop = FALSE]^2
    lambda2_rates <- 1 / hs_nu[hs_idx, , drop = FALSE] + sq_beta / (2 * matrix(hs_tau2, n_hs, p, byrow = TRUE))
    hs_lambda2[hs_idx, ] <- 1 / rgamma(n_hs * p, shape = 1, rate = lambda2_rates)
    hs_lambda2[hs_idx, ] <- pmin(pmax(hs_lambda2[hs_idx, ], 1e-10), 1e4)
  } else {
    # Vectorized Local Shrinkage (Pooled)
    hs_nu[hs_idx] <- 1 / rgamma(length(hs_idx), shape = 1, rate = 1 + 1 / hs_lambda2[hs_idx])
    ss_k <- rowSums(curr_beta[hs_idx, , drop = FALSE]^2)
    hs_lambda2[hs_idx] <- 1 / rgamma(length(hs_idx), shape = (p + 1) / 2, rate = 1 / hs_nu[hs_idx] + ss_k / (2 * hs_tau2))
    hs_lambda2[hs_idx] <- pmin(pmax(hs_lambda2[hs_idx], 1e-10), 1e4)
  }

  # 2. Update Global Shrinkage
  if (equation_specific) {
    if (length(hs_tau2) == 1) hs_tau2 <- rep(hs_tau2, p)
    if (length(hs_xi) == 1) hs_xi <- rep(hs_xi, p)

    hs_xi <- 1 / rgamma(p, shape = 1, rate = 1 + 1 / hs_tau2)
    ss_p <- colSums(curr_beta[hs_idx, , drop = FALSE]^2 / hs_lambda2[hs_idx, , drop = FALSE])
    hs_tau2 <- 1 / rgamma(p, shape = (length(hs_idx) + 1) / 2, rate = 1 / hs_xi + ss_p / 2)
    hs_tau2 <- pmin(pmax(hs_tau2, 1e-10), 1e4)
  } else {
    hs_xi <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / hs_tau2)
    ss_all_hs <- sum(curr_beta[hs_idx, , drop = FALSE]^2 / hs_lambda2[hs_idx])
    hs_tau2 <- 1 / rgamma(1, shape = (length(hs_idx) * p + 1) / 2, rate = 1 / hs_xi + ss_all_hs / 2)
    hs_tau2 <- pmin(pmax(hs_tau2, 1e-10), 1e4)
  }

  # 3. Conditional Slab Update (Finnish part)
  if (estimate_c2) {
    hs_zeta <- rig(1, 1 + 1 / hs_c2)
    shape_c <- (slab_df + length(hs_idx) * p) / 2
    rate_c <- (slab_df * slab_s2) / 2 + sum(curr_beta[hs_idx, ]^2) / 2
    hs_c2 <- rig(shape_c, rate_c)
  } else {
    hs_c2 <- slab_s2
  }

  if (equation_specific) {
    return(list(
      lambda2 = hs_lambda2, nu = hs_nu, tau2 = hs_tau2, xi = hs_xi,
      c2 = hs_c2, zeta = hs_zeta
    ))
  } else {
    return(list(
      lambda2 = as.vector(hs_lambda2), nu = as.vector(hs_nu),
      tau2 = hs_tau2, xi = hs_xi, c2 = hs_c2, zeta = hs_zeta
    ))
  }
}



# --- Robust Positive Definite Check ---
# Escalating ridge + Matrix::nearPD fallback to ensure TruncatedNormal doesn't crash
ensure_pd <- function(V, max_jitter = 1e-2) {
  # 1. Force absolute symmetry
  V <- 0.5 * (V + t(V))

  # 2. Fast escalating ridge using eigenvalues (more robust than chol)
  k <- nrow(V)
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  min_ev <- min(ev)

  if (min_ev > 1e-9) {
    # Matrix is already reasonably PD, but add a tiny floor for rtmvnorm safety
    return(V + diag(1e-11, k))
  }

  # NOTE: Global counter removed — it was not thread-safe for parallel chain
  # execution via future/mclapply and had no downstream consumer.

  jitter <- max(1e-10, abs(min_ev) * 1.5)
  is_pd <- FALSE

  cat(sprintf("\n[ensure_pd Warning] Matrix not PD (min_eigen = %g). Applying jitter...\n", min_ev))

  while (!is_pd && jitter <= max_jitter) {
    V_test <- V + diag(jitter, k)
    # Check PD via eigenvalues again
    ev_test <- eigen(V_test, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev_test) > 1e-9) {
      is_pd <- TRUE
      V <- V_test
    } else {
      jitter <- jitter * 10
    }
  }

  # 3. Nuclear fallback: Matrix::nearPD (handles complex singular cases)
  if (!is_pd) {
    cat("[ensure_pd Warning] Jitter limit reached. Falling back to Matrix::nearPD...\n")
    V <- as.matrix(nearPD(V, ensureSymmetry = TRUE, posd.tol = 1e-8)$mat)
  }

  # Final safety floor for TruncatedNormal
  return(V + diag(1e-11, k))
}

# --- R-native tree prediction for BART "slim trees" ---
# Walk the tree structure extracted from dbarts::getTrees()
predict_slim_bart_r <- function(X, tree_df) {
  n <- nrow(X)
  out <- numeric(n)

  # Group nodes by tree for faster iteration
  tree_list <- split(tree_df, tree_df$tree)

  for (df in tree_list) {
    # Each df is one tree
    vars <- df$var
    vals <- df$value
    n_nodes <- nrow(df)

    # Pre-calculate right children indices (Iterative approach to avoid recursion depth issues)
    right_children <- rep(-1, n_nodes)
    stack <- list()
    for (idx in 1:n_nodes) {
      if (vars[idx] != -1) {
        # Internal node: left child is idx + 1
        # We need to know where the right child is.
        # This is a bit tricky with preorder.
      }
    }

    # Simpler walker: since dbarts stores trees as binary, we can walk it.
    # Actually, the recursion in R is fine for BART trees (depth 4-6).
    right_children <- rep(-1, n_nodes)
    calc_right <- function(curr) {
      if (curr > n_nodes) {
        return(0)
      }
      if (vars[curr] == -1) {
        return(1)
      }
      l_size <- calc_right(curr + 1)
      r_idx <- curr + 1 + l_size
      right_children[curr] <<- r_idx
      r_size <- calc_right(r_idx)
      return(1 + l_size + r_size)
    }
    calc_right(1)

    # Walk tree for each observation
    for (i in 1:n) {
      curr <- 1
      while (vars[curr] != -1) {
        v <- vars[curr]
        if (X[i, v] <= vals[curr]) {
          curr <- curr + 1
        } else {
          curr <- right_children[curr]
        }
      }
      out[i] <- out[i] + vals[curr]
    }
  }
  return(out)
}

# --- Worker Refresh (Future/Parallel) ---
refresh_workers <- function() {
  cat("\n>>> Refreshing parallel workers to clear memory...\n")
  old_plan <- future::plan()
  future::plan(future::sequential)
  gc()
  future::plan(old_plan)
  cat("Workers refreshed.\n")
}

# --- Label Formatting for Reports ---
format_mnl_labels <- function(cov_names, sds) {
  # Map raw names to descriptive labels
  m <- c(
    intercept                = "Model Intercept",
    log_rent                 = "Agricultural Rent (log)",
    Slope_rad                = "Terrain Slope",
    Elevation                = "Elevation",
    log_Elevation            = "Elevation (log)",
    Aspect_cos_mean          = "North-South Exposure",
    Aspect_sin_mean          = "East-West Exposure",
    allPA_share              = "Protected Area Share",
    GHM_HI_2000              = "Human Modification (GHM)",
    log_RAI                  = "Road Accessibility (log)",
    RAI                      = "Road Accessibility",
    CISI                     = "Critical Infrastructure (CISI)",
    GDP_2000                 = "GDP",
    log1p_GDP_2000           = "GDP (log)",
    Pop_2000                 = "Population",
    log1p_Pop_2000           = "Population Density (log)"
  )

  # Neighborhood fallback map
  focal_map <- c(
    Arable = "Arable", Permanent = "Permanent Crops",
    Heterogeneous = "Hetero. Agri.", Pasture = "Pasture",
    Potential_Land = "Potent. Agri. Land", Urban = "Urban Area",
    Forest = "Forest", Sparse_Vegetation = "Sparse Veg.",
    Wetland = "Wetland", Inland_Waters = "Inland Water",
    Marine_Waters = "Marine Water",
    # LUM Detailed Classes
    Cropland_HI = "Cropland HI", Cropland_HIO = "Cropland HI Organic",
    Cropland_IR = "Cropland IR", Cropland_IRO = "Cropland IR Organic",
    Cropland_LI = "Cropland LI", Cropland_LIO = "Cropland LI Organic",
    Cropland_other = "Cropland Other", Cropland_other_O = "Cropland Other Organic",
    Forests_HI = "Forests HI", Forests_LI = "Forests LI",
    Forests_MI = "Forests MI", Forests_primary = "Forests Primary",
    Forests_unmanaged = "Forests Unmanaged",
    Forests_shortrotation = "Forests Short Rotation",
    Forests_SR = "Forests Short Rotation",
    Natural_other = "Natural Other", Natural_unmanaged = "Natural Unmanaged",
    Pasture_HI = "Pasture HI", Pasture_HIO = "Pasture HI Organic",
    Pasture_LI = "Pasture LI", Pasture_LIO = "Pasture LI Organic",
    Urban = "Urban Area",
    Waterbodies_inland = "Inland Water", Waterbodies_marine = "Marine Water",
    Grassland_HI = "Grassland HI", Grassland_LI = "Grassland LI",
    Settlements_HI = "Settlements HI", Settlements_LI = "Settlements LI",
    Wetlands_natural = "Wetlands Natural", Wetlands_peat = "Wetlands Peat"
  )

  new_labels <- character(length(cov_names))
  for (i in seq_along(cov_names)) {
    nm <- cov_names[i]
    base_name <- m[nm]

    # Handle Neighborhood
    if (is.na(base_name) && grepl("^focal_", nm)) {
      cls <- gsub("^focal_", "", nm)
      mapped <- focal_map[cls]
      # focal_map returns NA (not NULL) for unknown classes, which %||% does not
      # catch -> all unknown classes collapsed to "Neighbor: NA" (duplicate labels
      # crash factor(levels=...)). Fall back to the class name itself when unmapped.
      base_name <- paste0("Neighbor: ", if (is.na(mapped)) cls else unname(mapped))
    }

    # Handle Soil
    if (is.na(base_name)) {
      if (grepl("^OC_TOP", nm)) {
        base_name <- "Soil: Organic Carbon"
      } else if (grepl("^ROO", nm)) {
        base_name <- "Soil: Rooting Depth"
      } else if (grepl("^AWC_TOP", nm)) {
        base_name <- "Soil: Water Capacity"
      } else if (grepl("^VS", nm)) base_name <- "Soil: Stone Content"

      if (!is.na(base_name)) {
        # Add descriptive labels for soil levels (1=Low/Shallow, 5=High/Deep)
        lvl_map <- c(
          s0 = " (Level 0: Minimum)",
          s1 = " (Very Low/Shallow)",
          s2 = " (Low/Shallow)",
          s3 = " (Medium/Moderate)",
          s4 = " (High/Deep)",
          s5 = " (Very High/Deep)",
          s6 = " (Level 6: Maximum)"
        )

        for (l in seq_along(lvl_map)) {
          lvl_code <- names(lvl_map)[l] # e.g. "s1"
          lvl_num <- gsub("s", "", lvl_code) # e.g. "1"

          pat <- paste0("(_s|s|_)?", lvl_num, "(_|$)")
          if (grepl(pat, nm)) {
            base_name <- paste0(base_name, lvl_map[l])
            break
          }
        }

        if (grepl("_dev$", nm)) base_name <- paste0(base_name, " deviation")
      }
    }

    # Generic cleaning if still no label
    if (is.na(base_name)) {
      base_name <- nm
      base_name <- gsub("^log1p_", "log(1+", base_name)
      if (grepl("^log\\(1\\+", base_name)) base_name <- paste0(base_name, ")")
      base_name <- gsub("_dev", " deviation", base_name)
      if (base_name == "intercept") base_name <- "Intercept"
    }

    # Add SD scaling information
    s <- sds[nm]
    if (is.na(s) || s == 1) {
      new_labels[i] <- base_name
    } else {
      new_labels[i] <- sprintf("%s (sd=%.2f)", base_name, s)
    }
  }
  return(new_labels)
}

# --- Combine Multiple MCMC Chains ---
combine_chains <- function(res_list, keep_chains = FALSE) {
  n_chains <- length(res_list)
  if (n_chains == 1) {
    return(res_list[[1]])
  }

  if (keep_chains) {
    out <- res_list[[1]]
    
    # Identify all summary fields to combine
    fields <- c("postb_pooled", "postb_total", "post_log_lik", 
                "post_f", "post_sigma_re", "post_r", "post_log_lik_pointwise")

    for (f in fields) {
      if (is.null(out[[f]])) next
      chain_data <- lapply(res_list, function(res) res[[f]])
      if (is.null(chain_data[[1]])) next

      # 1. Get original dimensions (fallback to length if it's a vector)
      orig_dim <- dim(chain_data[[1]])
      if (is.null(orig_dim)) orig_dim <- length(chain_data[[1]])
      
      # 2. Append the chain dimension at the end
      new_dim <- c(orig_dim, n_chains)
      
      # 3. Create the new array (works universally for any N-dimensional structure)
      out[[f]] <- array(do.call(c, chain_data), dim = new_dim)
      
      # 4. Preserve and update dimnames
      orig_dimnames <- dimnames(chain_data[[1]])
      if (!is.null(orig_dimnames)) {
        new_dimnames <- orig_dimnames
        new_dimnames[[length(new_dimnames) + 1]] <- paste0("chain_", 1:n_chains)
        names(new_dimnames)[length(new_dimnames)] <- "chain"
        dimnames(out[[f]]) <- new_dimnames
      }
    }

    # Handle tree store separately: Keep them grouped by chain in a named list
    if (!is.null(out$tree_store)) {
      out$tree_store <- lapply(res_list, `[[`, "tree_store")
      names(out$tree_store) <- paste0("chain_", 1:n_chains)
    }

    # Horseshoe special handling
    if (!is.null(out$horseshoe)) {
      for (hf in names(out$horseshoe)) {
        if (grepl("^post_", hf)) {
          h_data <- lapply(res_list, function(res) res$horseshoe[[hf]])
          if (is.null(h_data[[1]])) {
            out$horseshoe[[hf]] <- NULL
            next
          }
          
          orig_dim <- dim(h_data[[1]])
          if (is.null(orig_dim)) orig_dim <- length(h_data[[1]])
          
          new_dim <- c(orig_dim, n_chains)
          out$horseshoe[[hf]] <- array(do.call(c, h_data), dim = new_dim)
        }
      }
    }

    # post_f_mean is a summary statistic. 
    # We leave this as the grand mean across all chains.
    if (!is.null(out$post_f_mean)) {
      f_means <- lapply(res_list, function(res) res$post_f_mean)
      out$post_f_mean <- Reduce("+", f_means) / n_chains
    }

    # Ensure the alias postb matches postb_total
    out$postb <- out$postb_total

    return(out)
  }

  # Classic collapsed combine_chains logic
  out <- res_list[[1]]
  # Identify all summary fields to combine
  fields <- c("postb_pooled", "postb_total", "post_log_lik", "post_f", "post_sigma_re", "post_r", "post_log_lik_pointwise")

  for (f in fields) {
    if (is.null(out[[f]])) next
    chain_data <- lapply(res_list, function(res) res[[f]])
    if (is.null(chain_data[[1]])) next

    d <- length(dim(chain_data[[1]]))
    if (d == 0) {
      out[[f]] <- do.call(c, chain_data)
    } else if (d == 2) {
      # Matrix: usually [iterations, parameters] or [parameters, iterations]
      # For log_lik_pointwise: [iterations, observations] -> combine along iterations (1)
      # For sigma_re: [parameters, iterations] -> combine along iterations (2)
      if (f == "post_log_lik_pointwise") {
        out[[f]] <- do.call(rbind, chain_data)
      } else {
        out[[f]] <- do.call(cbind, chain_data)
      }
    } else if (d == 3) {
      # Array: [covariates, categories, iterations]
      new_dim <- dim(chain_data[[1]])
      new_dim[3] <- sum(sapply(chain_data, function(x) dim(x)[3]))
      out[[f]] <- array(do.call(c, chain_data), dim = new_dim)
      dimnames(out[[f]]) <- dimnames(chain_data[[1]])
    } else if (d == 4) {
      # Array: [covariates, categories, groups, iterations]
      new_dim <- dim(chain_data[[1]])
      new_dim[4] <- sum(sapply(chain_data, function(x) dim(x)[4]))
      out[[f]] <- array(do.call(c, chain_data), dim = new_dim)
      dimnames(out[[f]]) <- dimnames(chain_data[[1]])
    }
  }

  # Handle tree store separately (list of lists)
  if (!is.null(out$tree_store)) {
    out$tree_store <- do.call(c, lapply(res_list, `[[`, "tree_store"))
  }

  # Horseshoe special handling
  if (!is.null(out$horseshoe)) {
    for (hf in names(out$horseshoe)) {
      if (grepl("^post_", hf)) {
        h_data <- lapply(res_list, function(res) res$horseshoe[[hf]])
        if (is.null(h_data[[1]])) {
          out$horseshoe[[hf]] <- NULL
          next
        }
        out$horseshoe[[hf]] <- do.call(cbind, h_data)
      }
    }
  }

  if (!is.null(out$post_f_mean)) {
    f_means <- lapply(res_list, function(res) res$post_f_mean)
    out$post_f_mean <- Reduce("+", f_means) / length(f_means)
  }

  # Ensure the alias postb matches postb_total
  out$postb <- out$postb_total

  return(out)
}

# --- Convergence Assessment ---
assess_convergence <- function(res_list, name, categories = NULL, baseline_idx = NULL,
                               type = c("pooled", "re")) {

  # Local fallback for %||% (base R >= 4.4 has it, but define defensively)
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  cat("\n=========================================================================\n")
  cat("Assessing convergence:", name, "\n")
  cat("=========================================================================\n")

  is_disk_backed <- is.character(res_list) && length(res_list) == 1 && dir.exists(res_list)
  dir_path <- NULL

  if (is_disk_backed) {
    dir_path <- res_list
    cat("    [Disk-backed mode detected. Loading global states...]\n")
    library(qs2)  # FIX #10: library() not require() — hard dep must error loudly

    # Auto-load metadata if categories/baseline were not provided
    meta_file <- file.path(dir_path, "model_metadata.qs")
    if (file.exists(meta_file)) {
      meta <- qs_read(meta_file)
      if (is.null(categories))   categories   <- c(meta$cat_names, meta$baseline_name)
      if (is.null(baseline_idx)) baseline_idx <- length(categories)
    }

    chain_files <- list.files(
      dir_path,
      pattern   = "posterior_batch_1_chain_[0-9]+\\.qs",
      full.names = TRUE
    )
    if (length(chain_files) == 0) stop("Could not find chain files in directory.")

    # FIX #6: replace per-iteration loop with simplify2array() for mu and sigma_re
    res_list <- lapply(chain_files, function(f) {
      batch_iters <- qs_read(f)

      # Stack mu: each element is [K, P] -> result is [K, P, n_iters]
      mu_arr <- simplify2array(lapply(batch_iters, `[[`, "mu"))

      # Stack sigma_re only where non-NULL; fill NA otherwise
      sigma_raw  <- lapply(batch_iters, `[[`, "sigma_re")
      has_sigma  <- !vapply(sigma_raw, is.null, logical(1))
      k_active   <- dim(mu_arr)[1]
      p_active   <- dim(mu_arr)[2]
      n_iters    <- dim(mu_arr)[3]
      sigma_arr  <- array(NA_real_, dim = c(k_active, p_active, n_iters))
      if (any(has_sigma)) {
        sigma_arr[, , has_sigma] <- simplify2array(sigma_raw[has_sigma])
      }

      rm(batch_iters, sigma_raw)  # free immediately
      list(postb_pooled = mu_arr, post_sigma_re = sigma_arr)
    })
  }

  if (is.null(categories) || is.null(baseline_idx)) {
    stop("categories and baseline_idx must be provided if not using a directory with model_metadata.qs")
  }

  # Category labels for the coefficient columns. The symmetric sampler (mnlogit_rcpp_sym)
  # returns ALL K categories (baseline included); the baseline-removed sampler returns K-1.
  # Derive the count + order from the array itself so both layouts work.
  .pp_cat <- if (!is.null(res_list[[1]]$postb_pooled)) dim(res_list[[1]]$postb_pooled)[2] else length(categories) - 1L
  .pp_nm  <- if (!is.null(res_list[[1]]$postb_pooled)) dimnames(res_list[[1]]$postb_pooled)[[2]] else NULL
  cat_names <- if (!is.null(.pp_nm) && length(.pp_nm) == .pp_cat) {
    .pp_nm                                                   # array carries its own category order
  } else if (!is.null(baseline_idx) && .pp_cat == length(categories) - 1L) {
    categories[-baseline_idx]                                # baseline-removed (K-1)
  } else if (!is.null(baseline_idx) && .pp_cat == length(categories)) {
    c(categories[-baseline_idx], categories[baseline_idx])   # symmetric (K), baseline last (matches heatplot)
  } else {
    categories
  }
  type      <- match.arg(type)

  # ---------------------------------------------------------------------------
  # Internal helper: flatten one chain's Rcpp result into a [n_iters x n_vars] matrix
  # ---------------------------------------------------------------------------
  extract_chain <- function(res) {
    k     <- dim(res$postb_pooled)[1]
    p     <- dim(res$postb_pooled)[2]
    n_ret <- dim(res$postb_pooled)[3]

    draws_beta <- aperm(res$postb_pooled, c(3, 1, 2))
    dim(draws_beta) <- c(n_ret, k * p)
    var_names_fixed <- dimnames(res$postb_pooled)[[1]] %||% res$var_names %||% paste0("V", seq_len(k))
    colnames(draws_beta) <- as.vector(outer(var_names_fixed, cat_names, paste, sep = "_"))
    out <- draws_beta

    # Random Effect Variances
    # NOTE: postb_total (individual pixel REs) is deliberately excluded here -
    # flattening tens of thousands of groups causes OOM. Pooled + sigma_re
    # diagnostics are sufficient for hierarchical mixing assessment.
    if (type == "re" && !is.null(res$post_sigma_re)) {
      if (length(dim(res$post_sigma_re)) == 3) {
        draws_sig      <- aperm(res$post_sigma_re, c(3, 1, 2))
        dim(draws_sig) <- c(n_ret, dim(res$post_sigma_re)[1] * dim(res$post_sigma_re)[2])
      } else {
        draws_sig <- t(res$post_sigma_re)   # legacy 2-D fallback
      }
      colnames(draws_sig) <- paste0("Sigma_RE_", seq_len(ncol(draws_sig)))
      out <- cbind(out, draws_sig)
    }

    # Horseshoe shrinkage parameters
    if (!is.null(res$horseshoe)) {
      if (!is.null(res$horseshoe$post_kappa_pooled)) {
        draws_kp            <- t(res$horseshoe$post_kappa_pooled)
        colnames(draws_kp)  <- paste0("Kappa_Pooled_", seq_len(ncol(draws_kp)))
        out <- cbind(out, draws_kp)
      }
      if (!is.null(res$horseshoe$post_kappa_re)) {
        draws_kr            <- t(res$horseshoe$post_kappa_re)
        colnames(draws_kr)  <- paste0("Kappa_RE_", seq_len(ncol(draws_kr)))
        out <- cbind(out, draws_kr)
      }
    }

    # BART / Log-Likelihood
    if (!is.null(res$post_f_mean) && !is.null(res$post_f)) {
      draws_f            <- t(apply(res$post_f, c(2, 3), mean))
      colnames(draws_f)  <- paste0("BART_avg_", cat_names)
      if (nrow(draws_f) == nrow(out)) out <- cbind(out, draws_f)
    }

    ll_total <- res$post_log_lik %||% res$log_lik
    if (!is.null(ll_total)) {
      if (is.matrix(ll_total) && ncol(ll_total) > 1) ll_total <- rowSums(ll_total)
      out <- cbind(out, matrix(ll_total, ncol = 1, dimnames = list(NULL, "Total_LogLik")))
    }
    out
  }

  mcmc_list <- lapply(res_list, extract_chain)
  n_iters   <- nrow(mcmc_list[[1]])
  n_chains  <- length(mcmc_list)
  n_vars    <- ncol(mcmc_list[[1]])

  # FIX #8: build draws array without a loop via simplify2array + aperm
  # simplify2array over a list of [n_iters x n_vars] matrices gives [n_iters, n_vars, n_chains]
  draws_arr <- aperm(simplify2array(mcmc_list), c(1, 3, 2))
  dimnames(draws_arr) <- list(
    Iteration = NULL,
    Chain     = NULL,
    Variable  = colnames(mcmc_list[[1]])
  )

  library(posterior)  # FIX #10
  draws_obj <- as_draws_array(draws_arr)

  cat(">>> Computing Bulk-ESS, Tail-ESS, and R-hat for fixed/global parameters...\n")
  full_stats <- summarise_draws(draws_obj)
  if (inherits(full_stats, "data.frame")) setDT(full_stats)

  # ---------------------------------------------------------------------------
  # Chunked / streamed Random Effects convergence
  # ---------------------------------------------------------------------------
  if (type == "re") {

    # CASE A: Memory-backed - postb_total is in RAM as a 4-D array
    if (!is_disk_backed &&
        !is.null(res_list[[1]]$postb_total) &&
        length(dim(res_list[[1]]$postb_total)) == 4) {

      cat(">>> Chunk-processing exact Random Effects (postb_total) convergence...\n")
      k_active <- dim(res_list[[1]]$postb_total)[1]
      p_len    <- dim(res_list[[1]]$postb_total)[2]
      G        <- dim(res_list[[1]]$postb_total)[3]

      var_names_fixed <- dimnames(res_list[[1]]$postb_pooled)[[1]] %||%
                         res_list[[1]]$var_names                     %||%
                         paste0("V", seq_len(k_active))
      var_names_re    <- var_names_fixed[seq_len(k_active)]
      G_names         <- res_list[[1]]$nuts0_names %||% paste0("G", seq_len(G))

      re_stats_list <- vector("list", p_len * k_active)
      list_idx <- 1L

      for (p_idx in seq_len(p_len)) {
        cn <- cat_names[p_idx]
        for (k_idx in seq_len(k_active)) {
          vn         <- var_names_re[k_idx]
          chunk_list <- lapply(res_list, function(r) t(r$postb_total[k_idx, p_idx, , ]))

          # FIX #8 (CASE A): replace chain loop with simplify2array
          # Each chunk_list[[ch]] is [n_iters x G]; simplify gives [n_iters, G, n_chains]
          chunk_arr <- aperm(simplify2array(chunk_list), c(1, 3, 2))
          dimnames(chunk_arr) <- list(
            Iteration = NULL, Chain = NULL,
            Variable  = paste(G_names, vn, cn, sep = "_")
          )

          chunk_stats <- summarise_draws(as_draws_array(chunk_arr))
          if (inherits(chunk_stats, "data.frame")) setDT(chunk_stats)
          re_stats_list[[list_idx]] <- chunk_stats
          list_idx <- list_idx + 1L
        }
      }

      full_stats <- rbind(full_stats, rbindlist(re_stats_list), fill = TRUE)

    # CASE B: Disk-backed - stream one batch x one chain at a time
    } else if (is_disk_backed) {

      cat(">>> Streaming Random Effects batches directly from disk...\n")

      batch_files_c1 <- list.files(
        dir_path,
        pattern    = "posterior_batch_[0-9]+_chain_1\\.qs",
        full.names = TRUE
      )
      n_batches <- length(batch_files_c1)
      # n_chains already derived from the loaded global states above

      if (n_batches > 0) {
        re_stats_list <- vector("list", n_batches)

        for (b in seq_len(n_batches)) {
          cat(sprintf("    Processing Batch %d / %d...\n", b, n_batches))  # FIX #13: \n not \r

          # Read batch dimensions from chain 1 only (avoid loading all chains at once)
          chain1_data  <- qs_read(file.path(dir_path,
                                            paste0("posterior_batch_", b, "_chain_1.qs")))
          b_iters      <- length(chain1_data)
          b_k_active   <- dim(chain1_data[[1]]$beta)[1]
          b_p_len      <- dim(chain1_data[[1]]$beta)[2]
          b_G          <- dim(chain1_data[[1]]$beta)[3]
          rm(chain1_data)

          # FIX #3: guard against category count mismatch
          if (b_p_len != length(cat_names)) {
            stop(sprintf(
              "Batch %d: beta dim[2] = %d but length(cat_names) = %d. Category mismatch.",
              b, b_p_len, length(cat_names)
            ))
          }

          b_var_names <- paste0("V", seq_len(b_k_active))
          b_G_names   <- paste0("Batch", b, "_G", seq_len(b_G))

          # Preallocate full [Iters, Chains, K*P*G] array
          batch_draws <- array(NA_real_,
                               dim = c(b_iters, n_chains, b_k_active * b_p_len * b_G))

          # FIX #9: load and process one chain at a time, free immediately
          # FIX #7: replace inner iteration loop with simplify2array
          for (ch in seq_len(n_chains)) {
            f_name     <- paste0("posterior_batch_", b, "_chain_", ch, ".qs")
            chain_data <- qs_read(file.path(dir_path, f_name))

            # simplify2array over iterations: each $beta is [K,P,G] -> [K, P, G, n_iters]
            raw        <- simplify2array(lapply(chain_data, `[[`, "beta"))
            dim(raw)   <- c(b_k_active * b_p_len * b_G, b_iters)  # flatten K*P*G, keep iters
            batch_draws[, ch, ] <- t(raw)

            rm(chain_data, raw)  # free this chain before loading the next
          }

          # Assign dimension names
          grid_names <- expand.grid(
            Var   = b_var_names,
            Cat   = cat_names,
            Group = b_G_names
          )
          dimnames(batch_draws) <- list(
            Iteration = NULL,
            Chain     = NULL,
            Variable  = paste(grid_names$Group, grid_names$Var, grid_names$Cat, sep = "_")
          )

          b_stats <- summarise_draws(as_draws_array(batch_draws))
          if (inherits(b_stats, "data.frame")) setDT(b_stats)
          re_stats_list[[b]] <- b_stats
          rm(batch_draws)
        }

        full_stats <- rbind(full_stats, rbindlist(re_stats_list), fill = TRUE)
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Block 1: Log-Likelihood & Model Fit
  # ---------------------------------------------------------------------------
  cat("\n[1] LOG-LIKELIHOOD & MODEL FIT:\n")
  ll_stats <- full_stats[variable == "Total_LogLik"]
  if (nrow(ll_stats) > 0) {
    cat(sprintf(
      "    Total LogLik Rhat: %.4f | ESS (Bulk): %.1f | ESS (Tail): %.1f\n",
      ll_stats$rhat, ll_stats$ess_bulk, ll_stats$ess_tail
    ))
  }

  # WAIC / LOO-CV
  ll_pw_list <- lapply(res_list, `[[`, "post_log_lik_pointwise")
  if (!is.null(ll_pw_list[[1]])) {
    cat("\n    Model Comparison (LOO-CV):\n")
    ll_combined <- abind::abind(ll_pw_list, along = 1)
    tryCatch({
      waic_res <- loo::waic(ll_combined)
      loo_res  <- loo::loo(ll_combined)
      cat(sprintf(
        "      WAIC:   %.2f (SE: %.2f) | p_waic: %.2f\n",
        waic_res$estimates["waic",     "Estimate"],
        waic_res$estimates["waic",     "SE"],
        waic_res$estimates["p_waic",   "Estimate"]
      ))
      cat(sprintf(
        "      LOO-CV: %.2f (SE: %.2f) | p_loo:  %.2f\n",
        loo_res$estimates["elpd_loo", "Estimate"] * -2,
        loo_res$estimates["elpd_loo", "SE"]       *  2,
        loo_res$estimates["p_loo",   "Estimate"]
      ))
    }, error = function(e) cat("      WAIC/LOO failed:", conditionMessage(e), "\n"))

  } else if (is_disk_backed) {
    # FIX #4: explicit notice instead of silent skip in disk-backed mode
    cat("\n    [WAIC/LOO skipped: post_log_lik_pointwise not in global chain state.",
        "Load from a dedicated pointwise batch if needed.]\n")
  }

  # ---------------------------------------------------------------------------
  # Block 2: Parameter Convergence (by Category)
  # ---------------------------------------------------------------------------
  cat("\n[2] PARAMETER CONVERGENCE (by Category):\n")
  diag_dir <- file.path("output/diagnostics",      tolower(gsub(" ", "_", name)))
  plot_dir <- file.path("output/plots/diagnostics", tolower(gsub(" ", "_", name)))
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  for (cn in cat_names) {
    cat(sprintf("\n    >>> CATEGORY: %s <<<\n", cn))
    pat        <- paste0("_", cn, "$")
    stats_cat  <- full_stats[grepl(pat, variable)]

    is_bart_cat  <- grepl("^BART_avg_", stats_cat$variable)
    is_re_cat    <- grepl("^G[0-9]+_",  stats_cat$variable)
    is_fixed_cat <- !is_bart_cat & !is_re_cat

    if (any(is_fixed_cat)) {
      cat("      Fixed Effects:\n")
      f_s <- stats_cat[is_fixed_cat]
      for (i in seq_len(nrow(f_s))) {
        var_name_clean <- sub(paste0("_", cn, "$"), "", f_s$variable[i])
        cat(sprintf(
          "        %-20s | Rhat: %.4f | ESS (Bulk): %.1f\n",
          var_name_clean, f_s$rhat[i], f_s$ess_bulk[i]
        ))
      }
    }

    if (any(is_re_cat)) {
      re_s <- stats_cat[is_re_cat]
      cat(sprintf("      Random Effects (Summary of %d regions):\n", nrow(re_s)))
      cat(sprintf(
        "        Rhat: Avg %.4f | Max %.4f | >1.1: %d\n",
        mean(re_s$rhat, na.rm = TRUE), max(re_s$rhat, na.rm = TRUE),
        sum(re_s$rhat > 1.1,  na.rm = TRUE)
      ))
      cat(sprintf(
        "        ESS:  Avg %.1f | Min %.1f | <100: %d\n",
        mean(re_s$ess_bulk, na.rm = TRUE), min(re_s$ess_bulk, na.rm = TRUE),
        sum(re_s$ess_bulk < 100, na.rm = TRUE)
      ))
    }

    # Robust safe_write_csv to handle OneDrive cloud storage locking timeouts
    safe_write_csv <- function(x, file) {
      for (i in 1:5) {
        res <- tryCatch({
          write.csv(x, file = file, row.names = FALSE)
          TRUE
        }, error = function(e) {
          Sys.sleep(1)
          FALSE
        })
        if (res) return(TRUE)
      }
      cat("      [Warning] Failed to write CSV after 5 attempts due to file locks. Skipping.\n")
      return(FALSE)
    }
    safe_write_csv(stats_cat, file.path(diag_dir, paste0("stats_", tolower(cn), ".csv")))
  }

  # ---------------------------------------------------------------------------
  # Block 3: Trace Plots & Worst Cases
  # ---------------------------------------------------------------------------
  tryCatch({
    library(bayesplot)  # FIX #10
    library(ggplot2)    # To enable labs() and theme() inside MCMC trace plots
    cat("\n[3] DIAGNOSTIC PLOTS:\n")

    bad_vars <- full_stats[order(rhat, decreasing = TRUE)]$variable

    # FIX #2: restrict candidate variables to those actually in draws_obj (fixed/global only)
    available_in_draws <- dimnames(draws_obj)[[3]]
    plot_vars <- unique(c(
      "Total_LogLik",
      head(bad_vars[!grepl("LogLik", bad_vars) & bad_vars %in% available_in_draws], 15)
    ))
    # Remove Total_LogLik if it's also not in draws_obj
    plot_vars <- plot_vars[plot_vars %in% available_in_draws]

    cat(sprintf("    Generating trace plots for %d worst parameters...\n", length(plot_vars)))

    for (p_name in plot_vars) {
      clean_p_name <- gsub("[^[:alnum:]_]", "_", p_name)
      pdf_file     <- file.path(plot_dir, paste0("trace_", clean_p_name, ".pdf"))

      # FIX #5: guard against zero-row lookup before building subtitle
      p_row       <- full_stats[variable == p_name]
      subtitle_str <- if (nrow(p_row) > 0) {
        sprintf("R-hat: %.3f | ESS (Bulk): %.1f", p_row$rhat, p_row$ess_bulk)
      } else {
        "Stats unavailable"
      }

      p_trace <- mcmc_trace(draws_obj, pars = p_name) +
        labs(title = paste("Trace Plot:", p_name), subtitle = subtitle_str) +
        theme_minimal()

      ggsave(pdf_file, p_trace, width = 10, height = 6)
    }

    # R-hat distribution plot across all parameters
    p_rhat <- ggplot(full_stats[is.finite(rhat)], aes(x = "Parameters", y = rhat)) +
      geom_jitter(alpha = 0.4, width = 0.2, color = "dodgerblue4") +
      geom_hline(yintercept = 1.1, linetype = "dashed", color = "red") +
      theme_minimal() +
      labs(title = paste("R-hat Distribution:", name), y = "R-hat")

    ggsave(
      file.path("output/plots/diagnostics",
                paste0("rhat_dist_", tolower(gsub(" ", "_", name)), ".png")),
      p_rhat, width = 8, height = 6
    )

  }, error = function(e) cat("    Diagnostic plotting failed:", conditionMessage(e), "\n"))

  return(invisible(full_stats))
}

# --- Covariate Ordering Helper ---
reorder_mnl_covariates <- function(cov_names, cat_names = NULL) {
  tier_patterns <- c(
    "^(OC_TOP|ROO|AWC_TOP|VS)",
    "^focal_",
    "^(allPA|GHM|RAI|CISI|GDP|Pop|log_RAI|log1p_Pop|log1p_GDP)",
    "^(Slope|Elevation|Aspect)",
    "^log_rent$",
    "^intercept$"
  )

  ranks <- rep(length(tier_patterns) + 1, length(cov_names))
  for (i in seq_along(tier_patterns)) {
    matches <- grepl(tier_patterns[i], cov_names, ignore.case = TRUE)
    ranks[matches & ranks > i] <- i
  }

  # Within the focal_ tier, sort to match category order (creates diagonal)
  focal_mask <- grepl("^focal_", cov_names)
  if (any(focal_mask) && !is.null(cat_names)) {
    # Extract the class name after "focal_" and match to cat_names order
    focal_classes <- gsub("^focal_", "", cov_names[focal_mask])
    # Assign sub-rank based on position in cat_names
    focal_order <- match(focal_classes, cat_names)
    # Unmatched focal classes go to the end
    focal_order[is.na(focal_order)] <- length(cat_names) + seq_len(sum(is.na(focal_order)))
    # Create a fractional rank so focal_ covariates sort within their tier
    focal_positions <- which(focal_mask)
    ranks[focal_positions] <- ranks[focal_positions] + focal_order / (max(focal_order) + 1)
  }

  return(cov_names[order(ranks)])
}

# --- Detect Sum-to-Constant Blocks ---
detect_constant_sum_blocks <- function(X, tol = 1e-4, mean_tol = 0.05,
                                        cv_tol = 0.05, robust = TRUE) {
  # Detect groups of columns whose row sums are (more or less) constant.
  # `tol` is the empty-row threshold; a block qualifies if the row-sum over valid rows is
  # near a constant ~1 within RELATIVE tolerances (mean within mean_tol of 1, robust
  # coefficient-of-variation < cv_tol). Robust (median/MAD) by default so a handful of
  # partial-coverage outlier rows (e.g. focal-neighbourhood border pixels: most sum ~1, a
  # few ~0.5) don't disqualify an otherwise-constant block the way plain mean/sd did.
  # Returns a named list of integer vectors (column indices).
  k <- ncol(X)
  cn <- colnames(X)
  if (is.null(cn)) cn <- paste0("V", seq_len(k))
  
  assigned <- rep(FALSE, k)
  blocks <- list()
  
  # Strategy: group by column name prefix, verify numerically
  get_prefix <- function(name) {
    if (grepl("^focal_", name)) return("focal")
    m <- regmatches(name, regexpr("^[A-Za-z_]+(?=_s?\\d)", name, perl = TRUE))
    if (length(m) > 0 && nchar(m) > 0) return(m)
    return(name)
  }
  
  prefixes <- sapply(cn, get_prefix)
  for (pfx in unique(prefixes)) {
    matching <- which(prefixes == pfx & !assigned)
    if (length(matching) < 2) next
    
    # We need to verify the block sums to a constant
    rs <- rowSums(X[, matching, drop = FALSE])
    
    # Exclude entirely empty rows (like border pixels missing focal data)
    valid_rows <- rs > tol
    if (sum(valid_rows) == 0) next
    
    valid_rs <- rs[valid_rows]

    # Centre + spread: robust (median / MAD) so rare partial-coverage outliers don't disqualify
    # an otherwise-constant block; non-robust falls back to mean / sd.
    if (length(valid_rs) == 1) {
      ctr <- valid_rs; spread <- 0
    } else if (robust) {
      ctr <- median(valid_rs); spread <- mad(valid_rs)
    } else {
      ctr <- mean(valid_rs); spread <- sd(valid_rs)
    }
    cv <- spread / max(abs(ctr), 1e-8)
    if (abs(ctr - 1) < mean_tol && cv < cv_tol) {
      blocks[[pfx]] <- matching
      assigned[matching] <- TRUE
    }
  }
  
  return(blocks)
}

# =============================================================================
# Deterministic RE-identifiability screen (shared across samplers).
# Flags (group x RE-covariate) cells whose random slope cannot be estimated, over
# the two X-side failure modes (the Y-side, A1, is handled per sampler via the
# separation / zero-prevalence gates):
#   A2  no within-group X variation : sd(X[group, r]) ~ 0      -> slope unidentified
#   A3  within-group X-X collinearity: |cor(r, r_keep)| > thresh -> drop lower-var one
# Returns a logical [G x q] mask (TRUE = pin that group's slope to 0). Samplers pin
# the flagged REs every sweep, deterministically.
# =============================================================================
screen_re_design <- function(X, group_idx, re_idx, sd_tol = 1e-8,
                             cor_thresh = 0.95, min_obs = 10L) {
  gl <- split(seq_len(nrow(X)), group_idx); G <- length(gl); q <- length(re_idx)
  re_mask <- matrix(FALSE, G, q, dimnames = list(NULL, colnames(X)[re_idx]))
  # Globally-constant RE columns (the intercept) carry a random *intercept*, which is
  # identified by the group mean even with no X variation -> EXEMPT from the A2 no-variation
  # rule (otherwise we'd pin every group's random intercept and force the FE to overfit).
  global_const <- apply(X[, re_idx, drop = FALSE], 2, function(c) { s <- sd(c); !is.finite(s) || s < sd_tol })
  for (gi in seq_len(G)) {
    R <- gl[[gi]]; Xg <- X[R, re_idx, drop = FALSE]; sds <- apply(Xg, 2, sd)
    re_mask[gi, (!is.finite(sds) | sds < sd_tol) & !global_const] <- TRUE   # A2: local no-variation (not the intercept)
    if (length(R) >= min_obs) {                                  # A3: collinearity
      cm <- suppressWarnings(abs(cor(Xg))); cm[is.na(cm)] <- 0; diag(cm) <- 0
      repeat {
        hit <- which(cm > cor_thresh, arr.ind = TRUE); hit <- hit[hit[, 1] < hit[, 2], , drop = FALSE]
        if (!nrow(hit)) break
        pr <- hit[1, ]; drop <- if (sds[pr[1]] <= sds[pr[2]]) pr[1] else pr[2]
        re_mask[gi, drop] <- TRUE; cm[drop, ] <- 0; cm[, drop] <- 0
      }
    }
  }
  re_mask
}

# Support-aware RE prior factor: per (group x covariate) participation ratio n_g / PR_{g,r},
# where PR = (sum d^2)^2 / sum d^4 on the WITHIN-GROUP-CENTERED column d = x - mean_g(x) =
# effective # observations informing that group's *slope* (deviation from the group mean).
# Centering matters: PR is scale- but NOT location-invariant, so a globally-centered sparse
# column would otherwise look near-constant (high PR). factor = (n_g / PR)^strength >= 1: ~1 for
# an evenly-spread covariate, large (-> hard pin) for one with little within-group spread. The
# graduated generalization of the A2 screen. Globally-constant columns (the random intercept,
# identified by the group mean) are EXEMPT (factor 1), mirroring screen_re_design(). Same group
# ordering as screen_re_design() so the [G x q] matrices align. strength 0 -> all 1.
re_support_design <- function(X, group_idx, re_idx, strength = 1, sd_tol = 1e-8) {
  gl <- split(seq_len(nrow(X)), group_idx); G <- length(gl); q <- length(re_idx)
  fac <- matrix(1, G, q, dimnames = list(NULL, colnames(X)[re_idx]))
  if (strength <= 0) return(fac)
  global_const <- apply(X[, re_idx, drop = FALSE], 2, function(c) { s <- sd(c); !is.finite(s) || s < sd_tol })
  for (gi in seq_len(G)) {
    Xg <- X[gl[[gi]], re_idx, drop = FALSE]; ng <- nrow(Xg)
    Xc <- sweep(Xg, 2, colMeans(Xg), "-")                 # within-group centering -> slope information
    s2 <- colSums(Xc^2); s4 <- colSums(Xc^4); pr <- s2^2 / pmax(s4, 1e-12)
    f <- (ng / pmax(pr, 1))^strength; f[global_const] <- 1
    fac[gi, ] <- f
  }
  fac
}

# Support-aware FIXED-effect prior factor: GLOBAL participation ratio n / PR per covariate,
# PR = (sum d^2)^2 / sum d^4 on the globally-centered column d = x - mean(x) (FE = global slope).
# factor = (n / PR)^strength >= 1, shrinking sparse covariates' fixed effects (mirrors the CLR
# FE support_factor). Globally-constant columns (intercept) exempt. strength 0 -> all 1. Returns
# a length-ncol(X) vector; callers apply it to the horseshoe / prior precision of the hs columns.
fe_support_factor <- function(X, strength = 1, sd_tol = 1e-8) {
  p <- ncol(X); f <- rep(1, p)
  if (strength <= 0) return(f)
  Xc <- sweep(X, 2, colMeans(X), "-")
  s2 <- colSums(Xc^2); s4 <- colSums(Xc^4); pr <- s2^2 / pmax(s4, 1e-12)
  f <- (nrow(X) / pmax(pr, 1))^strength
  gconst <- apply(X, 2, function(c) { s <- sd(c); !is.finite(s) || s < sd_tol })
  f[gconst | !is.finite(f)] <- 1
  f
}

# =============================================================================
# AUTOMATIC BART-covariate support screen (variable pruning for the BART block).
# A BART tree can only stably split where there is SUPPORT: a covariate present in
# a handful of pixels (a near-indicator share like soil class _s5, ~0.2% prevalence)
# lets trees carve tiny groups -> f-tail roughness / poor surface reproducibility,
# yet adds nothing OOS. The principled, automatic measure is the participation ratio
#   PR = (sum d^2)^2 / sum d^4   on the centered column d = x - mean(x)
# = effective # observations informing the covariate (a share nonzero in k pixels has
# PR ~= k; a well-spread covariate PR ~= n). Same metric as fe_support_factor /
# re_support_design. Auto-drops: near-constant (sd<tol), low-support (PR < min_pr),
# and near-duplicate (|cor|>cor_thresh, drop the lower-PR one). No covariate-name lists,
# no per-variable prevalence cutoffs to tune -- one support floor (default 0.5% of n).
# Returns kept/dropped column indices (into X) + the PR diagnostic table.
# =============================================================================
screen_bart_design <- function(X, bart_idx, min_pr = NULL, min_pr_frac = 0.005,
                               sd_tol = 1e-8, cor_thresh = 0.999) {
  Xb <- as.matrix(X[, bart_idx, drop = FALSE]); n <- nrow(Xb)
  if (is.null(min_pr)) min_pr <- max(20, min_pr_frac * n)
  Xc <- sweep(Xb, 2, colMeans(Xb), "-")
  s2 <- colSums(Xc^2); s4 <- colSums(Xc^4); pr <- s2^2 / pmax(s4, 1e-12)
  sds <- apply(Xb, 2, sd)
  drop <- !is.finite(sds) | sds < sd_tol | pr < min_pr
  keep0 <- which(!drop)                                       # collinearity among survivors
  if (length(keep0) > 1) {
    cm <- suppressWarnings(abs(cor(Xb[, keep0, drop = FALSE]))); cm[is.na(cm)] <- 0; diag(cm) <- 0
    repeat {
      hit <- which(cm > cor_thresh, arr.ind = TRUE); hit <- hit[hit[, 1] < hit[, 2], , drop = FALSE]
      if (!nrow(hit)) break
      pp <- keep0[hit[1, ]]; dcol <- hit[1, if (pr[pp[1]] <= pr[pp[2]]) 1 else 2]
      drop[keep0[dcol]] <- TRUE; cm[dcol, ] <- 0; cm[, dcol] <- 0
    }
  }
  list(keep_idx = bart_idx[!drop], drop_idx = bart_idx[drop], min_pr = min_pr,
       pr = setNames(round(pr, 1), colnames(X)[bart_idx]))
}

# =============================================================================
# Standard cross-sampler fit object. MNL / LNM / CLR all emit this so the
# OVERLAPPING parameters (beta coefficients) share dimensions + names:
#   sampler, B_mean [cov x cat], rhat [cov x cat], draws [cov x cat x ndraws x nchains],
#   cov_names, cat_names, Sigma_mean (NULL for MNL), baseline_class.
# `draws` must be [cov, cat, ndraws, nchains] (symmetric, all K categories).
# =============================================================================
build_fit_object <- function(sampler, draws, cov_names, cat_names,
                             baseline_class = NA, Sigma_mean = NULL,
                             dropped_blocks = NULL) {
  # Re-insert driver-level drop-one references: each dropped reference's per-class effect is the
  # zero-sum -sum(kept block columns), reconstructed per draw so B_mean + rhat include it. Appended
  # to the covariate dimension. Requires 4D draws [cov, cat, ndraws, nchains].
  if (!is.null(dropped_blocks) && length(dropped_blocks) > 0 && length(dim(draws)) == 4L) {
    d0 <- dim(draws)
    for (nm in names(dropped_blocks)) {
      b <- dropped_blocks[[nm]]
      kept <- match(b$kept_names, cov_names); kept <- kept[!is.na(kept)]
      if (length(kept) == 0L || b$ref_name %in% cov_names) next
      ref <- -apply(draws[kept, , , , drop = FALSE], c(2, 3, 4), sum)   # [cat, ndraws, nchains]
      draws <- abind::abind(draws, array(ref, c(1, d0[2], d0[3], d0[4])), along = 1)
      cov_names <- c(cov_names, b$ref_name)
    }
  }
  d <- dim(draws)
  dn <- vector("list", length(d)); dn[[1]] <- cov_names; dn[[2]] <- cat_names
  dimnames(draws) <- dn
  B_mean <- apply(draws, c(1, 2), mean)
  rhat <- if (length(d) >= 4 && d[4] > 1) {
    apply(draws, c(1, 2), function(M) {            # M is [ndraws x nchains]
      n <- nrow(M); W <- mean(apply(M, 2, var)); B <- n * var(colMeans(M))
      if (!is.finite(W) || W <= 0) NA_real_ else sqrt(((n - 1) / n * W + B / n) / W)
    })
  } else matrix(NA_real_, d[1], d[2])
  dimnames(B_mean) <- dimnames(rhat) <- list(cov_names, cat_names)
  list(sampler = sampler, B_mean = B_mean, rhat = rhat, draws = draws,
       cov_names = cov_names, cat_names = cat_names,
       Sigma_mean = Sigma_mean, baseline_class = baseline_class)
}
