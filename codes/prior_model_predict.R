# =============================================================================
# Phase 1 of the predict-wrapper (docs/predict_wrapper_design.md): reusable predictor on an
# ALREADY-ASSEMBLED design matrix X (focal etc. already in X). Reproduces the estimation-time
# numerics EXACTLY: BART f from the calibrated slim trees (reconstruct_bart_f_mean) + linear/RE
# (postb_total per group, appearance-order keying, unseen->pooled) + softmax. Standardization is
# NOT re-applied (betas are raw-scale; slim trees split on raw inputs).
# Requires reconstruct_bart_f_mean + predict_slim_bart_cpp (source codes/mnlogit_rcpp_sym.R).
# =============================================================================

# group_levels = unique(group_idx) IN THE SAMPLER'S APPEARANCE ORDER (postb_total[,,k] = group_levels[k]).
# class_map: OPTIONAL integer vector (length = #output categories, values in 1:J) mapping each output
#   category to a FITTED category -> nested classes inherit the parent's predictor ("same hierarchical
#   predictor value"). Children get eta_parent - log(#siblings) so a parent's total share is PRESERVED and
#   split equally by default. NULL/identity => current behaviour.
predict_shares <- function(fit, X, bart_cols, linear_cols, group_idx = NULL, group_levels = NULL,
                           type = c("mean", "predictive"), thin = 1L, class_map = NULL,
                           unseen_re = c("pooled", "marginal"), re_draws = 25L) {
  type <- match.arg(type); unseen_re <- match.arg(unseen_re)
  stopifnot(!is.null(fit$postb_total))
  J  <- dim(fit$postb_total)[2]
  # For an UNSEEN group the hierarchical predictive is beta_new ~ N(mu_pooled, Sigma_RE), NOT the pooled
  # POINT (which is overconfident -- validated: unseen-country log-score ~ EU-mean null). unseen_re="marginal"
  # marginalizes over the RE spread: Sigma_RE = the empirical between-group sd of the fitted REs (~0 for
  # pooled covariates, so only the RE covariates inflate). k x J sd matrix reused below.
  Xl <- as.matrix(X[, linear_cols, drop = FALSE])
  Xb <- as.matrix(X[, bart_cols,  drop = FALSE])
  # pp = the NON-BASELINE category positions. `seq_len(J-1)` is only correct when the baseline is
  # the LAST category; with any other baseline the BART columns were being placed in the wrong
  # slots. Take it from the fit when available.
  .bl <- fit$baseline
  bart_meta <- list(symmetric = TRUE, p_all = J,
                    pp = if (!is.null(.bl) && .bl >= 1 && .bl <= J) setdiff(seq_len(J), .bl) else seq_len(J - 1))
  has_re <- length(dim(fit$postb_total)) == 4L && !is.null(group_idx)
  if (has_re && is.null(group_levels)) stop("group_levels (appearance-order unique(group_idx)) required for RE prediction.")
  # broadcast fitted-category eta -> output categories via class_map, with -log(#siblings) to keep parent totals.
  expand <- function(eta) { if (is.null(class_map)) return(eta)
    K <- as.numeric(table(class_map)[as.character(class_map)])          # #siblings per output category
    sweep(eta[, class_map, drop = FALSE], 2, log(K), "-") }
  softmax <- function(E) { E <- expand(E); e <- exp(E - apply(E, 1, max)); e / rowSums(e) }

  if (type == "mean") {
    f <- reconstruct_bart_f_mean(fit$tree_store, Xb, bart_meta)          # n x J (posterior-mean BART f)
    # A LINEAR-ONLY fit has no tree_store, so reconstruct_bart_f_mean returns NULL and every
    # downstream `eta + ...` collapsed to numeric(0) -- softmax then failed with an opaque
    # "dim(X) must have a positive length". use_bart = FALSE is the PRODUCTION default, so this
    # path is the common one, not an edge case. Start from an explicit zero utility instead.
    eta <- if (is.null(f)) matrix(0, nrow(X), J) else f
    if (has_re) {
      B     <- apply(fit$postb_total, c(1, 2, 3), mean)                  # k x J x G
      Bpool <- apply(fit$postb_pooled, c(1, 2), mean)                    # k x J
      re_sd <- if (unseen_re == "marginal") apply(B, c(1, 2), sd) else NULL  # empirical Sigma_RE (k x J)
      k <- nrow(Bpool); marg_rows <- integer(0); marg_P <- NULL
      for (g in unique(group_idx)) {
        r <- which(group_idx == g); key <- match(g, group_levels)
        if (is.na(key) && unseen_re == "marginal") {
          acc <- matrix(0, length(r), J)                                 # MC-average softmax over beta_new = mu_pooled + N(0, Sigma_RE)
          for (mc in seq_len(re_draws)) {
            Bg <- Bpool + matrix(rnorm(k * J), k, J) * re_sd
            E  <- expand(eta[r, , drop = FALSE] + Xl[r, , drop = FALSE] %*% Bg)
            e  <- exp(E - apply(E, 1, max)); acc <- acc + e / rowSums(e)
          }
          marg_rows <- c(marg_rows, r); marg_P <- rbind(marg_P, acc / re_draws)
        } else {
          Bg <- if (is.na(key)) Bpool else B[, , key]
          eta[r, ] <- eta[r, ] + Xl[r, , drop = FALSE] %*% Bg            # marginal rows keep eta unchanged (result overwritten below)
        }
      }
      P <- softmax(eta)
      if (length(marg_rows)) P[marg_rows, ] <- marg_P
      return(P)
    } else {
      eta <- eta + Xl %*% apply(fit$postb_pooled, c(1, 2), mean)
    }
    return(softmax(eta))
  }

  # ---- posterior-predictive: average softmax over (thinned) draws ----
  has_bart <- !is.null(fit$tree_store) && length(fit$tree_store) > 0L
  nd <- if (has_bart) length(fit$tree_store) else dim(fit$postb_pooled)[3]
  ds <- seq(1L, nd, by = thin); Pbar <- matrix(0, nrow(X), J)
  for (d in ds) {
    eta <- if (!has_bart) matrix(0, nrow(X), J) else {
      g_ens <- vapply(seq_along(fit$tree_store[[d]]),
                      function(cc) as.numeric(predict_slim_bart_cpp(Xb, fit$tree_store[[d]][[cc]])),
                      numeric(nrow(X)))
      sweep(g_ens, 1, rowMeans(g_ens), "-")                              # CLR centering (symmetric)
    }
    if (has_re) {
      for (g in unique(group_idx)) {
        r <- which(group_idx == g); key <- match(g, group_levels)
        Bg <- if (is.na(key)) fit$postb_pooled[, , d] else fit$postb_total[, , key, d]
        eta[r, ] <- eta[r, ] + Xl[r, , drop = FALSE] %*% Bg
      }
    } else {
      eta <- eta + Xl %*% fit$postb_pooled[, , d]
    }
    Pbar <- Pbar + softmax(eta)
  }
  Pbar / length(ds)
}

# =============================================================================
# Phase 2: recipe + prior_model wrapper. build_recipe() runs at FIT (learns col order, group levels,
# res) and returns the assembled (X, Y, group_idx) + the recipe; apply_recipe() REPLAYS it at predict.
# The driver builds its training X via build_recipe -> fit & predict share one code path.
# Needs compute_focal_coord() (source experiments/focal/compute_focal_coord.R) + predict_shares().
# =============================================================================
suppressMessages(library(data.table))

# focal_<class> columns from the LU state (areas -> per-row shares -> coord queen-8 focal -> broadcast).
.add_focal <- function(raw, lu_classes, coord, slice, res, complete_only = TRUE, add_nodata = FALSE) {
  d <- as.data.table(copy(raw))
  pre <- grep("^focal_", names(d), value = TRUE); if (length(pre)) d[, (pre) := NULL]   # drop any stored focal_ before recompute
  sh <- d[, ..lu_classes]; rs <- rowSums(sh, na.rm = TRUE)
  sh <- as.data.table(lapply(sh, function(x) fifelse(rs > 0, x / rs, 0)))   # areas -> shares (driver L624-627)
  dn <- cbind(d[, c(coord, slice), with = FALSE], sh)
  fc <- compute_focal_coord(dn, classes = lu_classes, coord = coord, slice = slice, res = res, complete_only = complete_only)
  d2 <- merge(d, fc, by = c(coord, slice), all.x = TRUE, sort = FALSE)
  for (j in grep("^focal_", names(d2), value = TRUE)) set(d2, which(is.na(d2[[j]])), j, 0)
  # DERIVED column, not an LU class: the uncovered share of the neighbourhood. Kept OUT of
  # `lu_classes` so compute_focal_coord never looks for a "NODATA" layer, but it is what makes the
  # focal block sum to exactly 1 -- the precondition for the zero-sum const-sum reconstruction.
  # Mirrors the fit-time definition in run_prior_module_pixel_level_model.R.
  if (isTRUE(add_nodata)) {
    fc_cls <- paste0("focal_", lu_classes)
    set(d2, j = "focal_NODATA", value = pmax(0, 1 - rowSums(d2[, ..fc_cls])))
  }
  d2[]
}

.apply_transforms <- function(d, transforms) {
  for (tr in transforms) { fn <- match.fun(tr$fn); cols <- intersect(tr$cols, names(d))
    pfx <- if (!is.null(tr$prefix)) tr$prefix else ""          # prefix => write a RENAMED col (e.g. log1p_GDP); "" => in-place
    for (c in cols) set(d, j = paste0(pfx, c), value = fn(d[[c]])) }
  d
}

build_recipe <- function(raw, feature_cols, lu_classes, outcome_classes, coord = c("X","Y"),
                         slice = "out_year", group_col = "Grouping_Key", transforms = list(), res = NULL,
                         add_nodata = FALSE) {
  d <- as.data.table(copy(raw))
  if (is.null(res)) res <- min(diff(sort(unique(d[[coord[1]]]))))
  d <- .add_focal(d, lu_classes, coord, slice, res, add_nodata = add_nodata)
  focal_cols <- c(paste0("focal_", lu_classes), if (isTRUE(add_nodata)) "focal_NODATA")
  d <- .apply_transforms(d, transforms)
  col_order <- c("intercept", feature_cols, focal_cols)
  X <- cbind(intercept = 1, as.matrix(d[, c(feature_cols, focal_cols), with = FALSE])); X[!is.finite(X)] <- 0; X[,1] <- 1
  Y <- as.matrix(d[, ..outcome_classes])
  g <- as.integer(factor(d[[group_col]]))
  recipe <- list(feature_cols = feature_cols, lu_classes = lu_classes, outcome_classes = outcome_classes,
                 coord = coord, slice = slice, group_col = group_col, transforms = transforms, res = res,
                 col_order = col_order, focal_cols = focal_cols, add_nodata = isTRUE(add_nodata),
                 group_levels_factor = levels(factor(d[[group_col]])),   # for as.integer(factor) reproduction
                 group_levels_appear = unique(g))                        # sampler keying (appearance order)
  list(X = X, Y = Y, group_idx = g, recipe = recipe)
}

apply_recipe <- function(recipe, raw) {
  d <- as.data.table(copy(raw))
  d <- .add_focal(d, recipe$lu_classes, recipe$coord, recipe$slice, recipe$res,
                  add_nodata = isTRUE(recipe$add_nodata))
  d <- .apply_transforms(d, recipe$transforms)
  X <- cbind(intercept = 1, as.matrix(d[, c(recipe$feature_cols, recipe$focal_cols), with = FALSE])); X[!is.finite(X)] <- 0; X[,1] <- 1
  X <- X[, recipe$col_order, drop = FALSE]                              # exact training order
  g <- as.integer(factor(d[[recipe$group_col]], levels = recipe$group_levels_factor))
  list(X = X, group_idx = g)
}

# wrapper predict: raw X_t/Y_t -> exact-treatment design -> shares.
# class_map (optional): named/integer output-category -> fitted-category map for nested-class output
# (defaults to model$class_map, else identity).
predict_prior <- function(model, raw, type = c("mean","predictive"), thin = 1L, class_map = model$class_map) {
  a  <- apply_recipe(model$recipe, raw)
  bc <- match(model$bart_names,   colnames(a$X))
  lc <- match(model$linear_names, colnames(a$X))
  predict_shares(model$fit, a$X, bart_cols = bc, linear_cols = lc,
                 group_idx = a$group_idx, group_levels = model$recipe$group_levels_appear,
                 type = match.arg(type), thin = thin, class_map = class_map)
}

# =============================================================================
# Phase 3: dynamic projection. Iterate predict_prior feeding the predicted composition back as the LU
# state -> focal recomputes each step -> the neighbourhood co-evolves. DAMPED (auto-model criticality).
# Assumes outcome_classes align with lu_classes (predicted composition = the state that feeds focal).
# =============================================================================
project_prior <- function(model, raw, Y_init = NULL, n_steps = 25, damp = 0.5, tol = 1e-4,
                          type = "mean", verbose = TRUE) {
  d <- as.data.table(copy(raw)); lu <- model$recipe$lu_classes
  Y <- if (!is.null(Y_init)) as.matrix(Y_init) else as.matrix(d[, ..lu])
  Y <- Y / pmax(rowSums(Y), 1e-12)                                     # shares
  traj <- data.table(); s <- 0L; delta <- Inf
  repeat {
    s <- s + 1L
    d[, (lu) := as.data.table(Y)]                                      # current state -> focal source
    P <- predict_prior(model, d, type = type)                         # focal recomputed from Y
    Ynew  <- (1 - damp) * Y + damp * P
    delta <- max(abs(Ynew - Y)); mx <- max(colMeans(Ynew))
    traj  <- rbind(traj, data.table(step = s, max_dY = delta, max_class_share = mx))
    if (verbose) cat(sprintf("  step %2d | max|dY|=%.5f | max class share=%.3f\n", s, delta, mx))
    Y <- Ynew
    if (delta < tol || s >= n_steps) break
  }
  list(final = Y, trajectory = traj, steps = s, converged = delta < tol)
}

# =============================================================================
# F: package a fitted model + its recipe into a SELF-CONTAINED, shippable `prior_model` so predict/project
# run from RAW data through the SAME assembly the fit used (build_recipe) -> no train/predict drift.
# bart_names/linear_names = column names in recipe$col_order routed to BART vs linear (predict_prior matches
# them to the design). spatial_layout (optional) = the (coord, res) the focal used. Bundle is one .rds.
# =============================================================================
build_prior_model <- function(fit, recipe, bart_names, linear_names,
                              spatial_layout = NULL, class_names = NULL, class_map = NULL) {
  stopifnot(all(bart_names %in% recipe$col_order), all(linear_names %in% recipe$col_order))
  if (is.null(spatial_layout)) spatial_layout <- list(coord = recipe$coord, slice = recipe$slice, res = recipe$res)
  list(fit = fit, recipe = recipe, bart_names = bart_names, linear_names = linear_names,
       spatial_layout = spatial_layout,
       class_names = if (is.null(class_names)) recipe$outcome_classes else class_names,
       class_map = class_map)
}
save_prior_model <- function(model, path) { saveRDS(model, path); invisible(path) }

# Construct the recipe from a DRIVER that already assembled its design (dat_pixel + metadata), and VERIFY
# apply_recipe reproduces the driver's X_mat (self-consistency). Returns the recipe + per-block max|dX| so
# the driver can assert consistency before shipping. focal_cols in X_mat are terra(focal_df)-computed; the
# recipe recomputes them via compute_focal_coord -> report that drift separately (validated ~4e-4).
recipe_selfcheck <- function(recipe, raw, X_mat) {
  a <- apply_recipe(recipe, raw)
  common <- intersect(colnames(a$X), colnames(X_mat))
  d <- abs(a$X[, common, drop = FALSE] - X_mat[, common, drop = FALSE])
  foc <- grep("^focal_", common); non <- setdiff(seq_along(common), foc)
  df <- if (length(foc)) as.vector(d[, foc]) else 0
  list(max_dX_nonfocal = if (length(non)) max(d[, non]) else 0,
       max_dX_focal     = max(df), mean_dX_focal = mean(df),
       p99_dX_focal     = as.numeric(quantile(df, 0.99)), frac_focal_gt01 = mean(df > 0.01),
       ncol_match = length(common) == ncol(X_mat) && length(common) == ncol(a$X),
       cols_recipe = ncol(a$X), cols_driver = ncol(X_mat))
}

# =============================================================================
# COUNT model predictor: mu = r * exp(X.beta_g + f + offset) (NB(r,psi) param, log-link mean; NO softmax).
# Non-symmetric BART (per-column, f = mean_d predict_slim_bart_cpp, no CLR centering). Density = mu/exp(offset).
# type="mean" = plug-in (posterior-mean beta,r); "predictive" = posterior-predictive mean over draws.
# =============================================================================
.reconstruct_count_f_mean <- function(tree_store, Xb, p, fcap = Inf) {
  if (is.null(tree_store) || length(tree_store) == 0) return(NULL)
  n <- nrow(Xb); f <- matrix(0, n, p); nd <- 0L
  for (draw in tree_store) { if (is.null(draw)) next
    for (ip in seq_len(p)) f[, ip] <- f[, ip] +
      pmax(pmin(as.numeric(predict_slim_bart_cpp(Xb, draw[[ip]])), fcap), -fcap)  # match sampler cap (trees are uncapped)
    nd <- nd + 1L }
  if (nd == 0L) NULL else f / nd
}

predict_count <- function(fit, X, offset, linear_cols = seq_len(dim(fit$postb_total)[1]),
                          bart_cols = NULL, group_idx = NULL, group_levels = NULL,
                          type = c("mean","predictive"), thin = 1L, return = c("count","density"),
                          density_cap = NULL, cap_quantile = 1.0) {
  type <- match.arg(type); return <- match.arg(return)
  p  <- dim(fit$postb_total)[2]
  # DENSITY CAP: bound predicted density (mu/exp(offset) = animals per unit area) so no grid pixel
  # blows up to millions on the exp-link tail. NULL = off; a number (scalar or length-p) = absolute
  # ceiling; "auto" = the training density quantile `cap_quantile` (default the observed MAX, the
  # physical ceiling) from fit$train_density_q. Applied to DENSITY, per draw in the predictive path
  # (so one exploded draw can't drag the mean), then mu = density*exp(offset).
  dcap <- if (is.null(density_cap)) NULL
          else if (identical(density_cap, "auto")) {
            if (is.null(fit$train_density_q)) stop("density_cap='auto' needs fit$train_density_q (refit with the updated sampler)")
            rn <- sprintf("%g%%", as.numeric(cap_quantile) * 100)   # 1.0->"100%", 0.999->"99.9%"
            if (!(rn %in% rownames(fit$train_density_q)))
              stop("cap_quantile must be one of the stored quantiles: ", paste(rownames(fit$train_density_q), collapse=", "))
            fit$train_density_q[rn, ]
          } else density_cap
  if (!is.null(dcap)) dcap <- rep_len(as.numeric(dcap), p)
  # Mean param (default in current fits): beta models log(mu)=X.beta+f+offset (r-free)
  # -> mu = exp(eta+offset), NO r factor. Legacy natural param: mu = r*exp(eta+offset).
  mp <- isTRUE(fit$mean_param)
  Xl <- as.matrix(X[, linear_cols, drop = FALSE]); off <- matrix(offset, nrow(X), p)
  has_re   <- length(dim(fit$postb_total)) == 4L && !is.null(group_idx)
  has_bart <- !is.null(bart_cols) && !is.null(fit$tree_store) && length(fit$tree_store) > 0
  Xb <- if (has_bart) as.matrix(X[, bart_cols, drop = FALSE]) else NULL
  fc <- if (is.null(fit$bart_f_cap)) Inf else fit$bart_f_cap
  capmu <- function(m) { if (is.null(dcap)) return(m)
    for (ip in seq_len(p)) m[, ip] <- pmin(m[, ip] / exp(off[, ip]), dcap[ip]) * exp(off[, ip]); m }
  lin <- function(B, Bp) { eta <- matrix(0, nrow(X), p)
    if (has_re) for (g in unique(group_idx)) { r <- which(group_idx == g); key <- match(g, group_levels)
        eta[r, ] <- Xl[r, , drop = FALSE] %*% (if (is.na(key)) Bp else B[, , key]) }
    else eta <- Xl %*% Bp
    eta }

  if (type == "mean") {
    B  <- apply(fit$postb_total, c(1, 2, 3), mean); Bp <- apply(fit$postb_pooled, c(1, 2), mean)
    eta <- lin(B, Bp); if (has_bart) eta <- eta + .reconstruct_count_f_mean(fit$tree_store, Xb, p, fc)
    mu <- capmu((if (mp) 1 else mean(fit$post_r)) * exp(eta + off))
  } else {
    nd <- dim(fit$post_r)[2]; ds <- seq(1L, nd, by = thin); mu <- matrix(0, nrow(X), p)
    for (d in ds) { eta <- matrix(0, nrow(X), p)
      if (has_re) for (g in unique(group_idx)) { r <- which(group_idx == g); key <- match(g, group_levels)
          eta[r, ] <- Xl[r, , drop = FALSE] %*% (if (is.na(key)) fit$postb_pooled[, , d] else fit$postb_total[, , key, d]) }
      else eta <- Xl %*% fit$postb_pooled[, , d]
      if (has_bart) for (ip in seq_len(p)) eta[, ip] <- eta[, ip] +
        pmax(pmin(as.numeric(predict_slim_bart_cpp(Xb, fit$tree_store[[d]][[ip]])), fc), -fc)
      mu <- mu + capmu((if (mp) 1 else fit$post_r[, d]) * exp(eta + off)) }  # cap PER DRAW
    mu <- mu / length(ds)
  }
  if (return == "density") mu / exp(off) else mu
}

# combine_chains(keep_chains=TRUE) appends a "chain" dim: postb_pooled 3D->4D, postb_total 4D->5D,
# post_r 2D->3D, tree_store -> a per-chain named list. predict_count wants the NATIVE flat layout.
# Collapse the chain dim INTO the draws dim (draw fast, chain slow -- matches how combine concatenated
# memory, and how we concatenate the per-chain tree lists) so all chains' draws are pooled for predict.
.flatten_count_fit <- function(fit) {
  if (is.null(dim(fit$postb_pooled)) || length(dim(fit$postb_pooled)) <= 3) return(fit)  # native (no chain dim)
  collapse <- function(x) { d <- dim(x); if (length(d) < 2) return(x)
    dim(x) <- c(d[seq_len(length(d) - 2)], prod(d[(length(d) - 1):length(d)])); x }  # merge last 2 dims (column-major no-op)
  fit$postb_total  <- collapse(fit$postb_total)
  fit$postb_pooled <- collapse(fit$postb_pooled)
  if (!is.null(fit$post_r))    fit$post_r    <- collapse(fit$post_r)
  if (!is.null(fit$tree_store)) fit$tree_store <- do.call(c, unname(fit$tree_store))  # chain_1 draws, chain_2 draws, ...
  fit
}

# Grid predictor for the count model: apply a fitted (possibly multi-chain) count model to a NEW design
# (e.g. the 10km x NUTS3 downscale grid) with the density cap. group_country = the grid's RE-group label
# per row; fitted_group_levels = unique(group_idx_vec) from TRAINING (appearance order -- the count sampler
# keys REs by unique(group_idx), so the G dimension follows first-appearance, NOT sorted). Unseen groups
# fall back to the pooled mean. Returns counts (or density) with no pixel exceeding the density cap.
predict_count_grid <- function(fit, X, offset, group_country, fitted_group_levels,
                               train_group_idx = NULL, bart_cols = NULL, linear_cols = NULL,
                               type = c("mean","predictive"), density_cap = "auto", cap_quantile = 1.0,
                               return = c("count","density"), thin = 1L) {
  type <- match.arg(type); return <- match.arg(return)
  fit  <- .flatten_count_fit(fit)
  if (is.null(linear_cols)) linear_cols <- seq_len(dim(fit$postb_total)[1])
  # Map the grid's group labels to the TRAINING group codes. fitted_group_levels is the appearance-order
  # unique(train_group_idx); a grid label maps via its position among the training labels. If the caller
  # passes integer codes already aligned to training, they flow through; character labels are matched.
  grid_gi <- if (is.numeric(group_country)) as.integer(group_country)
             else match(as.character(group_country), as.character(fitted_group_levels))  # unseen -> NA -> pooled
  predict_count(fit, X, offset = offset, linear_cols = linear_cols, bart_cols = bart_cols,
                group_idx = grid_gi, group_levels = fitted_group_levels, type = type, thin = thin,
                return = return, density_cap = density_cap, cap_quantile = cap_quantile)
}

# =============================================================================
# F for the COUNT model (run_prior_module_count_model.R): one fit/predict code path, from raw admin data.
# NO focal (count uses own log1p_lu_area_* + a log-area OFFSET, not neighbour focal), so the recipe is just
# transforms + [intercept, time, features, year-dummies] assembly + offset + group keying. Mirrors
# build_recipe/apply_recipe but panel-aware (time col + year dummies) and offset-bearing.
# =============================================================================
apply_count_recipe <- function(recipe, raw) {
  d <- as.data.table(copy(raw))
  d <- .apply_transforms(d, recipe$transforms)                          # log1p skewed + lu_area (prefix log1p_)
  for (nm in names(recipe$year_map))                                    # year dummies from the time col
    set(d, j = nm, value = as.numeric(d[[recipe$time_col]] == recipe$year_map[[nm]]))
  M <- cbind(intercept = rep(1, nrow(d)))
  if (isTRUE(recipe$add_time)) M <- cbind(M, time = as.numeric(d[[recipe$time_col]]))
  M <- cbind(M, as.matrix(d[, recipe$feature_cols, with = FALSE])); M[!is.finite(M)] <- 0
  X <- M[, recipe$col_order, drop = FALSE]                              # exact training order (post covariate-selection)
  off <- log(pmax(as.numeric(d[[recipe$offset_col]]), recipe$offset_floor))
  g <- as.integer(factor(d[[recipe$group_col]], levels = recipe$group_levels_factor))
  list(X = X, offset = off, group_idx = g)
}

# selfcheck: apply_count_recipe reproduces the driver's X_mat + offset (asserted at DUMP time in the driver)
count_recipe_selfcheck <- function(recipe, raw, X_mat, offset_vec) {
  a <- apply_count_recipe(recipe, raw); common <- intersect(colnames(a$X), colnames(X_mat))
  list(ncol_match = ncol(a$X) == ncol(X_mat) && length(common) == ncol(X_mat),
       max_dX = if (length(common)) max(abs(a$X[, common, drop = FALSE] - X_mat[, common, drop = FALSE])) else NA,
       max_doffset = max(abs(a$offset - offset_vec)),
       group_match = all(a$group_idx == recipe$train_group_idx, na.rm = TRUE),
       cols_recipe = ncol(a$X), cols_driver = ncol(X_mat))
}

build_count_prior_model <- function(fit, recipe, linear_names, bart_names = NULL, class_name = NULL) {
  stopifnot(all(linear_names %in% recipe$col_order))
  list(fit = fit, recipe = recipe, linear_names = linear_names, bart_names = bart_names,
       class_name = class_name, kind = "count")
}

# raw admin data -> exact-treatment design -> count/density (density_cap defends the grid tail).
predict_count_prior <- function(model, raw, type = c("mean","predictive"), density_cap = "auto",
                                cap_quantile = 1.0, return = c("count","density"), thin = 1L) {
  a  <- apply_count_recipe(model$recipe, raw)
  lc <- match(model$linear_names, colnames(a$X))
  bc <- if (!is.null(model$bart_names)) { b <- match(model$bart_names, colnames(a$X)); b[!is.na(b)] } else NULL
  predict_count_grid(model$fit, a$X, offset = a$offset, group_country = a$group_idx,
                     fitted_group_levels = model$recipe$group_levels_appear,
                     bart_cols = if (length(bc)) bc else NULL, linear_cols = lc,
                     type = match.arg(type), density_cap = density_cap, cap_quantile = cap_quantile,
                     return = match.arg(return), thin = thin)
}

# =============================================================================
# population_averaged_effect(fit, group_idx, weights) — the SIZE-WEIGHTED estimand
# =============================================================================
# In a hierarchical model `mu` is the mean of the COUNTRY effect distribution: every country is one
# exchangeable draw, so Malta counts as much as Germany. For an EU-wide statement ("what does GHM_HI
# do to livestock density in Europe?") the decision-relevant quantity is the SIZE-WEIGHTED average
# marginal effect
#       theta_w = sum_g w_g * beta_g ,     w_g = herd share / area share / observation share
# which is a plain posterior functional -- computed per draw, so it carries honest uncertainty and
# needs NO refit. It differs from `mu` exactly when the effect is heterogeneous and correlated with
# size: on the BOV livestock fit GHM_HI had mu = 0.05 but a herd-weighted effect of +12.2, with 85%
# of countries positive (see also the `re_center` note in count_rcpp.R -- the two are complementary:
# re_center fixes WHICH parameter carries the average, this reports the average you actually want).
#
#   fit        : an mncount_rcpp / mnlogit_rcpp_sym fit with postb_total (RE)
#   group_idx  : the group vector the model was FIT with (raw ids)
#   weights    : per-observation weights (e.g. herd counts). NULL -> observation counts.
#   Returns a data.frame: covariate, mu (pooled), theta_w (weighted), sd, q025, q975, frac_pos.
population_averaged_effect <- function(fit, group_idx, weights = NULL, cov_names = NULL,
                                       target = 1L, probs = c(0.025, 0.975)) {
  stopifnot(!is.null(fit$postb_total))
  Bt <- fit$postb_total; Bp <- fit$postb_pooled
  if (length(dim(Bt)) == 5L) { d <- dim(Bt); dim(Bt) <- c(d[1], d[2], d[3], d[4] * d[5])   # drop chain dim
                               dp <- dim(Bp); dim(Bp) <- c(dp[1], dp[2], dp[3] * dp[4]) }
  # CRITICAL: postb_total[, , k] is the k-th group in APPEARANCE order (`unique(group_idx)`),
  # NOT group id k. Getting this wrong silently assigns each group another group's coefficients.
  lev <- unique(group_idx)
  G <- dim(Bt)[3]
  if (G != length(lev)) stop(sprintf("group mismatch: postb_total has %d slices, group_idx has %d distinct ids", G, length(lev)))
  w_g <- if (is.null(weights)) vapply(lev, function(g) sum(group_idx == g), 0) else
                               vapply(lev, function(g) sum(weights[group_idx == g], na.rm = TRUE), 0)
  w_g <- w_g / sum(w_g)
  nd <- dim(Bt)[4]; k <- dim(Bt)[1]
  th <- matrix(NA_real_, k, nd)                       # weighted effect per draw
  for (s in seq_len(nd)) th[, s] <- Bt[, target, , s] %*% w_g
  qs <- t(apply(th, 1, quantile, probs = probs, na.rm = TRUE))
  data.frame(
    covariate = cov_names %||% dimnames(Bt)[[1]] %||% paste0("V", seq_len(k)),
    mu        = rowMeans(Bp[, target, , drop = FALSE]),
    theta_w   = rowMeans(th),
    sd        = apply(th, 1, sd),
    q_lo      = qs[, 1], q_hi = qs[, 2],
    frac_pos  = rowMeans(apply(Bt[, target, , , drop = FALSE], c(1, 3), mean) > 0),
    stringsAsFactors = FALSE)
}
