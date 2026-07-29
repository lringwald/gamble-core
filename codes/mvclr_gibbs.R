# =============================================================================
# mvclr_gibbs.R — Multivariate CLR (logistic-normal) regression sampler  [v2]
# =============================================================================
# A lighter, conjugate alternative to the PG-MNL samplers. CLR-transforms the
# compositional response and fits a multivariate Gaussian regression:
#
#     z_i = CLR(Y_i)            (zero-sum across categories, K-vector)
#     z_i = X_i B + X_i^re b_{g(i)} + e_i,     e_i ~ N(0, Sigma)
#
# Features (parity goals with mnlogit_rcpp_sym, minus MNL-only machinery):
#   * Symmetric / zero-sum (CLR) coefficients — no privileged baseline, order-invariant
#   * CLR group horseshoe (regularized slab) on the pooled coefficients
#   * Hierarchical country random effects: INTERCEPT *or full SLOPES* via re_idx,
#     with per-covariate RE variances (automatic relevance: covariates that don't
#     vary by group are shrunk out)
#   * Explicit K x K residual covariance Sigma (Inverse-Wishart) — categories co-vary
#   * Conjugate Gibbs (no Polya-Gamma); RcppArmadillo accelerated matrix-normal draws
#
# NOT included (deliberately): separation/ spike-slab machinery (a Gaussian model has
# no separation pathology), CAR, and disk/recovery infra — see TODO at bottom.
# =============================================================================

.mvclr_reg_hs <- function(v, c2) c2 * v / (c2 + v)

# Try to compile the Rcpp core once; fall back to pure R if unavailable.
.mvclr_have_rcpp <- FALSE
try({
  suppressMessages({ library(Rcpp); library(RcppArmadillo) })
  Rcpp::sourceCpp("codes/mvclr_core.cpp")
  .mvclr_have_rcpp <- TRUE
}, silent = TRUE)

# Pure-R matrix-normal draw: B ~ MN(P^-1 RHS, P^-1, Sigma),  Schol = chol(Sigma) (upper)
.mn_draw_R <- function(P, RHS, Schol) {
  L <- chol(P)                                   # P = L'L (R chol is upper)
  M <- backsolve(L, backsolve(L, RHS, transpose = TRUE))   # P^-1 RHS
  Z <- matrix(rnorm(nrow(P) * ncol(RHS)), nrow(P), ncol(RHS))
  M + backsolve(L, Z) %*% Schol                  # backsolve(L,Z) = L^-1 Z = (P^-1)^{1/2}-type left factor
}

mvclr_gibbs <- function(Y, X, group_idx = NULL,
                        niter = 2000, nburn = 1000,
                        use_horseshoe = TRUE, horseshoe_idx = NULL,
                        use_re = TRUE, re_idx = NULL,    # covariate cols with random effects
                        use_bart = FALSE, bart_idx = NULL, n_trees_bart = 50,
                        A0 = 10, slab_s2 = 4,
                        iw_df = NULL, iw_scale = NULL,
                        re_a = 1, re_b = 1,
                        zero_eps = NULL, zero_prevalence_scaling = FALSE, standardize = TRUE, use_tempering = FALSE, tempering_T0 = 0.05,
                        init_state = NULL, save_to_disk = FALSE, disk_path = NULL,
                        use_re_spike_slab = FALSE, use_car = FALSE, car_adj = NULL, car_rho = 0.95,
                        use_rcpp = TRUE, re_asis = FALSE, re_screen = FALSE, support_prior_strength = 0,
                        init_jitter = 0, verbose = TRUE, ...) {
  rcpp <- use_rcpp && .mvclr_have_rcpp
  mn_draw <- if (rcpp) mn_draw_cpp else .mn_draw_R

  # ---- 1. CLR transform ------------------------------------------------------
  Y <- as.matrix(Y); n <- nrow(Y); K <- ncol(Y)
  cat_names <- colnames(Y) %||% paste0("cat", 1:K)
  Y <- Y / rowSums(Y)
  if (is.null(zero_eps)) zero_eps <- min(Y[Y > 0]) / 2
  if (zero_prevalence_scaling) {                          # prevalence-scaled eps (= driver STABILIZE_SPARSE_Y):
    spars <- colMeans(Y == 0)                             # deeply-zero classes keep more-extreme CLR -> zero-depth = signal
    Yp <- sweep(Y, 2, zero_eps * pmax(1 - spars, 0.02), "+")
  } else { Yp <- Y; Yp[Yp <= 0] <- zero_eps }
  Yp <- Yp / rowSums(Yp)
  Zclr <- log(Yp); Zclr <- Zclr - rowMeans(Zclr)
  V <- contr.helmert(K); V <- sweep(V, 2, sqrt(colSums(V^2)), "/")
  W <- Zclr %*% V; m <- K - 1L

  # ---- 2. Covariates ---------------------------------------------------------
  X <- as.matrix(X); p <- ncol(X)
  cov_names <- colnames(X) %||% paste0("V", 1:p)
  int_idx <- which(apply(X, 2, function(c) all(c == c[1])))
  Xmu <- rep(0, p); Xsd <- rep(1, p)
  if (standardize) for (j in setdiff(1:p, int_idx))
    if (length(unique(X[, j])) > 2 && sd(X[, j]) > 0) {
      Xmu[j] <- mean(X[, j]); Xsd[j] <- sd(X[, j]); X[, j] <- (X[, j] - Xmu[j]) / Xsd[j]
    }
  # Support-aware prior: shrink covariates whose information is concentrated in few high-leverage
  # observations more. n_eff = (sum x^2)^2 / sum x^4 (participation ratio): ~n for a dense covariate,
  # ~(nonzero count) for a sparse one. prior precision is scaled by (n/n_eff)^strength.
  support_factor <- rep(1, p)
  if (support_prior_strength > 0) {
    s2 <- colSums(X^2); s4 <- colSums(X^4); n_eff <- s2^2 / pmax(s4, 1e-12)
    support_factor <- (n / pmax(n_eff, 1))^support_prior_strength
  }

  # ---- 3. Groups (rows sorted by group for the RE block draw) ----------------
  if (use_re && !is.null(group_idx)) {
    gf0 <- as.integer(factor(group_idx)); grp_names <- levels(factor(group_idx))
    ord <- order(gf0); W <- W[ord, , drop = FALSE]; X <- X[ord, , drop = FALSE]; gf <- gf0[ord]
    G <- max(gf); gsize <- as.integer(table(gf)); gstart0 <- as.integer(cumsum(c(0, gsize[-G])))
    gidx <- split(seq_len(n), gf)
    if (is.null(re_idx)) re_idx <- if (length(int_idx)) int_idx else 1L   # default random intercept
    q <- length(re_idx); Xre <- X[, re_idx, drop = FALSE]
    # Deterministic RE-identifiability screen (A2 no-X-variation + A3 collinearity). Applied as an
    # IN-DRAW precision pin (huge prior precision on masked cells, per group) rather than a post-hoc
    # bre <- 0: post-hoc edits to bre inflate Sigma and run CLR away (Sigma scales the RE draw), but
    # an in-draw pin is a self-consistent draw. See memory clr-full-slopes-fragility.
    re_mask_mat <- if (re_screen) screen_re_design(X, gf, re_idx) else NULL
    if (re_screen && verbose) cat(sprintf("  RE screen pins %d of %d (group x covariate) slopes (in-draw)\n", sum(re_mask_mat), G * q))
    # Support-aware RE prior: scale each (group x covariate) RE prior precision by the
    # WITHIN-GROUP participation ratio n_g / PR_{g,r}. PR_{g,r} = (sum x^2)^2 / sum x^4 over
    # group g's rows = effective # observations informing that group's slope for covariate r.
    # A covariate that is constant within a group has PR=1 -> factor n_g -> pinned (the A2
    # limit, but graduated); a well-spread covariate has PR~n_g -> factor ~1 (untouched). This
    # is the smooth generalization of the binary A2 screen and shrinks information-sparse REs.
    re_support_mat <- NULL
    if (support_prior_strength > 0) {
      re_support_mat <- matrix(1, G, q)
      re_glob_const <- apply(Xre, 2, function(cc) { s <- sd(cc); !is.finite(s) || s < 1e-8 })
      for (gi in 1:G) {
        Xg <- Xre[gidx[[gi]], , drop = FALSE]; ng <- nrow(Xg)
        Xc <- sweep(Xg, 2, colMeans(Xg), "-")             # within-group centering -> slope information (PR is not location-invariant)
        s2g <- colSums(Xc^2); s4g <- colSums(Xc^4); ne_g <- s2g^2 / pmax(s4g, 1e-12)
        f <- (ng / pmax(ne_g, 1))^support_prior_strength; f[re_glob_const] <- 1   # exempt the random intercept
        re_support_mat[gi, ] <- f
      }
    }
  } else { use_re <- FALSE; G <- 0L; q <- 0L; re_mask_mat <- NULL; re_support_mat <- NULL }
  XtX <- crossprod(X)

  # ---- 4. Init & priors ------------------------------------------------------
  hs_idx <- setdiff(if (is.null(horseshoe_idx)) seq_len(p) else horseshoe_idx, int_idx)
  if (is.null(iw_df)) iw_df <- m + 2; if (is.null(iw_scale)) iw_scale <- diag(m)
  B <- matrix(0, p, m); bre <- matrix(0, max(G * q, 1L), m)   # rows blocked by group
  if (init_jitter > 0) B <- B + matrix(runif(p * m, -init_jitter, init_jitter), p, m)   # per-chain overdispersed init (honest Rhat)
  Sigma <- diag(m); SigInv <- diag(m); Schol <- chol(Sigma)
  lambda2 <- rep(1, p); nu <- rep(1, p); tau2 <- 1; xi <- 1; c2 <- slab_s2
  sigma_re2 <- rep(1, max(q, 1L)); prior_prec <- rep(1 / A0^2, p)
  if (!is.null(init_state)) {                                   # HOT-START
    for (nm in c("B", "bre", "sigma_re2", "lambda2", "tau2")) if (!is.null(init_state[[nm]])) assign(nm, init_state[[nm]])
    if (!is.null(init_state$Sigma)) { Sigma <- init_state$Sigma; SigInv <- chol2inv(chol(Sigma)); Schol <- chol(Sigma) }
  }
  delta_re <- rep(1L, max(G * q, 1L)); pi_re <- rep(0.5, max(q, 1L)); post_delta_sum <- rep(0, max(G * q, 1L))  # spike-slab RE
  do_car <- isTRUE(use_car) && !is.null(car_adj) && use_re; d_car <- if (do_car) rowSums(car_adj) else NULL     # CAR

  re_fitted <- function() {                     # n x m contribution of random effects
    if (!use_re) return(matrix(0, n, m))
    RF <- matrix(0, n, m)
    for (g in 1:G) { R <- gidx[[g]]; RF[R, ] <- Xre[R, , drop = FALSE] %*% bre[((g - 1) * q + 1):(g * q), , drop = FALSE] }
    RF
  }

  # ---- BART (per-category ensembles, backfit on the CLR residual; zero-sum) --
  f_w <- matrix(0, n, m)                                   # BART contribution in W-space
  if (use_bart) {
    suppressMessages(library(dbarts))
    if (is.null(bart_idx)) stop("use_bart = TRUE requires bart_idx (covariate columns for BART)")
    Xb <- X[, bart_idx, drop = FALSE]                     # already group-sorted & standardized
    Zc_sorted <- W %*% t(V)                                # CLR response (sorted), non-constant init
    bart_s <- lapply(1:K, function(k) dbarts::dbarts(Xb, Zc_sorted[, k],
      control = dbarts::dbartsControl(n.samples = 1, n.burn = 0, n.trees = n_trees_bart,
                                      keepTrees = TRUE, n.chains = 1, updateState = TRUE)))
    f_mat <- matrix(0, n, K)
  }

  nret <- niter - nburn
  postB <- array(0, c(p, K, nret), dimnames = list(cov_names, cat_names, NULL))
  postSigma <- array(0, c(K, K, nret)); post_sre <- matrix(0, max(q, 1L), nret)
  post_bre_sum <- if (use_re) matrix(0, G * q, K) else NULL    # accumulate RE (CLR space)
  post_f_sum <- if (use_bart) matrix(0, n, K) else NULL        # accumulate f (CLR, zero-sum)

  # ---- 5. Gibbs --------------------------------------------------------------
  # TEMPERING SETUP
  if (use_tempering) {
    cat(sprintf("Likelihood Tempering Enabled: linear from T=%.3f to 1.0 over first half of burnin (%d iters)\n", tempering_T0, floor(nburn/4)))
  }
  nburn_half <- max(1L, floor(nburn / 4))
  for (it in 1:niter) {
    temp_iter <- if (use_tempering && it <= nburn_half && nburn > 0) {
      tempering_T0 + (1.0 - tempering_T0) * (it / nburn_half)
    } else 1.0
    RF <- re_fitted()

    # (a) pooled B | Sigma, RE   (matrix-normal)
    if (use_horseshoe && length(hs_idx) > 0)
      prior_prec[hs_idx] <- support_factor[hs_idx] / pmax(.mvclr_reg_hs(tau2 * lambda2[hs_idx], c2), 1e-12)
    B <- mn_draw(XtX * temp_iter + diag(prior_prec, p), temp_iter * crossprod(X, W - RF - f_w), Schol)

    # (b) random effects | B, Sigma, sigma_re2
    if (use_re) {
      resid_re <- W - X %*% B - f_w
      # Ridge the RE precision so the per-group q x q system (crossprod(Xg) + Prc) stays PD
      # even with many random slopes / small groups (mirrors the MNL's half-Cauchy floor;
      # negligible when sigma_re2 is moderate, prevents chol failure when it blows up).
      Prc <- diag(1 / sigma_re2 + 1e-6, q)
      if (rcpp && !do_car && is.null(re_mask_mat) && is.null(re_support_mat)) {   # fast C++ path uses one shared Prc; per-group in-draw pin/support-scaling needs the R loop
        bre <- re_block_draw_cpp(Xre, resid_re, gstart0, gsize, Prc, Schol)
      } else {
        for (g in 1:G) { R <- gidx[[g]]; Xg <- Xre[R, , drop = FALSE]; blk <- ((g - 1) * q + 1):(g * q)
          base_prec <- if (do_car) diag(d_car[g] / sigma_re2 + 1e-6, q) else Prc
          if (!is.null(re_support_mat)) diag(base_prec) <- diag(base_prec) * re_support_mat[g, ]  # support-aware: scale RE precision by within-group participation ratio (A2 limit = pin)
          if (!is.null(re_mask_mat)) { mk <- re_mask_mat[g, ]; if (any(mk)) diag(base_prec)[mk] <- diag(base_prec)[mk] + 1e8 }  # in-draw pin: huge precision -> masked cells drawn ~0, consistently
          rhs <- temp_iter * crossprod(Xg, resid_re[R, , drop = FALSE])
          if (do_car) {                         # CAR: prior pulls each cov's RE toward its neighbours' mean
            nbr <- matrix(0, q, m); for (r in 1:q) nbr[r, ] <- car_rho * colSums(car_adj[g, ] * bre[seq(r, by = q, length.out = G), , drop = FALSE]) / sigma_re2[r]
            rhs <- rhs + nbr
          }
          bre[blk, ] <- .mn_draw_R(temp_iter * crossprod(Xg) + base_prec, rhs, Schol) }
      }
      if (use_re_spike_slab) {                  # spike-and-slab: include/prune each (group, covariate) RE block
        for (g in 1:G) for (r in 1:q) { row <- (g - 1) * q + r; R <- gidx[[g]]
          s <- sum(Xre[R, r]^2); tvec <- crossprod(Xre[R, r], resid_re[R, , drop = FALSE])   # data evidence (1 x m)
          sr <- max(sigma_re2[r], 1e-12)        # floor: avoid 1/sr = Inf -> NaN Bayes factor at high-dim RE
          logbf <- -0.5 * m * log(1 + sr * s) + 0.5 * sum((tvec %*% SigInv) * tvec) / (s + 1 / sr)
          bf <- exp(min(logbf, 700))
          p1 <- pi_re[r] * bf / (pi_re[r] * bf + (1 - pi_re[r]))
          if (!is.finite(p1)) p1 <- pi_re[r]    # guard: never let a bad BF produce NA inclusion
          delta_re[row] <- as.integer(runif(1) < p1)
          if (delta_re[row] == 0) bre[row, ] <- 0 }
        for (r in 1:q) { idxr <- seq(r, by = q, length.out = G); pi_re[r] <- rbeta(1, 1 + sum(delta_re[idxr]), 1 + G - sum(delta_re[idxr])) }
      }
      bre[!is.finite(bre)] <- 0                 # sanitize: a non-finite RE must never reach the Sigma draw
      for (r in 1:q) {                          # per-covariate RE variance (automatic relevance)
        ssr <- sum(bre[seq(r, by = q, length.out = G), ]^2)
        sigma_re2[r] <- 1 / rgamma(1, re_a + G * m / 2, re_b + ssr / 2)
      }
      # ---- RE-ASIS interweave: redraw each per-covariate scale in the non-centered ----
      # parameterization (b = sigma_r * btil, btil fixed) and MH-accept vs the IG prior.
      # Decouples sigma_r from the REs -> better joint mixing of variance and effects.
      if (re_asis) {
        RFc  <- re_fitted()
        ldig <- function(s2) -(re_a + 1) * log(s2) - re_b / s2     # IG log-kernel
        for (r in 1:q) {
          s_old <- sqrt(sigma_re2[r]); if (!is.finite(s_old) || s_old < 1e-8) next
          rows_r <- seq(r, by = q, length.out = G)
          btil <- bre[rows_r, , drop = FALSE] / s_old              # G x m standardized REs
          d <- matrix(0, n, m); contrib_r <- matrix(0, n, m)       # NCP design + current cov-r RE fit
          for (g in 1:G) { R <- gidx[[g]]; xr <- Xre[R, r]
            d[R, ] <- outer(xr, btil[g, ]); contrib_r[R, ] <- outer(xr, bre[rows_r[g], ]) }
          Rr <- resid_re - RFc + contrib_r                         # residual including ONLY cov-r's RE
          dS <- d %*% SigInv; A <- sum(dS * d); C <- sum(dS * Rr)  # NCP scalar regression on sigma_r
          if (!is.finite(A) || A <= 0) next
          s_prop <- rnorm(1, C / A, 1 / sqrt(A)); if (s_prop <= 1e-8) next
          la <- (ldig(s_prop^2) + log(s_prop)) - (ldig(sigma_re2[r]) + log(s_old))  # +Jacobian (sigma^2 -> sigma)
          if (is.finite(la) && log(runif(1)) < la) {
            f <- s_prop / s_old; bre[rows_r, ] <- bre[rows_r, ] * f
            RFc <- RFc + (f - 1) * contrib_r; sigma_re2[r] <- s_prop^2
          }
        }
      }
      RF <- re_fitted()
    }

    # (b2) BART | B, RE   — per-category backfit on the CLR residual, then zero-sum
    if (use_bart) {
      r_clr <- (W - X %*% B - RF) %*% t(V)                # n x K residual (zero-sum)
      for (k in 1:K) { bart_s[[k]]$setResponse(r_clr[, k]); f_mat[, k] <- bart_s[[k]]$run(0L, 1L)$train }
      f_mat <- f_mat - rowMeans(f_mat)                    # zero-sum across categories
      f_w <- f_mat %*% V                                  # back to W-space (offset)
    }

    # (c) Sigma | residuals   (Inverse-Wishart)
    resid <- W - X %*% B - RF - f_w
    resid[!is.finite(resid)] <- 0                                   # safety net: keep the scatter matrix finite
    Sc <- iw_scale + temp_iter * crossprod(resid); Sc <- (Sc + t(Sc)) / 2
    diag(Sc) <- diag(Sc) + 1e-6 * mean(diag(Sc)) + 1e-8            # RELATIVE ridge: caps condition number when a category collapses
    ScInv  <- tryCatch(chol2inv(chol(Sc)), error = function(e) ensure_pd(solve(Sc)))
    Wdraw  <- rWishart(1, iw_df + temp_iter * n, ScInv)[, , 1]; Wdraw <- (Wdraw + t(Wdraw)) / 2
    Sigma  <- tryCatch(chol2inv(chol(Wdraw)), error = function(e) ensure_pd(solve(Wdraw)))
    Sigma  <- (Sigma + t(Sigma)) / 2
    Schol  <- tryCatch(chol(Sigma), error = function(e) chol(ensure_pd(Sigma)))
    SigInv <- chol2inv(Schol)

    # (d) horseshoe (group / regularized) | B
    if (use_horseshoe && length(hs_idx) > 0) {
      ss <- rowSums(B[hs_idx, , drop = FALSE]^2)
      nu[hs_idx]      <- 1 / rgamma(length(hs_idx), 1, 1 + 1 / lambda2[hs_idx])
      lambda2[hs_idx] <- pmin(pmax(1 / rgamma(length(hs_idx), (m + 1) / 2, 1 / nu[hs_idx] + ss / (2 * tau2)), 1e-10), 1e6)
      xi   <- 1 / rgamma(1, 1, 1 + 1 / tau2)
      tau2 <- min(max(1 / rgamma(1, (length(hs_idx) * m + 1) / 2, 1 / xi + sum(ss / lambda2[hs_idx]) / 2), 1e-10), 1e6)
    }

    # (e) store (CLR / zero-sum K-space, raw X scale)
    if (it > nburn) {
      s <- it - nburn
      Bstd <- B %*% t(V); Bc <- Bstd / Xsd
      if (length(int_idx) > 0) { adj <- colSums((Xmu / Xsd) * Bstd)
        Bc[int_idx, ] <- sweep(Bc[int_idx, , drop = FALSE], 2, adj, "-") }
      postB[, , s] <- Bc; postSigma[, , s] <- V %*% Sigma %*% t(V)
      if (use_re) { post_sre[, s] <- sigma_re2; post_bre_sum <- post_bre_sum + bre %*% t(V) }
      if (use_bart) post_f_sum <- post_f_sum + f_mat
      if (use_re_spike_slab) post_delta_sum <- post_delta_sum + delta_re
    }
    if (verbose && it %% max(1, floor(niter / 10)) == 0) cat(sprintf("  mvclr iter %d/%d\n", it, niter))
  }

  result <- list(postB = postB, postSigma = postSigma,
       B_mean = apply(postB, c(1, 2), mean), Sigma_mean = apply(postSigma, c(1, 2), mean),
       post_sigma_re2 = if (use_re) post_sre else NULL,
       re_mean = if (use_re) post_bre_sum / nret else NULL,   # (G*q) x K, blocked by group (CLR)
       post_f_mean = if (use_bart) post_f_sum / nret else NULL, # n x K nonlinear effect (CLR, zero-sum)
       post_delta_mean = if (use_re_spike_slab) post_delta_sum / nret else NULL,
       final_state = list(B = B, bre = bre, Sigma = Sigma, sigma_re2 = sigma_re2, lambda2 = lambda2, tau2 = tau2),
       re_idx = if (use_re) re_idx else NULL, re_groups = if (use_re) grp_names else NULL,
       row_order = if (use_re && !is.null(group_idx)) ord else seq_len(n),  # n-row outputs are in this order
       cov_names = cov_names, cat_names = cat_names, K = K, p = p, m = m, n = n,
       used_rcpp = rcpp)
  if (isTRUE(save_to_disk) && !is.null(disk_path)) {           # disk-backed storage
    if (!dir.exists(disk_path)) dir.create(disk_path, recursive = TRUE)
    qs2::qs_save(result, file.path(disk_path, "mvclr_posterior.qs"))
    qs2::qs_save(result$final_state, file.path(disk_path, "mvclr_final_state.qs"))
  }
  result
}

# Recover a saved CLR posterior (round-trip with save_to_disk).
recover_mvclr_posterior <- function(disk_path) qs2::qs_read(file.path(disk_path, "mvclr_posterior.qs"))

`%||%` <- function(a, b) if (!is.null(a)) a else b
