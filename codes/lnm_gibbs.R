# =============================================================================
# lnm_gibbs.R — Logistic-Normal Multinomial sampler  [v2: Rcpp + hierarchy]
# =============================================================================
# PG-conjugate cousin of the Dirichlet-Multinomial. Multinomial on the LEVELS
# count N_i (= area; PG stays on the fast large-h path), with:
#   * per-pixel OVERDISPERSION  u_i ~ N(0, Sigma_od)  (cures the levels explosion;
#     Sigma_od = the identified DM-phi; u_i are PG augmentation, marginalized)
#   * hierarchical COUNTRY random effects (intercepts OR slopes via re_idx)
#   * group/CLR HORSESHOE on the fixed effects
#   * SYMMETRIC / zero-sum coefficient output (no privileged baseline)
#   * RcppArmadillo accelerated u / RE / beta draws (pure-R fallback)
#
#     eta_ij = X_i B_j + X^re_i b_{g(i),j} + u_ij ,  Y_i ~ Mult(N_i, softmax([eta_i,0]))
# Polson-Scott-Windle multinomial PG augmentation, baseline = last category.
# =============================================================================
suppressMessages(library(pg))
`%||%` <- function(a, b) if (!is.null(a)) a else b
.lnm_reg_hs <- function(v, c2) c2 * v / (c2 + v)

# Separation detector (port of mnlogit_rcpp_sym): for (covariate v, category j),
# the 5-95% support OVERLAP between v|"j present" and v|"j absent". Returns a
# continuous score [p x m]: 0 = disjoint (separated) .. 1 = full overlap (healthy).
# Drives the 3-tier prior (hard mask / soft graduated shrinkage / healthy).
.lnm_detect_separation <- function(X, Y, m, int_idx) {
  p <- ncol(X); ov <- matrix(1, p, m); cols <- setdiff(seq_len(p), int_idx)
  for (j in seq_len(m)) {
    is_j <- Y[, j] > 0; n1 <- sum(is_j); n0 <- sum(!is_j)
    if (n1 < 5 || n0 < 5) { ov[cols, j] <- 0; next }                      # absent / always-present
    for (v in cols) {
      v1 <- X[is_j, v]; va <- X[!is_j, v]
      q1l <- quantile(v1, .05, names = FALSE); q1h <- quantile(v1, .95, names = FALSE)
      qbl <- quantile(va, .05, names = FALSE); qbh <- quantile(va, .95, names = FALSE)
      inter <- min(q1h, qbh) - max(q1l, qbl); uni <- max(q1h, qbh) - min(q1l, qbl)
      ov[v, j] <- if (var(v1) < 1e-4) 0 else max(inter / (uni + 1e-9), 0)
    }
  }
  ov
}

.lnm_have_rcpp <- FALSE; .lnm_have_core <- FALSE
try({ suppressMessages({ library(Rcpp); library(RcppArmadillo) })
  Rcpp::sourceCpp("codes/lnm_core.cpp"); .lnm_have_rcpp <- TRUE }, silent = TRUE)
try({ Rcpp::sourceCpp("codes/lnm_gibbs_core.cpp"); .lnm_have_core <- TRUE }, silent = TRUE)
.lnm_have_symcore <- FALSE
try({ Rcpp::sourceCpp("codes/lnm_gibbs_core_sym.cpp"); .lnm_have_symcore <- TRUE }, silent = TRUE)

lnm_gibbs <- function(Y, X, group_idx = NULL,
                      niter = 1500, nburn = 750,
                      use_horseshoe = TRUE, horseshoe_idx = NULL,
                      use_re = TRUE, re_idx = NULL,
                      iw_df = NULL, iw_scale = NULL, A0 = 10, slab_s2 = 4,
                      re_a = 1, re_b = 1, standardize = TRUE, thin = 1, seed = 1,
                      const_sum_blocks = NULL, use_asis = FALSE, symmetric = FALSE, use_tempering = FALSE, tempering_T0 = 0.05,
                      use_overdispersion = TRUE, beta_clamp = NULL, separation_mask = FALSE,
                      separation_soft = FALSE, sep_overlap_hard = 0.0, sep_overlap_neutral = 0.5,
                      use_bart = FALSE, bart_idx = NULL, n_trees_bart = 50,
                      positive_idx = NULL, negative_idx = NULL, init_state = NULL,
                      save_to_disk = FALSE, disk_path = NULL, use_re_spike_slab = FALSE,
                      use_car = FALSE, car_adj = NULL, car_rho = 0.95,
                      use_rcpp = TRUE, use_cpp_core = TRUE, re_screen = FALSE,
                      support_prior_strength = 0, re_asis = FALSE, init_jitter = 0, verbose = FALSE, ...) {
  rcpp <- use_rcpp && .lnm_have_rcpp
  Y <- as.matrix(Y); n <- nrow(Y); K <- ncol(Y); m <- K - 1L
  cat_names <- colnames(Y) %||% paste0("cat", 1:K)
  N <- rowSums(Y); kappa <- Y[, 1:m, drop = FALSE] - N / 2

  X <- as.matrix(X); p <- ncol(X); cov_names <- colnames(X) %||% paste0("V", 1:p)

  # Sum-to-constant (compositional) blocks are handled DRIVER-LEVEL now (order-invariant drop-one
  # before dispatch + zero-sum reconstruction in build_fit_object), uniformly across all samplers --
  # so LNM no longer does its own Helmert-contrast carving (removed: order-DEPENDENT, and redundant).
  # `const_sum_blocks` is accepted for signature compatibility but ignored.
  const_bt <- NULL
  int_idx <- which(apply(X, 2, function(c) all(c == c[1])))
  Xmu <- rep(0, p); Xsd <- rep(1, p)
  if (standardize) for (j in setdiff(1:p, int_idx)) if (sd(X[, j]) > 1e-6) {   # floor: leave near-constant
    Xmu[j] <- mean(X[, j]); Xsd[j] <- sd(X[, j]); X[, j] <- (X[, j] - Xmu[j]) / Xsd[j] }   # contrasts (rare-class pairs) unstandardized -> horseshoe shrinks, no /tiny blow-up

  # group-sort rows for the RE block draw
  if (use_re && !is.null(group_idx)) {
    gf0 <- as.integer(factor(group_idx)); grp_names <- levels(factor(group_idx))
    ord <- order(gf0); X <- X[ord, ]; Y <- Y[ord, ]; N <- N[ord]; kappa <- kappa[ord, , drop = FALSE]; gf <- gf0[ord]
    G <- max(gf); gsize <- as.integer(table(gf)); gstart0 <- as.integer(cumsum(c(0, gsize[-G]))); gidx <- split(seq_len(n), gf)
    if (is.null(re_idx)) re_idx <- if (length(int_idx)) int_idx else 1L
    q <- length(re_idx); Xre <- X[, re_idx, drop = FALSE]
  } else { use_re <- FALSE; G <- 0L; q <- 0L; ord <- seq_len(n) }

  hs_idx <- setdiff(if (is.null(horseshoe_idx)) 1:p else horseshoe_idx, int_idx)
  if (is.null(iw_df)) iw_df <- m + 2; if (is.null(iw_scale)) iw_scale <- diag(m)

  # ---- FULLY SYMMETRIC (baseline-free) core ----------------------------------
  if (symmetric && .lnm_have_symcore) {
    Xre_a <- if (use_re) Xre else matrix(0, n, 0)
    gst_a <- if (use_re) as.integer(gstart0) else integer(0)
    gsz_a <- if (use_re) as.integer(gsize) else integer(0)
    # Deterministic RE-identifiability screen (A2 no-X-var + A3 collinearity); 0/1 per bre row
    # (g*q + r), pinned inside the C++ RE loop. Bounds the unidentified RE slopes.
    re_pin_a <- integer(max(as.integer(G) * q, 1L))
    re_support_a <- rep(1, max(as.integer(G) * q, 1L))   # within-group participation-ratio support factor (>=1), per bre row g*q+r
    if (use_re && re_screen && q >= 1) {
      re_pin_a <- as.integer(as.vector(t(screen_re_design(X, gf, re_idx))))
      if (verbose) cat(sprintf("  RE screen pinned %d of %d (group x covariate) slopes\n", sum(re_pin_a), G * q))
    }
    if (use_re && support_prior_strength > 0 && q >= 1) {
      re_support_a <- as.vector(t(re_support_design(X, gf, re_idx, support_prior_strength)))
      if (verbose) cat(sprintf("  RE support prior: max factor %.1f (participation ratio, strength %.2f)\n", max(re_support_a), support_prior_strength))
    }
    fe_support_a <- fe_support_factor(X, support_prior_strength)   # global participation-ratio FE shrinkage (length p)
    core <- lnm_gibbs_core_sym(X, Y, N, Xre_a, gst_a, gsz_a, as.integer(hs_idx - 1L),
                               niter, nburn, A0, slab_s2, K + 2, diag(K), re_a, re_b,
                               use_re, use_horseshoe && length(hs_idx) > 0, thin, seed, re_pin_a, re_support_a, fe_support_a, re_asis,
                               use_tempering, tempering_T0, init_jitter)
    btS <- function(Bs) { Bc <- Bs / Xsd                            # back-transform; already zero-sum
      if (length(int_idx) > 0) { adj <- colSums((Xmu / Xsd) * Bs); Bc[int_idx, ] <- sweep(Bc[int_idx, , drop = FALSE], 2, adj, "-") }
      Bc }
    Bm <- btS(matrix(core$B_sum / core$nret, p, K)); dimnames(Bm) <- list(cov_names, cat_names)
    nk <- ncol(core$Bdraws); Bd <- array(0, c(p, K, nk), dimnames = list(cov_names, cat_names, NULL))
    for (s in 1:nk) Bd[, , s] <- btS(matrix(core$Bdraws[, s], p, K))
    return(list(B_mean = Bm, B_interp = Bm, Sigma_mean = core$Sig_sum / core$nret,
                post_maxB = core$maxB, post_overdispersion = core$od, Bsym_draws = Bd,
                re_mean = if (use_re) core$bre_mean else NULL, re_groups = if (use_re) grp_names else NULL,
                row_order = ord, cov_names = cov_names, cat_names = cat_names, K = K, p = p, n = n,
                symmetric = TRUE, used_core = TRUE, used_rcpp = TRUE))
  }

  # ---- FULL C++ GIBBS CORE (whole loop in C++) -------------------------------
  if (use_cpp_core && .lnm_have_core) {
    Xre_a  <- if (use_re) Xre else matrix(0, n, 0)
    gst_a  <- if (use_re) as.integer(gstart0) else integer(0)
    gsz_a  <- if (use_re) as.integer(gsize) else integer(0)
    core <- lnm_gibbs_core(X, Y[, 1:m, drop = FALSE], N, Xre_a, gst_a, gsz_a,
                           as.integer(hs_idx - 1L), niter, nburn, A0, slab_s2, iw_df, iw_scale,
                           re_a, re_b, use_re, use_horseshoe && length(hs_idx) > 0, thin, seed, use_asis, use_tempering, tempering_T0)
    bt_sym <- function(Bs) {                             # baseline-coded std -> back-transformed symmetric
      Bf <- cbind(Bs, 0) / Xsd
      if (length(int_idx) > 0) { adj <- colSums((Xmu / Xsd) * cbind(Bs, 0)); Bf[int_idx, ] <- sweep(Bf[int_idx, , drop = FALSE], 2, adj, "-") }
      Bf - rowMeans(Bf)
    }
    Bsym <- bt_sym(core$B_sum / core$nret); dimnames(Bsym) <- list(cov_names, cat_names)
    nk <- ncol(core$Bdraws)
    Bsym_draws <- array(0, c(p, K, nk), dimnames = list(cov_names, cat_names, NULL))
    for (s in 1:nk) Bsym_draws[, , s] <- bt_sym(matrix(core$Bdraws[, s], p, m))
    # back-map contrast blocks -> interpretable per-category sum-to-zero effects
    backmap <- function(Bm) {
      if (is.null(const_bt)) return(Bm)
      bp <- unlist(lapply(const_bt, function(b) b$new_pos)); keep <- setdiff(seq_len(nrow(Bm)), bp)
      parts <- list(Bm[keep, , drop = FALSE]); rn <- list(rownames(Bm)[keep])
      for (nm in names(const_bt)) { b <- const_bt[[nm]]
        parts <- c(parts, list(b$C %*% Bm[b$new_pos, , drop = FALSE])); rn <- c(rn, list(b$orig_names)) }
      M <- do.call(rbind, parts); rownames(M) <- unlist(rn); colnames(M) <- cat_names; M
    }
    return(list(B_mean = Bsym, B_interp = backmap(Bsym), Sigma_mean = core$Sig_sum / core$nret,
                post_maxB = core$maxB, post_overdispersion = core$od,
                Bsym_draws = Bsym_draws, oddraws = core$oddraws,
                re_mean = if (use_re) core$bre_mean else NULL, tau_re2 = if (use_re) as.vector(core$tau_re2) else NULL,
                re_idx = if (use_re) re_idx else NULL,
                re_groups = if (use_re) grp_names else NULL, row_order = ord,
                cov_names = cov_names, cat_names = cat_names, K = K, p = p, m = m, n = n,
                used_rcpp = TRUE, used_core = TRUE))
  }

  B <- matrix(0, p, m); u <- matrix(0, n, m); Sigma <- diag(m); SigInv <- diag(m)
  bre <- matrix(0, max(G * q, 1L), m); tau_re2 <- rep(1, m)
  lambda2 <- rep(1, p); nu <- rep(1, p); tau2 <- 1; xi <- 1; c2 <- slab_s2; prior_prec <- rep(1 / A0^2, p)
  if (!is.null(init_state)) {                                    # HOT-START
    for (nm in c("B","u","bre","tau_re2","lambda2","tau2")) if (!is.null(init_state[[nm]])) assign(nm, init_state[[nm]])
    if (!is.null(init_state$Sigma)) { Sigma <- init_state$Sigma; SigInv <- chol2inv(chol(Sigma)) }
  }
  pos_c <- positive_idx; neg_c <- negative_idx                   # sign-constraint covariate indices
  delta_re <- matrix(1L, max(G, 1L), m); pi_re <- rep(0.5, m); post_delta_sum <- matrix(0, max(G, 1L), m)  # spike-slab RE
  do_car <- isTRUE(use_car) && !is.null(car_adj) && use_re                    # CAR spatial RE
  d_car <- if (do_car) rowSums(car_adj) else NULL
  post_bre_sum <- matrix(0, max(G * q, 1L), m)                                # RE posterior mean
  # BART setup (per-category dbarts on the PG working residual; centered for zero-sum)
  do_bart <- isTRUE(use_bart) && !is.null(bart_idx)
  f_mat <- matrix(0, n, m); post_f_sum <- matrix(0, n, m)
  if (do_bart) {
    if (!requireNamespace("dbarts", quietly = TRUE)) stop("use_bart needs the dbarts package")
    Xb <- X[, bart_idx, drop = FALSE]
    bart_s <- lapply(1:m, function(j) dbarts::dbarts(Xb, rnorm(n), weights = rep(1, n),
      control = dbarts::dbartsControl(n.samples = 1, n.burn = 0, n.trees = n_trees_bart, keepTrees = TRUE, n.chains = 1, updateState = TRUE)))
  }

  re_fitted <- function() {
    if (!use_re) return(matrix(0, n, m))
    RF <- matrix(0, n, m)
    for (g in 1:G) { R <- gidx[[g]]; RF[R, ] <- Xre[R, , drop = FALSE] %*% bre[((g - 1) * q + 1):(g * q), , drop = FALSE] }
    RF
  }
  # 3-tier separation prior: overlap -> per-cell precision MULTIPLIER (hard pin / soft graduated / healthy)
  sep_penalty <- NULL
  if (isTRUE(separation_mask) || isTRUE(separation_soft)) {
    ov <- .lnm_detect_separation(X, Y, m, int_idx); sep_penalty <- matrix(1, p, m)
    hard <- ov < sep_overlap_hard | ov <= 0
    soft <- !hard & ov < sep_overlap_neutral & isTRUE(separation_soft)
    sep_penalty[soft] <- 1 + 50 * (sep_overlap_neutral - ov[soft]) / (sep_overlap_neutral - sep_overlap_hard + 1e-9)
    sep_penalty[hard] <- 1e10                                     # hard: pin to ~0
    if (verbose) cat(sprintf("  separation: %d hard + %d soft of %d fixed cells (tier: %s)\n",
          sum(hard), sum(soft), length(ov), if (isTRUE(separation_soft)) "hard+soft" else "hard"))
  }
  re_sep_mask <- NULL
  if ((isTRUE(separation_mask) || isTRUE(separation_soft)) && use_re) {   # separated country x category REs (q=1)
    re_sep_mask <- matrix(FALSE, G, m)
    for (gg in 1:G) { R <- gidx[[gg]]; for (j in 1:m) if (sum(Y[R, j] > 0) < 5) re_sep_mask[gg, j] <- TRUE }
  }
  beta_draw <- function(r) if (rcpp && is.null(sep_penalty)) lnm_beta_draw_cpp(X, r, om, prior_prec) else {
    Bn <- matrix(0, p, m); for (j in 1:m) { pp <- if (is.null(sep_penalty)) prior_prec else prior_prec * sep_penalty[, j]
      P <- crossprod(X, X * om[, j]) + diag(pp, p)
      L <- chol(P); Bn[, j] <- backsolve(L, backsolve(L, crossprod(X, om[, j] * r[, j]), transpose = TRUE) + rnorm(p)) }; Bn }
  u_draw <- function(rb) if (rcpp) lnm_u_draw_cpp(om, rb, SigInv) else {
    un <- matrix(0, n, m); for (i in 1:n) { Li <- chol(diag(om[i, ], m) + SigInv)
      un[i, ] <- backsolve(Li, backsolve(Li, om[i, ] * rb[i, ], transpose = TRUE) + rnorm(m)) }; un }

  nret <- niter - nburn
  postB <- array(0, c(p, K, nret), dimnames = list(cov_names, cat_names, NULL))
  postSig <- array(0, c(m, m, nret)); post_od <- numeric(nret); post_maxB <- numeric(nret)

  for (it in 1:niter) {
    RF <- re_fitted()
    eta <- X %*% B + RF + u + f_mat
    expe <- exp(pmin(eta, 30)); denom_all <- 1 + rowSums(expe)
    om <- matrix(0, n, m); z <- matrix(0, n, m)
    for (j in 1:m) {
      cj <- log(pmax(denom_all - expe[, j], 1e-12))
      om[, j] <- pmax(suppressWarnings(rpg_hybrid(N, eta[, j] - cj)), 1e-9)
      z[, j] <- kappa[, j] / om[, j] + cj
    }
    if (use_horseshoe && length(hs_idx) > 0)
      prior_prec[hs_idx] <- 1 / pmax(.lnm_reg_hs(tau2 * lambda2[hs_idx], c2), 1e-12)

    B <- beta_draw(z - RF - u - f_mat)                           # fixed effects
    if (!is.null(beta_clamp)) B <- pmax(pmin(B, beta_clamp), -beta_clamp)   # crude separation cap
    if (!is.null(pos_c)) B[pos_c, ] <- abs(B[pos_c, ])           # sign constraints (reflect to enforce)
    if (!is.null(neg_c)) B[neg_c, ] <- -abs(B[neg_c, ])
    if (use_re) {                                                # country RE
      resid_re <- z - X %*% B - u - f_mat
      if (do_car) {                                              # CAR: smooth REs across neighbours (q=1 intercept)
        for (g in 1:G) { R <- gidx[[g]]
          for (j in 1:m) { sj <- sum(om[R, j]); tj <- sum(om[R, j] * resid_re[R, j])
            nbr <- sum(car_adj[g, ] * bre[, j]); prec <- sj + d_car[g] / tau_re2[j]
            bre[g, j] <- rnorm(1, (tj + (car_rho / tau_re2[j]) * nbr) / prec, 1 / sqrt(prec)) } }
      } else { tau_re2_d <- pmin(tau_re2, 1e6)   # floor RE precision (1/tau >= 1e-6) so per-group chol stays PD with many random slopes / small groups
        bre <- if (rcpp) lnm_re_draw_cpp(Xre, resid_re, om, gstart0, gsize, tau_re2_d) else {
        bn <- bre; for (g in 1:G) { R <- gidx[[g]]; Xg <- Xre[R, , drop = FALSE]
          for (j in 1:m) { P <- crossprod(Xg, Xg * om[R, j]) + diag(1 / tau_re2_d[j], q)
            L <- chol(P); bn[((g - 1) * q + 1):(g * q), j] <- backsolve(L, backsolve(L, crossprod(Xg, om[R, j] * resid_re[R, j]), transpose = TRUE) + rnorm(q)) } }; bn } }
      if (!is.null(re_sep_mask)) bre[re_sep_mask] <- 0                    # pin separated REs (q=1)
      if (use_re_spike_slab) {                                           # spike-and-slab RE inclusion (q=1 intercept)
        for (g in 1:G) { R <- gidx[[g]]
          for (j in 1:m) {
            sj <- sum(om[R, j]); tj <- sum(om[R, j] * resid_re[R, j]); t2 <- tau_re2[j]
            bf <- sqrt(1 / (1 + t2 * sj)) * exp(min(tj^2 * t2 / (2 * (1 + t2 * sj)), 700))
            p1 <- pi_re[j] * bf / (pi_re[j] * bf + (1 - pi_re[j]))
            delta_re[g, j] <- as.integer(runif(1) < p1)
            if (delta_re[g, j] == 0) bre[(g - 1) * q + 1, j] <- 0 } }
        for (j in 1:m) pi_re[j] <- rbeta(1, 1 + sum(delta_re[, j]), 1 + G - sum(delta_re[, j]))
      }
      for (j in 1:m) {                                                   # tau_re update (CAR quadratic form if spatial)
        rate_j <- if (do_car) re_b + 0.5 * (sum(d_car * bre[, j]^2) - car_rho * sum(bre[, j] * (car_adj %*% bre[, j]))) else re_b + sum(bre[, j]^2) / 2
        tau_re2[j] <- min(1 / rgamma(1, re_a + G * q / 2, max(rate_j, 1e-6)), 100) }
      RF <- re_fitted()
    }
    if (use_overdispersion) {                                    # overdispersion (skip at N~1: unidentified, re-creates ridge)
      u <- u_draw(z - X %*% B - RF - f_mat)
      Sigma  <- chol2inv(chol(rWishart(1, iw_df + n, chol2inv(chol(iw_scale + crossprod(u))))[, , 1]))
      SigInv <- chol2inv(chol(Sigma))
    }
    if (do_bart) {                                               # nonlinear terms (Bayesian backfitting on PG residual)
      rbart <- z - X %*% B - RF - u
      for (j in 1:m) { bart_s[[j]]$setWeights(om[, j]); bart_s[[j]]$setSigma(1)   # PG: var = 1/om exactly -> sigma fixed at 1
        bart_s[[j]]$setResponse(rbart[, j]); f_mat[, j] <- bart_s[[j]]$run(0L, 1L)$train }
      f_mat <- f_mat - rowMeans(f_mat)                           # zero-sum (symmetric)
    }

    if (use_horseshoe && length(hs_idx) > 0) {
      ss <- rowSums(B[hs_idx, , drop = FALSE]^2)
      nu[hs_idx]      <- 1 / rgamma(length(hs_idx), 1, 1 + 1 / lambda2[hs_idx])
      lambda2[hs_idx] <- pmin(pmax(1 / rgamma(length(hs_idx), (m + 1) / 2, 1 / nu[hs_idx] + ss / (2 * tau2)), 1e-10), 1e6)
      xi   <- 1 / rgamma(1, 1, 1 + 1 / tau2)
      tau2 <- min(max(1 / rgamma(1, (length(hs_idx) * m + 1) / 2, 1 / xi + sum(ss / lambda2[hs_idx]) / 2), 1e-10), 1e6)
    }

    if (it > nburn) {
      s <- it - nburn
      Bfull <- cbind(B, 0); Bc <- Bfull / Xsd
      if (length(int_idx) > 0) { adj <- colSums((Xmu / Xsd) * Bfull); Bc[int_idx, ] <- sweep(Bc[int_idx, , drop = FALSE], 2, adj, "-") }
      postB[, , s] <- Bc - rowMeans(Bc)                          # SYMMETRIC / zero-sum
      postSig[, , s] <- Sigma; post_od[s] <- mean(diag(Sigma)); post_maxB[s] <- max(abs(B))
      if (do_bart) post_f_sum <- post_f_sum + f_mat
      if (use_re_spike_slab) post_delta_sum <- post_delta_sum + delta_re
      if (use_re) post_bre_sum <- post_bre_sum + bre
    }
    if (verbose && it %% max(1, floor(niter / 10)) == 0)
      cat(sprintf("  lnm iter %d/%d  max|B|=%.2f  od=%.3f%s\n", it, niter, max(abs(B)), mean(diag(Sigma)), if (use_re) sprintf(" tau_re=%.2f", mean(tau_re2)) else ""))
  }
  result <- list(postB = postB, B_mean = apply(postB, c(1, 2), mean),
       postSigma = postSig, Sigma_mean = apply(postSig, c(1, 2), mean),
       post_overdispersion = post_od, post_maxB = post_maxB,
       post_f_mean = if (do_bart) post_f_sum / nret else NULL,
       post_delta_mean = if (use_re_spike_slab) post_delta_sum / nret else NULL,
       final_state = list(B = B, u = u, bre = bre, Sigma = Sigma, tau_re2 = tau_re2, lambda2 = lambda2, tau2 = tau2),
       re_mean = if (use_re) (post_bre_sum / nret) else NULL, re_idx = if (use_re) re_idx else NULL,
       re_groups = if (use_re) grp_names else NULL, row_order = ord,
       cov_names = cov_names, cat_names = cat_names, K = K, p = p, m = m, n = n, used_rcpp = rcpp)
  if (isTRUE(save_to_disk) && !is.null(disk_path)) {             # disk-backed storage
    if (!dir.exists(disk_path)) dir.create(disk_path, recursive = TRUE)
    qs2::qs_save(result, file.path(disk_path, "lnm_posterior.qs"))
    qs2::qs_save(result$final_state, file.path(disk_path, "lnm_final_state.qs"))
  }
  result
}

# Recover a saved LNM posterior (round-trip with save_to_disk).
recover_lnm_posterior <- function(disk_path) qs2::qs_read(file.path(disk_path, "lnm_posterior.qs"))
