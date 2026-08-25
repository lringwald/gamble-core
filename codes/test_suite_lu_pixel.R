# =============================================================================
# test_suite_lu_pixel.R
#
# Feature-specific validation for the SYMMETRIC additions in mnlogit_rcpp_sym.R
# / mnlogit_gibbs_core_sym.cpp — the parts that distinguish this sampler from
# the adaptive base:
#
#   Tier 3  Symmetric Horseshoe      (zero-sum coupling Msym, Helmert basis,
#                                      draw_beta_symhs_pooled, c_v * Msym)
#   Tier 2  Symmetric pooled RE var  (update_re_precision_hc_sym: one variance
#                                      shared symmetrically across categories)
#   Baseline invariance              (the motivation for symmetric coding)
#
# Layout:
#   PART A  Exact unit tests of the symmetric kernels (machine precision / MC)
#   PART B  Symmetric feature behaviour in the full sampler (end-to-end)
#
# Run:  Rscript codes/test_suite_lu_pixel.R
# Exit code is non-zero on any failure.
# =============================================================================

suppressMessages({ library(Rcpp); library(RcppArmadillo) })

`%||%` <- function(a, b) if (!is.null(a)) a else b
find_file <- function(fname) for (d in c("codes", ".")) {
  p <- file.path(d, fname); if (file.exists(p)) return(normalizePath(p))
} %||% stop("Cannot locate ", fname)
r_path <- find_file("mnlogit_rcpp_sym.R")

.tests_run <- 0L; .tests_failed <- 0L; .failures <- character(0)
check <- function(name, condition, detail = "") {
  .tests_run <<- .tests_run + 1L; ok <- isTRUE(condition)
  if (!ok) { .tests_failed <<- .tests_failed + 1L; .failures <<- c(.failures, paste(name, detail)) }
  cat(sprintf("  [%s] %-54s %s\n", if (ok) "PASS" else "FAIL", name, detail)); invisible(ok)
}
approx <- function(a, b, tol = 1e-10) max(abs(as.numeric(a) - as.numeric(b))) < tol
maxerr <- function(a, b) max(abs(as.numeric(a) - as.numeric(b)))

cat("================================================================\n")
cat(" test_suite_lu_pixel.R — symmetric-feature validation\n")
cat(sprintf(" BLAS: %s\n", extSoftVersion()["BLAS"]))
cat("================================================================\n\n")

cat(">>> Sourcing sampler (compiles symmetric C++ core)...\n")
suppressMessages(suppressWarnings(source(r_path)))   # defines draw_beta_symhs_pooled + compiles cpp
cat("    ready.\n\n")

# =============================================================================
# PART A — Exact unit tests of the symmetric kernels
# =============================================================================
cat("PART A — symmetric kernel unit tests\n")
set.seed(11)

## A1. Msym is the zero-sum quadratic form in baseline coordinates.
##     For baseline-coded beta (length p), beta' Msym beta must equal the sum of
##     squared coefficients after the FULL (p_all) vector is centered to zero-sum.
{
  ok <- TRUE; worst <- 0
  for (p_all in 3:7) {
    p <- p_all - 1
    Msym <- diag(1, p) - matrix(1 / p_all, p, p)        # exactly as in the sampler
    for (rep in 1:50) {
      beta <- rnorm(p)
      full <- c(beta, 0); centered <- full - mean(full)
      lhs <- as.numeric(t(beta) %*% Msym %*% beta)
      rhs <- sum(centered^2)
      worst <- max(worst, abs(lhs - rhs)); if (abs(lhs - rhs) > 1e-10) ok <- FALSE
    }
  }
  check("A1 Msym == zero-sum sum-of-squares", ok, sprintf("worst=%.1e", worst))
}

## A2. Helmert basis used for the symmetric HS projection is orthonormal and
##     spans the zero-sum subspace (columns sum to zero).
{
  ok <- TRUE; worst_orth <- 0; worst_zs <- 0
  for (p_all in 3:7) {
    Cmat <- contr.helmert(p_all)
    Cmat <- sweep(Cmat, 2, sqrt(colSums(Cmat^2)), "/")  # exactly as in the sampler
    worst_orth <- max(worst_orth, maxerr(crossprod(Cmat), diag(p_all - 1)))
    worst_zs <- max(worst_zs, max(abs(colSums(Cmat))))
  }
  check("A2 Helmert basis orthonormal", worst_orth < 1e-10, sprintf("err=%.1e", worst_orth))
  check("A2 Helmert basis zero-sum cols", worst_zs < 1e-10, sprintf("err=%.1e", worst_zs))
}

## A3. update_re_precision_hc_sym — symmetric pooling INVARIANTS.
##     Defining property: a single variance is shared symmetrically across all
##     p categories for each RE predictor (and across all members of a block).
{
  k <- 5; p_all <- 4; p <- p_all - 1; ng <- 6
  mu <- matrix(rnorm(k * p), k, p)
  beta_c <- array(rnorm(k * p * ng, 0, 0.5), c(k, p, ng)) + array(mu, c(k, p, ng))
  re_idx <- as.integer(0:(k - 1))
  re_mask <- array(1.0, c(k, p, ng)); y_mask <- matrix(1.0, p, ng)
  is_int <- as.integer(c(1, rep(0, k - 1)))
  prec_prev <- matrix(1.0, k, p); a_prev <- matrix(0.5, k, p)

  ## (a) non-block: prec constant across categories for each RE predictor
  out <- update_re_precision_hc_sym(beta_c, mu, re_idx, ng, re_mask, y_mask, is_int,
                                    prec_prev, a_prev, p_all, 1.0, NULL, NULL)
  row_spread <- max(apply(out$prec[2:k, , drop = FALSE], 1, function(r) diff(range(r))))
  check("A3a sym RE var constant across categories", row_spread < 1e-12,
        sprintf("max within-row spread=%.1e", row_spread))

  ## (b) block: all members of a const-sum block share one variance
  block_id <- rep(-1L, k); block_id[c(2, 3)] <- 0L      # cols 2,3 (1-based) form a block
  block_size <- as.integer(3)
  outb <- update_re_precision_hc_sym(beta_c, mu, re_idx, ng, re_mask, y_mask, is_int,
                                     prec_prev, a_prev, p_all, 1.0, block_id, block_size)
  block_spread <- diff(range(outb$prec[c(2, 3), ]))     # identical across both members & cats
  check("A3b block members share one variance", block_spread < 1e-12,
        sprintf("spread=%.1e", block_spread))

  ## (c) structurally masked predictor keeps its previous precision unchanged
  re_mask_m <- re_mask; re_mask_m[1, , ] <- 0.0
  prec_prev2 <- prec_prev; prec_prev2[1, ] <- 7.0
  outm <- update_re_precision_hc_sym(beta_c, mu, re_idx, ng, re_mask_m, y_mask, is_int,
                                     prec_prev2, a_prev, p_all, 1.0, NULL, NULL)
  check("A3c masked predictor precision preserved", approx(outm$prec[1, ], rep(7.0, p), 1e-12),
        sprintf("err=%.1e", maxerr(outm$prec[1, ], rep(7.0, p))))

  ## (d) variance pooling responds to data: smaller cross-group spread must yield
  ##     higher pooled precision (E[prec | beta_c] averaged over the Gamma draw).
  mean_prec <- function(sd_dev, M = 150) {
    set.seed(99)
    bc <- array(mu, c(k, p, ng)) + array(rnorm(k * p * ng, 0, sd_dev), c(k, p, ng))
    acc <- 0
    for (r in 1:M) acc <- acc + mean(update_re_precision_hc_sym(
      bc, mu, re_idx, ng, re_mask, y_mask, is_int, prec_prev, a_prev, p_all, 1.0, NULL, NULL)$prec[2:k, ])
    acc / M
  }
  p_small <- mean_prec(0.1); p_large <- mean_prec(1.5)
  check("A3d precision higher for smaller RE spread", p_small > p_large,
        sprintf("prec(sd=.1)=%.1f vs prec(sd=1.5)=%.2f", p_small, p_large))
}

## A4. draw_beta_symhs_pooled — the symmetric HS joint draw.
##     MC mean must equal P_joint^{-1} Pb_joint, where P_joint carries the
##     across-equation coupling c_v * Msym. Validates the coupled assembly.
{
  set.seed(22)
  n <- 200; k <- 4; p_all <- 4; p <- p_all - 1; baseline <- p_all
  pp <- as.integer((1:p_all)[-baseline])
  X <- matrix(rnorm(n * k), n, k); Xt <- t(X)
  omega <- matrix(runif(n * p, 0.5, 1.5), n, p)
  c_j_mat <- matrix(rnorm(n * p), n, p)
  kappa_w <- matrix(rnorm(n * p_all), n, p_all)
  prior_P <- diag(0.5, k); prior_Pb <- matrix(0, k, p)
  Msym <- diag(1, p) - matrix(1 / p_all, p, p)
  c_v <- c(0, 1.7, 0.9, 0)                               # couple predictors 2 and 3

  # Reference precision/RHS, replicating the function's assembly in plain R
  K <- k * p; P_joint <- matrix(0, K, K); Pb_joint <- numeric(K)
  for (ip in 1:p) {
    j <- pp[ip]; om <- omega[, ip]
    A <- prior_P + crossprod(X, X * om)
    rows <- ((ip - 1) * k + 1):(ip * k)
    P_joint[rows, rows] <- A
    Pb_joint[rows] <- prior_Pb[, ip] + crossprod(X, kappa_w[, j] + om * c_j_mat[, ip])
  }
  for (v in which(c_v > 0)) { idx <- v + (0:(p - 1)) * k; P_joint[idx, idx] <- P_joint[idx, idx] + c_v[v] * Msym }
  analytic <- solve(P_joint, Pb_joint)

  N <- 4000; acc <- numeric(K)
  for (r in 1:N) acc <- acc + as.numeric(draw_beta_symhs_pooled(
    X, Xt, kappa_w, omega, c_j_mat, prior_P, prior_Pb, pp,
    matrix(0, n, p), FALSE, c_v, Msym, NULL, NULL))
  emp <- acc / N
  tol <- 5 * max(sqrt(diag(solve(P_joint)))) / sqrt(N)
  check("A4 draw_beta_symhs_pooled MC mean", maxerr(emp, analytic) < tol,
        sprintf("err=%.4f tol=%.4f", maxerr(emp, analytic), tol))
}

# =============================================================================
# PART B — Symmetric feature behaviour in the full sampler
# =============================================================================
cat("\nPART B — symmetric features end-to-end\n")

# Synthetic MNL with a genuinely NULL predictor (true effect 0 in every category)
make_data_sym <- function(n = 1600, p_all = 3, ng = 12, seed = 7, re_sd = 0.35) {
  set.seed(seed)
  X <- cbind(intercept = 1, x1 = rnorm(n), x2 = rnorm(n), xnull = rnorm(n))
  k <- ncol(X); group_idx <- sample(1:ng, n, replace = TRUE)
  beta_true <- matrix(0, k, p_all, dimnames = list(colnames(X), NULL))
  beta_true["x1", 1:(p_all - 1)] <- c(1.1, -0.7)[1:(p_all - 1)]
  beta_true["x2", 1:(p_all - 1)] <- c(-0.9, 0.8)[1:(p_all - 1)]
  beta_true["intercept", 1:(p_all - 1)] <- c(0.2, -0.3)[1:(p_all - 1)]
  # xnull: 0 everywhere
  re <- array(rnorm(k * p_all * ng, 0, re_sd), c(k, p_all, ng)); re[, p_all, ] <- 0
  Y <- matrix(0, n, p_all)
  for (i in 1:n) {
    eta <- sapply(1:p_all, function(j) sum(X[i, ] * (beta_true[, j] + re[, j, group_idx[i]])))
    pr <- exp(eta - max(eta)); pr <- pr / sum(pr)
    Y[i, ] <- rmultinom(1, sample(25:55, 1), pr)
  }
  colnames(Y) <- paste0("cat", 1:p_all)
  list(X = X, Y = Y, group_idx = group_idx, beta_true = beta_true)
}

fit_sym <- function(d, baseline, seed = 1, niter = 500, nburn = 250) {
  set.seed(seed)
  suppressWarnings(suppressMessages(
    mnlogit_rcpp_sym(X = d$X, Y = d$Y, baseline = baseline, group_idx = d$group_idx,
                     niter = niter, nburn = nburn, use_re = TRUE, use_ncp = TRUE,
                     use_horseshoe = TRUE, symmetric_hs = TRUE, horseshoe_idx = 2:ncol(d$X),
                     use_spike_slab = TRUE, use_car = FALSE, use_bart = FALSE,
                     standardize = TRUE, chain_id = 1L)
  ))
}

d <- make_data_sym()

## B1. Symmetric pooled RE variance: in a real fit, the returned per-predictor RE
##     sd must be (near-)constant across categories for each RE predictor — the
##     end-to-end manifestation of the A3a invariant.
{
  cat("  -- fit: symmetric_hs, baseline = p_all\n")
  fit <- fit_sym(d, baseline = ncol(d$Y))
  ran <- is.list(fit) && !is.null(fit$sigma_beta_pooled)
  check("B1 sym RE fit runs & finite", ran && all(is.finite(fit$post_log_lik)))
  if (ran) {
    sb <- fit$sigma_beta_pooled                          # k x p (internal, baseline-coded)
    re_rows <- 1:nrow(sb)
    spread <- apply(sb[re_rows, , drop = FALSE], 1, function(r) diff(range(r)))
    check("B1 RE sd constant across categories", max(spread) < 1e-6,
          sprintf("max within-row spread=%.1e", max(spread)))
  }
  assign("fit_base_phigh", if (ran) fit else NULL, inherits = TRUE)
}

## B2. BASELINE INVARIANCE — the defining motivation for symmetric coding.
##     Refitting with a different baseline must leave the zero-sum posterior
##     means ~unchanged (up to MC error), because postb_pooled is stored on the
##     zero-sum scale and the symmetric prior treats categories exchangeably.
{
  cat("  -- fit: symmetric_hs, baseline = 1\n")
  fit_b1 <- fit_sym(d, baseline = 1L, seed = 2)
  ok <- !is.null(fit_base_phigh) && is.list(fit_b1)
  check("B2 alt-baseline fit runs & finite", ok && all(is.finite(fit_b1$post_log_lik)))
  if (ok) {
    m_a <- apply(fit_base_phigh$postb_pooled, c(1, 2), mean)   # zero-sum, cat order = Y columns
    m_b <- apply(fit_b1$postb_pooled, c(1, 2), mean)
    err <- max(abs(m_a - m_b))
    check("B2 zero-sum means baseline-invariant", err < 0.30,
          sprintf("max|Δ| across baselines = %.3f (tol 0.30)", err))
  }
}

## B3. Symmetric HORSESHOE shrinkage — a truly null predictor's zero-sum contrasts
##     must shrink markedly relative to a genuinely active predictor.
{
  if (!is.null(fit_base_phigh)) {
    m <- apply(fit_base_phigh$postb_pooled, c(1, 2), mean)     # zero-sum k x p_all
    mag <- function(name) max(abs(m[name, ]))
    null_mag <- mag("xnull"); active_mag <- mag("x1")
    check("B3 symHS shrinks null predictor", null_mag < 0.5 * active_mag,
          sprintf("|xnull|=%.3f vs |x1|=%.3f", null_mag, active_mag))
  }
}


# =============================================================================
cat("\nPART C — 2026-08-13 fixes and new functionality\n")
# WHY PART C EXISTS. PART B never caught the NULL-index bug because `fit_sym` passes an EXPLICIT
# `horseshoe_idx` (so it is never NULL) and its design has no prefix-grouped columns (so no const-sum
# reference is ever dropped and the post-drop remap never runs). Both conditions are required to
# reproduce the failure, so the suite was blind to it. C1/C2 close that gap.

# design WITH a constant-sum block, so `cs_drop` is non-empty and the remap at the drop block fires
make_data_cs <- function(n = 1200, p_all = 3, ng = 10, seed = 11) {
  set.seed(seed)
  w <- matrix(rgamma(n * 3, 1), n, 3); w <- w / rowSums(w)     # sums to EXACTLY 1 -> const-sum block
  X <- cbind(intercept = 1, x1 = rnorm(n), soil_s1 = w[,1], soil_s2 = w[,2], soil_s3 = w[,3])
  k <- ncol(X); gi <- sample(1:ng, n, replace = TRUE)
  bt <- matrix(0, k, p_all); bt[2, 1:(p_all-1)] <- c(1.0, -0.6)[1:(p_all-1)]
  Y <- matrix(0, n, p_all)
  for (i in 1:n) { eta <- sapply(1:p_all, function(j) sum(X[i,] * bt[,j]))
    pr <- exp(eta - max(eta)); pr <- pr/sum(pr); Y[i,] <- rmultinom(1, 40, pr) }
  colnames(Y) <- paste0("cat", 1:p_all)
  list(X = X, Y = Y, group_idx = gi)
}
dcs <- make_data_cs()
fit_cs <- function(..., seed = 3) { set.seed(seed)
  suppressWarnings(mnlogit_rcpp_sym(X = dcs$X, Y = dcs$Y, baseline = 1, group_idx = dcs$group_idx,
    niter = 260, nburn = 120, use_re = TRUE, use_ncp = TRUE, use_horseshoe = TRUE,
    const_sum_blocks = "auto", re_idx = c(1, 2), chain_id = 1L, ...)) }

## C1. REGRESSION for the NULL -> integer(0) bug. With horseshoe_idx = NULL AND a const-sum drop,
##     the post-drop remap used to turn NULL into integer(0), emptying hs_idx (k_hs = 0). The
##     signature is a NEGATIVE global scale: "tau0_pooled=-0.0331 (p0=-0.5)".
{
  out <- capture.output(fit_cs(), type = "output")
  ln  <- grep("Horseshoe Calibration", out, value = TRUE)
  ok_found <- length(ln) > 0
  tau0 <- if (ok_found) as.numeric(sub(".*tau0_pooled=([-0-9.eE]+).*", "\\1", ln[1])) else NA
  p0   <- if (ok_found) as.numeric(sub(".*\\(p0=([-0-9.eE]+)\\).*", "\\1", ln[1])) else NA
  check("C1 horseshoe calibration emitted", ok_found)
  check("C1 tau0 POSITIVE (was -0.0331 when hs_idx empty)", isTRUE(tau0 > 0), sprintf("tau0=%.5f", tau0))
  check("C1 p0 >= 1 (was -0.5 when k_hs = 0)", isTRUE(p0 >= 1), sprintf("p0=%.2f", p0))
}

## C2. The calibrated horseshoe must actually REACH the draw: changing p0 must change the fit.
##     Before the fix these were bit-identical because hs_idx was empty.
{
  # NB pick a p0 that DIFFERS from the default. This design has 4 active covariates after the
  # const-sum drop, so the default is max(1, round(4*0.3)) = 1 -- comparing against p0_mu = 1
  # compares a value with itself and trivially "passes" as identical (my first version did).
  a <- fit_cs(p0_mu = NULL); b <- fit_cs(p0_mu = 3)
  Ba <- apply(a$postb_pooled, c(1,2), mean); Bb <- apply(b$postb_pooled, c(1,2), mean)
  check("C2 p0_mu changes the posterior (HS reaches the draw)", !identical(Ba, Bb),
        sprintf("max|diff|=%.4g", maxerr(Ba, Bb)))
}

## C3. gibbs_step_re_ncp(return_lik) — opt-in likelihood return must be EXACT and INERT when off.
{
  set.seed(2); k <- 4; p <- 2; p_all <- 3; ng <- 5; n <- 200
  Xk <- cbind(1, matrix(rnorm(n*(k-1)), n, k-1)); gi <- sample(ng, n, TRUE)
  idx <- lapply(1:ng, function(m) which(gi == m))
  args <- list(Xk, t(Xk), matrix(rnorm(n*p_all), n, p_all), matrix(runif(n*p,.5,1.5), n, p),
               matrix(0, n, p), diag(1, k), matrix(0, k, p), matrix(0, k, p),
               matrix(0, k, p), matrix(.5, k, p), array(0, c(k,p,ng)), as.integer(1:p),
               idx, lapply(idx, function(i) Xk[i,,drop=FALSE]),
               lapply(idx, function(i) t(Xk[i,,drop=FALSE])), ng, matrix(0, n, p_all), FALSE,
               as.integer(0:(k-1)), matrix(1, k, ng), matrix(1, p, ng), as.integer(c(1, rep(0,k-1))),
               matrix(1, k, ng), matrix(1, k, p), FALSE)
  set.seed(9); r0 <- do.call(gibbs_step_re_ncp, c(args, list(FALSE)))
  set.seed(9); r1 <- do.call(gibbs_step_re_ncp, c(args, list(TRUE)))
  check("C3 return_lik=TRUE leaves the draw unchanged", identical(r0$mu, r1$mu),
        sprintf("max|diff|=%.3g", maxerr(r0$mu, r1$mu)))
  check("C3 P_lik returned with correct shape", !is.null(r1$P_lik) && all(dim(r1$P_lik) == c(k,k,p)))
  check("C3 Pb_lik returned with correct shape", !is.null(r1$Pb_lik) && all(dim(r1$Pb_lik) == c(k,p)))
  check("C3 P_lik symmetric (X'Omega X)", isTRUE(approx(r1$P_lik[,,1], t(r1$P_lik[,,1]), 1e-10)))
  check("C3 P_lik absent when return_lik=FALSE", is.null(r0$P_lik))
}

## C4. Regularised (Finnish) slab in update_re_precision_hc — the NON-symmetric variant. The cap
##     lived only in the _sym variant, so `re_regularize`/`collapse_slab_c2` were inert on the path
##     the sampler actually uses (proved by bit-identical fits at slab_c2 = 100 vs 4).
{
  set.seed(4); kk <- 5; pp_ <- 3; G <- 20; TRUE_SD <- 3
  bc <- array(rnorm(kk*pp_*G, 0, TRUE_SD), c(kk,pp_,G)); mu0 <- matrix(0, kk, pp_)
  chain <- function(c2, reg, nit = 40) { pr <- matrix(1, kk, pp_); au <- matrix(.5, kk, pp_)
    for (i in seq_len(nit)) { r <- update_re_precision_hc(beta_c = bc, mu_pooled = mu0,
        re_idx = as.integer(0:(kk-1)), n_groups = G, re_mask = array(1, c(kk,pp_,G)),
        y_mask = matrix(1, pp_, G), is_intercept = as.integer(c(1, rep(0, kk-1))),
        prec_prev = pr, a_aux_prev = au, re_scale_A = 1.0, re_regularize = reg, slab_c2 = c2)
      pr <- r$prec; au <- r$a_aux }
    median(r$sigma) }
  s_loose <- chain(100, TRUE); s_tight <- chain(1, TRUE); s_off <- chain(1, FALSE)
  check("C4 loose cap recovers the true RE sd", abs(s_loose - TRUE_SD) < 1.2,
        sprintf("sigma=%.3f (true %.1f, cap 10)", s_loose, TRUE_SD))
  check("C4 tight cap clips sigma at sqrt(c2)", s_tight <= 1.02, sprintf("sigma=%.3f (cap 1)", s_tight))
  check("C4 re_regularize=FALSE ignores slab_c2", s_off > 1.5, sprintf("sigma=%.3f (uncapped)", s_off))
}

## C5. Parameter-count-aware slab: collapse_slab_c2 = "auto" must give c2 = R^2 / K.
{
  out <- capture.output(fit_cs(collapse_slab_c2 = "auto", re_slab_range = 3), type = "output")
  ln <- grep("Auto slab", out, value = TRUE)
  got <- if (length(ln)) as.numeric(sub(".*c2 = R\\^2/K = ([0-9.]+).*", "\\1", ln[1])) else NA
  check("C5 auto slab emitted", length(ln) > 0)
  check("C5 auto slab c2 == R^2/K", isTRUE(abs(got - 9/2) < 1e-6), sprintf("c2=%.4f expected %.4f", got, 9/2))
}

## C6. Hierarchical horseshoe on the RE scales (re_hs_global) — must be OFF by default and inert.
{
  a <- fit_cs(); b <- fit_cs(re_hs_global = FALSE)
  check("C6 re_hs_global=FALSE is the default (bit-identical)",
        identical(apply(a$postb_pooled, c(1,2), mean), apply(b$postb_pooled, c(1,2), mean)))
  cc <- fit_cs(re_hs_global = TRUE, re_hs_tau0 = 0.1)
  check("C6 re_hs_global=TRUE runs and returns a tau trace",
        !is.null(cc$post_re_tau) && all(is.finite(cc$post_re_tau)) && all(cc$post_re_tau > 0))
}

## C7. Symmetric coupling now covers the NON-RE covariates too (complement block). Its defining
##     property is BASELINE-INVARIANCE of the zero-sum posterior; the diagonal penalty is not
##     invariant because it acts in baseline-removed coordinates.
##
##     `hs_kernel_live = !sym` IS LOAD-BEARING (added 2026-08-21). The symmetric penalty travels
##     through c_v and is always applied, but the DIAGONAL penalty travels through hs_prec_kernel,
##     which is bound once before the Gibbs loop to a still-zero hs_prec_mat and never rebinds -- so
##     without this flag the "diagonal" arm carries NO FE horseshoe at all and C7 silently compares
##     symmetric-HS-plus-ridge against RIDGE ALONE. The 1.80x it used to report for diagonal was the
##     ridge's baseline dependence, not the diagonal horseshoe's. Both arms now apply their own
##     penalty, which is what the test claims to compare.
{
  zs <- function(bl, sym, seed) { set.seed(seed)
    f <- suppressWarnings(mnlogit_rcpp_sym(X = dcs$X, Y = dcs$Y, baseline = bl,
      group_idx = dcs$group_idx, niter = 300, nburn = 140, use_re = TRUE, use_ncp = TRUE,
      use_horseshoe = TRUE, symmetric_hs = sym, hs_kernel_live = !sym,
      const_sum_blocks = "auto",
      re_idx = c(1, 2), chain_id = 1L))
    apply(f$postb_pooled, c(1,2), mean) }
  ratio <- function(sym) { a1 <- zs(1, sym, 11); a2 <- zs(1, sym, 22); cr <- zs(3, sym, 11)
    mean(abs(a1 - cr)) / max(mean(abs(a1 - a2)), 1e-12) }
  r_diag <- ratio(FALSE); r_sym <- ratio(TRUE)
  # The old comparative assertion ("symmetric is MORE baseline-invariant than diagonal") is GONE.
  # It was never measurable here: with both horseshoes live it reverses (symmetric 1.06x vs diagonal
  # 0.67x), and both sit inside MC noise because on a 1200-row, 3-category synthetic the likelihood
  # swamps the prior. It is a PRIOR property, so C7a below asserts it EXACTLY instead; what remains
  # end-to-end is the claim that actually holds and is worth regression-testing.
  check("C7 symmetric cross-baseline shift is within MC noise", r_sym < 1.6,
        sprintf("symmetric=%.2fx (diagonal=%.2fx, both live)", r_sym, r_diag))
}

## C7a. BASELINE INVARIANCE OF THE PENALTY ITSELF — exact, no MCMC, no noise.
##      For a zero-sum coefficient vector z (sum z = 0), baseline coding with baseline j gives
##      b_k = z_k - z_j. Then  b'Msym b = ||b||^2 - (sum b)^2/p_all = sum z^2  for EVERY j, while the
##      diagonal penalty ||b||^2 = sum z^2 + p_all * z_j^2 moves with the arbitrary baseline choice.
##      That identity is the entire reason symmetric_hs exists, so test it directly.
{
  set.seed(4)
  p_all <- 6
  z <- rnorm(p_all); z <- z - mean(z)                    # a zero-sum coefficient vector
  Msym <- diag(1, p_all - 1) - matrix(1 / p_all, p_all - 1, p_all - 1)
  q_sym <- q_diag <- numeric(p_all)
  for (j in seq_len(p_all)) {
    b <- (z - z[j])[-j]                                  # baseline-removed coords for baseline j
    q_sym[j]  <- as.numeric(t(b) %*% Msym %*% b)
    q_diag[j] <- sum(b^2)
  }
  check("C7a symmetric penalty is EXACTLY baseline-invariant",
        max(abs(q_sym - q_sym[1])) < 1e-10,
        sprintf("spread=%.2e (value %.6f == sum z^2 %.6f)", max(abs(q_sym - q_sym[1])),
                q_sym[1], sum(z^2)))
  check("C7a symmetric quadratic form equals the zero-sum norm",
        abs(q_sym[1] - sum(z^2)) < 1e-10, sprintf("diff=%.2e", abs(q_sym[1] - sum(z^2))))
  check("C7a diagonal penalty is NOT baseline-invariant (the motivation)",
        diff(range(q_diag)) > 1e-6,
        sprintf("range %.4f-%.4f (%.2fx)", min(q_diag), max(q_diag), max(q_diag)/min(q_diag)))
}


## C8. FULL BAYES ON THE RE SLAB (estimate_slab_c2). c2 must be (a) traced, (b) IDENTIFIED -- two very
##     different starting values must converge to the same neighbourhood, which is exactly the check
##     the learned global tau (re_hs_global) FAILED -- and (c) inert when off.
{
  c2run <- function(est, start, seed = 5) { set.seed(seed)
    suppressWarnings(mnlogit_rcpp_sym(X = dcs$X, Y = dcs$Y, baseline = 1, group_idx = dcs$group_idx,
      niter = 400, nburn = 180, use_re = TRUE, use_ncp = TRUE, use_horseshoe = TRUE,
      const_sum_blocks = "auto", re_idx = c(1, 2), re_regularize = TRUE,
      collapse_slab_c2 = start, estimate_slab_c2 = est, slab_df_re = 10, slab_s2_re = 4,
      chain_id = 1L)) }
  off <- c2run(FALSE, 4)
  check("C8 estimate_slab_c2=FALSE leaves no c2 trace", is.null(off$post_slab_c2))
  a <- c2run(TRUE, 4); b <- c2run(TRUE, 100)
  ok_tr <- !is.null(a$post_slab_c2) && all(is.finite(a$post_slab_c2)) && all(a$post_slab_c2 > 0)
  check("C8 c2 traced, finite and positive", ok_tr)
  ma <- median(a$post_slab_c2); mb <- median(b$post_slab_c2)
  # identified => the ratio of the two posterior medians is near 1 despite a 25x gap in starting value
  check("C8 c2 IDENTIFIED (starts 4 and 100 converge)", max(ma, mb) / min(ma, mb) < 2.5,
        sprintf("median from 4 = %.2f vs from 100 = %.2f (ratio %.2f)", ma, mb, max(ma,mb)/min(ma,mb)))
  check("C8 c2 does not run off to the funnel regime", mb < 25,
        sprintf("median from start 100 = %.2f (start was 100)", mb))
}

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n================================================================\n")
cat(sprintf(" RESULT: %d/%d tests passed", .tests_run - .tests_failed, .tests_run))
if (.tests_failed > 0) {
  cat(sprintf("  (%d FAILED)\n", .tests_failed)); cat(" Failures:\n")
  for (f in .failures) cat("   -", f, "\n")
} else cat("  — ALL PASS\n")
cat("================================================================\n")
quit(status = if (.tests_failed > 0) 1L else 0L)
