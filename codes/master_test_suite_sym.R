# =============================================================================
# master_test_suite_sym.R
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
# Run:  Rscript codes/master_test_suite_sym.R
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
cat(" master_test_suite_sym.R — symmetric-feature validation\n")
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
