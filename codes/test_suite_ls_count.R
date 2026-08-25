# =============================================================================
# test_suite_ls_count.R
#
# Feature-specific validation for the COUNT model additions in count_rcpp.R
# / count_gibbs_core.cpp.
#
# Layout:
#   PART A  Exact unit tests of the count kernels (machine precision / MC)
#   PART B  Count feature behaviour in the full sampler (end-to-end)
#
# Run:  Rscript codes/test_suite_ls_count.R
# Exit code is non-zero on any failure.
# =============================================================================

suppressMessages({ library(Rcpp); library(RcppArmadillo); library(MASS) })

`%||%` <- function(a, b) if (!is.null(a)) a else b
find_file <- function(fname) for (d in c("codes", ".")) {
  p <- file.path(d, fname); if (file.exists(p)) return(normalizePath(p))
} %||% stop("Cannot locate ", fname)

.tests_run <- 0L; .tests_failed <- 0L; .failures <- character(0)
check <- function(name, condition, detail = "") {
  .tests_run <<- .tests_run + 1L; ok <- isTRUE(condition)
  if (!ok) { .tests_failed <<- .tests_failed + 1L; .failures <<- c(.failures, paste(name, detail)) }
  cat(sprintf("  [%s] %-54s %s\n", if (ok) "PASS" else "FAIL", name, detail)); invisible(ok)
}
approx <- function(a, b, tol = 1e-10) max(abs(as.numeric(a) - as.numeric(b))) < tol
maxerr <- function(a, b) max(abs(as.numeric(a) - as.numeric(b)))

cat("================================================================\n")
cat(" test_suite_ls_count.R — count-feature validation\n")
cat(sprintf(" BLAS: %s\n", extSoftVersion()["BLAS"]))
cat("================================================================\n\n")

cat(">>> Sourcing sampler (compiles count C++ core)...\n")
suppressMessages(suppressWarnings(source(find_file("count_rcpp.R"))))
cat("    ready.\n\n")

# =============================================================================
# PART A — Exact unit tests of the count kernels
# =============================================================================
cat("PART A — count kernel unit tests\n")
set.seed(11)

## A1. .count_re_precision_pooled — pooled RE precision
##     A single variance is shared across all p count columns for each RE predictor.
{
  k <- 5; p <- 4; ng <- 6
  mu <- matrix(rnorm(k * p), k, p)
  beta_c <- array(rnorm(k * p * ng, 0, 0.5), c(k, p, ng)) + array(mu, c(k, p, ng))
  re_idx <- as.integer(1:k)
  re_mask <- array(1.0, c(k, p, ng)); y_mask <- matrix(1.0, p, ng)
  is_int <- as.integer(c(1, rep(0, k - 1)))
  prec_prev <- matrix(1.0, k, p)
  re_support <- matrix(1L, k, ng)

  out <- .count_re_precision_pooled(beta_c, mu, re_idx, ng, re_mask, y_mask, is_int,
                                    prec_prev, 1.0, FALSE, 4.0, re_support, p)
  row_spread <- max(apply(out$prec[2:k, , drop = FALSE], 1, function(r) diff(range(r))))
  check("A1 pooled RE var constant across columns", row_spread < 1e-12,
        sprintf("max within-row spread=%.1e", row_spread))
}

## A2. update_country_shrinkage_hc_cpp — Country-Gatekeeper with Regularized Slab Cap
##     Verifies the `slab_c2_country = 4` capping logic binds tau_m to <= sqrt(4) = 2.0.
{
  # Synthetic setup
  k <- 2; p <- 3; ng <- 10
  beta_c <- array(rnorm(k*p*ng, 0, 8.0), c(k, p, ng)) # Large RE to force tau expansion
  mu_pooled <- matrix(0, k, p)
  sigma_mat <- matrix(0.1, k, p) # small sigma -> very large standardized deviation (8.0/0.1 = 80)
  tau_country <- rep(1.0, ng)
  nu_country <- rep(1.0, ng)
  re_mask <- array(1.0, c(k, p, ng))
  y_mask <- matrix(1.0, p, ng)
  is_int <- as.integer(rep(0, k))
  re_support <- matrix(1L, k, ng)

  out_capped <- update_country_shrinkage_hc_cpp(
    beta_c, mu_pooled, sigma_mat, as.integer(0), 
    re_mask, y_mask, is_int, tau_country, nu_country, 
    1.0, 4.0, re_support, FALSE
  )
  out_uncapped <- update_country_shrinkage_hc_cpp(
    beta_c, mu_pooled, sigma_mat, as.integer(0), 
    re_mask, y_mask, is_int, tau_country, nu_country, 
    1.0, -1.0, re_support, FALSE # slab_c2 <= 0 means no cap
  )
  max_capped <- max(out_capped$tau_country)
  max_uncapped <- max(out_uncapped$tau_country)
  
  check("A2 Country shrinkage cap is enforced (<= 2.0)", max_capped <= 2.01, sprintf("max_capped=%.3f", max_capped))
  check("A2 Uncapped shrinkage is much larger", max_uncapped > 4.0, sprintf("max_uncapped=%.3f", max_uncapped))
}

# =============================================================================
# PART B — Count feature behaviour in the full sampler (end-to-end)
# =============================================================================
cat("\nPART B — count features end-to-end\n")

# Synthetic Negative Binomial dataset
make_data_count <- function(n = 600, p = 2, ng = 6, seed = 7) {
  set.seed(seed)
  X <- cbind(intercept = 1, x1 = rnorm(n))
  k <- ncol(X); group_idx <- sample(1:ng, n, replace = TRUE)
  beta_true <- matrix(0, k, p, dimnames = list(colnames(X), NULL))
  beta_true["intercept", ] <- c(0.5, -0.2)
  beta_true["x1", ] <- c(1.2, -0.8)
  
  # Country Shrinkage simulation: 
  # ng=1,2 have near zero RE (tau_m small). ng=5,6 have large RE (tau_m large)
  re_sd <- c(0.01, 0.01, 0.5, 0.5, 1.5, 1.5) 
  re <- array(0, c(k, p, ng), dimnames=list(colnames(X), NULL, NULL))
  for(g in 1:ng) {
    re["intercept", , g] <- rnorm(p, 0, 0.2)
    re["x1", , g] <- rnorm(p, 0, re_sd[g])
  }
  
  offset_mat <- matrix(log(10), n, p) # 10x exposure
  r_true <- c(2.5, 5.0) # overdispersion params
  
  Y <- matrix(0, n, p)
  for (i in 1:n) {
    eta <- sapply(1:p, function(j) sum(X[i, ] * beta_true[, j]) + re["intercept", j, group_idx[i]] + re["x1", j, group_idx[i]] * X[i, "x1"] + offset_mat[i, j])
    # rnM(n, mu, theta) where variance = mu + mu^2/theta.
    # We use theta = r. E[Y] = exp(eta).
    mu <- exp(eta)
    for(j in 1:p) Y[i, j] <- rnbinom(1, size = r_true[j], mu = mu[j])
  }
  colnames(Y) <- paste0("count", 1:p)
  list(X = X, Y = Y, group_idx = group_idx, offset_mat = offset_mat, beta_true = beta_true, r_true = r_true)
}

fit_count <- function(d, seed = 1, niter = 300, nburn = 100, use_shrinkage=FALSE, offset=NULL) {
  set.seed(seed)
  suppressWarnings(suppressMessages(
    mncount_rcpp(X = d$X, Y = d$Y, family = "negbin", group_idx = d$group_idx, offset = offset,
               niter = niter, nburn = nburn, use_re = TRUE, use_horseshoe = FALSE, 
               re_idx = c(1, 2), chain_id = 1L, re_prec_pooled = TRUE,
               use_country_shrinkage = use_shrinkage, slab_c2_country = 4.0)
  ))
}

d <- make_data_count()

## B1. Offset Invariance
##     Adding a +log(10) offset shifts the fitted intercept down by exactly log(10)
{
  cat("  -- fit: offset=0 vs offset=log(10)\n")
  fit_base <- fit_count(d, offset = matrix(0, nrow(d$Y), ncol(d$Y)))
  fit_off <- fit_count(d, offset = d$offset_mat)
  
  m_base <- apply(fit_base$postb_pooled, c(1, 2), mean)
  m_off <- apply(fit_off$postb_pooled, c(1, 2), mean)
  
  shift_intercept <- m_base["intercept",] - m_off["intercept",]
  shift_slope <- m_base["x1",] - m_off["x1",]
  
  check("B1 intercept exactly shifted by offset", maxerr(shift_intercept, rep(log(10), ncol(d$Y))) < 0.35,
        sprintf("diff=%.3f", maxerr(shift_intercept, rep(log(10), ncol(d$Y)))))
  check("B1 slope is invariant to offset", maxerr(shift_slope, rep(0, ncol(d$Y))) < 0.30,
        sprintf("diff=%.3f", maxerr(shift_slope, rep(0, ncol(d$Y)))))
}

## B2. Dispersion Identification (r parameter)
##     The NB dispersion `r` should be positive and finite.
{
  cat("  -- fit: NB dispersion\n")
  r_est <- apply(fit_off$post_r, 2, mean)
  check("B2 r (dispersion) traced and positive", all(r_est > 0) && all(is.finite(r_est)),
        sprintf("Est r: [%.2f, %.2f]", r_est[1], r_est[2]))
}

## B3. Country Shrinkage acts on random effects
##     Verifies the shrinkage parameters are traced and strictly bounded by the slab cap.
{
  cat("  -- fit: use_country_shrinkage=TRUE\n")
  fit_cs <- fit_count(d, offset = d$offset_mat, use_shrinkage = TRUE)
  
  if (!is.null(fit_cs$post_tau_country)) {
    tau_est <- apply(fit_cs$post_tau_country, 1, mean)
    
    check("B3 tau_country traced and positive", all(tau_est > 0),
          sprintf("min_tau=%.3f", min(tau_est)))
    check("B3 tau_country does not exceed cap", max(tau_est) <= 2.1,
          sprintf("max_tau=%.3f", max(tau_est)))
  } else {
    check("B3 Country Shrinkage tau trace missing", FALSE, "tau_country trace missing")
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
