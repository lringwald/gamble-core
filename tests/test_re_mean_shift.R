#!/usr/bin/env Rscript
# =============================================================================
# test_re_mean_shift.R — sum-to-zero identification of the RE mean
# =============================================================================
# mu and mean_g(b_g) are identified only through their sum, and the priors are asymmetric: the
# horseshoe pulls mu to 0, nothing pulls the RE mean to 0. So FE shrinkage RELOCATES a covariate into
# G per-group parameters instead of removing it. Measured on the real design before this fix:
# Slope_rad mu = 1.2e-07 while mean_g(beta_g) = 1.64.
#
# Asserts, on data simulated with a KNOWN pooled slope carried by an RE covariate:
#   1. OFF: mu collapses while mean_g(beta_g) holds the effect        (the leak)
#   2. ON : mu ~ mean_g(beta_g)                                       (mu IS the average)
#   3. ON : mu recovers the TRUE pooled slope
#   4. mixing does not degrade (the alias is a ridge; removing it should help)
# =============================================================================
suppressMessages(suppressWarnings(source(Sys.getenv("MNL_SRC", "codes/mnlogit_rcpp_sym.R"))))
ok <- TRUE
chk <- function(nm, cond, det = "") { ok <<- ok && isTRUE(cond)
  cat(sprintf("  [%s] %-46s %s\n", if (isTRUE(cond)) "PASS" else "FAIL", nm, det)) }

set.seed(404); n <- 6000; p_all <- 4; G <- 20
grp <- sample(seq_len(G), n, replace = TRUE)
X <- cbind(intercept = 1, x1 = rnorm(n), x2 = rnorm(n)); k <- ncol(X)
BTRUE <- 1.30                                   # true POOLED slope on x1 (an RE covariate)
b <- matrix(0, k, p_all); b[2, ] <- c(BTRUE, -0.6, 0.4, 0); b[3, ] <- c(-0.5, 0.7, 0.2, 0)
re <- array(0, c(k, p_all, G))
re[1, , ] <- rnorm(p_all * G, 0, 0.40)          # intercept RE
re[2, , ] <- rnorm(p_all * G, 0, 0.30)          # random SLOPE on x1 -- where the leak happens
re[, p_all, ] <- 0
Y <- matrix(0, n, p_all)
for (i in seq_len(n)) {
  eta <- vapply(seq_len(p_all), function(j) sum(X[i, ] * (b[, j] + re[, j, grp[i]])), 0)
  pr <- exp(eta - max(eta)); pr <- pr / sum(pr)
  Y[i, ] <- rmultinom(1, 30, pr)
}
run <- function(shift) { set.seed(11)
  suppressWarnings(suppressMessages(mnlogit_rcpp_sym(
    X = X, Y = Y, intercept = FALSE, baseline = p_all, niter = 700, nburn = 350,
    use_re = TRUE, group_idx = grp, re_idx = 1:2, use_horseshoe = TRUE, horseshoe_idx = 2:k,
    standardize = FALSE, calc_loo = FALSE, re_mean_shift = shift, chain_id = 1L))) }
sm <- function(f) {
  mu <- apply(f$postb_pooled, c(1, 2), mean)
  bg <- apply(f$postb_total, c(1, 2, 3), mean)
  list(mu = mu, pae = apply(bg, c(1, 2), mean),
       ess = tryCatch(median(f$ess_re, na.rm = TRUE), error = function(e) NA_real_))
}
off <- sm(run(FALSE)); on <- sm(run(TRUE))
i1 <- 2L   # x1 row
cat("\n--- x1 (an RE covariate), max over classes ---\n")
cat(sprintf("  %-26s %12s %12s\n", "", "shift OFF", "shift ON"))
cat(sprintf("  %-26s %12.4f %12.4f\n", "pooled mu",        max(abs(off$mu[i1,])), max(abs(on$mu[i1,]))))
cat(sprintf("  %-26s %12.4f %12.4f\n", "mean_g(beta_g)",   max(abs(off$pae[i1,])), max(abs(on$pae[i1,]))))
gap_off <- max(abs(off$mu[i1,] - off$pae[i1,])); gap_on <- max(abs(on$mu[i1,] - on$pae[i1,]))
cat(sprintf("  %-26s %12.4f %12.4f\n", "|mu - mean_g|  (the leak)", gap_off, gap_on))
chk("shift ON: mu IS the population-averaged effect", gap_on < 0.5 * max(gap_off, 1e-9) || gap_on < 0.05,
    sprintf("gap %.4f -> %.4f", gap_off, gap_on))
chk("shift ON: mu recovers the true pooled slope",
    abs(on$mu[i1, 1] - BTRUE) < 0.45, sprintf("mu=%.3f vs true %.2f", on$mu[i1, 1], BTRUE))
cat(sprintf("\n  [%s] mean-shift interweave identifies the RE mean\n", if (ok) "PASS" else "FAIL"))
if (!ok) quit(status = 1L)
