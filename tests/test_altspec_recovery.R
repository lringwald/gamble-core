# ---------------------------------------------------------------------------
# alt_spec_Z RECOVERY TEST. Simulate V_ij = x_i'beta_j + delta * z_ij with a KNOWN delta and check
# the sampler returns it. Before 2026-08-26 delta was disconnected from the likelihood (the additive
# utility channel was gated on use_bart), so it random-walked from its prior -- the failure this
# test exists to catch. Reports the posterior against truth AND against the prior, because a
# disconnected delta looks like the PRIOR, not like noise.
# ---------------------------------------------------------------------------
suppressMessages(suppressWarnings(source(Sys.getenv("MNL_SRC", "codes/mnlogit_rcpp_sym.R"))))
DELTA_TRUE <- as.numeric(Sys.getenv("DELTA_TRUE", "1.5"))
alt_prior_sd <- 10
set.seed(2026)
n <- 4000; p_all <- 4; k <- 3
X <- cbind(intercept = 1, x1 = rnorm(n), x2 = rnorm(n))
beta <- matrix(0, k, p_all)
beta[2, ] <- c( 0.9, -0.6, 0.4, 0)
beta[3, ] <- c(-0.5,  0.8, 0.2, 0)
beta[1, ] <- c( 0.2, -0.3, 0.1, 0)
# alternative-specific attribute: a share vector per pixel (exactly the shape a t-1 LU state has)
Zr <- matrix(rgamma(n * p_all, 1.2), n, p_all); Z <- Zr / rowSums(Zr)
U <- X %*% beta + DELTA_TRUE * Z
P <- exp(U - apply(U, 1, max)); P <- P / rowSums(P)
Y <- t(apply(P, 1, function(pr) rmultinom(1, 40, pr)))
colnames(Y) <- paste0("c", 1:p_all)
fit <- suppressWarnings(suppressMessages(mnlogit_rcpp_sym(
  X = X, Y = Y, intercept = FALSE, baseline = p_all, niter = 1200, nburn = 600,
  use_re = FALSE, use_horseshoe = FALSE, standardize = FALSE,
  alt_spec_Z = list(z = list(Z = Z, coef = "shared", scale = "none")), alt_spec_prior_sd = alt_prior_sd, alt_spec_allow_unvalidated = TRUE,
  calc_loo = FALSE, chain_id = 1L)))
d <- as.numeric(fit$post_delta)
cat(sprintf("\n=============== alt_spec delta RECOVERY ===============\n"))
cat(sprintf("  true delta        : %.3f\n", DELTA_TRUE))
cat(sprintf("  posterior mean    : %.3f\n", mean(d)))
cat(sprintf("  posterior 95%% CI  : [%.3f, %.3f]\n", quantile(d, .025), quantile(d, .975)))
cat(sprintf("  posterior sd      : %.4f   (prior sd %.0f -> a DISCONNECTED delta shows ~the prior)\n", sd(d), alt_prior_sd))
cat(sprintf("  covers truth      : %s\n", if (quantile(d,.025) <= DELTA_TRUE && DELTA_TRUE <= quantile(d,.975)) "YES" else "NO"))
bm <- apply(fit$postb_pooled, c(1,2), mean)
cat(sprintf("\n  beta[x1] recovered: %s   (truth %s)\n",
    paste(sprintf("%+.2f", bm[2, ]), collapse=" "), paste(sprintf("%+.2f", beta[2, ] - mean(beta[2, ])), collapse=" ")))
# PASS/FAIL. The discriminating statistic is the posterior sd against the PRIOR sd: a delta that is
# disconnected from the likelihood still COVERS the truth (its posterior is the prior, which is wide
# and centred at 0), so a coverage check alone passes the broken sampler. Measured pre-fix:
# posterior sd 9.894 vs prior 10, CI [-23.1, 16.9]. Post-fix: sd 0.088, CI [1.35, 1.70].
# STREAMED-PATH REGRESSION. post_delta[, s] used to live inside `if (!save_posterior_to_disk)`, so
# every STREAM=TRUE run (the production default) returned the zero matrix it was initialised with:
# delta drawn correctly each sweep, then not recorded. Assert it is recorded on BOTH paths.
fit_s <- suppressWarnings(suppressMessages(mnlogit_rcpp_sym(
  X = X, Y = Y, intercept = FALSE, baseline = p_all, niter = 400, nburn = 200,
  use_re = FALSE, use_horseshoe = FALSE, standardize = FALSE,
  alt_spec_Z = list(z = list(Z = Z, coef = "shared", scale = "none")), alt_spec_prior_sd = alt_prior_sd, alt_spec_allow_unvalidated = TRUE,
  save_posterior_to_disk = TRUE, disk_path = tempfile("altspec_stream"),
  calc_loo = FALSE, chain_id = 1L)))
ds <- as.numeric(fit_s$post_delta)
cat(sprintf("  streamed path     : mean %+.3f  sd %.4f  %s\n", mean(ds), sd(ds),
            if (sd(ds) > 1e-8) "recorded" else "ALL ZERO -- not recorded"))
streamed_ok <- sd(ds) > 1e-8

informative <- sd(d) < 0.25 * alt_prior_sd
covers <- quantile(d, .025) <= DELTA_TRUE && DELTA_TRUE <= quantile(d, .975)
ok <- informative && covers && streamed_ok
cat(sprintf("\n  [%s] alt_spec delta: data-informed, covers truth, recorded on BOTH paths\n", if (ok) "PASS" else "FAIL"))
if (!ok) quit(status = 1L)
