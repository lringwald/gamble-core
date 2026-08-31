#!/usr/bin/env Rscript
# =============================================================================
# test_altspec_sampler.R — the alternative-specific (conditional-logit) block, at SAMPLER level
# =============================================================================
# Consolidates what were three files (recovery / multiblock / scaling). Covers, in order:
#   1. a single shared delta recovers, and is RECORDED on both the in-RAM and STREAMED paths
#   2. two blocks (spatial + temporal) recover SIMULTANEOUSLY and stay separately identified
#   3. per_class recovers (delta_b pinned at 0)
#   4. symmetric recovers, is zero-sum, and is BASELINE-INVARIANT
#   5. scale="sd" makes deltas comparable across regressors of different spread
#
# EVERY assertion gates on posterior sd vs PRIOR sd, not on coverage. A parameter disconnected from
# the likelihood still COVERS the truth, because its posterior IS the wide prior -- that is exactly
# how this block sat broken behind a validation flag. Coverage alone passes a broken sampler.
# =============================================================================
suppressMessages(suppressWarnings(source(Sys.getenv("MNL_SRC", "codes/mnlogit_rcpp_sym.R"))))

# ---------------------------------------------------------------------------
# alt_spec_Z RECOVERY TEST. Simulate V_ij = x_i'beta_j + delta * z_ij with a KNOWN delta and check
# the sampler returns it. Before 2026-08-26 delta was disconnected from the likelihood (the additive
# utility channel was gated on use_bart), so it random-walked from its prior -- the failure this
# test exists to catch. Reports the posterior against truth AND against the prior, because a
# disconnected delta looks like the PRIOR, not like noise.
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# MULTI-BLOCK alternative-specific recovery: a SPATIAL Y-lag block and a TEMPORAL (t-1 LU state)
# block, each with its own delta, plus a per_class block. This is the structure the transition
# feature needs: separate effects, separate coefficients, drawn jointly (they are correlated).
# ---------------------------------------------------------------------------
set.seed(4242)
n <- 5000; p_all <- 4; k <- 3
prior_sd <- 10
D_SPAT <- 1.20     # spatial Y-lag coefficient
D_TEMP <- -0.80    # temporal (t-1 state) coefficient -- opposite sign on purpose
X <- cbind(intercept = 1, x1 = rnorm(n), x2 = rnorm(n))
beta <- matrix(0, k, p_all); beta[2, ] <- c(.7,-.5,.3,0); beta[3, ] <- c(-.4,.6,.1,0)
mk <- function() { M <- matrix(rgamma(n * p_all, 1.3), n, p_all); M / rowSums(M) }
Zs <- mk(); Zt_ <- mk()                       # correlated-in-shape share matrices
U <- X %*% beta + D_SPAT * Zs + D_TEMP * Zt_
P <- exp(U - apply(U, 1, max)); P <- P / rowSums(P)
Y <- t(apply(P, 1, function(pr) rmultinom(1, 40, pr)))
fit <- suppressWarnings(suppressMessages(mnlogit_rcpp_sym(
  X = X, Y = Y, intercept = FALSE, baseline = p_all, niter = 1200, nburn = 600,
  use_re = FALSE, use_horseshoe = FALSE, standardize = FALSE, calc_loo = FALSE,
  alt_spec_Z = list(spatial  = list(Z = Zs,  coef = "shared", scale = "none"),
                    temporal = list(Z = Zt_, coef = "shared", scale = "none")),
  alt_spec_prior_sd = prior_sd, alt_spec_allow_unvalidated = TRUE, chain_id = 1L)))
pd <- fit$post_delta
cat("\n========== MULTI-BLOCK alt_spec RECOVERY ==========\n")
truth <- c(spatial = D_SPAT, temporal = D_TEMP)
ok <- TRUE
for (nm in rownames(pd)) {
  d <- pd[nm, ]; lo <- quantile(d, .025); hi <- quantile(d, .975)
  cov <- lo <= truth[[nm]] && truth[[nm]] <= hi
  inf <- sd(d) < 0.25 * prior_sd
  ok <- ok && cov && inf
  cat(sprintf("  %-9s true %+.2f | post %+.3f [%+.3f, %+.3f] sd %.4f | covers %s\n",
              nm, truth[[nm]], mean(d), lo, hi, sd(d), if (cov) "YES" else "NO"))
}
cat(sprintf("  correlation(spatial, temporal) draws: %+.3f  (drawn JOINTLY, so the ridge is handled)\n",
            cor(pd["spatial", ], pd["temporal", ])))
cat(sprintf("\n  [%s] both blocks recovered and separately identified\n", if (ok) "PASS" else "FAIL"))

# ---- per_class block (the transition "degrees": 1 delta = scalar inertia, p deltas = per-class) ----
# The baseline handling DIFFERS by coef type and getting it wrong silently attenuates the effect:
#   shared    -> regressor is (Z_ij - Z_ib), because delta*(Z_ij) - delta*(Z_ib) factorises.
#   per_class -> delta_j Z_ij - delta_b Z_ib does NOT factorise; delta_b is pinned at 0 and the
#                regressor is the RAW Z_ij. Using the differenced column here read
#                true +1.50/-1.00/+0.50 as +0.78/-0.31/+0.35.
set.seed(77)
dk <- c(1.5, -1.0, 0.5, 0)
Zp <- mk()
Up <- X %*% beta + sweep(Zp, 2, dk, "*")
Pp <- exp(Up - apply(Up, 1, max)); Pp <- Pp / rowSums(Pp)
Yp <- t(apply(Pp, 1, function(pr) rmultinom(1, 40, pr)))
fp <- suppressWarnings(suppressMessages(mnlogit_rcpp_sym(
  X = X, Y = Yp, intercept = FALSE, baseline = p_all, niter = 1200, nburn = 600,
  use_re = FALSE, use_horseshoe = FALSE, standardize = FALSE, calc_loo = FALSE,
  alt_spec_Z = list(lag = list(Z = Zp, coef = "per_class", scale = "none")),
  alt_spec_prior_sd = prior_sd, alt_spec_allow_unvalidated = TRUE, chain_id = 1L)))
pdp <- fp$post_delta; ok2 <- TRUE
cat("\n---------- per_class block ----------\n")
for (i in seq_len(nrow(pdp))) {
  d <- pdp[i, ]; lo <- quantile(d, .025); hi <- quantile(d, .975)
  cv <- lo <= dk[i] && dk[i] <= hi; ok2 <- ok2 && cv
  cat(sprintf("  %-8s true %+.2f | post %+.3f [%+.3f, %+.3f] covers %s\n",
              rownames(pdp)[i], dk[i], mean(d), lo, hi, if (cv) "YES" else "NO"))
}
cat(sprintf("\n  [%s] per_class block recovered\n", if (ok2) "PASS" else "FAIL"))

# ---- symmetric block: zero-sum constraint instead of a pinned baseline ----------------------------
# per_class pins delta_b = 0, so every coefficient is relative to an ARBITRARY reference class --
# and `baseline` is which.max(colSums(Y)), i.e. data-dependent and different per node. `symmetric`
# imposes sum_j delta_j = 0 instead, making the coefficients baseline-INVARIANT and each one read as
# that class's effect relative to the AVERAGE. Design per equation is
#   eta_ij - eta_ib = sum_k delta_k [1{k=j} Z_ij + Z_ib]
set.seed(31)
dsym <- c(1.2, -0.9, 0.4, -0.7); dsym <- dsym - mean(dsym)      # zero-sum truth
Zy <- mk()
Uy <- X %*% beta + sweep(Zy, 2, dsym, "*")
Py <- exp(Uy - apply(Uy, 1, max)); Py <- Py / rowSums(Py)
Yy <- t(apply(Py, 1, function(pr) rmultinom(1, 40, pr)))
run_sym <- function(bl) { set.seed(9)
  fy <- suppressWarnings(suppressMessages(mnlogit_rcpp_sym(
    X = X, Y = Yy, intercept = FALSE, baseline = bl, niter = 900, nburn = 450,
    use_re = FALSE, use_horseshoe = FALSE, standardize = FALSE, calc_loo = FALSE,
    alt_spec_Z = list(lag = list(Z = Zy, coef = "symmetric", scale = "none")),
    alt_spec_prior_sd = prior_sd, alt_spec_allow_unvalidated = TRUE, chain_id = 1L)))
  d <- fy$post_delta; full <- rbind(d, -colSums(d))
  ppv <- setdiff(seq_len(p_all), bl); out <- matrix(0, p_all, ncol(d)); out[c(ppv, bl), ] <- full; out
}
sa <- run_sym(p_all); sb <- run_sym(2L)
ok3 <- abs(sum(rowMeans(sa))) < 1e-8 && max(abs(rowMeans(sa) - dsym)) < 0.15 &&
       max(abs(rowMeans(sa) - rowMeans(sb))) < 0.15
cat("\n---------- symmetric block ----------\n")
for (j in seq_len(p_all))
  cat(sprintf("  cls%-3d true %+.2f | bl=last %+.3f | bl=2 %+.3f\n", j, dsym[j], mean(sa[j,]), mean(sb[j,])))
cat(sprintf("  zero-sum: %+.1e | baseline shift: %.4f\n", sum(rowMeans(sa)),
            max(abs(rowMeans(sa) - rowMeans(sb)))))
cat(sprintf("\n  [%s] symmetric block: recovers truth, zero-sum, baseline-invariant\n",
            if (ok3) "PASS" else "FAIL"))

# ---- scale="sd" (the DEFAULT): delta is per-sd, so classes with different regressor spread become
# comparable. Four columns simulated with the SAME standardised effect but 5x different spread give
# raw deltas 6.8/6.7/31.7/31.0 -- a pure units artifact -- and per-sd deltas all near 0.800.
set.seed(88); nn <- 6000; pa <- 5
Xs <- cbind(intercept = 1, x1 = rnorm(nn), x2 = rnorm(nn))
bs <- matrix(0, 3, pa); bs[2, ] <- c(.6,-.4,.3,.2,0); bs[3, ] <- c(-.3,.5,.1,-.2,0)
rw <- cbind(rgamma(nn,6), rgamma(nn,6), rgamma(nn,0.15), rgamma(nn,0.15), rgamma(nn,3))
Zc <- rw / rowSums(rw); sdc <- apply(Zc, 2, sd)
dt <- 0.8 / sdc; dt[pa] <- 0
Uc <- Xs %*% bs + sweep(Zc, 2, dt, "*")
Pc <- exp(Uc - apply(Uc,1,max)); Pc <- Pc / rowSums(Pc)
Yc <- t(apply(Pc, 1, function(pr) rmultinom(1, 40, pr)))
set.seed(4)
fc <- suppressWarnings(suppressMessages(mnlogit_rcpp_sym(
  X = Xs, Y = Yc, intercept = FALSE, baseline = pa, niter = 800, nburn = 400,
  use_re = FALSE, use_horseshoe = FALSE, standardize = FALSE, calc_loo = FALSE,
  alt_spec_Z = list(lag = list(Z = Zc, coef = "per_class", scale = "sd")),
  alt_spec_prior_sd = prior_sd, alt_spec_allow_unvalidated = TRUE, chain_id = 1L)))
dm <- rowMeans(fc$post_delta)
ok4 <- max(abs(dm - 0.8)) < 0.12
cat("\n---------- scale=\"sd\" comparability ----------\n")
cat(sprintf("  true standardised effect 0.800 for all; fitted: %s  (spread %.3f)\n",
            paste(sprintf("%.3f", dm), collapse = " "), diff(range(dm))))
cat(sprintf("\n  [%s] scale=\"sd\" makes deltas comparable across regressor scales\n",
            if (ok4) "PASS" else "FAIL"))
if (!ok || !ok2 || !ok3 || !ok4) quit(status = 1L)
