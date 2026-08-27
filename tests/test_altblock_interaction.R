# ---------------------------------------------------------------------------
# COVARIATE-MODULATED INERTIA:  delta_i = delta_0 + gamma * w_i
# The interaction w_m * z_ij is itself alternative-specific, so it is just another alt-spec block.
#   gamma < 0 -> the covariate REDUCES persistence -> makes a transition more likely
#   gamma > 0 -> it locks the pixel in
# This is how a STATIC variable informs the EVOLUTION equation rather than only the initial state.
# ---------------------------------------------------------------------------
suppressMessages(suppressWarnings({source("codes/mnlogit_rcpp_sym.R"); source("codes/nested_cut.R")}))
set.seed(19); n <- 7000; p_all <- 4; cats <- paste0("c", 1:p_all)
w <- scale(rnorm(n))[, 1]
Zr <- matrix(rgamma(n * p_all, 1.3), n, p_all); Z <- Zr / rowSums(Zr)
colnames(Z) <- paste0("prev_", cats)
X <- cbind(intercept = 1, x1 = rnorm(n), w = w, Z)
b <- matrix(0, ncol(X), p_all); b[2, ] <- c(.6,-.4,.3,0); b[3, ] <- c(.2,-.1,.15,0)
D0 <- 1.2; GAM <- -0.9
U <- X %*% b + sweep(Z, 1, (D0 + GAM * w), "*")
P <- exp(U - apply(U, 1, max)); P <- P / rowSums(P)
Y <- t(apply(P, 1, function(pr) rmultinom(1, 40, pr))); colnames(Y) <- cats
Sys.setenv(NCUT_ALT_BLOCKS = "inertia:prev_:shared:none,modul:prev_:shared:none:w")
f <- suppressWarnings(suppressMessages(nested_cut_fit(
  X = X, Y = Y, tree = setNames(as.list(cats), cats),
  use_iv = FALSE, use_re = FALSE, niter = 900, nburn = 450, thin = 1L, M = 1L)))
Sys.unsetenv("NCUT_ALT_BLOCKS")
d <- f$root$delta_draws[[1]]; ok <- TRUE
cat("\n=== delta_i = delta_0 + gamma * w_i ===\n")
for (nm in rownames(d)) {
  tv <- if (nm == "inertia") D0 else GAM
  lo <- quantile(d[nm, ], .025); hi <- quantile(d[nm, ], .975)
  cv <- lo <= tv && tv <= hi; ok <- ok && cv
  cat(sprintf("  %-8s true %+.2f | post %+.3f [%+.3f, %+.3f]  covers %s\n",
              nm, tv, mean(d[nm, ]), lo, hi, if (cv) "YES" else "NO"))
}
cat(sprintf("\n  [%s] covariate-modulated inertia recovered\n", if (ok) "PASS" else "FAIL"))
if (!ok) quit(status = 1L)
