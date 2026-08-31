#!/usr/bin/env Rscript
# =============================================================================
# test_altspec_nested.R — the block routed through nested_cut
# =============================================================================
# Consolidates what were three files (nested wiring / IV path / interactions). Covers:
#   1. NCUT_ALT_BLOCKS builds an n x K_node matrix aligned to the node's ALTERNATIVES, at leaf AND
#      nest, via one resolver; source columns leave the design; scale factor recorded for replay
#   2. rebuilding a block on new X reproduces the fitted block EXACTLY (the IV precondition)
#   3. the nested/IV path runs, and delta genuinely ENTERS the inclusive value
#   4. covariate-modulated inertia: delta_i = delta_0 + gamma * w_i
# =============================================================================
suppressMessages(suppressWarnings({source("codes/mnlogit_rcpp_sym.R"); source("codes/nested_cut.R")}))

# ---------------------------------------------------------------------------
# NATIVE ALT-SPEC BLOCKS through nested_cut: does NCUT_ALT_BLOCKS build an n x K_node matrix
# aligned to the node's ALTERNATIVES, remove the source columns from the design, and reach the
# sampler? Covers BOTH node types via one resolver (.ncut_alt_fine_of): at a leaf an alternative is
# a fine class, at a nest it is the set of fine classes under that child.
# ---------------------------------------------------------------------------
ok <- TRUE
chk <- function(nm, cond, det = "") { ok <<- ok && isTRUE(cond)
  cat(sprintf("  [%s] %-52s %s\n", if (isTRUE(cond)) "PASS" else "FAIL", nm, det)) }

cats <- c("Wheat", "Maize", "Forest", "Urban")
set.seed(5); n <- 1200
Fsh <- matrix(rgamma(n * 4, 1.2), n, 4); Fsh <- Fsh / rowSums(Fsh)
colnames(Fsh) <- paste0("focal_", cats)
X <- cbind(intercept = 1, x1 = rnorm(n), x2 = rnorm(n), Fsh)
cat("\n--- unit: block construction ---\n")
# pin scale explicitly: nested_cut scales blocks by default (so the factor can be REPLAYED
# at IV/predict time), which would change these raw-value assertions.
Sys.setenv(NCUT_ALT_BLOCKS = "spatial:focal_:shared:none")
sp <- .ncut_alt_specs()
chk("spec parsed", length(sp) == 1 && sp[[1]]$name == "spatial" && sp[[1]]$coef == "shared" && sp[[1]]$scale == "none")

# LEAF: alternatives are the fine classes
ab_leaf <- .ncut_alt_blocks(X, cats, .ncut_alt_fine_of(cats), sp)
chk("leaf block is n x K", nrow(ab_leaf$blocks$spatial$Z) == n && ncol(ab_leaf$blocks$spatial$Z) == 4L)
chk("leaf block columns are the alternatives", identical(colnames(ab_leaf$blocks$spatial$Z), cats))
chk("leaf consumed the 4 focal cols", length(ab_leaf$drop) == 4L)
chk("leaf block rows sum to 1", max(abs(rowSums(ab_leaf$blocks$spatial$Z) - 1)) < 1e-10)

# NEST: alternatives are children; a child's column must be the SUM over its fine classes
node <- list(Crop = c("Wheat", "Maize"), Nature = c("Forest"), Built = c("Urban"))
ab_nest <- .ncut_alt_blocks(X, names(node), .ncut_alt_fine_of(node), sp)
Zn <- ab_nest$blocks$spatial$Z
chk("nest block is n x n_children", nrow(Zn) == n && ncol(Zn) == 3L)
chk("nest Crop col == focal_Wheat + focal_Maize (renormalised)",
    max(abs(Zn[, "Crop"] - (Fsh[, "focal_Wheat"] + Fsh[, "focal_Maize"]))) < 1e-10)

# the DEFAULT (scale="sd") standardises and records the factor so it can be replayed
sp_sd <- { Sys.setenv(NCUT_ALT_BLOCKS = "spatial:focal_:shared"); .ncut_alt_specs() }
ab_sd <- .ncut_alt_blocks(X, cats, .ncut_alt_fine_of(cats), sp_sd)
chk("default scale='sd' standardises the block", abs(sd(ab_sd$blocks$spatial$Z) - 1) < 0.05,
    sprintf("sd=%.3f", sd(ab_sd$blocks$spatial$Z)))
chk("scale factor is recorded for replay", !is.null(ab_sd$blocks$spatial$def$z_scale))
chk("block def carries what a rebuild needs",
    all(c("prefix","coef","alt_names") %in% names(ab_sd$blocks$spatial$def)))
Zr <- .ncut_alt_rebuild(ab_sd$blocks$spatial$def, X, .ncut_alt_fine_of(cats))
chk("rebuild on the same X reproduces the fitted block", max(abs(Zr - ab_sd$blocks$spatial$Z)) < 1e-12,
    sprintf("max|diff|=%.2e", max(abs(Zr - ab_sd$blocks$spatial$Z))))
Sys.setenv(NCUT_ALT_BLOCKS = "spatial:focal_:shared:none")

cat("\n--- integration: through nested_cut_fit (flat tree) ---\n")
b <- matrix(0, ncol(X), length(cats))
b[2, ] <- c(.8, -.5, .3, 0); b[3, ] <- c(-.4, .6, .1, 0)
U <- X %*% b + 1.1 * Fsh                       # true spatial-lag delta = 1.1
P <- exp(U - apply(U, 1, max)); P <- P / rowSums(P)
Y <- t(apply(P, 1, function(pr) rmultinom(1, 30, pr))); colnames(Y) <- cats
tree <- setNames(as.list(cats), cats)
fit <- suppressWarnings(suppressMessages(nested_cut_fit(
  X = X, Y = Y, tree = tree, use_iv = FALSE, use_re = FALSE,
  niter = 400, nburn = 200, thin = 1L, M = 1L)))
chk("nested_cut_fit runs with the block on", is.list(fit) && !is.null(fit$root),
    if (is.list(fit)) paste("top:", paste(names(fit), collapse=",")) else "not a list")
Sys.unsetenv("NCUT_ALT_BLOCKS")
cat(sprintf("\n  [%s] native alt-spec blocks wire through nested_cut\n", if (ok) "PASS" else "FAIL"))
if (!ok) quit(status = 1L)

# ---------------------------------------------------------------------------
# ALT-SPEC BLOCKS ON THE NESTED / IV PATH -- the configuration that used to crash with
#   Error in Xn %*% Bpool : non-conformable arguments   (.ncut_node_iv -> .ncut_eta)
# because the inclusive value rebuilt a child's utilities from the FULL X while the child was
# fitted on its reduced design, and without delta. Both are now replayed from the stored
# x_cols/alt_cols + alt_def, with the FIT-TIME scale factor.
# ---------------------------------------------------------------------------
ok <- TRUE
chk <- function(nm, cond, det = "") { ok <<- ok && isTRUE(cond)
  cat(sprintf("  [%s] %-50s %s\n", if (isTRUE(cond)) "PASS" else "FAIL", nm, det)) }

set.seed(23); n <- 2500
fine <- c("W","M","F1","F2","U")
Fsh <- matrix(rgamma(n*length(fine), 1.2), n, length(fine)); Fsh <- Fsh/rowSums(Fsh)
colnames(Fsh) <- paste0("prev_", fine)
X <- cbind(intercept = 1, x1 = rnorm(n), x2 = rnorm(n), Fsh)
b <- matrix(0, ncol(X), length(fine)); b[2, ] <- c(.7,-.5,.3,.1,0); b[3, ] <- c(-.3,.4,.2,-.1,0)
U <- X %*% b + 1.3 * Fsh
P <- exp(U - apply(U,1,max)); P <- P/rowSums(P)
Y <- t(apply(P, 1, function(pr) rmultinom(1, 30, pr))); colnames(Y) <- fine
tree <- list(Crop = c("W","M"), Forest = c("F1","F2"), Urban = "U")   # a REAL nest: IV children

Sys.setenv(NCUT_ALT_BLOCKS = "temporal:prev_:shared")
f <- try(suppressWarnings(suppressMessages(nested_cut_fit(
  X = X, Y = Y, tree = tree, use_iv = TRUE, use_re = FALSE,
  niter = 400, nburn = 200, thin = 1L, M = 2L))), silent = TRUE)
Sys.unsetenv("NCUT_ALT_BLOCKS")
chk("nested fit WITH IV and blocks completes", !inherits(f, "try-error"),
    if (inherits(f, "try-error")) sub("\n.*", "", conditionMessage(attr(f, "condition"))) else "")
if (!inherits(f, "try-error")) {
  r <- f$root
  chk("root is a nest with IV children", r$type == "nest" && length(r$iv_children) > 0)
  chk("root kept its pre-xmap design width", !is.null(r$alt_cols))
  chk("root carries block definitions", !is.null(r$alt_def))
  kids <- r$children
  chk("children carry x_cols + alt_def",
      all(vapply(kids, function(k) !is.null(k$x_cols) || k$type %in% c("singleton","even"), TRUE)))
  # the real check: recompute a child's inclusive value on fresh X of the SAME shape
  iv <- try(.ncut_node_iv(kids[[r$iv_children[1]]], X, 1L), silent = TRUE)
  chk("child inclusive value recomputes on X", !inherits(iv, "try-error") &&
      length(iv) == nrow(X) && all(is.finite(iv)),
      if (inherits(iv, "try-error")) sub("\n.*", "", conditionMessage(attr(iv, "condition"))) else
        sprintf("range [%.2f, %.2f]", min(iv), max(iv)))
  # and that delta actually enters it: zeroing delta must CHANGE the IV
  k1 <- kids[[r$iv_children[1]]]
  if (!is.null(k1$delta_draws)) {
    # zero delta WITHOUT changing its structure: lapply over a bare matrix iterates ELEMENTS and
    # silently turns a [n_delta x draws] matrix into a list of scalars.
    k0 <- k1
    k0$delta_draws <- if (is.list(k1$delta_draws)) lapply(k1$delta_draws, function(d) d * 0)
                      else k1$delta_draws * 0
    iv0 <- .ncut_node_iv(k0, X, 1L)
    chk("delta CONTRIBUTES to the inclusive value", max(abs(iv - iv0)) > 1e-6,
        sprintf("max|IV - IV(delta=0)| = %.4f", max(abs(iv - iv0))))
  }
}
cat(sprintf("\n  [%s] alt-spec blocks work on the nested/IV path\n", if (ok) "PASS" else "FAIL"))
if (!ok) quit(status = 1L)

# ---------------------------------------------------------------------------
# COVARIATE-MODULATED INERTIA:  delta_i = delta_0 + gamma * w_i
# The interaction w_m * z_ij is itself alternative-specific, so it is just another alt-spec block.
#   gamma < 0 -> the covariate REDUCES persistence -> makes a transition more likely
#   gamma > 0 -> it locks the pixel in
# This is how a STATIC variable informs the EVOLUTION equation rather than only the initial state.
# ---------------------------------------------------------------------------
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
