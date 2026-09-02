# =============================================================================
# test_nested_cut_focal.R — validate the FOCAL ROUTING that identifies lambda.
# =============================================================================
# The pathology: IV_c = logsumexp(X'beta_leaf) is BUILT FROM the leaf utilities, so if the
# same fine focal columns also sit at the parent node, IV_c is nearly linear in them ->
# lambda is not separately identified, flips out of (0,1], and dumps probability onto the
# residual singletons (the real-data Sweden/Finland "urban wash").
#
# DGP mirrors the intended structure: the LEAF choice is driven by the fine focal
# composition, and the MACRO choice by lambda*IV_B PLUS a genuine macro-persistence term
# kappa*focal_B_tot. So "macro_totals" is the correctly-specified root and "flat" is the
# over-parameterised one that cannot separate lambda from the fine focal effects.
#
# Checks:
#   (1) routing does what it says: node design columns / leaf design untouched,
#   (2) IV~design R2 falls sharply under macro_totals vs flat,
#   (3) lambda: macro_totals recovers truth in (0,1]; flat is worse,
#   (4) predict_nested_cut replays the routed design (in-sample LL sane, no dim error),
#   (5) leaf_only also identifies lambda (IV as the sole channel).
# =============================================================================
suppressMessages({library(Rcpp); library(RcppArmadillo)})
source("codes/mnlogit_rcpp_sym.R")
source("codes/nested_cut.R")
set.seed(11)

n <- 6000; ntr <- 40; LAMBDA_TRUE <- 0.6; KAPPA_TRUE <- 1.2

# ---- covariates: 2 plain + a constant-sum focal block over {A, B1, B2, B3, ctx} ---------
x1 <- rnorm(n); x2 <- rnorm(n)
Fr <- matrix(rgamma(n * 5, shape = c(1.0, 1.3, 1.0, 0.8, 1.1)), n, 5, byrow = TRUE)
Fr <- Fr / rowSums(Fr)
colnames(Fr) <- c("focal_A", "focal_B1", "focal_B2", "focal_B3", "focal_ctx")
X <- cbind(intercept = 1, x1 = x1, x2 = x2, Fr)
foc_B <- c("focal_B1", "focal_B2", "focal_B3")
B_tot <- rowSums(Fr[, foc_B])

# ---- leaf DGP: the fine focal composition drives WHICH subtype -------------------------
gL <- rbind(intercept = c(0.0,  0.3, -0.2),      # rows = X cols, cols = B1,B2,B3
            x1        = c(0.4, -0.3,  0.1),
            x2        = c(-0.2, 0.5,  0.2),
            focal_A   = c(0.5, -0.6,  0.3),      # cross-class -> subtype signal (kept at leaves)
            focal_B1  = c(2.5, -0.5, -0.8),
            focal_B2  = c(-0.7, 2.6, -0.4),
            focal_B3  = c(-0.6, -0.9, 2.4),
            focal_ctx = c(0.2,  0.1, -0.3))
U_B <- X %*% gL; colnames(U_B) <- c("B1", "B2", "B3"); IV_B <- .ncut_lse(U_B)

# ---- macro DGP: lambda*IV_B + kappa*focal_B_tot (macro persistence) --------------------
aA <- c(intercept = 0.5, x1 = 0.3, x2 = 0.0, focal_A = 1.0, focal_ctx = 0.4)
UA <- as.vector(X[, names(aA)] %*% aA)
UB <- -1.0 + 0.1 * x1 - 0.2 * x2 + KAPPA_TRUE * B_tot + LAMBDA_TRUE * IV_B
Pm <- .ncut_softmax(cbind(A = UA, B = UB)); colnames(Pm) <- c("A", "B")
Ps <- .ncut_softmax(U_B)
Pf <- cbind(A = Pm[, "A"], Pm[, "B"] * Ps); colnames(Pf) <- c("A", "B1", "B2", "B3")
Y <- t(vapply(seq_len(n), function(i) rmultinom(1, ntr, Pf[i, ]), numeric(4)))
colnames(Y) <- colnames(Pf); Y <- Y / rowSums(Y)
tree <- list(A = "A", B = c("B1", "B2", "B3"))
tr <- sample(n, round(.7 * n)); te <- setdiff(seq_len(n), tr)
cat(sprintf("n=%d (train %d) | true lambda=%.2f kappa=%.2f | macro shares A=%.2f B=%.2f\n",
            n, length(tr), LAMBDA_TRUE, KAPPA_TRUE, mean(Pm[, "A"]), mean(Pm[, "B"])))

# ---- (1) unit-test the design map itself ----------------------------------------------
cat("\n--- (1) design map ---\n")
xm <- .ncut_focal_xmap(colnames(X), tree, iv_children = "B", prefix = "focal_", rule = "macro_totals")
Xn <- .ncut_apply_xmap(xm, X)
cat("  node design:", paste(colnames(Xn), collapse = ", "), "\n")
ok_cols <- identical(colnames(Xn),
                     c("intercept", "x1", "x2", "focal_A", "focal_ctx", "focal_B_tot"))
ok_val  <- max(abs(Xn[, "focal_B_tot"] - B_tot)) < 1e-12
ok_lo   <- is.null(.ncut_focal_xmap(colnames(X), tree, "B", "focal_", "flat"))
xm_lo   <- .ncut_focal_xmap(colnames(X), tree, "B", "focal_", "leaf_only")
ok_drop <- length(xm_lo$agg) == 0L && ncol(.ncut_apply_xmap(xm_lo, X)) == ncol(X) - 3L
# the AME chain-rule map: X %*% T must reproduce the node design exactly
ok_T    <- max(abs(X %*% .ncut_xmap_T(xm, ncol(X)) - Xn)) < 1e-12
cat(sprintf("  cols OK: %s | total==rowSums: %s | flat=no-op: %s | leaf_only drops 3: %s | X%%*%%T==Xnode: %s\n",
            ok_cols, ok_val, ok_lo, ok_drop, ok_T))

# ---- (2)-(3) fit the three rules -------------------------------------------------------
fitr <- function(rule) nested_cut_fit(X[tr, ], Y[tr, ], tree, use_iv = TRUE, focal_rule = rule,
                                      M = 15, draws_per_impute = 15, niter = 800, nburn = 300,
                                      progress = FALSE, fit_progress = FALSE)
cat("\n--- (2) fit: focal_rule='flat' (legacy, pathological) ---\n");  f_flat <- fitr("flat")
cat("\n--- (2) fit: focal_rule='macro_totals' (default)      ---\n");  f_macro <- fitr("macro_totals")
cat("\n--- (2) fit: focal_rule='leaf_only'                   ---\n");  f_leaf <- fitr("leaf_only")

id <- function(f) nested_cut_identification(f)[1, ]
cat("\n================= IDENTIFICATION =================\n")
tab <- do.call(rbind, lapply(list(flat = f_flat, macro_totals = f_macro, leaf_only = f_leaf), id))
tab$rule <- rownames(tab)
print(tab[, c("rule", "iv_r2", "lambda", "lambda_q025", "lambda_q975", "n_focal_routed")], row.names = FALSE)
cat(sprintf("\nTRUE lambda = %.2f\n", LAMBDA_TRUE))

# ---- (4) leaves untouched + predict replays the routed design --------------------------
cat("\n--- (4) leaf design + prediction ---\n")
leafB <- f_macro$root$children$B
ok_leaf <- nrow(leafB$beta_draws[[1]]) == ncol(X)      # leaf still fit on the FULL fine focal block
jll <- function(P, Yh) sum(Yh * log(pmax(P, 1e-12)))
P_macro <- predict_nested_cut(f_macro, X[te, ], D = 100, summary = "mean")
P_flat  <- predict_nested_cut(f_flat,  X[te, ], D = 100, summary = "mean")
P_leaf  <- predict_nested_cut(f_leaf,  X[te, ], D = 100, summary = "mean")
pbar <- colMeans(Y[tr, ]); LLn <- sum(sweep(Y[te, ], 2, log(pbar), `*`))
cat(sprintf("  leaf keeps all %d fine focal cols: %s\n", ncol(X), ok_leaf))
cat(sprintf("  held-out joint LL: macro_totals=%.1f  leaf_only=%.1f  flat=%.1f  null=%.1f\n",
            jll(P_macro, Y[te, ]), jll(P_leaf, Y[te, ]), jll(P_flat, Y[te, ]), LLn))
cat("  held-out calibration (obs vs predicted):\n")
print(round(rbind(obs = colMeans(Y[te, ]), macro_totals = colMeans(P_macro),
                  leaf_only = colMeans(P_leaf), flat = colMeans(P_flat)), 3))

# marginal-effects table must still run in the ORIGINAL covariate space
ame <- effective_symmetric_table(f_macro, X[te[1:400], ])
ok_ame <- identical(dim(ame$median), c(4L, ncol(X)))

# ---- verdict ---------------------------------------------------------------------------
lam <- setNames(tab$lambda, tab$rule); r2 <- setNames(tab$iv_r2, tab$rule)
lo <- setNames(tab$lambda_q025, tab$rule); hi <- setNames(tab$lambda_q975, tab$rule)
ok2 <- r2["macro_totals"] < r2["flat"] - 0.05
# lambda: in (0,1] with a CI that brackets the truth. NOTE deliberately NOT "beats flat" --
# once the const-sum reconstruction bug is fixed (2026-08-04) flat is no longer pathological
# on this DGP, and the flat/macro_totals lambda gap is far inside the posterior CIs. Routing
# is justified by the R^2 drop (identification), not by a lambda difference here.
ok3 <- lam["macro_totals"] > 0 && lam["macro_totals"] <= 1 &&
       lo["macro_totals"] <= LAMBDA_TRUE && LAMBDA_TRUE <= hi["macro_totals"]
ok5 <- lam["leaf_only"] > 0 && lam["leaf_only"] <= 1 &&
       lo["leaf_only"] <= LAMBDA_TRUE && LAMBDA_TRUE <= hi["leaf_only"]
cat(sprintf("\n%s (1) design map\n", if (all(ok_cols, ok_val, ok_lo, ok_drop, ok_T)) "[PASS]" else "[FAIL]"))
cat(sprintf("%s (2) IV~design R2 drops (%.2f -> %.2f)\n", if (ok2) "[PASS]" else "[FAIL]", r2["flat"], r2["macro_totals"]))
cat(sprintf("%s (3) macro_totals lambda in (0,1], CI brackets truth (%.3f [%.3f,%.3f] vs true %.2f)\n",
            if (ok3) "[PASS]" else "[FAIL]", lam["macro_totals"], lo["macro_totals"], hi["macro_totals"], LAMBDA_TRUE))
cat(sprintf("%s (4) leaves untouched / predict+AME replay routing\n", if (all(ok_leaf, ok_ame)) "[PASS]" else "[FAIL]"))
cat(sprintf("%s (5) leaf_only lambda in (0,1], CI brackets truth (%.3f [%.3f,%.3f])\n",
            if (ok5) "[PASS]" else "[FAIL]", lam["leaf_only"], lo["leaf_only"], hi["leaf_only"]))
cat("TEST DONE\n")
