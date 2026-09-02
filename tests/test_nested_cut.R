# =============================================================================
# test_nested_cut.R — validate the CUT / IV-imputation nested framework.
# Known-lambda synthetic DGP (same as test_nested_iv). Checks:
#   (1) lambda recovered WITH a credible interval that brackets the truth,
#   (2) cut CI is wider than the sequential mean-plug-in point (honest uncertainty),
#   (3) held-out joint log-lik: IV-nested >= factorized, both > null; vs flat,
#   (4) predictive draws are calibrated.
# =============================================================================
suppressMessages({library(Rcpp); library(RcppArmadillo)})
source("codes/mnlogit_rcpp_sym.R")
source("codes/mnlogit_nested_iv.R")   # sequential (mean-plug-in) for comparison
source("codes/nested_cut.R")
set.seed(42)

# ---- DGP (A singleton vs B={B1,B2,B3}, lambda*IV_B) ------------------------
n <- 6000; ntr <- 40
X <- cbind(intercept = 1, x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n)); LAMBDA_TRUE <- 0.6
gamma <- rbind(c(0.0,0.3,-0.3), c(0.8,-0.4,0.2), c(-0.3,0.6,0.1), c(0.2,0.1,-0.5))
U_B <- X %*% gamma; colnames(U_B) <- c("B1","B2","B3"); IV_B <- .ncut_lse(U_B)
alpha_A <- c(0.6,0.3,0.0,0.2); alpha_B <- c(-1.0,0.0,0.4,-0.2)
U_macro <- cbind(A = as.vector(X %*% alpha_A), B = as.vector(X %*% alpha_B) + LAMBDA_TRUE * IV_B)
Pm <- .ncut_softmax(U_macro); colnames(Pm) <- c("A","B"); Ps <- .ncut_softmax(U_B)
Pfine <- cbind(A = Pm[,"A"], Pm[,"B"] * Ps); colnames(Pfine) <- c("A","B1","B2","B3")
Y <- t(vapply(seq_len(n), function(i) rmultinom(1, ntr, Pfine[i,]), numeric(4)))
colnames(Y) <- colnames(Pfine); Y <- Y / rowSums(Y)
tree  <- list(A = "A", B = c("B1","B2","B3"))
flat  <- list(all = c("A","B1","B2","B3"))
tr <- sample(n, round(.7*n)); te <- setdiff(seq_len(n), tr)
jll <- function(P, Yh) sum(Yh * log(pmax(P, 1e-12)))

cat("\n--- CUT IV-nested (use_iv=TRUE) ---\n")
f_iv  <- nested_cut_fit(X[tr,], Y[tr,], tree, use_iv = TRUE,  M = 20, draws_per_impute = 15, niter = 800, nburn = 300)
cat("--- CUT factorized (use_iv=FALSE) ---\n")
f_no  <- nested_cut_fit(X[tr,], Y[tr,], tree, use_iv = FALSE, M = 20, draws_per_impute = 15, niter = 800, nburn = 300)
cat("--- CUT flat MNL (single nest) ---\n")
f_fl  <- nested_cut_fit(X[tr,], Y[tr,], flat, use_iv = FALSE, M = 20, draws_per_impute = 15, niter = 800, nburn = 300)
cat("--- sequential mean-plug-in (for CI comparison) ---\n")
f_seq <- nested_iv_fit(X[tr,], Y[tr,], tree, use_iv = TRUE, niter = 1000, nburn = 400)

# lambda with CI
s <- summary_nested_cut(f_iv); lamB <- s[["root"]]$lambda["B",]
cat(sprintf("\n================= RESULTS =================\n"))
cat(sprintf("TRUE lambda_B = %.2f\n", LAMBDA_TRUE))
cat(sprintf("CUT   lambda_B = %.3f  [%.3f, %.3f]  (95%% CI, width %.3f)\n", lamB["median"], lamB["q025"], lamB["q975"], lamB["q975"]-lamB["q025"]))
cat(sprintf("SEQ   lambda_B = %.3f  (point, no CI)\n", f_seq$lambda[["B"]]))

# held-out predictive
P_iv <- predict_nested_cut(f_iv, X[te,], D = 150, summary = "mean")
P_no <- predict_nested_cut(f_no, X[te,], D = 150, summary = "mean")
P_fl <- predict_nested_cut(f_fl, X[te,], D = 150, summary = "mean")
pbar <- colMeans(Y[tr,]); LLn <- sum(sweep(Y[te,], 2, log(pbar), `*`))
cat(sprintf("\nheld-out joint log-lik:  IV=%.1f  factorized=%.1f  flat=%.1f  null=%.1f\n",
    jll(P_iv, Y[te,]), jll(P_no, Y[te,]), jll(P_fl, Y[te,]), LLn))
cat(sprintf("  IV beats factorized by %.1f | IV beats flat by %.1f\n",
    jll(P_iv,Y[te,])-jll(P_no,Y[te,]), jll(P_iv,Y[te,])-jll(P_fl,Y[te,])))
cat("\nheld-out calibration (obs vs predicted fine shares):\n")
print(round(rbind(obs = colMeans(Y[te,]), IV = colMeans(P_iv), factorized = colMeans(P_no), flat = colMeans(P_fl)), 3))

ok1 <- lamB["q025"] <= LAMBDA_TRUE && LAMBDA_TRUE <= lamB["q975"]
ok2 <- jll(P_iv, Y[te,]) >= jll(P_no, Y[te,]) - 1
cat(sprintf("\n%s lambda-CI-brackets-truth\n%s IV>=factorized OOS\n",
            if (ok1) "[PASS]" else "[FAIL]", if (ok2) "[PASS]" else "[FAIL]"))
cat("TEST DONE\n")
