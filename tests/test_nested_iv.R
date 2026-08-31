# =============================================================================
# test_nested_iv.R — synthetic validation of the IV-nested framework
# DGP with a KNOWN nest structure and a KNOWN lambda (inclusive-value coupling).
# Checks: (1) lambda is recovered, (2) IV coef is positive & in ~(0,1],
#         (3) the IV model beats the no-IV factorized model on HELD-OUT joint log-lik.
# =============================================================================
suppressMessages({library(Rcpp); library(RcppArmadillo)})
source("codes/mnlogit_rcpp_sym.R")
source("codes/mnlogit_nested_iv.R")
set.seed(42)

# ---- DGP -------------------------------------------------------------------
n <- 6000; ntr <- 40                      # pixels, multinomial trials/pixel
X <- cbind(intercept = 1, x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
k <- ncol(X)
LAMBDA_TRUE <- 0.6

# sub-model B: 3 crop-like classes, utilities depend on X
gamma <- rbind(c(0.0, 0.3, -0.3),        # intercepts B1,B2,B3
               c(0.8, -0.4, 0.2),        # x1
               c(-0.3, 0.6, 0.1),        # x2
               c(0.2, 0.1, -0.5))        # x3
U_B <- X %*% gamma; colnames(U_B) <- c("B1","B2","B3")   # n x 3
IV_B <- .inclusive_value(U_B)             # true inclusive value (logsum)

# macro: nest A (singleton) vs nest B (the 3-class crop nest), with lambda*IV_B
alpha_A <- c(0.6, 0.3, 0.0, 0.2)
alpha_B <- c(-1.0, 0.0, 0.4, -0.2)
U_macro <- cbind(A = as.vector(X %*% alpha_A), B = as.vector(X %*% alpha_B) + LAMBDA_TRUE * IV_B)
Pm <- .softmax_rows(U_macro); colnames(Pm) <- c("A","B")   # n x 2 nest probs
Ps <- .softmax_rows(U_B)                  # n x 3 within-B probs

# joint fine probs: A, B1, B2, B3
Pfine <- cbind(A = Pm[, "A"], Pm[, "B"] * Ps)
colnames(Pfine) <- c("A", "B1", "B2", "B3")
cat(sprintf("realized nest shares: A=%.2f B=%.2f | fine: %s\n",
    mean(Pm[,"A"]), mean(Pm[,"B"]), paste(sprintf("%s=%.2f", colnames(Pfine), colMeans(Pfine)), collapse=" ")))

# draw multinomial counts -> shares
Y <- t(vapply(seq_len(n), function(i) rmultinom(1, ntr, Pfine[i, ]), numeric(4)))
colnames(Y) <- colnames(Pfine); Y <- Y / rowSums(Y)

nest_structure <- list(A = "A", B = c("B1", "B2", "B3"))

# ---- train / test split ----------------------------------------------------
tr <- sample(n, round(0.7 * n)); te <- setdiff(seq_len(n), tr)
jll <- function(P, Yh) sum(Yh * log(pmax(P, 1e-12)))   # joint multinomial log-lik (shares)

cat("\n--- fitting IV nested (use_iv=TRUE) ---\n")
fit_iv <- nested_iv_fit(X[tr,], Y[tr,], nest_structure, use_iv = TRUE, niter = 1500, nburn = 500)
cat("--- fitting factorized nested (use_iv=FALSE) ---\n")
fit_no <- nested_iv_fit(X[tr,], Y[tr,], nest_structure, use_iv = FALSE, niter = 1500, nburn = 500)

P_iv <- predict_nested_iv(fit_iv, X[te,])
P_no <- predict_nested_iv(fit_no, X[te,])

lam_B <- fit_iv$root$lambda[["B"]]   # lambda is on the ROOT NODE; fit_iv$lambda is NULL and
                                     # silently degrades sprintf() to character(0) (prints nothing)
stopifnot(is.numeric(lam_B), length(lam_B) == 1L, is.finite(lam_B))
cat("\n================= RESULTS =================\n")
cat(sprintf("TRUE lambda_B = %.2f  |  ESTIMATED lambda_B = %.3f\n", LAMBDA_TRUE, lam_B))
cat(sprintf("held-out joint log-lik:  IV = %.1f   no-IV = %.1f   (IV better by %.1f)\n",
    jll(P_iv, Y[te,]), jll(P_no, Y[te,]), jll(P_iv, Y[te,]) - jll(P_no, Y[te,])))
# null (global shares) held-out ll for scale
pbar <- colMeans(Y[tr,]); LLn <- sum(sweep(Y[te,], 2, log(pbar), `*`))
cat(sprintf("null (global-share) held-out log-lik = %.1f\n", LLn))
cat(sprintf("held-out per-pixel calibration (obs vs pred fine shares):\n"))
print(round(rbind(obs = colMeans(Y[te,]), pred_IV = colMeans(P_iv), pred_noIV = colMeans(P_no)), 3))
ok_lambda <- abs(lam_B - LAMBDA_TRUE) < 0.2
ok_better <- jll(P_iv, Y[te,]) > jll(P_no, Y[te,])
stopifnot(length(ok_lambda) == 1L, length(ok_better) == 1L)   # guard against vacuous assertions
cat(sprintf("\nPASS lambda-recovery: %s | PASS IV-beats-noIV: %s\n", ok_lambda, ok_better))
# use_iv=FALSE must PREDICT (it appended IV columns the fit never saw before the iv_children fix)
ok_nofit <- length(fit_no$root$iv_children) == 0L && all(is.finite(P_no)) &&
            max(abs(rowSums(P_no) - 1)) < 1e-8
cat(sprintf("PASS no-IV predicts (iv_children empty, rows sum to 1): %s\n", ok_nofit))
if (!all(ok_lambda, ok_better, ok_nofit)) { cat("TEST DONE\n"); quit(status = 1) }
cat("TEST DONE\n")
