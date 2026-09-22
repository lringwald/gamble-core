# Mundlak device: recover gamma, and show that OMITTING it biases beta.
# Truth:  eta_ic[j] = mu_j + beta_j x_ic + gamma_j xbar_c + u_cj   (all zero-sum over j)
# x_ic = xbar_c + within noise, so omitting the xbar term forces the country RE to absorb it
# and part of it leaks into beta -- the classic RE bias Mundlak exists to fix.
suppressMessages({library(Rcpp); library(RcppArmadillo)})
sys.source("codes/mnlogit_rcpp_sym.R", environment()); source("codes/mundlak.R")
source("tests/helper_recovery.R")
zs <- function(v) v - mean(v)
set.seed(21)
G <- 40; npc <- 150; n <- G*npc; K <- 3
gidx  <- rep(seq_len(G), each = npc)
xbar  <- rnorm(G, 0, 1)                        # between-country level
x     <- xbar[gidx] + rnorm(n, 0, 1)           # 50/50 within/between
mu_t  <- zs(c( 0.30, -0.10, -0.20))
be_t  <- zs(c( 1.00, -0.60, -0.40))            # WITHIN effect (the estimand)
ga_t  <- zs(c(-1.20,  0.90,  0.30))            # group-mean loading, opposite sign to beta
u     <- t(apply(matrix(rnorm(G*K, 0, 0.35), G, K), 1, zs))
eta   <- outer(rep(1,n), mu_t) + outer(x, be_t) + outer(xbar[gidx], ga_t) + u[gidx, ]
P     <- exp(eta - apply(eta,1,max)); P <- P/rowSums(P)
Y     <- t(apply(P, 1, function(p) rmultinom(1,1,p))); colnames(Y) <- paste0("c",1:K)
X     <- cbind(intercept = 1, x = x)

md <- build_mundlak_design(X, gidx, mean_cols = "x", re_cols = "intercept", verbose = TRUE)
fit <- function(Xf, ri) { set.seed(5); capture.output(f <- mnlogit_rcpp_sym(
  X = Xf, Y = Y, intercept = FALSE, baseline = which.max(colSums(Y)), symmetric = TRUE,
  niter = 1500, nburn = 600, thin = 2L, use_re = TRUE, group_idx = gidx, re_idx = ri,
  use_horseshoe = FALSE)); f }
f_no <- fit(X,      which(colnames(X) == "intercept"))
f_md <- fit(md$X,   which(colnames(md$X) == "intercept"))
cf <- function(f, nm) apply(f$postb_pooled[nm, , , drop=FALSE], 2, mean)

b_no <- cf(f_no, "x"); b_md <- cf(f_md, "x"); g_md <- cf(f_md, md$names[1])
cat("\n=== beta on x (WITHIN effect) ===\n")
cat(sprintf("%-10s %8s %8s %8s\n", "class", "truth", "no-MDL", "with-MDL"))
for (j in 1:K) cat(sprintf("%-10s %8.3f %8.3f %8.3f\n", colnames(Y)[j], be_t[j], b_no[j], b_md[j]))
cat(sprintf("\nmean |error| : no-MDL %.3f   with-MDL %.3f\n",
    mean(abs(b_no-be_t)), mean(abs(b_md-be_t))))
cat("\n=== gamma (group-mean loading) ===\n")
for (j in 1:K) cat(sprintf("%-10s truth %7.3f   est %7.3f\n", colnames(Y)[j], ga_t[j], g_md[j]))
tv_record_vec("Mundlak gamma", ga_t, g_md, names = colnames(Y),
              note = "group-mean loading; the device is an interpretation fix")
cat(sprintf("mean |error| : %.3f\n", mean(abs(g_md-ga_t))))

ok_g <- mean(abs(g_md - ga_t)) < 0.35
ok_b <- mean(abs(b_md - be_t)) < mean(abs(b_no - be_t))
ok_s <- all(sign(g_md) == sign(ga_t))
cat(sprintf("\n[%s] gamma recovered (mean|err| %.3f < 0.35)\n", if(ok_g)"PASS" else "FAIL", mean(abs(g_md-ga_t))))
cat(sprintf("[%s] Mundlak REDUCES beta bias (%.3f -> %.3f)\n", if(ok_b)"PASS" else "FAIL",
    mean(abs(b_no-be_t)), mean(abs(b_md-be_t))))
cat(sprintf("[%s] gamma signs correct\n", if(ok_s)"PASS" else "FAIL"))
if (!all(ok_g, ok_b, ok_s)) quit(status = 1)
cat("TEST DONE\n")
