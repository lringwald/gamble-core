# predict_gamble(): one entry point, correct on every variant it dispatches to.
suppressMessages({library(Rcpp); library(RcppArmadillo)})
source("codes/mnl_aux_func.R"); source("codes/mnlogit_rcpp_sym.R")
source("codes/prior_model_predict.R"); source("codes/nested_cut.R")
source("codes/mundlak.R"); source("codes/predict_gamble.R")
np <- 0L; nf <- 0L
ok <- function(c_, m) { if (isTRUE(c_)) { np <<- np+1L; cat(sprintf("[PASS] %s\n", m)) }
                        else { nf <<- nf+1L; cat(sprintf("[FAIL] %s\n", m)) } }
sm <- function(U) { m <- apply(U,1,max); E <- exp(U-m); E/rowSums(E) }

set.seed(11); n <- 1500; G <- 6; K <- 3; k <- 4
X <- cbind(intercept=1, matrix(rnorm(n*(k-1)),n,k-1,dimnames=list(NULL,paste0("x",1:(k-1)))))
g <- sample(G, n, TRUE)
eta <- X %*% matrix(rnorm(k*K,0,.6),k,K); P0 <- sm(eta)
Y <- t(apply(P0,1,function(p) rmultinom(1,1,p))); colnames(Y) <- paste0("c",1:K)

# ---- flat MNL, no RE ----
f1 <- mnlogit_rcpp_sym(X=X, Y=Y, intercept=FALSE, baseline=K, symmetric=TRUE,
                       niter=400, nburn=150, thin=2L, calc_loo=FALSE)
P1 <- predict_gamble(f1, X)
ok(is.matrix(P1) && all(dim(P1) == c(n, K)), "flat MNL: returns n x J")
ok(max(abs(rowSums(P1) - 1)) < 1e-10, "flat MNL: rows sum to 1")
ok(all(P1 >= 0 & P1 <= 1), "flat MNL: probabilities in [0,1]")
# predicted composition should track the observed one
ok(max(abs(colMeans(P1) - colMeans(Y))) < 0.05, "flat MNL: mean predicted shares track observed")

# ---- flat MNL with RE ----
f2 <- mnlogit_rcpp_sym(X=X, Y=Y, intercept=FALSE, baseline=K, symmetric=TRUE, niter=400,
                       nburn=150, thin=2L, use_re=TRUE, group_idx=g, re_idx=1L, calc_loo=FALSE)
P2 <- predict_gamble(f2, X, group_idx = g, group_levels = unique(g))
ok(max(abs(rowSums(P2) - 1)) < 1e-10, "flat MNL + RE: rows sum to 1")

# ---- alternative-specific delta: predict_shares IGNORES it, predict_gamble must not ----
Z <- matrix(rnorm(n*K), n, K, dimnames=list(NULL, colnames(Y)))
etaD <- X %*% matrix(rnorm(k*K,0,.5),k,K) + 1.2*Z
YD <- t(apply(sm(etaD),1,function(p) rmultinom(1,1,p))); colnames(YD) <- colnames(Y)
f3 <- mnlogit_rcpp_sym(X=X, Y=YD, intercept=FALSE, baseline=K, symmetric=TRUE, niter=400, nburn=150,
                       thin=2L, alt_spec_Z=list(zz=list(Z=Z, coef="shared")),
                       alt_spec_allow_unvalidated=TRUE, calc_loo=FALSE)
P_no <- predict_shares(f3, X, bart_cols=NULL, linear_cols=seq_len(k))       # delta DROPPED
P_yes <- predict_gamble(f3, X, alt_spec_Z = list(list(Z = Z, coef = "shared")))
ok(max(abs(rowSums(P_yes) - 1)) < 1e-10, "alt-spec: rows sum to 1")
ok(mean(abs(P_yes - P_no)) > 1e-3, "alt-spec: delta CHANGES the prediction (predict_shares drops it)")
ll <- function(P) sum(YD * log(pmax(P, 1e-12)))
ok(ll(P_yes) > ll(P_no), sprintf("alt-spec: including delta fits better (%.1f vs %.1f)", ll(P_yes), ll(P_no)))

# ---- Mundlak columns rebuilt from the FIT's stored scaling ----
md <- build_mundlak_design(X, g, mean_cols="x1", re_cols="intercept", verbose=FALSE)
f4 <- mnlogit_rcpp_sym(X=md$X, Y=Y, intercept=FALSE, baseline=K, symmetric=TRUE, niter=400,
                       nburn=150, thin=2L, use_re=TRUE, group_idx=g, re_idx=1L, calc_loo=FALSE)
f4$mundlak_def <- md
P4 <- predict_gamble(f4, X, group_idx = g, group_levels = unique(g))   # X WITHOUT the MDL cols
ok(max(abs(rowSums(P4) - 1)) < 1e-10, "mundlak: rebuilt from bare X, rows sum to 1")
P4b <- predict_gamble(f4, md$X, group_idx = g, group_levels = unique(g))  # already augmented
ok(max(abs(P4 - P4b)) < 1e-12, "mundlak: rebuild == passing the augmented design")

# ---- wrong design must ERROR, not silently mispredict ----
e <- try(predict_gamble(f1, X[, 1:2, drop=FALSE]), silent=TRUE)
ok(inherits(e, "try-error"), "wrong column count errors instead of mispredicting")

# ---- disk round trip: recovered fit predicts identically to the streamed return value ----
dd <- file.path(tempdir(), "pg_disk"); unlink(dd, recursive=TRUE); dir.create(dd, showWarnings=FALSE)
f5 <- mnlogit_rcpp_sym(X=X, Y=Y, intercept=FALSE, baseline=K, symmetric=TRUE, niter=400, nburn=150,
                       thin=2L, save_posterior_to_disk=TRUE, disk_path=dd, calc_loo=FALSE)
r5 <- recover_mnlogit_posterior(dd)
ok(!is.null(r5$baseline), "disk: baseline recovered")
Pr <- predict_gamble(r5, X)
ok(max(abs(rowSums(Pr) - 1)) < 1e-10, "disk: recovered fit predicts")

# ---- DISPATCH: nested/count route to their own predictors (checked without a full fit) ----
# The nested paths are exercised end-to-end elsewhere (predict_nested_cut on a real saved fit
# returns 200x27 with rowSums 1); here we only assert that predict_gamble ROUTES correctly.
env <- environment()
called <- NULL
assign("predict_nested_cut", function(fit, X_new, ...) { called <<- "nested_cut"; matrix(1/3, nrow(X_new), 3) }, envir = env)
assign("predict_nested_iv",  function(fit, X_new, ...) { called <<- "nested_iv";  matrix(1/3, nrow(X_new), 3) }, envir = env)
stub_nc <- structure(list(), class = "nested_cut"); stub_iv <- structure(list(), class = "nested_iv")
invisible(predict_gamble(stub_nc, X, group_idx = g)); ok(identical(called, "nested_cut"), "dispatch: class nested_cut -> predict_nested_cut")
invisible(predict_gamble(stub_iv, X));                ok(identical(called, "nested_iv"),  "dispatch: class nested_iv  -> predict_nested_iv")
e2 <- try(predict_gamble(list(post_r = 1), X), silent = TRUE)
ok(inherits(e2, "try-error") && grepl("offset", conditionMessage(attr(e2, "condition"))),
   "dispatch: count fit without offset errors clearly")

cat(sprintf("\nRESULT: %d/%d checks passed\n", np, np + nf))
if (nf > 0L) quit(status = 1)
