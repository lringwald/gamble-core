#!/usr/bin/env Rscript
# =============================================================================
# test_nsample.R — nsample states the KEPT draws; niter stays the total
# =============================================================================
# niter is a total, so asking for N draws means passing N + nburn. Getting that backwards is not a
# small error: a caller that passed a total where it meant a sample ran "0 sweeps (burn 4000)", and
# the generic symptom (niter <= nburn) only surfaces once the design is loaded. nsample removes the
# arithmetic; these checks pin that it removes it EXACTLY -- N in, N kept -- and that the old
# calling convention still works, since every existing script uses it.
# =============================================================================
suppressMessages(source("codes/mnl_aux_func.R"))
suppressMessages(source("codes/mnlogit_rcpp_sym.R"))
np <- 0L; nf <- 0L
ok <- function(c_, m) { if (isTRUE(c_)) { np <<- np+1L; cat(sprintf("[PASS] %s\n", m)) }
                        else { nf <<- nf+1L; cat(sprintf("[FAIL] %s\n", m)) } }
set.seed(1); n <- 300; J <- 3
X <- cbind(intercept = 1, a = rnorm(n), b = rnorm(n))
Y <- t(apply(matrix(runif(n*J), n), 1, function(r) r/sum(r)))
fit <- function(...) tryCatch(mnlogit_rcpp_sym(X = X, Y = Y, intercept = FALSE, thin = 1,
                                               calc_loo = FALSE, ...),
                              error = function(e) conditionMessage(e))
kept <- function(r) if (is.character(r)) NA_integer_ else dim(r$postb_pooled)[3]

ok(identical(kept(fit(nsample = 200, nburn = 100)), 200L),
   "nsample=200 keeps exactly 200 draws")
ok(identical(kept(fit(nsample = 150)), 150L),
   "nsample honours the default nburn without the caller computing a total")
ok(identical(kept(fit(niter = 300, nburn = 100)), 200L),
   "the classic niter/nburn convention is unchanged (300-100 = 200 kept)")

e1 <- fit(nsample = 200, niter = 300, nburn = 100)
ok(is.character(e1) && grepl("not both", e1),
   "passing BOTH nsample and niter is refused, not silently resolved")

e2 <- fit(niter = 100, nburn = 200)
ok(is.character(e2) && grepl("niter must be > nburn", e2) && grepl("nsample", e2),
   "the niter<=nburn error still fires AND points at nsample as the fix")

e3 <- fit(nsample = 0, nburn = 10)
ok(is.character(e3) && grepl("nsample must be", e3), "nsample < 1 is rejected")

# the other two samplers must accept the argument, or the convention is only half adopted
for (f in c("codes/mnlogit_rcpp.R", "codes/count_rcpp.R")) {
  src <- paste(readLines(f, warn = FALSE), collapse = "\n")
  ok(grepl("nsample\\s*=\\s*NULL", src), sprintf("%s takes nsample", basename(f)))
}
cat(sprintf("RESULT: %d/%d checks passed\n", np, np + nf))
if (nf > 0L) quit(status = 1)
