# Fast unit tests for nested_cut helpers -- no model fitting, runs in seconds.
# Exists because the M=1 IV bug below needed no fit to catch, yet went undetected: the slow
# end-to-end tests use M > 1, and every real run that DID use M=1 died inside the sampler with a
# generic "X contains NA/Inf" that named neither the column nor the node.
source("codes/nested_cut.R")
np <- 0L; nf <- 0L
ok <- function(cond, msg) { if (isTRUE(cond)) { np <<- np + 1L; cat(sprintf("[PASS] %s\n", msg)) }
                            else { nf <<- nf + 1L; cat(sprintf("[FAIL] %s\n", msg)) } }

mk <- function(M, seed = 1) { set.seed(seed); lapply(seq_len(M), function(i)
  matrix(rnorm(6), 3, 2, dimnames = list(NULL, c("IV_A","IV_B")))) }

# --- M = 1: no spread to carry. Fc is identically zero and the 1/sqrt(M-1) scaling is 1/0, so
# every perturbation field used to come out mu +/- sqrt(3)*NaN. Field 1 was clean, fields 2+ were
# NaN in every cell, and the nest fit got an all-NaN IV column.
r1 <- .ncut_sigma_iv(mk(1), 1L)
ok(length(r1$fields) == 1L, "M=1 yields exactly ONE field (no spread to carry)")
ok(all(vapply(r1$fields, function(f) all(is.finite(f)), TRUE)), "M=1 fields are finite (was NaN)")
ok(abs(sum(r1$weights) - 1) < 1e-12, "M=1 weights sum to 1")
ok(identical(colnames(r1$fields[[1]]), c("IV_A","IV_B")), "M=1 keeps IV column names")
ok(max(abs(r1$fields[[1]] - mk(1)[[1]])) < 1e-12, "M=1 mean field IS the single draw")

# --- M >= 2 unchanged: 1 mean + 2 per retained direction, weights normalised
r3 <- .ncut_sigma_iv(mk(3), 1L)
ok(length(r3$fields) == 3L, "M=3, Q=1 yields 3 fields (mean + 2 perturbations)")
ok(all(vapply(r3$fields, function(f) all(is.finite(f)), TRUE)), "M=3 fields are finite")
ok(abs(sum(r3$weights) - 1) < 1e-12, "M=3 weights sum to 1")
ok(max(abs(r3$fields[[1]] - Reduce(`+`, mk(3)) / 3)) < 1e-12, "M=3 field 1 is the mean field")
# the two perturbations straddle the mean symmetrically
ok(max(abs((r3$fields[[2]] + r3$fields[[3]]) / 2 - r3$fields[[1]])) < 1e-10,
   "M=3 perturbations are symmetric about the mean")

# --- gaussianity gate must not report off all-zero scores at M=1
g1 <- .ncut_iv_gaussianity(mk(1), 1L)
ok(isTRUE(g1$degenerate), "M=1 gaussianity is flagged degenerate, not silently 'gaussian'")

cat(sprintf("\nRESULT: %d/%d unit checks passed\n", np, np + nf))
if (nf > 0L) quit(status = 1)
