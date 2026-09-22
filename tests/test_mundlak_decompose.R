#!/usr/bin/env Rscript
# =============================================================================
# test_mundlak_decompose.R — the country coefficient is mu + gamma'xbar_c + u_c
# =============================================================================
# Under the Mundlak expansion the sampler's b_g is mu + u_c only: the contextual part lives in
# separate MDL_* columns whose VALUE is the country's own mean. Reporting b_g as "the country
# slope" therefore describes the country as if its national context sat at the sample mean.
#
# That is not a rounding difference. The fixture below is built so a country has a NEGATIVE
# residual and the LARGEST total -- naive reporting ranks it last, the truth ranks it first. This
# test pins the decomposition, and pins that inversion, because the inversion is the reason the
# decomposition has to exist.
# =============================================================================
`%||%` <- function(a, b) if (!is.null(a)) a else b
source("codes/mundlak.R")
source("tests/helper_recovery.R")
np <- 0L; nf <- 0L
ok <- function(c_, m) { if (isTRUE(c_)) { np <<- np+1L; cat(sprintf("[PASS] %s\n", m)) }
                        else { nf <<- nf+1L; cat(sprintf("[FAIL] %s\n", m)) } }

set.seed(1); G <- 4; D <- 400
mean_cols <- c("log1p_GDP", "log1p_Pop")
xbar <- rbind(NL = c(2.0, 2.0), PT = c(-1.0, -0.5), DE = c(0.5, 0.8), RO = c(-1.5, -1.3))
colnames(xbar) <- mean_cols
def <- list(xbar = xbar, center = c(0, 0), scale = c(1, 1), mean_cols = mean_cols,
            re_cols = "intercept", names = paste0("MDL_", mean_cols), group_levels = seq_len(G))
cov_names <- c("intercept", paste0("MDL_", mean_cols))
MU <- 0.10; GAMMA <- c(0.15, 0.05); U <- c(NL = -0.05, PT = 0.02, DE = 0.01, RO = -0.03)
pp <- array(0, c(3, 1, D), dimnames = list(cov_names, NULL, NULL))
pp["intercept", 1, ]     <- rnorm(D, MU, 1e-3)
pp["MDL_log1p_GDP", 1, ] <- rnorm(D, GAMMA[1], 1e-3)
pp["MDL_log1p_Pop", 1, ] <- rnorm(D, GAMMA[2], 1e-3)
pt <- array(0, c(3, 1, G, D))
for (g in seq_len(G)) pt[1, 1, g, ] <- pp["intercept", 1, ] + U[g]

d <- mundlak_decompose(pp, pt, def, cov_names, "intercept", group_labels = rownames(xbar))
true_ctx <- as.numeric(xbar %*% GAMMA)

ok(max(abs(d$contextual - true_ctx)) < 5e-3, "contextual term equals gamma' xbar_c")
ok(max(abs(d$u - U[d$group])) < 5e-3,        "residual u_c is recovered")
ok(max(abs(d$total - (d$mu + d$contextual + d$u))) < 1e-9,
   "total = mu + contextual + u, exactly")
ok(all(d$total_lo <= d$total & d$total <= d$total_hi),
   "the interval brackets the reported total")

naive <- d$mu + d$u
ok(which.max(d$total) != which.max(naive),
   sprintf("naive mu+u ranks '%s' top; the true total ranks '%s' top -- the inversion",
           d$group[which.max(naive)], d$group[which.max(d$total)]))
ok(d$group[which.max(d$total)] == "NL" && d$u[d$group == "NL"] < 0,
   "the country with the largest total has a NEGATIVE residual")

ve <- mundlak_variance_explained(d)
ok(ve$explained > 0.9, sprintf("context explains %.1f%% of the between-country spread", 100*ve$explained))

for (i in seq_len(nrow(d)))
  tv_record("Mundlak country total", d$group[i], true_ctx[i] + MU + U[d$group[i]], d$total[i],
            d$total_lo[i], d$total_hi[i], "mu + gamma'xbar_c + u_c, not mu + u_c")

cat(sprintf("RESULT: %d/%d checks passed\n", np, np + nf))
if (nf > 0L) quit(status = 1)
