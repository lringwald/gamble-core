#!/usr/bin/env Rscript
# =============================================================================
# recover_bart_fit.R — report a BART fit whose own post-fit stage never ran
# =============================================================================
#   Rscript postprocess/recover_bart_fit.R <batch dir> <design dump> [thin]
#
# WHY THIS EXISTS, AND WHY IT IS NOT engine.R.
# engine.R reconstructs the design from dat_pixel_FULL -- the PRE-FILTER dump, "including dummy
# years" -- so it lands on a different row and group set than the fit used (65,313 / 29 against
# 64,472 / 30 on the 2026-09-19 GLOBIOM run) and stops. Worse, it has no BART support at all: it
# builds utilities from the linear coefficients only, and on a fit whose terrain sits entirely in
# BART that silently deletes 22.9% of the model's gain over null. A report that confident and that
# wrong is worse than no report.
#
# The right inputs are a design REBUILT THROUGH THE DRIVER -- same code path, therefore same
# filtering, verified row-for-row -- and predict_shares(), which evaluates the stored slim trees.
#
#   DRIVER_USE_BART=TRUE DRIVER_BART_COLS=topo DRIVER_ADD_COORDS=TRUE \
#   DRIVER_DUMP_INPUTS=TRUE DRIVER_DUMP_EXIT=TRUE DRIVER_DUMP_PATH=<dump> \
#   Rscript drivers/run_lu_pixel_model.R        # ~8 min, no refit
#
# The design MUST match the fit: bart_idx indexes the full design, so a dump built without
# ADD_COORDS puts lon/lat nowhere and the indices point at the wrong columns.
# =============================================================================
suppressMessages({library(qs2); library(data.table)})
suppressMessages(source("codes/mnl_aux_func.R")); suppressMessages(source("codes/mnlogit_rcpp_sym.R"))
suppressMessages(source("codes/prior_model_predict.R"))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 2) stop("usage: recover_bart_fit.R <batch dir> <design dump> [thin]")
D <- a[1]; DUMP <- a[2]; THIN <- if (length(a) > 2) as.integer(a[3]) else 30L

m   <- qs_read(file.path(D, "model_metadata.qs"))
inp <- readRDS(DUMP)
X <- as.matrix(inp$X_mat); X[!is.finite(X)] <- 0
Y <- as.matrix(inp$Y_pixel); Y <- Y / rowSums(Y)
g <- as.integer(inp$group_idx_vec)
bi <- as.integer(m$bart_idx); li <- setdiff(seq_len(ncol(X)), bi)

# FAIL LOUDLY ON A MISMATCHED DESIGN. Every number below is meaningless if the dump is not the one
# the fit used, and the failure would otherwise be silent -- plausible predictions from the wrong
# columns. cov_names holds the LINEAR columns only, which is what makes this check exact.
if (!identical(colnames(X)[li], m$cov_names))
  stop("design mismatch: linear columns differ from the fit's cov_names.\n",
       "  Rebuild the dump with the SAME driver settings (ADD_COORDS especially).")
if (!identical(as.integer(m$group_levels), unique(g)))
  stop("design mismatch: group appearance order differs from the fit's group_levels.")
cat(sprintf(">>> design verified: %d rows x %d cols | %d linear | %d BART: %s\n",
            nrow(X), ncol(X), length(li), length(bi), paste(colnames(X)[bi], collapse = ", ")))

chs <- sort(unique(as.integer(sub(".*_chain_(\\d+)\\.qs$", "\\1",
        list.files(D, "^posterior_batch_.*_chain_\\d+\\.qs$")))))
ll0 <- sum(Y * log(pmax(matrix(colMeans(Y), nrow(Y), ncol(Y), byrow = TRUE), 1e-12)))
res <- list()
for (ci in chs) {
  f <- recover_mnlogit_posterior(D, chain_id = ci)
  P <- predict_shares(f, X, bart_cols = bi, linear_cols = li, group_idx = g,
                      group_levels = as.integer(m$group_levels), type = "mean", thin = THIN)
  ll <- sum(Y * log(pmax(P, 1e-12)))
  # The linear-only arm is not a curiosity: it is the check that the trees were actually evaluated.
  B  <- apply(f$postb_total, c(1,2,3), mean); gpos <- match(g, as.integer(m$group_levels))
  Xl <- X[, li, drop = FALSE]
  U  <- t(vapply(seq_len(nrow(Xl)), function(i) as.numeric(Xl[i, ] %*% B[, , gpos[i]]), numeric(ncol(Y))))
  U  <- U - apply(U, 1, max); Pl <- exp(U); Pl <- Pl / rowSums(Pl)
  llL <- sum(Y * log(pmax(Pl, 1e-12)))
  cat(sprintf("chain %d: McFadden %.4f (linear-only %.4f) | BART adds %.0f nats\n",
              ci, 1 - ll/ll0, 1 - llL/ll0, ll - llL))
  res[[as.character(ci)]] <- list(P = P, ll = ll, ll_linear = llL, ll0 = ll0)
}
out <- file.path(dirname(D), paste0(basename(D), "_recovered_predictions.rds"))
saveRDS(list(per_chain = res, ll0 = ll0, cats = m$cat_names, bart_cols = colnames(X)[bi]), out)
cat(sprintf("\nwrote %s\n", out))
