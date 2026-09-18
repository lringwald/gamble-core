#!/usr/bin/env Rscript
# =============================================================================
# Extract ONE country's rotated parameters for the downscaling engine
# =============================================================================
#   Rscript .../extract_country.R <run dir> <fitted group> [CAPRI code]
#   e.g.    Rscript .../extract_country.R results/.../prior_2026-09-10_2145 Germany DE
#
# Emits covariate x from_class x to_class for that country: median contrast, 95% interval, AND a
# per-row Rhat computed on the CONTRAST ITSELF (not borrowed from the two levels), so a consumer can
# filter to the rows that actually converged rather than trusting the file wholesale.
#
# RELABELLING a GLOB_country group to a CAPRI code is only sound where the two geometries agree.
# Measured by 1km cell: Germany/DE is 99.1% and 97.2% either way (Jaccard 0.96) -- fine. It is NOT
# fine everywhere: Switzerland maps to its best CAPRI code at 0.41 purity, Serbia-Monte 0.70,
# Luxembourg 0.93. The script measures the pair it is asked for and refuses below a threshold.
# =============================================================================
suppressMessages({library(Rcpp); library(RcppArmadillo); library(data.table); library(matrixStats)})
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 2) stop("usage: extract_country.R <run dir> <fitted group> [CAPRI code]")
RUN <- a[1]; if (basename(RUN) == "posterior") RUN <- dirname(RUN)
GRP <- a[2]; CODE <- if (length(a) >= 3) a[3] else GRP
MINPUR <- as.numeric(Sys.getenv("EX_MIN_PURITY", "0.95"))
D <- file.path(RUN, "posterior")
source("codes/mnl_aux_func.R"); source("codes/mnlogit_rcpp_sym.R")
inp <- readRDS(Sys.getenv("BM_INPUT", "output/pixel_model_inputs_BMLEH_Los1_CAPRI.rds"))
gi <- match(GRP, inp$re_group_names)
if (is.na(gi)) stop(GRP, " is not a fitted group. Have: ", paste(inp$re_group_names, collapse = ", "))

# ---- does the relabel hold for THIS pair? -----------------------------------------------------
if (!identical(GRP, CODE)) {
  ok <- tryCatch({
    suppressMessages(library(arrow))
    gw <- Sys.getenv("GAMBLE_GRIDWORK_DIR", "../LAMASUS_gridwork/output")
    gf <- sort(list.files(gw, "^one_kmID_master_mapping_.*\\.parquet$", full.names = TRUE))
    g <- as.data.table(arrow::read_parquet(gf[length(gf)]))
    g[, `:=`(cap2 = substr(as.character(CAPRI_NUTS), 1, 2), glob = as.character(GLOB_country))]
    n_g <- g[glob == GRP, .N]; n_c <- g[cap2 == CODE, .N]; n_b <- g[glob == GRP & cap2 == CODE, .N]
    cat(sprintf("relabel %s -> %s : P(%s|%s) %.4f | P(%s|%s) %.4f | Jaccard %.4f\n",
                GRP, CODE, CODE, GRP, n_b/n_g, GRP, CODE, n_b/n_c, n_b/(n_g + n_c - n_b)))
    min(n_b/n_g, n_b/n_c)
  }, error = function(e) { cat("  (grid unavailable; relabel unchecked)\n"); NA_real_ })
  if (!is.na(ok) && ok < MINPUR)
    stop(sprintf("relabel %s -> %s is only %.3f pure, below EX_MIN_PURITY=%.2f. The two geometries do not agree well enough for this to be a rename.", GRP, CODE, ok, MINPUR))
}

chs <- sort(unique(as.integer(sub(".*_chain_(\\d+)\\.qs$", "\\1",
        list.files(D, "^posterior_batch_.*_chain_\\d+\\.qs$")))))
F <- lapply(chs, function(ci) recover_mnlogit_posterior(D, chain_id = ci))
G <- lapply(F, function(f) f$postb_total[, , gi, ])       # k x p x draws, per chain
K <- dim(G[[1]])[1]; P <- dim(G[[1]])[2]
cn <- F[[1]]$cat_names; vn <- F[[1]]$var_names
nd <- vapply(G, function(x) dim(x)[3], 1L)
cat(sprintf("%s: %d chains x %d draws | %d covariates x %d classes\n", GRP, length(G), nd[1], K, P))

rh <- function(M) { n <- nrow(M); W <- mean(apply(M, 2, var)); if (!is.finite(W) || W <= 0) return(NA_real_)
  sqrt(((n - 1)/n * W + n * var(colMeans(M))/n)/W) }
out <- vector("list", P * (P - 1L)); z <- 1L
pb <- utils::txtProgressBar(min = 0, max = P * (P - 1L), style = 3)
for (fi in seq_len(P)) for (ti in seq_len(P)) {
  if (fi == ti) next
  per_chain <- lapply(G, function(x) x[, ti, ] - x[, fi, ])        # k x draws
  allm <- do.call(cbind, per_chain)
  md <- matrixStats::rowMedians(allm); qq <- matrixStats::rowQuantiles(allm, probs = c(0.025, 0.975))
  rr <- vapply(seq_len(K), function(i) rh(do.call(cbind, lapply(per_chain, function(m) m[i, ]))), 1)
  out[[z]] <- data.table(ks = vn, from_class = cn[fi], to_class = cn[ti], group = CODE,
                         value_median = md, value_q025 = qq[, 1], value_q975 = qq[, 2], rhat = rr)
  z <- z + 1L; utils::setTxtProgressBar(pb, z - 1L)
}
close(pb)
res <- rbindlist(out[seq_len(z - 1L)])
res[, converged := is.finite(rhat) & rhat <= 1.05]
f1 <- file.path(RUN, sprintf("beta_rotated_%s.rds", CODE)); f2 <- sub("\\.rds$", ".csv.gz", f1)
saveRDS(res, f1); fwrite(res, f2)
cat(sprintf("\nrows %s | converged (rhat<=1.05) %.1f%%\n", format(nrow(res), big.mark = ","), 100*mean(res$converged)))
cat("\nconvergence of the EXPORTED contrasts, by destination class (worst first):\n")
s <- res[, .(rows = .N, pct_ok = round(100*mean(converged), 1), rhat_med = round(median(rhat, na.rm = TRUE), 3)), by = to_class][order(pct_ok)]
print(head(s, 8)); cat("\nbest:\n"); print(head(s[order(-pct_ok)], 5))
cat(sprintf("\nwrote %s\n      %s\n", f1, f2))
