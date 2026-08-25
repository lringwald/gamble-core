#!/usr/bin/env Rscript
# =============================================================================
# score_nested_cut.R — score one or more saved nested_cut fits on a common design
# =============================================================================
# Usage:  Rscript codes/score_nested_cut.R <fit1.rds> [fit2.rds ...]
#   env: SCORE_N  (subsample rows for speed; 0 = all, default 0)
#        SCORE_D  (posterior draws per fit; default = fit$M)
#
# THIS IS IN-SAMPLE. Every fit here trained on ALL rows of its design, so these numbers measure
# reproduction, not generalisation, and they STRUCTURALLY FAVOUR the richer random-effects block.
# Use them to check that a fit reproduces its own data (and how it handles rare classes), NOT to
# choose between RE configurations -- that needs a train/test refit (experiments/mixing/
# re_idx_tradeoff.R does it properly on one node).
#
# WHY predict_nested_cut AND NOT THE COEFFICIENT TABLE. The likelihood only ever sees
# beta_g = mu + b_g, and softmax is nonlinear, so neither mu nor mean_g(beta_g) reproduces a
# group-heterogeneous MNL. group_idx MUST be passed: it is matched against the fit's group_levels
# (appearance order), and unseen groups fall back to pooled.
# =============================================================================
suppressMessages({library(qs2); library(data.table)})
source("codes/mnlogit_rcpp_sym.R"); source("codes/nested_cut.R"); source("codes/nest_trees.R")

fits <- commandArgs(trailingOnly = TRUE)
if (!length(fits)) stop("pass at least one nested_cut_*.rds")
SCORE_N <- as.integer(Sys.getenv("SCORE_N", "0"))
SCORE_D <- as.integer(Sys.getenv("SCORE_D", "0"))

softmax_null <- function(Y) { s <- colSums(Y) / sum(Y); matrix(s, nrow(Y), ncol(Y), byrow = TRUE) }
jll  <- function(P, Y) sum(Y * log(pmax(P, 1e-12)))

rows <- list(); comp <- list()
for (f in fits) {
  z <- readRDS(f); nm <- basename(f)
  # --- the design this fit RECORDS, not the newest one on disk (staleness guard) ---
  dsn <- z$inputs
  if (is.null(dsn) || is.na(dsn) || !file.exists(dsn)) stop(nm, ": recorded design not found: ", dsn)
  if (!is.null(z$inputs_mtime) && !is.na(z$inputs_mtime) &&
      !identical(format(file.mtime(dsn)), z$inputs_mtime))
    warning(nm, ": design MTIME CHANGED since the fit (", z$inputs_mtime, " -> ",
            format(file.mtime(dsn)), ") -- scoring against a possibly different design", call. = FALSE)
  inp <- readRDS(dsn)
  X <- as.matrix(inp$X_mat); if (is.null(colnames(X))) colnames(X) <- inp$col_names
  X[!is.finite(X)] <- 0
  Y <- as.matrix(inp$Y_pixel); Y <- Y / rowSums(Y)
  # KEY ON THE INTEGER group_idx, NOT the country NAMES. The fit's group_levels are the sampler's
  # unique(group_idx) in APPEARANCE ORDER, stored as characters of the INTEGER ids ("20","24",...),
  # whereas the saved fit's `group` field holds human labels ("Portugal"). Matching labels against
  # integer levels returns all-NA and predict_nested_cut then silently falls back to POOLED
  # coefficients -- which do not reproduce a group-heterogeneous MNL at all. This produced a fully
  # bogus A/B/C comparison once already (2026-08-18).
  g <- inp$group_idx_vec
  if (nrow(X) != length(g)) stop(nm, ": design has ", nrow(X), " rows but group_idx_vec has ", length(g))
  idx <- seq_len(nrow(X))
  if (SCORE_N > 0 && SCORE_N < nrow(X)) { set.seed(42); idx <- sort(sample(nrow(X), SCORE_N)) }
  Xs <- X[idx, , drop = FALSE]; Ys <- Y[idx, z$cats, drop = FALSE]; gs <- g[idx]

  # FAIL LOUDLY on a keying mismatch: a silent pooled fallback looks like a plausible score and is not.
  lev <- .ncut_group_levels(z$fit$root)
  if (!is.null(lev)) {
    nmiss <- sum(is.na(match(as.character(gs), lev)))
    if (nmiss == length(gs)) stop(nm, ": NONE of the ", length(gs), " rows matched the fit's group_levels ",
      "(first levels: ", paste(head(lev, 4), collapse = ","), ") -- wrong group key, scoring would be pooled-only")
    if (nmiss > 0) cat(sprintf("  note: %d/%d rows in groups unseen by the fit -> pooled for those\n", nmiss, length(gs)))
  }

  D <- if (SCORE_D > 0) SCORE_D else z$fit$M
  t0 <- Sys.time()
  P <- predict_nested_cut(z$fit, Xs, D = D, summary = "mean", group_idx = gs)
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  ll <- jll(P, Ys); ll0 <- jll(softmax_null(Ys), Ys)
  rows[[nm]] <- data.frame(fit = nm,
    re = z$re_cols %||% NA, sym_hs = z$symmetric_hs %||% NA, n = nrow(Xs), D = D,
    jll = round(ll, 1), null = round(ll0, 1), mcfadden = round(1 - ll / ll0, 4),
    mad = round(mean(abs(P - Ys)), 5), mins = round(mins, 1), stringsAsFactors = FALSE)
  comp[[nm]] <- data.frame(class = z$cats,
    obs = 100 * colSums(Ys) / sum(Ys), prd = 100 * colSums(P) / sum(P), stringsAsFactors = FALSE)
  cat(sprintf("scored %-52s jll %.1f  McF %.4f  [%.1f min]\n", nm, ll, 1 - ll / ll0, mins))
}

cat("\n================ IN-SAMPLE SCORE (see header: favours richer RE) ================\n")
print(do.call(rbind, rows), row.names = FALSE)

cat("\n================ COMPOSITION: observed vs predicted (% of area) ================\n")
base <- comp[[1]][order(-comp[[1]]$obs), "class"]
tab <- data.frame(class = base, obs = round(comp[[1]]$obs[match(base, comp[[1]]$class)], 3))
for (nm in names(comp)) tab[[paste0("ratio_", substr(nm, nchar(nm) - 7, nchar(nm) - 4))]] <-
  round(comp[[nm]]$prd[match(base, comp[[nm]]$class)] / comp[[nm]]$obs[match(base, comp[[nm]]$class)], 3)
print(tab, row.names = FALSE)
cat("\nratio = predicted/observed share. <1 = under-predicted. Watch the RAREST classes:\n",
    "a systematic ratio far below 1 on small classes is the composition bias to chase.\n")
