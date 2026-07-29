#!/usr/bin/env Rscript
# =============================================================================
# run_nested_cut.R — fit the CUT nested MNL on a real branch (AGMIP/GLOBIOM/BIOCLIMA)
# =============================================================================
# Flexible: the classification tree is built by codes/nest_trees.R (any scheme + subnodes).
# Usage:  Rscript run_nested_cut.R <BRANCH> <arg2> [M] [NITER] [USE_RE] [IVMODE] [STREAM] [N_CHAINS] [N_CORES]
#   BRANCH   AGMIP | GLOBIOM | BIOCLIMA
#   arg2     integer -> subsample n pixels (0 = all), design reconstructed from the newest saved
#                       dat_pixel_FULL + metadata for the branch (auto-discovered).
#            *.rds   -> FROM_INPUTS: read that pixel_model_inputs.rds dump directly (recommended).
#   M imputations   NITER sweeps   USE_RE TRUE/FALSE   IVMODE draws|moments|auto (default auto)
#   STREAM TRUE/FALSE (disk-stream sub-fits, default TRUE)
#   N_CHAINS int (default 1; >1 e.g. 4 -> per-node Rhat/ESS convergence check in fit$convergence)
#   N_CORES  int (default min(N_CHAINS, ncpu-1); chains fitted in parallel via fork on macOS/Linux)
suppressMessages({library(qs2); library(data.table)})
`%||%` <- function(a,b) if (is.null(a)) b else a
source("codes/mnlogit_rcpp_sym.R"); source("codes/nested_cut.R"); source("codes/nest_trees.R")
args <- commandArgs(trailingOnly = TRUE)
BR   <- toupper(args[1] %||% "AGMIP")
# arg2 is EITHER a subsample size (integer) OR a path to a pixel_model_inputs.rds dump
# (from run_prior_module_pixel_level_model.R with DRIVER_DUMP_INPUTS=TRUE) -> FROM_INPUTS mode:
# fits straight off X_mat/Y_pixel/group_idx, no dat_pixel_FULL / model_metadata.qs needed.
A2   <- args[2] %||% "0"
FROM_INPUTS <- if (nzchar(A2) && grepl("\\.rds$", A2) && file.exists(A2)) A2 else NULL
SUB  <- if (is.null(FROM_INPUTS)) as.integer(A2) else 0L
M    <- as.integer(args[3] %||% "25")
NIT  <- as.integer(args[4] %||% "1000")
USE_RE <- as.logical(args[5] %||% "TRUE")     # RE by region (production); FALSE = pooled (fast smoke)
IVMODE <- args[6] %||% "auto"                 # "draws" (exact, M fits) | "moments" (sigma-point, ~3 fits) | "auto" (moments where IV ~Gaussian, else draws; recommended)
STREAM <- as.logical(args[7] %||% "TRUE")     # stream sub-fits to disk (safer RAM at scale); FALSE = in-RAM
NCHAINS <- as.integer(args[8] %||% "1")       # >1 (e.g. 4) -> per-node Rhat/ESS convergence check (fit$convergence)
NCORES  <- as.integer(args[9] %||% as.character(max(1L, min(NCHAINS, parallel::detectCores() - 1L))))  # chains in parallel (fork; macOS/Linux). 1 = serial
STORE   <- Sys.getenv("NCUT_STORE_DIR", sprintf("output/ncut_store_%s", BR))  # persist node draws -> kill-safe/resumable; env NCUT_STORE_DIR=none to disable
if (identical(tolower(STORE), "none") || !nzchar(STORE)) STORE <- NULL

RECON_MIN_GROUP_PIXELS <- 100L   # reconstruction fallback: drop RE regions with fewer pixels (RE stability)
socio <- c("log1p_GDP","log1p_Pop","GHM_HI","log1p_RAI","CISI","allPA_share","Slope_rad","Elevation")
if (!is.null(FROM_INPUTS)) {
  # ---- read the assembled inputs dump directly (no metadata / dat reconstruction) ----
  cat(sprintf(">>> FROM_INPUTS: %s\n", FROM_INPUTS))
  inp <- readRDS(FROM_INPUTS)
  X <- as.matrix(inp$X_mat); if (is.null(colnames(X)) && !is.null(inp$col_names)) colnames(X) <- inp$col_names
  X[!is.finite(X)] <- 0
  Y <- as.matrix(inp$Y_pixel); Y <- Y / rowSums(Y); cats <- colnames(Y)
  group_idx <- as.integer(inp$group_idx_vec)
  coordX <- inp$coord_X; coordY <- inp$coord_Y
  grp_lab <- if (!is.null(inp$re_group_names)) inp$re_group_names[group_idx] else as.character(group_idx)
} else {
  # ---- reconstruct the design from the NEWEST saved dat_pixel_FULL + metadata for this branch ----
  dat_f <- Sys.glob(sprintf("output/dat_pixel_FULL_*%s*.rds", BR))
  dir_f <- Sys.glob(sprintf("output/saved_model_outputs/*%s*", BR))
  if (!length(dat_f) || !length(dir_f))
    stop(sprintf("No saved dat_pixel_FULL / metadata for %s under output/. Pass a pixel_model_inputs.rds dump as arg2, or run the data-prep first (docs/DEVELOPMENT.md).", BR))
  dat_f <- dat_f[which.max(file.mtime(dat_f))]; dir_f <- dir_f[which.max(file.mtime(dir_f))]
  cat(sprintf(">>> RECONSTRUCT  dat: %s | meta: %s\n", dat_f, dir_f))
  meta <- qs_read(file.path(dir_f, "model_metadata.qs")); covs <- meta$cov_names; cats <- meta$cat_names
  covx <- setdiff(covs, "intercept")
  d <- as.data.table(readRDS(dat_f))
  for (v in c("RAI","Pop","GDP")) { lv <- paste0("log1p_",v); if (lv %in% covs && v %in% names(d)) d[[lv]] <- log1p(d[[v]]) }
  miss <- setdiff(covx, names(d)); if (length(miss)) stop("dat is missing design columns: ", paste(head(miss, 8), collapse=", "))
  Ym <- as.matrix(d[, ..cats]); d <- d[rowSums(Ym) > 0]; Ym <- as.matrix(d[, ..cats]); d <- d[rowSums(Ym) >= 1]
  gt <- table(d$Grouping_Key); d <- d[Grouping_Key %in% names(gt)[gt >= RECON_MIN_GROUP_PIXELS]]  # data-driven, not a hardcoded list
  X <- cbind(intercept = 1, as.matrix(d[, ..covx])); X[!is.finite(X)] <- 0
  Yraw <- as.matrix(d[, ..cats]); Y <- Yraw / rowSums(Yraw)
  group_idx <- as.integer(as.factor(d$Grouping_Key))
  coordX <- d$X; coordY <- d$Y; grp_lab <- as.character(d$Grouping_Key)
}
# common: optional subsample + socio RE index
if (SUB > 0 && SUB < nrow(X)) { set.seed(1); sel <- sort(sample(nrow(X), SUB))
  X <- X[sel, , drop = FALSE]; Y <- Y[sel, , drop = FALSE]; group_idx <- group_idx[sel]
  coordX <- coordX[sel]; coordY <- coordY[sel]; grp_lab <- grp_lab[sel] }
re_idx <- sort(intersect(c(which(colnames(X) == "intercept"), which(colnames(X) %in% socio)), seq_len(ncol(X))))

cat(sprintf("\n[%s] n=%d  cov=%d  classes=%d  RE-groups=%d  RE-covs=%d  use_re=%s  M=%d niter=%d\n",
            BR, nrow(X), ncol(X), length(cats), length(unique(group_idx)), length(re_idx), USE_RE, M, NIT))
cat(">>> nesting tree:\n"); tree <- build_nest_tree(BR, cats)

t0 <- Sys.time()
fit <- nested_cut_fit(X, Y, tree, use_iv = TRUE, use_re = USE_RE, group_idx = group_idx, re_idx = re_idx,
                      M = M, draws_per_impute = 15, niter = NIT, nburn = round(NIT/3), thin = 2L, min_pixels = 200,
                      iv_mode = IVMODE, moment_rank = 2L, stream_disk = STREAM, n_chains = NCHAINS, n_cores = NCORES, store_dir = STORE)
cat(sprintf(">>> fit done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))
if (NCHAINS > 1L && !is.null(fit$convergence)) {
  cat(sprintf("\n>>> convergence (%d chains) per node — Rhat / ESS:\n", NCHAINS)); print(fit$convergence) }

# ---- report lambdas (dissimilarity) per internal nest, with 95% CI ----------
s <- summary_nested_cut(fit)
cat("\n>>> Inclusive-value dissimilarity  lambda  (median [95% CI]) per nest:\n")
for (nm in names(s)) if (!is.null(s[[nm]]$lambda)) {
  L <- s[[nm]]$lambda
  for (r in rownames(L)) cat(sprintf("    %-28s %-14s  %.3f  [%.3f, %.3f]\n", nm, r, L[r,"median"], L[r,"q025"], L[r,"q975"]))
}
ts <- format(Sys.time(), "%Y-%m-%d_%H%M")
saveRDS(list(fit = fit, tree = tree, cats = cats, re_idx = re_idx, branch = BR,
             coordX = coordX, coordY = coordY, group = grp_lab, obs = Y, n_chains = NCHAINS),
        sprintf("output/nested_cut_%s_%s.rds", BR, ts))
cat(sprintf(">>> saved output/nested_cut_%s_%s.rds\n", BR, ts))
