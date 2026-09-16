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
# Fix macOS Accelerate/OpenMP fork bug
Sys.setenv(OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)

suppressMessages({library(qs2); library(data.table)})
`%||%` <- function(a,b) if (is.null(a) || length(a)==0 || is.na(a)) b else a
source("codes/mnlogit_rcpp_sym.R"); source("codes/nested_cut.R"); source("codes/nest_trees.R")
# Accept arguments from an interactive session too: define NCUT_ARGS (character vector) before
# source()ing this file, e.g. from run/nested.R in RStudio. Falls back to the command line.
args <- if (exists("NCUT_ARGS", inherits = TRUE)) as.character(get("NCUT_ARGS", inherits = TRUE)) else commandArgs(trailingOnly = TRUE)
BR   <- toupper(args[1] %||% "AGMIP")
# arg2 is EITHER a subsample size (integer) OR a path to a pixel_model_inputs.rds dump
# (from run_prior_module_pixel_level_model.R with DRIVER_DUMP_INPUTS=TRUE) -> FROM_INPUTS mode:
# fits straight off X_mat/Y_pixel/group_idx, no dat_pixel_FULL / model_metadata.qs needed.
A2   <- args[2] %||% "0"
FROM_INPUTS <- if (nzchar(A2) && grepl("\\.rds$", A2) && file.exists(A2)) A2 else NULL
SUB  <- if (is.null(FROM_INPUTS)) as.integer(A2) else as.integer(Sys.getenv("NCUT_SUB", "0"))
M    <- as.integer(args[3] %||% "25")
NIT  <- as.integer(args[4] %||% "1000")
USE_RE <- as.logical(args[5] %||% "TRUE")     # RE by region (production); FALSE = pooled (fast smoke)
IVMODE <- args[6] %||% "auto"                 # "draws" (exact, M fits) | "moments" (sigma-point, ~3 fits) | "auto" (moments where IV ~Gaussian, else draws; recommended)
STREAM <- as.logical(args[7] %||% "TRUE")     # stream sub-fits to disk (safer RAM at scale); FALSE = in-RAM
NCHAINS <- as.integer(args[8] %||% "1")       # >1 (e.g. 4) -> per-node Rhat/ESS convergence check (fit$convergence)
NCORES  <- as.integer(args[9] %||% as.character(max(1L, min(NCHAINS, parallel::detectCores() - 1L))))  # chains in parallel (fork; macOS/Linux). 1 = serial
STORE   <- Sys.getenv("NCUT_STORE_DIR", sprintf("output/ncut_store_%s", BR))  # persist node draws -> kill-safe/resumable; env NCUT_STORE_DIR=none to disable
# focal routing at IV nodes (lambda identification; see .ncut_focal_xmap):
#   macro_totals (default) = fine focal of each IV child -> one focal_<child>_tot beside the IVs
#   leaf_only              = drop them (IV is the sole channel -> cleanest lambda)
#   flat                   = legacy no routing (lambda unidentified; comparison only)
FOCAL_RULE <- Sys.getenv("NCUT_FOCAL_RULE", "macro_totals")
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
  # CURATED TREE: the driver dumps `class_nest` (leaf + ancestor-path columns) when the mapping
  # carries a `<CLASS_COL>_nest` column. It is the ONLY route that can restructure the hierarchy
  # without renaming classes, so it must WIN over the coded prefix scheme -- otherwise a curated
  # run silently fits the prefix tree (e.g. Wetlands_natural as its own 1-member nest instead of
  # sitting under Natural). NCUT_TREE=<scheme> forces the coded builder for comparison runs.
  CLASS_NEST <- inp$class_nest
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

# --- Mundlak (NCUT_MUNDLAK) -------------------------------------------------------------------
# Normally the pixel driver adds these at design-assembly time (DRIVER_MUNDLAK) and they arrive
# inside the dumped X_mat, so nothing is needed here. This switch covers an OLDER dump, or the
# reconstruct branch. Guarded on the columns not already being present: applying it twice would
# be a silent duplicate design (build_mundlak_design() also refuses the name collision).
if (isTRUE(as.logical(Sys.getenv("NCUT_MUNDLAK", "FALSE")))) {
  if (any(grepl("^MDL_", colnames(X)))) {
    cat(">>> NCUT_MUNDLAK: columns already present in the design (from DRIVER_MUNDLAK); not re-adding.\n")
  } else {
    source("codes/mundlak.R")
    .mc <- intersect(trimws(strsplit(Sys.getenv("NCUT_MUNDLAK_COLS",
             "CISI,log1p_Pop,log1p_GDP,GHM_HI"), ",")[[1]]), colnames(X))
    # see run_lu_pixel_model.R: intercept-only can displace the correlation into the slopes
    .rc <- trimws(strsplit(Sys.getenv("NCUT_MUNDLAK_RE", "intercept,log1p_Pop"), ",")[[1]])
    .rc <- .rc[.rc == "intercept" | .rc %in% colnames(X)]
    if (!length(.mc) || !length(.rc)) {
      warning("NCUT_MUNDLAK set but no usable NCUT_MUNDLAK_COLS / NCUT_MUNDLAK_RE in X; skipping.")
    } else {
      .md <- build_mundlak_design(X, group_idx, mean_cols = .mc, re_cols = .rc)
      X <- .md$X
      cat(sprintf(">>> Mundlak ON: +%d column(s) [means: %s | RE rows: %s]\n",
                  length(.md$names), paste(.mc, collapse = ","), paste(.rc, collapse = ",")))
    }
  }
}

# common: optional subsample + socio RE index
if (SUB > 0 && SUB < nrow(X)) { set.seed(1); sel <- sort(sample(nrow(X), SUB))
  X <- X[sel, , drop = FALSE]; Y <- Y[sel, , drop = FALSE]; group_idx <- group_idx[sel]
  coordX <- coordX[sel]; coordY <- coordY[sel]; grp_lab <- grp_lab[sel]
  # A subsample can starve the rarest classes (Forests_SR is 0.08% of area) -> report it, so a smoke
  # whose lambdas are driven by an empty leaf is visible rather than mistaken for a real result.
  .sh <- colSums(Y) / sum(Y); .dead <- names(.sh)[.sh <= 0]
  cat(sprintf(">>> SUBSAMPLE n=%d (smoke; NOT a production fit) | rarest class %s %.3f%%%s\n",
              SUB, names(which.min(.sh[.sh > 0])), 100 * min(.sh[.sh > 0]),
              if (length(.dead)) sprintf(" | ZERO-AREA now: %s", paste(.dead, collapse = ",")) else "")) }
# RE BLOCK. Default stays intercept+socio for continuity with earlier fits, but `intercept` is the
# configuration to prefer for new fits. RE-MEASURED 2026-08-17 on the 27-class design (Forests, n=6000,
# 4 chains x 2000-4000 iter, country-stratified 75/25): random slopes are NO LONGER harmful -- the old
# "worse than no RE at all" result (McFadden 0.188 vs 0.331) was a PRIOR artifact, cured by the RE slab
# + support prior + a working FE horseshoe. Slopes now TIE the intercept on held-out (-1060.8 vs
# -1058.7). They just do not earn their keep: 1404 extra per-group parameters (9 x 6 x 26) buy 0 nats
# and cost convergence (RE Rhat 1.52 vs 1.09).
RE_COLS <- Sys.getenv("NCUT_RE_COLS", "intercept+socio")
re_idx <- switch(RE_COLS,
  "intercept"       = which(colnames(X) == "intercept"),
  "intercept+socio" = sort(intersect(c(which(colnames(X) == "intercept"), which(colnames(X) %in% socio)), seq_len(ncol(X)))),
  stop("NCUT_RE_COLS must be 'intercept' or 'intercept+socio', got: ", RE_COLS))

# FACTORIZED vs IV-NESTED. use_iv=FALSE fits each level independently (no logsum coupling, no lambda),
# which is the honest model when lambda is not identified -- and it is far cheaper, since the root is
# fit ONCE instead of once per imputation. NOTE: with IV off there is no focal routing anywhere
# (xmap is gated on use_iv, nested_cut.R:520), so every node sees the FULL focal block; that block IS
# constant-sum once focal_NODATA exists, which puts the drop-one + zero-sum reconstruction path in
# charge of the focal columns. FOCAL_RULE is inert here.
USE_IV <- as.logical(Sys.getenv("NCUT_USE_IV", "TRUE"))
cat(sprintf("\n[%s] n=%d  cov=%d  classes=%d  RE-groups=%d  RE-covs=%d(%s)  use_re=%s  use_iv=%s  M=%d niter=%d  focal_rule=%s\n",
            BR, nrow(X), ncol(X), length(cats), length(unique(group_idx)), length(re_idx), RE_COLS,
            USE_RE, USE_IV, M, NIT, if (USE_IV) FOCAL_RULE else "n/a (factorized)"))
NCUT_TREE <- Sys.getenv("NCUT_TREE", "")      # non-empty -> force this coded scheme, ignore the curated column
if (!exists("CLASS_NEST")) CLASS_NEST <- NULL
# NCUT_TREE=flat -> LEVEL (flat) MNL: one softmax over every fine class, no nesting, no IV.
# nested_cut treats a single-node tree as exactly that (see codes/nested_cut.R:5).
if (identical(tolower(NCUT_TREE), "flat")) {
  cat(sprintf(">>> nesting tree: FLAT / LEVEL model -- one softmax over %d classes\n", length(cats)))
  # Each class is a LEAF CHILD OF THE ROOT -> the root node is the single softmax over all classes.
  # NOT list(all = cats): that wraps them in an intermediate node, so nested_cut fits the real
  # softmax at root/all and then tries to fit the root itself as a choice among ONE alternative
  # (p_all = 1 -> empty mu_pooled -> "invalid 'type' (list)" in the separation scan).
  tree <- setNames(as.list(cats), cats)
  USE_IV <- FALSE                     # no internal nodes -> nothing for an inclusive value to carry
} else if (!is.null(CLASS_NEST) && !nzchar(NCUT_TREE)) {
  .cc <- names(CLASS_NEST)[1]; .nc <- names(CLASS_NEST)[2]
  cat(sprintf(">>> nesting tree: CURATED from dump (%s -> %s, %d mapped leaves)\n",
              .cc, .nc, nrow(CLASS_NEST)))
  tree <- nest_tree_from_mapping(as.data.frame(CLASS_NEST), cats, class_col = .cc, nest_col = .nc)
  if (!nest_tree_check(tree, cats)) stop("curated nest tree does not cover the observed classes -- see MISSING/EXTRA above")
} else {
  SCHEME <- if (nzchar(NCUT_TREE)) NCUT_TREE else BR
  cat(sprintf(">>> nesting tree: CODED scheme %s%s\n", SCHEME,
              if (is.null(CLASS_NEST)) " (no curated class_nest in the dump)" else " (forced via NCUT_TREE)"))
  tree <- build_nest_tree(SCHEME, cats)
}

t0 <- Sys.time()
fit <- nested_cut_fit(X, Y, tree, use_iv = USE_IV, use_re = USE_RE, group_idx = group_idx, re_idx = re_idx,
                      M = M, draws_per_impute = 15, niter = NIT, nburn = round(NIT/3), thin = 2L, min_pixels = 200,
                      iv_mode = IVMODE, moment_rank = 2L, stream_disk = STREAM, n_chains = NCHAINS, n_cores = NCORES,
                      focal_rule = FOCAL_RULE, store_dir = STORE)
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
# ---- lambda IDENTIFICATION: R^2 of each IV on its own node design (high => lambda not
# separately identified from the direct effects; that is what the focal routing controls) ----
idt <- nested_cut_identification(fit)
if (!is.null(idt)) { cat("\n>>> IV identification (focal_rule = ", FOCAL_RULE, "):\n", sep = "")
  print(idt[, c("node","iv_child","iv_r2","lambda","lambda_q025","lambda_q975","n_focal_routed")], row.names = FALSE)
  if (any(idt$iv_r2 > 0.9, na.rm = TRUE))
    cat("    WARNING: IV~design R2 > 0.9 at some node -> lambda is weakly identified there.\n") }

ts <- format(Sys.time(), "%Y-%m-%d_%H%M")
out_f <- sprintf("output/nested_cut_%s%s_%s.rds", BR, if (USE_IV) "" else "_factorized", ts)
saveRDS(list(fit = fit, tree = tree, cats = cats, re_idx = re_idx, branch = BR, identification = idt,
             focal_rule = if (USE_IV) FOCAL_RULE else NA_character_, use_iv = USE_IV, re_cols = RE_COLS,
             # PRIOR SWITCHES that change the fit but are set by ENV, so nothing else records them.
             # Without this, two fits differing only in symmetric_hs are indistinguishable once the
             # shell history is gone (exactly what happened to the 2026-08-18 B/C pair).
             # DEFAULTS MUST TRACK .ncut_fit_block (both flipped to TRUE on 2026-09-16). This is a
             # RECORD of what was fitted: if the default here disagrees with the one the sampler
             # actually used, an unset run is recorded as diagonal while being fitted symmetric --
             # a provenance field that lies, which is worse than no field at all.
             symmetric_hs = isTRUE(as.logical(Sys.getenv("NCUT_SYM_HS", "TRUE"))),
             hs_no_intercept = isTRUE(as.logical(Sys.getenv("NCUT_HS_NO_INTERCEPT", "TRUE"))),
             # PROVENANCE: which design and which tree this fit actually used. A silently reused/stale
             # design dump has cost a full multi-hour run before -- record it so the fit can be audited.
             inputs = FROM_INPUTS %||% NA_character_,
             inputs_mtime = if (!is.null(FROM_INPUTS)) format(file.mtime(FROM_INPUTS)) else NA_character_,
             tree_source = if (!is.null(CLASS_NEST) && !nzchar(NCUT_TREE)) "curated:class_nest" else paste0("coded:", if (nzchar(NCUT_TREE)) NCUT_TREE else BR),
             coordX = coordX, coordY = coordY, group = grp_lab, obs = Y, n_chains = NCHAINS,
             # The INTEGER key prediction needs: group_levels are unique(group_idx) in appearance
             # order stored as characters of these ids. `group` above is human labels and will NOT
             # match -- passing it to predict_nested_cut silently degrades to pooled coefficients.
             group_idx = group_idx),
        out_f)
cat(sprintf(">>> saved %s\n", out_f))
