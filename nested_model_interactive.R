# =============================================================================
# nested_model_interactive.R  —  run the CUT nested land-use model from RStudio
# =============================================================================
# HOW TO USE
#   1. Open this file in RStudio; set the working directory to the repo ROOT
#      (Session > Set Working Directory > To Source File Location, then up if needed;
#       the guard below checks it).
#   2. Edit the CONTROL PANEL. Then either Source the whole file, or run it
#      section-by-section (Ctrl/Cmd-Alt-B runs to here; or step block by block).
#   3. The fitted object lives in `fit`; the INSPECT block at the bottom has ready
#      calls for lambdas, predictions, the effective coefficient table, and maps.
#
# Prerequisite: a design to fit on. Options (CONTROL PANEL):
#   INPUTS = "auto"  -> newest output/pixel_model_inputs_<BRANCH>_*.rds
#   INPUTS = "<path>"-> that dump
#   BUILD_INPUTS=TRUE-> build the design, then use it. BUILD_FROM = "reconstruct" (fast, from an existing
#                       dat_pixel_FULL + metadata) OR "dataprep" (run the raster data-prep to build
#                       dat_pixel_FULL from scratch — needs the gridwork data ~5GB; see docs/DEVELOPMENT.md §1).
# =============================================================================

## ============================ CONTROL PANEL ================================ ##
BRANCH       <- "GLOBIOM"     # "AGMIP" | "GLOBIOM" | "BIOCLIMA"
PIXEL_RES    <- 10L           # grid resolution in km (10 ~ 64k pixels | 5 ~ 180k); selects which dat_pixel_FULL to reconstruct from
PIXEL_INTERSECT <- "NUTS3"    # region the pixels were split by in data-prep (matches the dat filename "by_<X>"): "NUTS3" | "CAPRI_NUTS" | "native"
INPUTS       <- "auto"        # "auto" = newest output/pixel_model_inputs_<BRANCH>_<INTERSECT>_<RES>km_*.rds | a dump path | NULL = reconstruct on the fly (not saved)
BUILD_INPUTS <- NULL         # TRUE = (re)build the design, then use it (how = BUILD_FROM)
BUILD_FROM   <- "reconstruct" # "reconstruct" = fast, from an EXISTING dat_pixel_FULL + metadata (no gridwork)
                              # "dataprep"    = run the FULL raster data-prep to BUILD dat_pixel_FULL from scratch (needs gridwork ~5GB)
RECON_SEARCH <- NULL          # extra dir to ALSO search for saved dat_pixel_FULL + metadata when reconstructing
                              # (e.g. "../LAMASUS_high_res_prior_models/output"); the dump is always saved to ./output. NULL = only ./output

PRESET      <- "custom"        # "smoke" (fast) | "production" (full) | "custom" (use the values below verbatim)

# --- knobs (a PRESET overrides M/NITER/USE_RE/SUBSAMPLE unless PRESET = "custom") ---
SUBSAMPLE   <- 0              # 0 = all pixels; else n pixels for a quick test
USE_RE      <- TRUE           # TRUE = country random effects (production) | FALSE = pooled (fast)
M           <- 30             # inclusive-value imputations
NITER       <- 6000           # Gibbs sweeps (speed config: 6000 usually enough; check fit$convergence Rhat)
IVMODE      <- "moments"      # speed config: 1+2*MOMENT_RANK fits per nest (5) instead of up to M (30). "auto"|"draws" = slower/exact
MOMENT_RANK <- 2L             # sigma-point directions for moments/auto
DRAWS_PER_IMPUTE <- 15L
THIN        <- 5L
MIN_PIXELS  <- 200L           # nests sparser than this degrade to an even split
STREAM_DISK <- TRUE           # stream each sub-fit to disk (safer RAM at scale, niter-independent) | FALSE = faster in-RAM for small runs
STORE_DIR   <- sprintf("output/ncut_store_%s_%s_%dkm", BRANCH, PIXEL_INTERSECT, PIXEL_RES)  # persist node draws here -> KILL-SAFE / RESUMABLE: a re-source with the same config skips finished nodes. NULL = off. (Change a model knob -> use a fresh path or delete this folder.)
N_CHAINS    <- 4L             # >1 (e.g. 4) runs a multi-chain convergence check -> per-node Rhat/ESS in fit$convergence. Cost ~ N_CHAINS x the leaf + per-nest-conv fits (imputations stay single-chain)
N_CORES     <- max(1L, min(N_CHAINS, parallel::detectCores() - 1L))  # chains fitted in parallel (fork; macOS/Linux only, Windows -> serial). 1 = serial. Cuts the N_CHAINS wall-clock ~ N_CORES x
SAVE        <- TRUE           # write output/nested_cut_<BRANCH>_<ts>.rds
SEED        <- 1L
## ========================================================================== ##

# ---- preset expansion -------------------------------------------------------
if (PRESET == "smoke")      { SUBSAMPLE <- 8000; USE_RE <- FALSE; M <- 15; NITER <- 800 }
if (PRESET == "production") { SUBSAMPLE <- 0;    USE_RE <- TRUE;  M <- 30; NITER <- 12000; N_CHAINS <- 4L }

# ---- environment + dependencies --------------------------------------------
if (!dir.exists("codes")) stop("Set the working directory to the gamble-core repo ROOT (no 'codes/' here).")
suppressMessages({ library(qs2); library(data.table) })
`%||%` <- function(a, b) if (is.null(a)) b else a
source("codes/mnlogit_rcpp_sym.R"); source("codes/nested_cut.R"); source("codes/nest_trees.R")

# ---- inputs: resolve to `inp` (a pixel_model_inputs list) -------------------
RECON_MIN_GROUP_PIXELS <- 100L   # reconstruction: drop RE regions with fewer pixels than this (RE stability)
socio <- c("log1p_GDP","log1p_Pop","GHM_HI","log1p_RAI","CISI","allPA_share","Slope_rad","Elevation")

# Reconstruct a pixel_model_inputs list from the NEWEST saved dat_pixel_FULL + metadata for `branch`.
# When save=TRUE, writes a DATED, branch-named dump carrying provenance (source_dat / source_meta / built),
# so every dump is traceable to the flat run it came from.
reconstruct_inputs <- function(branch, res = PIXEL_RES, intersect = PIXEL_INTERSECT,
                               min_group_pixels = RECON_MIN_GROUP_PIXELS, save = TRUE,
                               search = c("output", RECON_SEARCH)) {
  search <- unique(search[!is.na(search) & nzchar(search)])
  itag <- if (nzchar(intersect) && !intersect %in% c("native","NULL")) sprintf("by_%s*", intersect) else ""
  dat_f <- unlist(lapply(search, function(s) Sys.glob(file.path(s, sprintf("dat_pixel_FULL_%dkm*%s%s*.rds", res, itag, branch)))))
  dir_f <- unlist(lapply(search, function(s) Sys.glob(file.path(s, sprintf("saved_model_outputs/%dkm*%s%s*", res, itag, branch)))))
  if (!length(dat_f) || !length(dir_f))
    stop(sprintf("reconstruct_inputs(%s, %dkm, %s): no saved dat_pixel_FULL / metadata found in {%s}. Run the data-prep (docs/DEVELOPMENT.md §1), set RECON_SEARCH / PIXEL_RES / PIXEL_INTERSECT, or point INPUTS at a dump.",
                 branch, res, intersect, paste(search, collapse = ", ")))
  dat_f <- dat_f[which.max(file.mtime(dat_f))]; dir_f <- dir_f[which.max(file.mtime(dir_f))]
  message(">>> reconstruct [", branch, "]  dat: ", basename(dat_f), " | meta: ", basename(dir_f))
  meta <- qs_read(file.path(dir_f, "model_metadata.qs")); covs <- meta$cov_names; cts <- meta$cat_names; covx <- setdiff(covs, "intercept")
  d <- as.data.table(readRDS(dat_f))
  for (v in c("RAI","Pop","GDP")) { lv <- paste0("log1p_",v); if (lv %in% covs && v %in% names(d)) d[[lv]] <- log1p(d[[v]]) }
  miss <- setdiff(covx, names(d)); if (length(miss)) stop("dat is missing design columns: ", paste(head(miss, 8), collapse=", "))
  Ym <- as.matrix(d[, ..cts]); d <- d[rowSums(Ym) > 0]; Ym <- as.matrix(d[, ..cts]); d <- d[rowSums(Ym) >= 1]
  gt <- table(d$Grouping_Key); d <- d[Grouping_Key %in% names(gt)[gt >= min_group_pixels]]  # data-driven, not hardcoded
  X <- cbind(intercept = 1, as.matrix(d[, ..covx])); X[!is.finite(X)] <- 0
  Y <- as.matrix(d[, ..cts]); Y <- Y / rowSums(Y); gf <- as.factor(d$Grouping_Key)
  inp <- list(X_mat = X, Y_pixel = Y, group_idx_vec = as.integer(gf), re_group_names = levels(gf),
              col_names = colnames(X), coord_X = d$X, coord_Y = d$Y,
              branch = branch, pixel_res = res, pixel_intersect = intersect,
              source_dat = dat_f, source_meta = dir_f, built = format(Sys.time(), "%Y-%m-%d_%H%M"))
  ilab <- if (nzchar(intersect)) intersect else "native"
  if (save) { p <- sprintf("output/pixel_model_inputs_%s_%s_%dkm_%s.rds", branch, ilab, res, inp$built); saveRDS(inp, p); inp$path <- p
    message(">>> saved dated dump -> ", p) }
  inp
}

# Per-branch env for the raster data-prep (selects the class scheme + LUM lineage).
# VERIFY against the DRIVER_* block at the top of run_prior_module_pixel_level_model.R.
BRANCH_DATAPREP_ENV <- list(
  AGMIP    = c("DRIVER_CLASS_COLS=AgMIP_label", "DRIVER_CROP_SPLIT=TRUE"),
  GLOBIOM  = c("DRIVER_CLASS_COLS=GLOBIOM_UNFCCC,GLOBIOM_mngmt", "DRIVER_LUM_LINEAGE=globiom", "DRIVER_INCLUDE_ORGANIC=TRUE"),
  BIOCLIMA = c("DRIVER_LUM_LINEAGE=bioclima")
)

# Build dat_pixel_FULL FROM SCRATCH: shells out to the pixel driver (data-prep only, exits before the fit),
# which reads the gridwork rasters and writes dat_pixel_FULL_* + output/pixel_model_inputs.rds. Then loads it.
run_dataprep <- function(branch, res = PIXEL_RES, intersect = PIXEL_INTERSECT) {
  env <- c(BRANCH_DATAPREP_ENV[[branch]] %||% character(0),
           sprintf("DRIVER_PIXEL_RES=%d", res),
           sprintf("DRIVER_PIXEL_INTERSECT=%s", if (identical(intersect, "native")) "" else intersect),
           "DRIVER_RUN_MODE=production", "DRIVER_DUMP_INPUTS=TRUE", "DRIVER_DUMP_EXIT=TRUE")
  message(">>> DATA-PREP [", branch, " ", res, "km ", intersect, "]: running run_prior_module_pixel_level_model.R",
          " (reads gridwork ~5GB; several minutes)\n    env: ", paste(env, collapse = " "))
  st <- system2("Rscript", "run_prior_module_pixel_level_model.R", env = env)
  if (!identical(as.integer(st), 0L)) stop("data-prep failed (exit ", st, "). Check console output, GRIDWORK data, and the branch env config.")
  g <- "output/pixel_model_inputs.rds"; if (!file.exists(g)) stop("data-prep finished but did not write ", g)
  inp <- readRDS(g)
  inp$branch <- branch; inp$pixel_res <- res; inp$pixel_intersect <- intersect
  inp$built <- format(Sys.time(), "%Y-%m-%d_%H%M"); inp$source_dat <- "data-prep (fresh from gridwork)"; inp$source_meta <- NA
  ilab <- if (nzchar(intersect)) intersect else "native"
  p <- sprintf("output/pixel_model_inputs_%s_%s_%dkm_%s.rds", branch, ilab, res, inp$built); saveRDS(inp, p); inp$path <- p
  message(">>> saved dated dump -> ", p, "  (a fresh dat_pixel_FULL_* was also written by the data-prep)")
  inp
}

if (isTRUE(BUILD_INPUTS)) {
  inp <- if (identical(BUILD_FROM, "dataprep")) run_dataprep(BRANCH) else reconstruct_inputs(BRANCH, save = TRUE)
} else if (is.null(INPUTS)) {
  inp <- reconstruct_inputs(BRANCH, save = FALSE)                         # reconstruct on the fly (not persisted)
} else if (identical(INPUTS, "auto")) {
  cand <- Sys.glob(sprintf("output/pixel_model_inputs_%s_%s_%dkm_*.rds", BRANCH, if (nzchar(PIXEL_INTERSECT)) PIXEL_INTERSECT else "native", PIXEL_RES))
  if (length(cand)) { f <- cand[which.max(file.mtime(cand))]; message(">>> inputs: newest dump  ", basename(f)); inp <- readRDS(f) }
  else { message(">>> no saved dump for ", BRANCH, " — reconstructing on the fly (set BUILD_INPUTS<-TRUE to persist a dated dump)"); inp <- reconstruct_inputs(BRANCH, save = FALSE) }
} else if (file.exists(INPUTS)) {
  message(">>> inputs: ", INPUTS); inp <- readRDS(INPUTS)
} else stop("INPUTS='", INPUTS, "' not found. Use \"auto\", a dump path, NULL, or set BUILD_INPUTS<-TRUE.")

# ---- unpack + optional subsample -------------------------------------------
X <- as.matrix(inp$X_mat); if (is.null(colnames(X)) && !is.null(inp$col_names)) colnames(X) <- inp$col_names; X[!is.finite(X)] <- 0
Y <- as.matrix(inp$Y_pixel); Y <- Y / rowSums(Y); cats <- colnames(Y)
group_idx <- as.integer(inp$group_idx_vec); coordX <- inp$coord_X; coordY <- inp$coord_Y
grp_lab <- if (!is.null(inp$re_group_names)) inp$re_group_names[group_idx] else as.character(group_idx)
inputs_src <- inp$path %||% (if (is.character(INPUTS)) INPUTS else "on-the-fly reconstruction")
if (SUBSAMPLE > 0 && SUBSAMPLE < nrow(X)) { set.seed(SEED); sel <- sort(sample(nrow(X), SUBSAMPLE))
  X <- X[sel,,drop=FALSE]; Y <- Y[sel,,drop=FALSE]; group_idx <- group_idx[sel]
  coordX <- coordX[sel]; coordY <- coordY[sel]; grp_lab <- grp_lab[sel] }
re_idx <- sort(intersect(c(which(colnames(X)=="intercept"), which(colnames(X) %in% socio)), seq_len(ncol(X))))

# ---- build + inspect the nesting tree (edit `tree` here to override) --------
tree <- build_nest_tree(BRANCH, cats)   # prints the tree + coverage; assign your own list to override
cat(sprintf("\n[%s] n=%d  cov=%d  classes=%d  RE-groups=%d  RE-covs=%d  use_re=%s | M=%d niter=%d iv_mode=%s\n",
            BRANCH, nrow(X), ncol(X), length(cats), length(unique(group_idx)), length(re_idx), USE_RE, M, NITER, IVMODE))

# ---- FIT --------------------------------------------------------------------
t0 <- Sys.time()
fit <- nested_cut_fit(X, Y, tree, use_iv = TRUE, use_re = USE_RE, group_idx = group_idx, re_idx = re_idx,
                      M = M, draws_per_impute = DRAWS_PER_IMPUTE, niter = NITER, nburn = round(NITER/3),
                      thin = THIN, min_pixels = MIN_PIXELS, iv_mode = IVMODE, moment_rank = MOMENT_RANK,
                      stream_disk = STREAM_DISK, n_chains = N_CHAINS, n_cores = N_CORES, store_dir = STORE_DIR)
cat(sprintf(">>> fit done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))

if (N_CHAINS > 1L && !is.null(fit$convergence)) {
  cat(sprintf("\n>>> convergence (%d chains) per node — Rhat / ESS:\n", N_CHAINS)); print(fit$convergence) }

# ---- lambda per nest + mode used -------------------------------------------
s <- summary_nested_cut(fit)
cat("\n>>> lambda (median [95% CI]) per internal nest:\n")
for (nm in names(s)) if (!is.null(s[[nm]]$lambda)) { L <- s[[nm]]$lambda
  for (r in rownames(L)) cat(sprintf("    %-26s %-16s %.3f  [%.3f, %.3f]\n", nm, r, L[r,"median"], L[r,"q025"], L[r,"q975"])) }
cat("\n>>> iv_mode used per node:\n"); print(nested_cut_modes(fit))

if (SAVE) { ts <- format(Sys.time(), "%Y-%m-%d_%H%M")
  saveRDS(list(fit=fit, tree=tree, cats=cats, re_idx=re_idx, branch=BRANCH,
               coordX=coordX, coordY=coordY, group=grp_lab, obs=Y,
               inputs_source=inputs_src, iv_mode=IVMODE, use_re=USE_RE, M=M, niter=NITER, n_chains=N_CHAINS),  # provenance
          f <- sprintf("output/nested_cut_%s_%s.rds", BRANCH, ts))
  cat(sprintf(">>> saved %s  (inputs: %s)\n", f, inputs_src)) }

# ============================ INSPECT (interactive) ========================= #
# Run these lines as needed after the fit.
#
# ## per-nest zero-sum coefficient matrices + lambda
#   str(s, max.level = 2)
#   s[["root"]]$coef
#
# ## predicted fine-class shares
#   pred_mean <- predict_nested_cut(fit, X, D = 200, summary = "mean")   # [n x J], O(n*J) memory (safe)
#   round(rbind(obs = colMeans(Y), pred = colMeans(pred_mean)), 4)
#   # full posterior draws materialise [n x J x D] in RAM (~2.5 GB at n=64k, D=200) —
#   # only for uncertainty, and prefer a pixel subset:  P <- predict_nested_cut(fit, X[idx,], D = 200)
#
# ## fine-class effective symmetric table (flat-comparable marginal effects, median + CI)
#   eff <- effective_symmetric_table(fit, X)
#   round(eff$median, 3)
#
# ## quick grid map of one predicted class (needs ggplot2)
#   library(ggplot2); j <- which(cats == "Cropland_HI")
#   ggplot(data.frame(x=coordX, y=coordY, share=pred_mean[,j]), aes(x,y,fill=share)) +
#     geom_tile() + coord_equal() + scale_fill_viridis_c() + theme_void() + ggtitle(cats[j])
