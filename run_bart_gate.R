#!/usr/bin/env Rscript
# =============================================================================
# run_bart_gate.R — full-sample GLOBIOM flat MNL, WITH and WITHOUT BART on topography
# =============================================================================
# Fits both arms on the same data, split and seeds, then scores them per class on a common
# held-out set. Arms differ in ONE thing: whether the six topographic covariates
# (slope, elevation, aspect cos/sin, lon, lat) enter as linear terms or as a tree ensemble.
#
#   Rscript run_bart_gate.R
#
# Resumable: each arm is written as soon as it finishes and is SKIPPED on a rerun, so a crash
# or a kill costs you at most the arm in flight. Re-running after both arms exist just redoes
# the (cheap) scoring.
#
# Env knobs (all optional):
#   BG_NITER   8000   total sweeps          BG_NBURN 2000   burn-in
#   BG_THIN       8   storage thinning      BG_CHAINS   1   chains per arm
#   BG_NPIX       0   0 = all 64,178        BG_TESTFRAC 0.2 held-out fraction
#   BG_SPLIT random   random | country      BG_ARMS  linear,bart
#   BG_OUT  output/bart_gate                BG_INPUT output/pixel_model_inputs_GLOBIOM_init2000.rds
#
# RUNTIME, measured at 20k pixels x 1000 sweeps (31.5 min linear, 42.6 min BART, 1 chain):
#   full sample, 1 chain per arm, both arms —  4000 sweeps ~16 h | 8000 ~32 h | 16000 ~63 h
# The default 8000 leaves 6000 retained draws. That is ample for a HELD-OUT PREDICTION
# comparison, which is what this script measures. It is NOT the 16000/4000 regime validated for
# RE-VARIANCE inference (median RE-var ESS ~173 there vs ~85 here) — raise BG_NITER if you intend
# to read sigma_v or its intervals off this fit.
# =============================================================================
Sys.setenv(OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)
t_start <- Sys.time()
`%||%` <- function(a,b) if (!is.null(a)) a else b
suppressMessages({library(Rcpp); library(RcppArmadillo)})

NITER <- as.integer(Sys.getenv("BG_NITER","8000")); NBURN <- as.integer(Sys.getenv("BG_NBURN","2000"))
THIN  <- as.integer(Sys.getenv("BG_THIN","8"));     NCH   <- as.integer(Sys.getenv("BG_CHAINS","1"))
NPIX  <- as.integer(Sys.getenv("BG_NPIX","0"));     TESTF <- as.numeric(Sys.getenv("BG_TESTFRAC","0.2"))
SPLIT <- Sys.getenv("BG_SPLIT","random");           ARMS  <- trimws(strsplit(Sys.getenv("BG_ARMS","linear,bart"),",")[[1]])
OUT   <- Sys.getenv("BG_OUT","output/bart_gate")
INPUT <- Sys.getenv("BG_INPUT","output/pixel_model_inputs_GLOBIOM_init2000.rds")

# ---- preflight: fail NOW, not 6 hours in ------------------------------------------------------
if (!file.exists(INPUT)) stop("input not found: ", INPUT,
  "\n  produce it with: DRIVER_DUMP_INPUTS=TRUE DRIVER_DUMP_EXIT=TRUE Rscript run_lu_pixel_model.R")
for (p in c("sf","dbarts","qs2")) if (!requireNamespace(p, quietly=TRUE)) stop("missing package: ", p)
for (f in c("codes/mnl_aux_func.R","codes/mnlogit_rcpp_sym.R")) if (!file.exists(f)) stop("run from the repo root; missing ", f)
if (NBURN >= NITER) stop("BG_NBURN must be < BG_NITER")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
source("codes/mnl_aux_func.R"); source("codes/mnlogit_rcpp_sym.R")

# ---- design ------------------------------------------------------------------------------------
inp <- readRDS(INPUT)
X <- as.matrix(inp$X_mat); X[!is.finite(X)] <- 0
pts <- sf::st_transform(sf::st_as_sf(data.frame(x=inp$coord_X, y=inp$coord_Y), coords=c("x","y"), crs=3035), 4326)
cc <- sf::st_coordinates(pts)
X <- cbind(X, lon=cc[,1], lat=cc[,2])          # RAW degrees: BART splits on raw inputs
Y <- as.matrix(inp$Y_pixel); g <- inp$group_idx_vec; w <- inp$weights_pixel
if (NPIX > 0 && NPIX < nrow(X)) { set.seed(1); s <- sort(sample(nrow(X), NPIX))
  X <- X[s,]; Y <- Y[s,]; g <- g[s]; w <- w[s] }
keep <- colSums(Y) > 0; Y <- Y[, keep]; Y <- Y/rowSums(Y); cats <- colnames(Y); J <- ncol(Y)

set.seed(20260904)
te <- if (identical(SPLIT,"country")) {
  # leave-country-out: the honest spatial test. Random CV is ~20x optimistic here because the
  # country RE carries most of the log-score skill; both arms get the same handicap either way.
  hold <- sample(unique(g), max(1, round(length(unique(g)) * TESTF)))
  which(g %in% hold)
} else sort(sample(nrow(X), round(nrow(X) * TESTF)))
tr <- setdiff(seq_len(nrow(X)), te)

icpt <- which(colnames(X) == "intercept")
topo <- c("Slope_rad","Elevation","Aspect_cos_mean","Aspect_sin_mean","lon","lat")
tcol <- which(colnames(X) %in% topo)
if (length(tcol) != length(topo)) stop("missing topo column(s): ", paste(setdiff(topo, colnames(X)), collapse=", "))
re_i <- which(colnames(X) %in% c("intercept","log1p_GDP","log1p_Pop","GHM_HI","CISI"))
bl   <- which.max(colSums(Y))

eta_h <- (31.5*(nrow(X)/20000)*(NITER/1000)/60) + (42.6*(nrow(X)/20000)*(NITER/1000)/60)
cat(sprintf("\n%s\nGLOBIOM flat MNL — BART gate\n%s\n", strrep("=",70), strrep("=",70)))
cat(sprintf("  pixels %d (train %d / test %d, %s split)\n", nrow(X), length(tr), length(te), SPLIT))
cat(sprintf("  classes %d | covariates %d | RE groups %d | RE covs %d\n", J, ncol(X), length(unique(g)), length(re_i)))
cat(sprintf("  BART covariates: %s\n", paste(colnames(X)[tcol], collapse=", ")))
cat(sprintf("  %d sweeps (burn %d, thin %d) x %d chain(s) per arm\n", NITER, NBURN, THIN, NCH))
cat(sprintf("  arms: %s | out: %s\n", paste(ARMS, collapse=", "), OUT))
cat(sprintf("  ROUGH ETA for both arms: %.1f h (started %s)\n\n", eta_h, format(t_start, "%H:%M")))

fit_arm <- function(use_b, chain_id = 1L, dpath = NULL) {
  set.seed(99 + chain_id)          # BG_CHAINS was read and PRINTED but never used, so every run
  mnlogit_rcpp_sym(                # was single-chain regardless -- and no Rhat was possible.
    chain_id = chain_id,
    # Stream the FULL posterior rather than a curated summary: any diagnostic can then be computed
    # later without refitting, and it reuses the tested recover_mnlogit_posterior() path instead of
    # a bespoke extraction that might quietly omit a block. ~0.3 GB per chain at these settings.
    save_posterior_to_disk = !is.null(dpath), disk_path = dpath,
    X = X[tr,], Y = Y[tr,], intercept = FALSE, baseline = bl,
    symmetric = TRUE, symmetric_hs = TRUE,
    niter = NITER, nburn = NBURN, thin = as.integer(THIN), y_weight = w[tr],
    use_re = TRUE, group_idx = g[tr], re_idx = re_i,
    use_horseshoe = TRUE, horseshoe_idx = setdiff(seq_len(ncol(X)), icpt),
    fe_support_strength = 0,                      # (n/PR)^2 on a live HS kernel costs ~1274 nats
    re_regularize = TRUE, estimate_slab_c2 = TRUE,
    use_bart = use_b,
    bart_idx    = if (use_b) tcol else NULL,
    linear_idx  = if (use_b) setdiff(seq_len(ncol(X)), tcol) else seq_len(ncol(X)),
    n_trees_bart = 50, bart_symmetric = TRUE, bart_base = 0.90, bart_power = 3.0, bart_k = 2.0,
    bart_k_prevalence = TRUE,                     # per-class leaf shrinkage; see docs section 3b
    store_bart_trees = TRUE, do_slim_trees = TRUE,
    calc_loo = FALSE)
}
score <- function(f, use_b) {
  lc <- if (use_b) setdiff(seq_len(ncol(X)), tcol) else seq_len(ncol(X))
  B  <- apply(f$postb_total, c(1,2,3), mean)
  gl <- unique(g[tr]); gi <- match(g[te], gl); rows <- which(!is.na(gi))
  if (!length(rows)) stop("no test row belongs to a fitted group — with BG_SPLIT=country the RE ",
                          "cannot transfer; score the pooled mean instead or use the random split.")
  Xl <- X[te, lc, drop=FALSE]
  U <- t(vapply(rows, function(i) as.numeric(Xl[i,] %*% B[,,gi[i]]), numeric(J)))
  if (use_b) {
    .sym <- if (!is.null(f$bart_symmetric)) f$bart_symmetric else f$bart$symmetric
    bm <- list(symmetric = isTRUE(.sym), p_all = J,
               pp = f$bart$pp %||% setdiff(seq_len(J), f$baseline))
    ts <- f$tree_store[seq(1, length(f$tree_store), by = max(1, length(f$tree_store) %/% 60))]
    U <- U + reconstruct_bart_f_mean(ts, X[te,,drop=FALSE][rows, tcol, drop=FALSE], bm)
  }
  U <- U - apply(U,1,max); P <- exp(U); list(P = P/rowSums(P), rows = rows)
}

for (nm in ARMS) {
  fp <- file.path(OUT, sprintf("arm_%s.rds", nm))
  if (file.exists(fp)) { cat(sprintf(">>> %-6s already present, skipping (delete %s to refit)\n", nm, fp)); next }
  ub <- identical(nm, "bart"); t0 <- Sys.time()
  dp <- file.path(OUT, sprintf("posterior_%s", nm)); dir.create(dp, recursive=TRUE, showWarnings=FALSE)
  cat(sprintf(">>> %-6s fitting (%d chain%s), posterior -> %s\n", nm, NCH, if (NCH>1) "s" else "", dp)); flush.console()
  for (ci in seq_len(NCH)) invisible(fit_arm(ub, ci, dp))
  # streamed runs return no postb_total in RAM -- read chain 1 back for scoring, which also
  # exercises the recovery path the later diagnostics will use
  f <- recover_mnlogit_posterior(dp, chain_id = 1L); sc <- score(f, ub)
  # Per-chain posterior traces for EVERY structural block, so convergence can be checked without
  # refitting. Full fit objects are far too large to keep at this size; these are the summaries
  # each block is actually judged on.
  saveRDS(list(P = sc$P, rows = sc$rows, arm = nm, minutes = as.numeric(difftime(Sys.time(), t0, units="mins")),
               niter = NITER, nburn = NBURN, thin = THIN, nchains = NCH, cats = cats, te = te,
               posterior_dir = dp), fp)
  cat(sprintf(">>> %-6s done in %.1f h -> %s\n", nm, as.numeric(difftime(Sys.time(), t0, units="hours")), fp))
  rm(f); gc()
}

# ---- scoring ------------------------------------------------------------------------------------
have <- ARMS[file.exists(file.path(OUT, sprintf("arm_%s.rds", ARMS)))]
if (length(have) < 2) { cat(sprintf("\nonly %d arm(s) present; rerun to complete the comparison.\n", length(have))); quit(save="no") }
A <- lapply(have, function(nm) readRDS(file.path(OUT, sprintf("arm_%s.rds", nm)))); names(A) <- have
rows <- Reduce(intersect, lapply(A, `[[`, "rows")); Yt <- Y[te,][match(rows, A[[1]]$rows),,drop=FALSE]
LL <- lapply(A, function(a) colSums(Yt * log(pmax(a$P[match(rows, a$rows),,drop=FALSE], 1e-12))))
obs <- colMeans(Yt); rare <- obs < 0.02
nullLL <- sum(Yt * log(pmax(matrix(obs, nrow(Yt), J, byrow=TRUE), 1e-12)))
cat(sprintf("\n%s\nHELD-OUT (%d pixels, %d classes)\n%s\n", strrep("=",70), nrow(Yt), J, strrep("=",70)))
cat(sprintf("%-8s %14s %10s %12s %12s\n", "arm", "log-lik", "McFadden",
            sprintf("rare(n=%d)", sum(rare)), sprintf("common(n=%d)", sum(!rare))))
for (nm in have) cat(sprintf("%-8s %14.1f %10.4f %12.1f %12.1f\n", nm, sum(LL[[nm]]),
    1 - sum(LL[[nm]])/nullLL, sum(LL[[nm]][rare]), sum(LL[[nm]][!rare])))
if (all(c("linear","bart") %in% have)) {
  d <- LL$bart - LL$linear
  cat(sprintf("\nBART - linear: total %+.1f | rare(n=%d) %+.1f | common(n=%d) %+.1f | improved %d/%d classes\n",
      sum(d), sum(rare), sum(d[rare]), sum(!rare), sum(d[!rare]), sum(d>0), J))
  cat(sprintf("\n%-24s %9s %12s %11s %11s\n","class","obs","dLL","pred_lin","pred_bart"))
  Pl <- A$linear$P[match(rows, A$linear$rows),,drop=FALSE]; Pb <- A$bart$P[match(rows, A$bart$rows),,drop=FALSE]
  for (j in order(-d)) cat(sprintf("%-24s %9.4f %12.1f %11.4f %11.4f\n", cats[j], obs[j], d[j], mean(Pl[,j]), mean(Pb[,j])))
  saveRDS(list(LL=LL, obs=obs, cats=cats, d=d, n_test=nrow(Yt)), file.path(OUT,"comparison.rds"))
  cat(sprintf("\nsaved %s\n", file.path(OUT,"comparison.rds")))
}
cat(sprintf("\ntotal elapsed %.1f h\n", as.numeric(difftime(Sys.time(), t_start, units="hours"))))
