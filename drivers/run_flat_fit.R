#!/usr/bin/env Rscript
# =============================================================================
# drivers/run_flat_fit.R — fit the FLAT multinomial pixel land-use model
# =============================================================================
# Standard flat multinomial logit on any design dump (GLOBIOM, BMLEH, etc.).
# Features:
#   - Direct fitting from an assembled design dump (pixel_model_inputs*.rds)
#   - Hybrid spatial hold-out evaluation (codes/spatial_split.R: random + spatial block)
#   - Modern validated sampler configuration (symmetric HS, ASIS, full Bayes c2)
#   - Real-time multicore unbuffered progress logging (console + heartbeat files)
#   - Kill-safe resumable segmented sampling (final_state_chain_*.qs)
#   - In-pipeline scoring (held-out log-likelihood, McFadden R2, optimism)
#   - Automated post-estimation convergence diagnostics (diagnose_posterior.R)
# =============================================================================

Sys.setenv(OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)
t0 <- Sys.time()

suppressMessages({
  library(data.table)
  library(Rcpp)
  library(RcppArmadillo)
  library(qs2)
  library(future)
  library(future.apply)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
SMOKE <- "--smoke" %in% args || identical(Sys.getenv("RUN_MODE"), "smoke")

# ---- 1. Resolve Design Dump --------------------------------------------------
INPUT <- Sys.getenv("DESIGN_PATH", Sys.getenv("BM_INPUT", Sys.getenv("GB_INPUT", "")))
if (!nzchar(INPUT) || identical(INPUT, "auto") || !file.exists(INPUT)) {
  # Candidate search
  cands <- c(
    list.files("output/designs", pattern = "^pixel_model_inputs.*\\.rds$", full.names = TRUE),
    list.files("output", pattern = "^pixel_model_inputs.*\\.rds$", full.names = TRUE),
    if (dir.exists("/mnt/wdrv/gamble-core/output/designs"))
      list.files("/mnt/wdrv/gamble-core/output/designs", pattern = "^pixel_model_inputs.*\\.rds$", full.names = TRUE)
  )
  cands <- unique(cands[file.exists(cands)])
  if (length(cands) > 0) {
    # Pick newest
    cands <- cands[order(file.mtime(cands), decreasing = TRUE)]
    INPUT <- cands[1]
    cat(sprintf(">>> Auto-detected design dump: %s\n", INPUT))
  } else {
    stop("No design dump found. Provide DESIGN_PATH or run TASK=flat_design first.")
  }
}

if (!file.exists(INPUT)) stop("Design dump not found: ", INPUT)

for (f in c("codes/mnl_aux_func.R", "codes/mnlogit_rcpp_sym.R", "codes/spatial_split.R")) {
  if (!file.exists(f)) stop("Run from gamble-core repo root; missing ", f)
}
source("codes/mnl_aux_func.R")
source("codes/mnlogit_rcpp_sym.R")
source("codes/spatial_split.R")

# ---- 2. MCMC & Model Settings ------------------------------------------------
.get_env_val <- function(keys, default = "") {
  for (k in keys) {
    v <- Sys.getenv(k, "")
    if (nzchar(v) && !identical(v, "auto") && !identical(v, "default")) return(v)
  }
  default
}
.get_env_int <- function(keys, default) {
  v <- .get_env_val(keys)
  if (nzchar(v)) {
    iv <- suppressWarnings(as.integer(v))
    if (!is.na(iv)) return(iv)
  }
  as.integer(default)
}
.get_env_num <- function(keys, default) {
  v <- .get_env_val(keys)
  if (nzchar(v)) {
    nv <- suppressWarnings(as.numeric(v))
    if (!is.na(nv)) return(nv)
  }
  as.numeric(default)
}

NITER  <- .get_env_int(c("NITER", "DRIVER_NITER", "GB_NITER", "BM_NITER"), if (SMOKE) 600L else 6000L)
NBURN  <- .get_env_int(c("NBURN", "DRIVER_NBURN", "GB_NBURN", "BM_NBURN"), if (SMOKE) 200L else 2000L)
THIN   <- .get_env_int(c("THIN", "DRIVER_THIN", "GB_THIN", "BM_THIN"), 1L)
NCH    <- .get_env_int(c("N_CHAINS", "DRIVER_NCHAINS", "GB_CHAINS", "BM_CHAINS"), if (SMOKE) 1L else 4L)
NPIX   <- .get_env_int(c("SUBSAMPLE", "GB_NPIX", "BM_NPIX"), if (SMOKE) 5000L else 0L)
TESTF  <- .get_env_num(c("TESTFRAC", "GB_TESTFRAC", "BM_TESTFRAC"), 0.2)
NCORES <- .get_env_int(c("N_CORES", "DRIVER_NCORES", "GB_CORES", "BM_CORES"), max(1L, min(NCH, parallel::detectCores())))

RE_ASIS     <- isTRUE(as.logical(.get_env_val(c("RE_ASIS", "DRIVER_RE_ASIS", "GB_RE_ASIS"), "TRUE")))
RE_SUPPORT  <- .get_env_num(c("RE_SUPPORT", "GB_RE_SUPPORT"), 1.0)
SLAB_EST    <- isTRUE(as.logical(.get_env_val(c("SLAB_C2", "GB_SLAB_C2"), "TRUE")))
SLAB_VAL    <- {
  .v <- Sys.getenv("SLAB_C2_VAL", Sys.getenv("GB_SLAB_C2_VAL", "4"))
  if (identical(.v, "auto") || !nzchar(.v)) "auto" else as.numeric(.v)
}
USE_BART    <- isTRUE(as.logical(.get_env_val(c("USE_BART", "DRIVER_USE_BART", "GB_BART"), "FALSE")))

# Output path
RUN_TAG <- Sys.getenv("RUN_ID", Sys.getenv("RUN_TAG", ""))
default_tag <- format(Sys.time(), "flat_%Y-%m-%d_%H%M")
if (nzchar(RUN_TAG) && RUN_TAG != "auto" && RUN_TAG != "none") {
  default_tag <- paste0(default_tag, "_", RUN_TAG)
}

OUT_BASE <- Sys.getenv("OUTPUT_DIR", "results/gamble_model")
if (identical(OUT_BASE, "output")) OUT_BASE <- "output/results"
OUT <- Sys.getenv("GB_OUT", Sys.getenv("BM_OUT", file.path(OUT_BASE, default_tag)))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Resume directory
RESUME <- Sys.getenv("RESUME", Sys.getenv("GB_RESUME", Sys.getenv("BM_RESUME", "")))
RESUME_DIR <- ""
if (nzchar(RESUME)) {
  cand <- c(RESUME, file.path(RESUME, "posterior"))
  hit <- cand[vapply(cand, function(d) length(list.files(d, "^final_state_chain_.*\\.qs$")) > 0, logical(1))]
  if (length(hit)) {
    RESUME_DIR <- hit[1]
    if (NBURN == 2000L && !nzchar(Sys.getenv("NBURN"))) NBURN <- 0L
  }
}

SCORE_ONLY <- isTRUE(as.logical(.get_env_val(c("SCORE_ONLY", "GB_SCORE_ONLY"), "FALSE")))
PROG_SEC   <- .get_env_num(c("PROGRESS_SEC", "GB_PROGRESS_SEC"), 30.0)

# ---- 3. Load Design Dump -----------------------------------------------------
cat(sprintf(">>> Loading design dump: %s\n", INPUT))
inp <- readRDS(INPUT)
X <- as.matrix(inp$X_mat); X[!is.finite(X)] <- 0
Y <- as.matrix(inp$Y_pixel); g <- inp$group_idx_vec
w <- if (!is.null(inp$weights_pixel)) inp$weights_pixel else rep(1, nrow(X))

cx <- inp$coord_X; cy <- inp$coord_Y

# Subsample if requested
if (NPIX > 0 && NPIX < nrow(X)) {
  set.seed(1)
  s <- sort(sample(nrow(X), NPIX))
  X <- X[s, , drop = FALSE]
  Y <- Y[s, , drop = FALSE]
  g <- g[s]
  w <- w[s]
  if (!is.null(cx)) { cx <- cx[s]; cy <- cy[s] }
}

keep <- colSums(Y) > 0
Y <- Y[, keep, drop = FALSE]
Y <- Y / rowSums(Y)
cats <- colnames(Y)
J <- ncol(Y)
bl <- which.max(colSums(Y))

# Intercept & RE variables
icpt <- which(colnames(X) == "intercept")
if (!length(icpt)) stop("Design has no `intercept` column; the sampler needs it.")

re_mode <- Sys.getenv("RE_BLOCK", "intercept+socio")
if (identical(re_mode, "intercept")) {
  re_i <- icpt
} else {
  .re_spec <- Sys.getenv("DRIVER_RE_VARS", Sys.getenv("GB_RE_VARS",
                "log1p_GDP,log1p_Pop,GHM_HI,GHM_TI,CISI,allPA_share"))
  re_candidates <- trimws(strsplit(.re_spec, ",")[[1]])
  re_i <- sort(unique(c(icpt, which(colnames(X) %in% re_candidates))))
}

# ---- 4. Train / Held-out Spatial Split ----------------------------------------
SPLIT_MODE  <- Sys.getenv("SPLIT", Sys.getenv("GB_SPLIT", "hybrid"))
BLOCK_MULT  <- as.integer(Sys.getenv("BLOCK_MULT", "5"))
BUFFER_CELL <- as.integer(Sys.getenv("BUFFER", "1"))
BLOCK_SHARE <- as.numeric(Sys.getenv("BLOCK_SHARE", "0.5"))

sp <- make_spatial_split(cx, cy, g, test_frac = TESTF, mode = SPLIT_MODE,
                         block_mult = BLOCK_MULT, focal_reach = 1L,
                         buffer_cells = BUFFER_CELL, block_share = BLOCK_SHARE,
                         seed = 20260916L)
tr <- sp$train; te <- sp$test

cat(sprintf("\n%s\nFLAT MNL MODEL | BART %s\n%s\n",
            strrep("=", 70), if (USE_BART) "ON" else "off", strrep("=", 70)))
cat(sprintf("  design %s\n  %d pixels | %d classes | %d covariates | %d RE groups\n",
            INPUT, nrow(X), J, ncol(X), length(unique(g))))
print_spatial_split(sp)
cat(sprintf("  %d sweeps (burn %d, thin %d) x %d chain(s) | baseline %s\n  out %s\n",
            NITER, NBURN, THIN, NCH, cats[bl], OUT))
cat(sprintf("  RE on %d covariate(s): %s\n", length(re_i), paste(colnames(X)[re_i], collapse = ", ")))
cat(sprintf("  RE-scale ASIS: %s | re_support_strength: %g\n", if (RE_ASIS) "ON" else "off", RE_SUPPORT))
cat(sprintf("  RE slab: %s\n", if (SLAB_EST) "c2 SAMPLED (standard)"
            else sprintf("c2 FIXED at %s", format(SLAB_VAL))))

if (nzchar(RESUME_DIR)) {
  .have <- length(list.files(RESUME_DIR, "^final_state_chain_.*\\.qs$"))
  cat(sprintf("  HOT START from %s (%d chain state(s))\n", RESUME_DIR, .have))
} else {
  cat("  cold start\n")
}

pv <- sort(colMeans(Y))
cat(sprintf("\n  rarest 5 classes: %s\n", paste(sprintf("%s %.4f%%", names(pv)[1:min(5, J)],
    100 * pv[1:min(5, J)]), collapse = " | ")))
cat(sprintf("  classes below 0.1%% of area: %d of %d\n\n", sum(pv < 0.001), J))

# ---- 5. Segment Lineage & State Setup ----------------------------------------
SEG_ID <- as.integer(Sys.getenv("SEG_ID", if (nzchar(RESUME_DIR)) "0" else "1"))
if (SEG_ID == 0L) {
  pl <- file.path(dirname(RESUME_DIR), "lineage.txt")
  if (!file.exists(pl)) pl <- file.path(RESUME_DIR, "lineage.txt")
  SEG_ID <- 2L
  if (file.exists(pl)) {
    .pn <- suppressWarnings(as.integer(sub("^segment\\s+", "",
             grep("^segment", readLines(pl), value = TRUE)[1])))
    if (!is.na(.pn)) SEG_ID <- .pn + 1L
  }
}

dp <- file.path(OUT, "posterior")
dir.create(dp, recursive = TRUE, showWarnings = FALSE)

.clean_husk <- function() {
  if (dir.exists(dp) && !length(list.files(dp, "\\.qs$"))) unlink(OUT, recursive = TRUE)
}
on.exit(.clean_husk(), add = TRUE)

if (!SCORE_ONLY) {
  writeLines(c(sprintf("segment %d", SEG_ID), sprintf("design  %s", INPUT),
               sprintf("parent  %s", if (nzchar(RESUME_DIR)) RESUME_DIR else "(cold start)"),
               sprintf("sweeps  %d (burn %d, thin %d) x %d chains", NITER, NBURN, THIN, NCH),
               sprintf("built   %s", format(Sys.time()))),
             file.path(OUT, "lineage.txt"))

  saveRDS(list(
    segment = SEG_ID, design = INPUT,
    parent  = if (nzchar(RESUME_DIR)) RESUME_DIR else NA_character_,
    built   = Sys.time(), sampler = "mnlogit_rcpp_sym",
    mcmc = list(niter = NITER, nburn = NBURN, thin = THIN, chains = NCH, cores = NCORES),
    data = list(pixels_total = nrow(X), pixels_train = length(tr), test_frac = TESTF,
                subsample_cap = NPIX, classes = J, covariates = ncol(X),
                re_groups = length(unique(g)), baseline = cats[bl],
                re_idx = colnames(X)[re_i], group_col = inp$re_group_col),
    split = sp$meta,
    switches = list(symmetric = TRUE, symmetric_hs = TRUE, intercept = FALSE,
                    use_re = TRUE, re_regularize = TRUE, re_asis = RE_ASIS,
                    use_horseshoe = TRUE, hs_no_intercept = TRUE,
                    fe_support_strength = 0, re_support_strength = RE_SUPPORT,
                    estimate_slab_c2 = SLAB_EST, collapse_slab_c2 = SLAB_VAL,
                    use_bart = USE_BART, calc_loo = FALSE)),
    file.path(OUT, "run_config.rds"))
}

# ---- 6. Sampling Function ----------------------------------------------------
fit_one <- function(ci, dpath) {
  set.seed(99 + ci + 1000L * SEG_ID)
  ist <- NULL
  if (nzchar(RESUME_DIR)) {
    f <- file.path(RESUME_DIR, sprintf("final_state_chain_%d.qs", ci))
    if (file.exists(f)) ist <- qs2::qs_read(f)
  }
  mnlogit_rcpp_sym(
    init_state = ist,
    X = X[tr, , drop = FALSE], Y = Y[tr, , drop = FALSE],
    intercept = FALSE, baseline = bl,
    symmetric = TRUE, symmetric_hs = TRUE,
    niter = NITER, nburn = NBURN, thin = as.integer(THIN), y_weight = w[tr],
    use_re = TRUE, group_idx = g[tr], re_idx = re_i,
    use_horseshoe = TRUE, horseshoe_idx = setdiff(seq_len(ncol(X)), icpt),
    fe_support_strength = 0, re_support_strength = RE_SUPPORT,
    re_regularize = TRUE, estimate_slab_c2 = SLAB_EST, collapse_slab_c2 = SLAB_VAL,
    re_asis = RE_ASIS,
    use_bart = USE_BART,
    save_posterior_to_disk = TRUE, disk_path = dpath, chain_id = ci, calc_loo = FALSE,
    progress_cb = make_progress_cb(sprintf("chain %d", ci), PROG_SEC, NITER)
  )
}

# ---- 7. Execute Sampling -----------------------------------------------------
if (SCORE_ONLY) {
  cat("  (skipping sampling, scoring existing runs)\n")
} else if (NCH > 1L && NCORES > 1L) {
  plan_type <- if (.Platform$OS.type != "windows" && future::supportsMulticore()) future::multicore else future::multisession
  future::plan(plan_type, workers = min(NCH, NCORES))
  cat(sprintf(">>> Parallel plan: %s with %d workers\n",
              if (.Platform$OS.type != "windows" && future::supportsMulticore()) "multicore (shared unbuffered stderr)" else "multisession",
              min(NCH, NCORES)))
  invisible(future.apply::future_lapply(seq_len(NCH), function(ci) {
    source("codes/mnl_aux_func.R"); source("codes/mnlogit_rcpp_sym.R")
    fit_one(ci, dp); NULL
  }, future.seed = TRUE, future.stdout = NA))
  future::plan(future::sequential)
} else {
  for (ci in seq_len(NCH)) invisible(fit_one(ci, dp))
}

# ---- 8. Recover & Score Held-Out Set -----------------------------------------
.ancestors <- function(d) {
  out <- character(0)
  repeat {
    lf <- file.path(d, "lineage.txt"); if (!file.exists(lf)) break
    pl <- grep("^parent", readLines(lf), value = TRUE); if (!length(pl)) break
    par <- trimws(sub("^parent\\s+", "", pl[1]))
    if (!nzchar(par) || par == "(cold start)") break
    pdir <- if (basename(par) == "posterior") dirname(par) else par
    if (!dir.exists(pdir) || pdir %in% out) break
    out <- c(out, pdir); d <- pdir
  }
  file.path(out, "posterior")
}

.pool <- Sys.getenv("GB_SCORE_POOL", "lineage")
.dirs <- if (identical(.pool, "self")) dp else unique(c(dp, .ancestors(OUT)))
.dirs <- .dirs[dir.exists(.dirs)]

.acc <- NULL; .nd <- 0L; .used <- 0L
for (.d in .dirs) {
  .chs <- sort(unique(as.integer(sub(".*_chain_(\\d+)\\.qs$", "\\1",
            list.files(.d, "^posterior_batch_.*_chain_\\d+\\.qs$")))))
  for (.ci in .chs) {
    .fi <- tryCatch(recover_mnlogit_posterior(.d, chain_id = .ci), error = function(e) NULL)
    if (is.null(.fi) || is.null(.fi$postb_total)) next
    .acc <- if (is.null(.acc)) apply(.fi$postb_total, c(1, 2, 3), sum)
            else .acc + apply(.fi$postb_total, c(1, 2, 3), sum)
    .nd <- .nd + dim(.fi$postb_total)[4]; .used <- .used + 1L
  }
}

if (!is.null(.acc)) {
  cat(sprintf("\nScoring on %d chain-segment(s) across %d segment(s), %d draws pooled\n",
              .used, length(.dirs), .nd))
  B <- .acc / .nd
  gl <- unique(g[tr])

  .score_idx <- function(idx, label) {
    if (!length(idx)) return(NULL)
    gi <- match(g[idx], gl); ok <- which(!is.na(gi))
    if (!length(ok)) return(NULL)
    ii <- idx[ok]; gg <- gi[ok]; Xi <- X[ii, , drop = FALSE]
    U <- t(vapply(seq_along(ii), function(k) as.numeric(Xi[k, ] %*% B[, , gg[k]]), numeric(J)))
    U <- U - apply(U, 1, max); P <- exp(U); P <- P / rowSums(P)
    Yt <- Y[ii, , drop = FALSE]
    ll <- sum(Yt * log(pmax(P, 1e-12)))
    nl <- sum(Yt * log(pmax(matrix(colMeans(Yt), nrow(Yt), J, byrow = TRUE), 1e-12)))
    list(label = label, idx = ii, P = P, Yt = Yt, ll = ll, nl = nl, mcfadden = 1 - ll / nl, n = nrow(Yt))
  }

  S_all <- .score_idx(te, "ALL held-out")
  S_rnd <- .score_idx(sp$test_random, "random (interpolation)")
  S_blk <- .score_idx(sp$test_block,  "block (spatial transfer)")

  cat(sprintf("\n%-26s %10s %12s %10s\n", "HELD-OUT", "n", "log-lik", "McFadden"))
  for (S in Filter(Negate(is.null), list(S_all, S_rnd, S_blk))) {
    cat(sprintf("%-26s %10d %12.1f %10.4f\n", S$label, S$n, S$ll, S$mcfadden))
  }

  if (!is.null(S_rnd) && !is.null(S_blk)) {
    cat(sprintf("\n  OPTIMISM of the random split: %+.4f McFadden (%.4f random vs %.4f block)\n",
                S_rnd$mcfadden - S_blk$mcfadden, S_rnd$mcfadden, S_blk$mcfadden))
  }

  if (!is.null(S_all)) {
    P <- S_all$P; Yt <- S_all$Yt
    per <- data.table(class = cats, obs = colMeans(Yt), pred = colMeans(P),
                      dLL = colSums(Yt * log(pmax(P, 1e-12))))
    saveRDS(list(P = P, rows = S_all$idx, per_class = per, cats = cats, ll = S_all$ll,
                 mcfadden = S_all$mcfadden,
                 score_all = S_all[c("label", "ll", "nl", "mcfadden", "n")],
                 score_random = if (is.null(S_rnd)) NULL else S_rnd[c("label", "ll", "nl", "mcfadden", "n")],
                 score_block  = if (is.null(S_blk)) NULL else S_blk[c("label", "ll", "nl", "mcfadden", "n")],
                 optimism_mcfadden = if (!is.null(S_rnd) && !is.null(S_blk)) S_rnd$mcfadden - S_blk$mcfadden else NA_real_,
                 split = sp$meta,
                 niter = NITER, nburn = NBURN, thin = THIN, nchains = NCH, input = INPUT, posterior_dir = dp,
                 scored_dirs = .dirs, draws_pooled = .nd, segment = SEG_ID),
            file.path(OUT, "prior_fit.rds"))
    fwrite(per, file.path(OUT, "per_class_fit.csv"))
    cat(sprintf("\nSaved fit summary -> %s\n", file.path(OUT, "prior_fit.rds")))
  }
}

# ---- 9. Convergence Diagnostics ---------------------------------------------
if (isTRUE(as.logical(Sys.getenv("DIAGNOSE", "TRUE"))) && file.exists("postprocess/diagnose_posterior.R")) {
  cat(sprintf("\n%s\nPOSTERIOR CONVERGENCE DIAGNOSTICS\n%s\n", strrep("=", 70), strrep("=", 70)))
  tryCatch({
    system2("Rscript", c("postprocess/diagnose_posterior.R", dp))
  }, error = function(e) {
    cat(sprintf("WARNING: diagnose_posterior.R failed: %s\n", conditionMessage(e)))
  })
}

cat(sprintf("\n>>> Flat model fitting finished successfully in %.2f hours\n",
            as.numeric(difftime(Sys.time(), t0, units = "hours"))))
