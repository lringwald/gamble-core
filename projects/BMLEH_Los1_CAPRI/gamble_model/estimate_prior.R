#!/usr/bin/env Rscript
# =============================================================================
# BMLEH_Los1_CAPRI — estimate the prior module
# =============================================================================
# FLAT multinomial logit on the BMLEH target classification. Moderate scope, NO BART:
# the linear-vs-BART gate on GLOBIOM found BART worth +249 nats held-out but at ~1.7x the fit
# time, concentrated in abundant terrain-driven classes -- not the right first cut for a target
# with 33 crop classes, many of them rare. Establish the linear baseline first.
#
#   Rscript projects/BMLEH_Los1_CAPRI/gamble_model/estimate_prior.R [--smoke]
#
# Needs a design dump carrying the BMLEH target classes. Produce one with the pixel driver:
#   DRIVER_DUMP_INPUTS=TRUE DRIVER_DUMP_EXIT=TRUE DRIVER_TARGET_RULES=<target_rules.R> \
#   DRIVER_TARGET_SHARES=<crop_shares_nuts2.csv> Rscript run_lu_pixel_model.R
#
# SEGMENTED SAMPLING. Every run writes final_state_chain_<id>.qs, the sampler's raw internal state
# (beta, mu, prec_beta, z_c, sigma_re, a_re, horseshoe, and the estimated scalars). Point BM_RESUME
# at a finished run to continue those chains instead of starting cold:
#
#   Rscript .../estimate_prior.R                                  # segment 1, burns in
#   BM_RESUME=results/.../prior_2026-09-10_1200 Rscript ...       # segment 2, continues
#
# A resumed segment defaults to NBURN 0 -- the chain is already burnt in, so discarding again would
# throw away good draws. Segments are separate directories; pool them for diagnostics rather than
# concatenating files. Keep the DESIGN fixed across segments: the state is indexed by covariate and
# class, so resuming onto a different dump is meaningless. The seed differs per segment by design,
# otherwise a resumed segment would replay the previous one's random numbers.
#
# Env: BM_INPUT (design dump)  BM_NITER 6000  BM_NBURN 2000  BM_THIN 8  BM_CHAINS 4
#      BM_NPIX 0 (all)  BM_TESTFRAC 0.2  BM_OUT results/gamble_model/BMLEH_Los1_CAPRI
#      BM_RESUME (previous run dir; enables hot start)  BM_SEG (segment label)
# =============================================================================
Sys.setenv(OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)
t0 <- Sys.time()
suppressMessages({library(Rcpp); library(RcppArmadillo); library(data.table)})
args  <- commandArgs(trailingOnly = TRUE); SMOKE <- "--smoke" %in% args
RESUME <- Sys.getenv("BM_RESUME", "")
# Score an existing run without sampling: the draws are already on disk, so a fit produced by an
# older version of this script (which scored chain 1 only) can be re-scored without a refit.
SCORE_ONLY <- isTRUE(as.logical(Sys.getenv("BM_SCORE_ONLY", "FALSE")))
# The RE slab scale. estimate_slab_c2=TRUE samples it, which is the default but is where this model's
# convergence fails: c2 and the RE precision are aliased (only tau_eff = tau_raw + 1/c2 reaches the
# likelihood), chains settle at different points on that ridge and stay. BM_SLAB_C2=FALSE fixes c2 at
# BM_SLAB_C2_VAL instead -- the control that decides whether the ridge is what holds Rhat up.
# RE-scale ASIS interweaving. The RE-variance block is the worst-mixing part of this sampler --
# measured sigma_re is still drifting upward after 5000 sweeps (p ~ 1e-40), which is what holds the
# joint Rhat at 1.7 rather than any mode or ridge. re_asis alternates the centred and non-centred
# parameterisation of the RE scale, which is the standard remedy for exactly that.
RE_ASIS <- isTRUE(as.logical(Sys.getenv("BM_RE_ASIS", "FALSE")))
SLAB_EST <- isTRUE(as.logical(Sys.getenv("BM_SLAB_C2", "TRUE")))
SLAB_VAL <- { .v <- Sys.getenv("BM_SLAB_C2_VAL", "100")
              if (identical(.v, "auto")) "auto" else as.numeric(.v) }
# A resumed chain is already burnt in; burning again discards good draws. Explicit BM_NBURN wins.
NITER <- as.integer(Sys.getenv("BM_NITER", if (SMOKE) "600"  else "6000"))
NBURN <- as.integer(Sys.getenv("BM_NBURN", if (nzchar(RESUME)) "0" else if (SMOKE) "200" else "2000"))
# thin 1 (was 8, matched to GLOBIOM 2026-09-17): thinning never improves ESS, it only discards
# draws, and storage is not the constraint (~60 MB per 500 draws x 4 chains).
THIN  <- as.integer(Sys.getenv("BM_THIN",  "1"))
NCH   <- as.integer(Sys.getenv("BM_CHAINS", if (SMOKE) "1" else "4"))
NPIX  <- as.integer(Sys.getenv("BM_NPIX",  if (SMOKE) "5000" else "0"))
TESTF <- as.numeric(Sys.getenv("BM_TESTFRAC", "0.2"))
NCORES<- as.integer(Sys.getenv("BM_CORES", as.character(max(1L, parallel::detectCores() - 2L))))
# Per-chain progress heartbeat, in seconds. 0 silences it. One throttled line per chain, because
# four txtProgressBars cannot share a terminal line.
PROG_SEC <- as.numeric(Sys.getenv("BM_PROGRESS_SEC", "30"))
run_tag <- Sys.getenv("RUN_ID", Sys.getenv("RUN_TAG", ""))
default_tag <- format(Sys.time(), "prior_%Y-%m-%d_%H%M")
if (nzchar(run_tag) && run_tag != "auto" && run_tag != "none") {
  default_tag <- paste0(default_tag, "_", run_tag)
}
OUT   <- Sys.getenv("BM_OUT", file.path("results/gamble_model/BMLEH_Los1_CAPRI", default_tag))
INPUT <- Sys.getenv("BM_INPUT", "output/designs/pixel_model_inputs_BMLEH_Los1_CAPRI.rds")
# SPLIT (ported from GLOBIOM 2026-09-17). Hybrid: half the hold-out is country-stratified random
# (interpolation), half is whole 150 km spatial blocks with a 1-cell buffer (transfer). Both are
# scored from the same fit, so the gap between them measures the optimism of a random split rather
# than assuming it. See codes/spatial_split.R for why the buffer, not the block size, does the work.
SPLIT_MODE  <- Sys.getenv("BM_SPLIT", "hybrid")               # hybrid | random | country | block
BLOCK_MULT  <- as.integer(Sys.getenv("BM_BLOCK_MULT", "5"))   # block side = MULT x 3 cells
BUFFER_CELL <- as.integer(Sys.getenv("BM_BUFFER", "1"))
BLOCK_SHARE <- as.numeric(Sys.getenv("BM_BLOCK_SHARE", "0.5"))

# Resolve the resume directory to the folder that actually holds final_state_chain_*.qs
RESUME_DIR <- ""
if (nzchar(RESUME)) {
  cand <- c(RESUME, file.path(RESUME, "posterior"))
  # resuming needs the chain STATE; scoring needs only the DRAWS -- a run killed mid-flight has
  # batches but no final state, and is still perfectly scorable
  want <- if (SCORE_ONLY) "^posterior_batch_[0-9]+_chain_.*\\.qs$" else "^final_state_chain_.*\\.qs$"
  hit  <- cand[vapply(cand, function(d) length(list.files(d, want)) > 0, logical(1))]
  if (!length(hit)) stop(sprintf("BM_RESUME has no %s under %s",
    if (SCORE_ONLY) "posterior_batch_*.qs to score" else "final_state_chain_*.qs to resume from", RESUME))
  RESUME_DIR <- hit[1]
}
if (!file.exists(INPUT)) stop("design dump not found: ", INPUT,
  "\n  build it with the pixel driver (see the header), or point BM_INPUT at an existing dump.")
for (f in c("codes/mnl_aux_func.R","codes/mnlogit_rcpp_sym.R","codes/spatial_split.R"))
  if (!file.exists(f)) stop("run from the gamble-core root; missing ", f)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source("codes/mnl_aux_func.R"); source("codes/mnlogit_rcpp_sym.R"); source("codes/spatial_split.R")

inp <- readRDS(INPUT)
X <- as.matrix(inp$X_mat); X[!is.finite(X)] <- 0
Y <- as.matrix(inp$Y_pixel); g <- inp$group_idx_vec
w <- if (!is.null(inp$weights_pixel)) inp$weights_pixel else rep(1, nrow(X))
# Coordinates ride along for the spatial split. The dump stores them aligned 1:1 with X_mat rows, so
# they MUST be subsampled with the same index -- subsetting rows but not coordinates would misalign
# the tiling against the data silently, and every block would hold out the wrong pixels.
cx <- inp$coord_X; cy <- inp$coord_Y
if (NPIX > 0 && NPIX < nrow(X)) { set.seed(1); s <- sort(sample(nrow(X), NPIX))
  X <- X[s,]; Y <- Y[s,]; g <- g[s]; w <- w[s]
  if (!is.null(cx)) { cx <- cx[s]; cy <- cy[s] } }
keep <- colSums(Y) > 0; Y <- Y[, keep, drop = FALSE]; Y <- Y / rowSums(Y); cats <- colnames(Y)

# ---- does this design actually belong to THIS project? ----------------------------------------
# A design dump is just an .rds; nothing in it says which class column built it. A dump produced
# before BMLEH_Los1_label existed was fitted for 12 hours without complaint, silently missing the
# energy class and the grassland split. Check the class set against the project's own markers, and
# the RE labels against the CAPRI keys the downscaler needs, before spending the compute.
.want_cls <- c("Cropland_permanent_energy", "Grassland_intensive", "Grassland_extensive",
               "Cropland_arable_other_crop", "Artificial")
.miss_cls <- setdiff(.want_cls, cats)
.capri    <- all(grepl("^[A-Z]{2}$", inp$re_group_names))
if (length(.miss_cls) || !.capri) {
  .msg <- c(
    if (length(.miss_cls)) sprintf("  class set is NOT the BMLEH target: missing %s
    (the dump predates BMLEH_Los1_label -- re-run `run.sh design`)",
                                   paste(.miss_cls, collapse = ", ")),
    if (!.capri) sprintf("  RE groups are NOT CAPRI codes: %s ...
    (re-run `run.sh design`, which sets DRIVER_RE_GROUP_COL=CAPRI_NUTS; these are not relabellable afterwards)",
                         paste(head(inp$re_group_names, 3), collapse = ", ")))
  if (isTRUE(as.logical(Sys.getenv("BM_ALLOW_STALE_DESIGN", "FALSE")))) {
    cat(sprintf("\n!! STALE DESIGN, proceeding because BM_ALLOW_STALE_DESIGN=TRUE:\n%s\n\n", paste(.msg, collapse = "\n")))
  } else {
    stop(sprintf("this design does not match the BMLEH_Los1_CAPRI project:\n%s\n  set BM_ALLOW_STALE_DESIGN=TRUE to fit it anyway (diagnostics only).",
                 paste(.msg, collapse = "\n")))
  }
}
sp <- make_spatial_split(cx, cy, g, test_frac = TESTF, mode = SPLIT_MODE,
                         block_mult = BLOCK_MULT, focal_reach = 1L,
                         buffer_cells = BUFFER_CELL, block_share = BLOCK_SHARE,
                         seed = 20260909L)
tr <- sp$train; te <- sp$test
icpt  <- which(colnames(X) == "intercept")
if (!length(icpt)) stop("design has no `intercept` column; the sampler and the horseshoe both need it.")
# RE covariates (matched to GLOBIOM 2026-09-17). A covariate gets a COUNTRY slope only if its EFFECT
# is plausibly mediated by national policy or development -- not merely because it varies in space.
# This was hardcoded to 5 and omitted GHM_TI and allPA_share: GHM_HI/GHM_TI are two threat groups of
# the SAME index, so shrinking one by country and pooling the other is arbitrary, and allPA_share
# binds through national ENFORCEMENT, which is the clearest policy-mediated effect in the design.
# Terrain, climate and soil stay POOLED (biophysical: the relationship is not nationally mediated).
.re_spec <- Sys.getenv("BM_RE_VARS", "log1p_GDP,log1p_Pop,GHM_HI,GHM_TI,CISI,allPA_share")
re_i  <- sort(unique(c(icpt, which(colnames(X) %in% trimws(strsplit(.re_spec, ",")[[1]])))))
bl    <- which.max(colSums(Y)); J <- ncol(Y)

cat(sprintf("\n%s\nBMLEH_Los1_CAPRI prior | FLAT MNL, no BART\n%s\n", strrep("=",70), strrep("=",70)))
cat(sprintf("  design %s\n  %d pixels | %d classes | %d covariates | %d RE groups\n",
    INPUT, nrow(X), J, ncol(X), length(unique(g))))
print_spatial_split(sp)
cat(sprintf("  %d sweeps (burn %d, thin %d) x %d chain(s) | baseline %s\n  out %s\n",
    NITER, NBURN, THIN, NCH, cats[bl], OUT))
cat(sprintf("  RE on %d covariate(s): %s\n", length(re_i), paste(colnames(X)[re_i], collapse = ", ")))
cat(sprintf("  RE-scale ASIS: %s | re_support_strength: 1\n", if (RE_ASIS) "ON" else "off"))
cat(sprintf("  RE slab: %s\n", if (SLAB_EST) "c2 SAMPLED (estimate_slab_c2=TRUE)"
            else sprintf("c2 FIXED at %s (control for the aliasing ridge)", format(SLAB_VAL))))
if (nzchar(RESUME_DIR)) {
  .have <- length(list.files(RESUME_DIR, "^final_state_chain_.*\\.qs$"))
  cat(sprintf("  HOT START from %s (%d chain state(s))\n", RESUME_DIR, .have))
  if (.have < NCH) cat(sprintf("  NOTE: only %d of %d chains have a saved state; the rest start cold\n", .have, NCH))
} else cat("  cold start\n")
cat("\n")
# rare classes are the point of this target -- report them before fitting, not after
pv <- sort(colMeans(Y)); cat(sprintf("  rarest 5 classes: %s\n", paste(sprintf("%s %.4f%%", names(pv)[1:min(5,J)],
    100*pv[1:min(5,J)]), collapse=" | ")))
cat(sprintf("  classes below 0.1%% of area: %d of %d\n\n", sum(pv < 0.001), J))

# Segment lineage has to be known BEFORE fit_one() runs: SEG_ID offsets the per-chain seed, and a
# seed that is computed after the fit would leave every segment replaying the same random stream.
SEG_ID <- as.integer(Sys.getenv("BM_SEG", if (nzchar(RESUME_DIR)) "0" else "1"))
if (SEG_ID == 0L) {  # derive from the parent's lineage when not given
  pl <- file.path(dirname(RESUME_DIR), "lineage.txt")
  if (!file.exists(pl)) pl <- file.path(RESUME_DIR, "lineage.txt")
  SEG_ID <- 2L
  if (file.exists(pl)) {
    .pn <- suppressWarnings(as.integer(sub("^segment\\s+", "",
             grep("^segment", readLines(pl), value = TRUE)[1])))
    if (!is.na(.pn)) SEG_ID <- .pn + 1L
  }
}

fit_one <- function(ci, dpath) {
  # Segment-dependent seed: reusing the cold seed would make a resumed segment replay the previous
  # one's random numbers, which looks like sampling but adds no information.
  set.seed(99 + ci + 1000L * SEG_ID)
  ist <- NULL
  if (nzchar(RESUME_DIR)) {
    f <- file.path(RESUME_DIR, sprintf("final_state_chain_%d.qs", ci))
    if (file.exists(f)) ist <- qs2::qs_read(f)
  }
  mnlogit_rcpp_sym(
    init_state = ist,
    X = X[tr,], Y = Y[tr,], intercept = FALSE, baseline = bl,
    symmetric = TRUE, symmetric_hs = TRUE,
    niter = NITER, nburn = NBURN, thin = as.integer(THIN), y_weight = w[tr],
    use_re = TRUE, group_idx = g[tr], re_idx = re_i,
    use_horseshoe = TRUE, horseshoe_idx = setdiff(seq_len(ncol(X)), icpt),
    fe_support_strength = 0,           # (n/PR)^2 on a live HS kernel costs ~1274 nats
    # re_support_strength = 1 EXPLICITLY (2026-09-17). This used to be omitted, so it fell through
    # support_prior_strength to the sampler default 0 -- i.e. the RE sparse-group gate was OFF here
    # while GLOBIOM ran it at 1. That was an inheritance, not a decision; see sec. 8.1b.
    re_support_strength = 1,
    re_regularize = TRUE, estimate_slab_c2 = SLAB_EST, collapse_slab_c2 = SLAB_VAL,
    re_asis = RE_ASIS,
    use_bart = FALSE,                  # deliberate: linear baseline first
    save_posterior_to_disk = TRUE, disk_path = dpath, chain_id = ci, calc_loo = FALSE,
    # One throttled line per chain (see make_progress_cb in codes/mnl_aux_func.R). Without it the
    # sampler cat()s every iteration and future_lapply swallows the lot until the worker exits.
    progress_cb = make_progress_cb(sprintf("chain %d", ci), PROG_SEC, NITER))
}
if (SCORE_ONLY) {
  if (!nzchar(RESUME_DIR)) stop("BM_SCORE_ONLY needs BM_RESUME pointing at the run to score")
  OUT <- dirname(RESUME_DIR)
  cat(sprintf("  SCORE ONLY: no sampling, pooling what is already in %s\n\n", OUT))
}
dp <- file.path(OUT, "posterior"); dir.create(dp, recursive = TRUE, showWarnings = FALSE)
# Remove the output directory again if the run dies before writing anything. An empty husk is not
# harmless: it sorts newest by name, so `run.sh more` would pick it over the finished segment beside
# it and fail with "no final_state_chain_*.qs to resume from".
.clean_husk <- function() {
  if (dir.exists(dp) && !length(list.files(dp, "\\.qs$"))) unlink(OUT, recursive = TRUE)
}
on.exit(.clean_husk(), add = TRUE)
if (!SCORE_ONLY) writeLines(c(sprintf("segment %d", SEG_ID), sprintf("design  %s", INPUT),
             sprintf("parent  %s", if (nzchar(RESUME_DIR)) RESUME_DIR else "(cold start)"),
             sprintf("sweeps  %d (burn %d, thin %d) x %d chains", NITER, NBURN, THIN, NCH),
             sprintf("built   %s", format(Sys.time()))),
           file.path(OUT, "lineage.txt"))

# Record the SETTINGS, not only the lineage. lineage.txt says which design and how many sweeps; it
# says nothing about whether ASIS was on, whether c2 was estimated, or which structural switches the
# sampler ran with -- and those are exactly what a reader has to know before quoting a number from
# this run. Written machine-readable so the report renders it, rather than being reconstructed from
# the shell history of whoever happened to launch it.
if (!SCORE_ONLY) saveRDS(list(
  segment = SEG_ID, design = INPUT,
  parent  = if (nzchar(RESUME_DIR)) RESUME_DIR else NA_character_,
  built   = Sys.time(), sampler = "mnlogit_rcpp_sym",
  mcmc = list(niter = NITER, nburn = NBURN, thin = THIN, chains = NCH, cores = NCORES),
  # The SPLIT is part of what produced any number quoted from this run: a random-split McFadden and a
  # block-split McFadden are not the same quantity. Realised geometry, not the requested one.
  split = sp$meta,
  data = list(pixels_total = nrow(X), pixels_train = length(tr), test_frac = TESTF,
              subsample_cap = NPIX, classes = J, covariates = ncol(X),
              re_groups = length(unique(g)), baseline = cats[bl],
              re_idx = colnames(X)[re_i], group_col = inp$re_group_col),
  # every non-default argument the sampler is actually called with, including the ones that are
  # hardcoded here rather than env-driven -- a switch being constant is not a reason to hide it.
  switches = list(symmetric = TRUE, symmetric_hs = TRUE, intercept = FALSE,
                  use_re = TRUE, re_regularize = TRUE, re_asis = RE_ASIS,
                  use_horseshoe = TRUE, fe_support_strength = 0,
                  estimate_slab_c2 = SLAB_EST, collapse_slab_c2 = SLAB_VAL,
                  use_bart = FALSE, calc_loo = FALSE)),
  file.path(OUT, "run_config.rds"))
if (SCORE_ONLY) {
  cat("  (skipping the sampler)\n")
} else if (NCH > 1L && NCORES > 1L) {
  suppressMessages(library(future.apply)); future::plan(future::multisession, workers = min(NCH, NCORES))
  # future.stdout = NA: do NOT capture worker output. The default buffers everything a worker prints
  # and relays it only on completion, so per-chain progress arrives hours late, all at once.
  invisible(future.apply::future_lapply(seq_len(NCH), function(ci) {
    source("codes/mnl_aux_func.R"); source("codes/mnlogit_rcpp_sym.R"); fit_one(ci, dp); NULL },
    future.seed = TRUE, future.stdout = NA)); future::plan(future::sequential)
} else if (!SCORE_ONLY) for (ci in seq_len(NCH)) invisible(fit_one(ci, dp))

# Score on EVERY chain and EVERY segment, not chain 1 of the current run. Two reasons: with 4 chains
# a single-chain score discards three quarters of the draws, and under segmented sampling the
# posterior IS the union of the segments -- scoring only the newest one would silently throw away
# all the earlier sampling the resume exists to preserve.
.ancestors <- function(d) {          # walk parent links recorded in lineage.txt
  out <- character(0)
  repeat {
    lf <- file.path(d, "lineage.txt"); if (!file.exists(lf)) break
    pl <- grep("^parent", readLines(lf), value = TRUE); if (!length(pl)) break
    par <- trimws(sub("^parent\\s+", "", pl[1])); if (!nzchar(par) || par == "(cold start)") break
    pdir <- if (basename(par) == "posterior") dirname(par) else par
    if (!dir.exists(pdir) || pdir %in% out) break
    out <- c(out, pdir); d <- pdir
  }
  file.path(out, "posterior")
}
# BM_SCORE_POOL: "lineage" (default) pools this segment and every ancestor -- the right thing when
# the segments are one continued chain. "self" scores THIS segment only, for when the earlier ones are
# being treated as burn-in and discarded.
.pool <- Sys.getenv("BM_SCORE_POOL", "lineage")
.dirs <- if (identical(.pool, "self")) dp else unique(c(dp, .ancestors(OUT)))
.dirs <- .dirs[dir.exists(.dirs)]
.acc <- NULL; .nd <- 0L; .used <- 0L
for (.d in .dirs) {
  .chs <- sort(unique(as.integer(sub(".*_chain_(\\d+)\\.qs$", "\\1",
            list.files(.d, "^posterior_batch_.*_chain_\\d+\\.qs$")))))
  for (.ci in .chs) {
    .fi <- tryCatch(recover_mnlogit_posterior(.d, chain_id = .ci), error = function(e) NULL)
    if (is.null(.fi) || is.null(.fi$postb_total)) next
    .acc <- if (is.null(.acc)) apply(.fi$postb_total, c(1,2,3), sum)
            else .acc + apply(.fi$postb_total, c(1,2,3), sum)
    .nd <- .nd + dim(.fi$postb_total)[4]; .used <- .used + 1L
  }
}
if (is.null(.acc)) stop("no recoverable posterior batches under: ", paste(.dirs, collapse = ", "))
cat(sprintf("\nscoring on %d chain-segment(s) across %d segment(s), %d draws pooled\n",
            .used, length(.dirs), .nd))
B <- .acc / .nd; gl <- unique(g[tr])
# Score an arbitrary index set (ported from GLOBIOM). The null is computed WITHIN the set, so a block
# score is measured against the composition of the blocks themselves, not of the whole map.
.score_idx <- function(idx, label) {
  if (!length(idx)) return(NULL)
  gi <- match(g[idx], gl); ok <- which(!is.na(gi))
  if (!length(ok)) {
    cat(sprintf("  %-26s UNSCORABLE: every group is unseen in training\n", label)); return(NULL)
  }
  if (length(ok) < length(idx))
    cat(sprintf("  %-26s note: %d of %d rows dropped (group unseen in training)\n",
                label, length(idx) - length(ok), length(idx)))
  ii <- idx[ok]; gg <- gi[ok]; Xi <- X[ii, , drop = FALSE]
  U <- t(vapply(seq_along(ii), function(k) as.numeric(Xi[k,] %*% B[,,gg[k]]), numeric(J)))
  U <- U - apply(U,1,max); P <- exp(U); P <- P/rowSums(P)
  Yt <- Y[ii,,drop=FALSE]
  ll <- sum(Yt*log(pmax(P,1e-12)))
  nl <- sum(Yt*log(pmax(matrix(colMeans(Yt), nrow(Yt), J, byrow=TRUE),1e-12)))
  list(label=label, idx=ii, P=P, Yt=Yt, ll=ll, nl=nl, mcfadden=1-ll/nl, n=nrow(Yt))
}
cat("\n")
S_all <- .score_idx(te, "ALL held-out")
S_rnd <- .score_idx(sp$test_random, "random (interpolation)")
S_blk <- .score_idx(sp$test_block,  "block (spatial transfer)")
cat(sprintf("\n%-26s %10s %12s %10s\n", "HELD-OUT", "n", "log-lik", "McFadden"))
for (S in Filter(Negate(is.null), list(S_all, S_rnd, S_blk)))
  cat(sprintf("%-26s %10d %12.1f %10.4f\n", S$label, S$n, S$ll, S$mcfadden))
# The gap is the quantity of interest: how much a random split flatters this model.
if (!is.null(S_rnd) && !is.null(S_blk))
  cat(sprintf("\n  OPTIMISM of the random split: %+.4f McFadden (%.4f random vs %.4f block)\n",
              S_rnd$mcfadden - S_blk$mcfadden, S_rnd$mcfadden, S_blk$mcfadden))
P <- S_all$P; Yt <- S_all$Yt; rows <- S_all$idx
per <- data.table(class = cats, obs = colMeans(Yt), pred = colMeans(P), dLL = colSums(Yt*log(pmax(P,1e-12))))
cat("\nworst-fitting classes (largest negative contribution per unit area):\n")
print(head(per[obs > 0][order(dLL/pmax(obs,1e-9))][, .(class, obs=round(obs,5), pred=round(pred,5), dLL=round(dLL,1))], 6))
saveRDS(list(P=P, rows=rows, per_class=per, cats=cats, ll=S_all$ll, mcfadden=S_all$mcfadden,
             score_all=S_all[c("label","ll","nl","mcfadden","n")],
             score_random=if (is.null(S_rnd)) NULL else S_rnd[c("label","ll","nl","mcfadden","n")],
             score_block =if (is.null(S_blk)) NULL else S_blk[c("label","ll","nl","mcfadden","n")],
             optimism_mcfadden = if (!is.null(S_rnd) && !is.null(S_blk)) S_rnd$mcfadden - S_blk$mcfadden else NA_real_,
             split=sp$meta,
             niter=NITER, nburn=NBURN, thin=THIN, nchains=NCH, input=INPUT, posterior_dir=dp,
             scored_dirs=.dirs, draws_pooled=.nd, segment=SEG_ID),
        file.path(OUT, "prior_fit.rds"))
fwrite(per, file.path(OUT, "per_class_fit.csv"))
cat(sprintf("\nsaved %s\n", file.path(OUT, "prior_fit.rds")))
cat(sprintf("diagnose with: Rscript diagnose_posterior.R %s\n", dp))
cat(sprintf("elapsed %.1f h\n", as.numeric(difftime(Sys.time(), t0, units="hours"))))
