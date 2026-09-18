#!/usr/bin/env Rscript
# ============================================================================
# Full-fledged multi-chain convergence check on the REAL GLOBIOM pixel design.
#   - 64,178 pixels x 68 covariates x 26 country REs (subsample via NPIX)
#   - N chains from overdispersed inits -> split-Rhat + bulk/tail ESS
#   - Arms: production baseline, and re_mean_shift ON (the RE-mixing claim)
# Env: NPIX (default 20000), NCHAIN (4), NITER (2000), NBURN (800), ARMS
# ============================================================================
suppressMessages({library(Rcpp); library(RcppArmadillo); library(posterior); library(future.apply)})
sys.source("codes/mnlogit_rcpp_sym.R", environment())

NPIX   <- as.integer(Sys.getenv("NPIX",   "20000"))
NCHAIN <- as.integer(Sys.getenv("NCHAIN", "4"))
NITER  <- as.integer(Sys.getenv("NITER",  "2000"))
NBURN  <- as.integer(Sys.getenv("NBURN",  "800"))
ARMS   <- strsplit(Sys.getenv("ARMS", "baseline,mean_shift"), ",")[[1]]
# symmetric_hs is a SEPARATE argument from symmetric: `symmetric` sets the zero-sum coding of the
# response, `symmetric_hs` selects the symmetric (zero-sum) SHRINKAGE prior and, via re_prec_sym,
# the _sym RE-variance updater -- the only one that accepts dof_mean_pinned. Passing only
# `symmetric = TRUE` silently measures the DIAGONAL horseshoe.
SYM_HS <- isTRUE(as.logical(Sys.getenv("SYM_HS", "TRUE")))
THIN    <- as.integer(Sys.getenv("THIN", "2"))       # storage thinning; postb_total is k*p*G*draws
WORKERS <- as.integer(Sys.getenv("WORKERS", "4"))    # raise on a many-core box, or run arms concurrently
TAG     <- Sys.getenv("TAG", "")                     # label so concurrent runs are distinguishable

inp <- readRDS("output/designs/pixel_model_inputs_GLOBIOM_init2000.rds")
set.seed(1)
idx <- if (NPIX < nrow(inp$X_mat)) sort(sample(nrow(inp$X_mat), NPIX)) else seq_len(nrow(inp$X_mat))
X <- as.matrix(inp$X_mat)[idx, , drop = FALSE]
Y <- as.matrix(inp$Y_pixel)[idx, , drop = FALSE]
g <- inp$group_idx_vec[idx]
w <- if (!is.null(inp$weights_pixel)) inp$weights_pixel[idx] else NULL
# drop classes with no mass in the subsample (they carry no information and break the baseline pick)
keepY <- colSums(Y) > 0; Y <- Y[, keepY, drop = FALSE]; Y <- Y / rowSums(Y)
base_i <- which.max(colSums(Y))
re_idx <- which(colnames(X) %in% c("intercept", "log1p_GDP", "log1p_Pop", "GHM_HI", "CISI"))

cat(sprintf("\n%sREAL design: %d pixels x %d cov, %d classes, %d groups | %d chains x %d iter (burn %d, thin %d)\n",
            if (nzchar(TAG)) paste0("[", TAG, "] ") else "",
            nrow(X), ncol(X), ncol(Y), length(unique(g)), NCHAIN, NITER, NBURN, THIN))
cat(sprintf("RE covariates: %s\n", paste(colnames(X)[re_idx], collapse = ", ")))
cat(sprintf("shrinkage: symmetric_hs = %s  (=> re_prec_sym = %s, the _sym RE-variance updater)\n\n",
            SYM_HS, SYM_HS))

plan(multisession, workers = min(NCHAIN, WORKERS))
fit_arm <- function(arm) {
  future_lapply(seq_len(NCHAIN), function(ch) {
    # sourceCpp pointers do NOT survive into multisession workers ("NULL value passed as
    # symbol address") -- the worker must compile its own. Same as run_lu_pixel_model.R:1732.
    source("codes/mnl_aux_func.R")
    source("codes/mnlogit_rcpp_sym.R")
    set.seed(1000 + ch)
    mnlogit_rcpp_sym(
      X = X, Y = Y, intercept = FALSE, baseline = base_i, symmetric = TRUE,
      symmetric_hs = SYM_HS,
      niter = NITER, nburn = NBURN, thin = as.integer(THIN), y_weight = w,
      use_re = TRUE, group_idx = g, re_idx = re_idx,
      use_horseshoe = TRUE, re_regularize = TRUE, estimate_slab_c2 = TRUE,
      re_mean_shift     = arm %in% c("mean_shift", "mean_shift_nodof"),
      # isolation arm: the shift WITHOUT the degrees-of-freedom change
      re_mean_shift_dof = if (identical(arm, "mean_shift_nodof")) FALSE else NULL,
      chain_id = ch, calc_loo = FALSE, save_posterior_to_disk = FALSE)
  }, future.seed = TRUE)
}

# ---- extractors: each returns an (iter x param) matrix for ONE chain -------------------
# Checking only mu and sigma_re is not enough. Every defect found in this sampler lived in the
# SHRINKAGE AUXILIARIES (frozen HS kernel, fe_support inflation, kappa/sigma aliasing), and those
# corrupt tau/lambda/c2 long before mu looks wrong. log_lik is the standard global sentinel.
as_mat <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.vector(x)) return(matrix(x, ncol = 1))
  if (length(dim(x)) == 2L) return(t(x))                       # param x draw -> draw x param
  if (length(dim(x)) == 3L) return(matrix(aperm(x, c(3,1,2)), nrow = dim(x)[3]))
  if (length(dim(x)) == 4L) return(matrix(aperm(x, c(4,1,2,3)), nrow = dim(x)[4]))
  NULL
}
# postb_total is k x p x G x draws -> subsample cells so the table stays readable
sub_cols <- function(M, cap = 400) {
  if (is.null(M) || ncol(M) <= cap) return(M)
  set.seed(7); M[, sort(sample(ncol(M), cap)), drop = FALSE]
}
BLOCKS <- list(
  "log_lik"    = function(f) as_mat(f$post_log_lik),           # global sentinel
  "mu (FE)"    = function(f) as_mat(f$postb_pooled),           # estimand
  "b_g (RE)"   = function(f) sub_cols(as_mat(f$postb_total)),  # the RE deviations themselves
  "sigma_re"   = function(f) as_mat(f$post_sigma_re),          # estimand
  "lambda (HS)"= function(f) sub_cols(as_mat(f$horseshoe$post_kappa_pooled)), # NB: nested under $horseshoe
  "hs_c2"      = function(f) as_mat(f$horseshoe$post_c2),      # horseshoe slab
  "re_tau"     = function(f) as_mat(f$post_re_tau),            # global RE shrinkage scale
  "slab_c2"    = function(f) as_mat(f$post_slab_c2),           # estimated slab (full Bayes)
  "kappa"      = function(f) as_mat(f$post_kappa)              # joint gate (NULL unless enabled)
)

summarise_arm <- function(chains, label) {
  cat(sprintf("[%s]\n", label))
  cat(sprintf("  %-12s %5s | %-28s | %-22s\n", "block", "n", "Rhat med/max  >1.01  >1.05", "ESS_bulk med/min"))
  out <- list()
  for (nm in names(BLOCKS)) {
    lst <- lapply(chains, BLOCKS[[nm]])
    if (any(vapply(lst, is.null, TRUE))) {
      # Say WHY. "not produced" read as a sampler gap once when it was a wrong extractor path.
      why <- switch(nm, "re_tau" = "re_hs_global = FALSE (global RE horseshoe off)",
                        "kappa"  = "joint_fe_re_shrink = FALSE",
                        "absent from the fit object -- check the extractor path")
      cat(sprintf("  %-12s   -- %s\n", nm, why)); next }
    n <- min(vapply(lst, nrow, 1L)); pdim <- ncol(lst[[1]])
    if (n < 4L || pdim < 1L) { cat(sprintf("  %-12s   -- too few draws\n", nm)); next }
    A <- array(unlist(lapply(lst, function(m) m[seq_len(n), , drop = FALSE])), dim = c(n, pdim, length(lst)))
    d <- posterior::as_draws_array(aperm(A, c(1, 3, 2)))
    su <- posterior::summarise_draws(d, "rhat", "ess_bulk")
    su <- su[is.finite(su$rhat) & is.finite(su$ess_bulk), ]
    if (!nrow(su)) { cat(sprintf("  %-12s   -- degenerate (constant)\n", nm)); next }
    cat(sprintf("  %-12s %5d | %.3f / %.3f  %5.1f%% %5.1f%% | %8.0f / %6.0f\n",
        nm, nrow(su), median(su$rhat), max(su$rhat),
        100*mean(su$rhat > 1.01), 100*mean(su$rhat > 1.05),
        median(su$ess_bulk), min(su$ess_bulk)))
    out[[nm]] <- su
  }
  invisible(out)
}

res <- list()
for (arm in ARMS) {
  t0 <- Sys.time()
  ch <- fit_arm(arm)
  cat(sprintf("\n=== arm '%s' (%.1f min) ===\n", arm, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  res[[arm]] <- summarise_arm(ch, arm)
}
if (length(ARMS) > 1) {
  cat("\n=== mixing comparison (the re_mean_shift claim), per block ===\n")
  for (blk in names(BLOCKS)) {
    a <- res[[ARMS[1]]][[blk]]; b <- res[[ARMS[2]]][[blk]]
    if (is.null(a) || is.null(b)) next
    cat(sprintf("  %-12s ESS_bulk med %6.0f -> %6.0f (%+5.0f%%) | max Rhat %.3f -> %.3f\n",
      blk, median(a$ess_bulk), median(b$ess_bulk),
      100*(median(b$ess_bulk)/median(a$ess_bulk) - 1), max(a$rhat), max(b$rhat)))
  }
}
cat("\nCONVERGENCE RUN DONE\n")
