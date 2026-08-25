# =============================================================================
# mnlogit_rcpp_sym.R  —  Adaptive sampler with Phases 1, 2, and 4
#
# Based on mnlogit_rcpp_adaptive.R with the following additions:
#   Tier 3 Symmetric Horseshoe
#   Tier 2 Symmetric Pooled RE Variances
#
# Usage:
#   Rcpp::sourceCpp("mnlogit_gibbs_core_sym.cpp")   # one-time compile
#   source("mnlogit_rcpp_sym.R")
#   fit <- mnlogit_rcpp_sym(X, Y, ...)
# =============================================================================

library(progressr)
library(MASS)
library(TruncatedNormal)
library(dbarts)
library(loo)
library(matrixStats)
library(Matrix)
library(pg)
library(BayesLogit)
library(Rcpp)
library(qs2)
library(RcppArmadillo)

aux_path <- "codes/mnl_aux_func.R"
if (!file.exists(aux_path)) aux_path <- "mnl_aux_func.R"
source(aux_path)

# --- Compile C++ (idempotent — skips if already loaded) ---
# The exists() guard alone can find stale function objects from a previous
# R session whose compiled pointers are now NULL.  We also do a test call
# to detect that case and force recompilation.
.needs_cpp_compile <- TRUE
if (exists("chol_sample_precision_cpp", mode = "function")) {
  .needs_cpp_compile <- tryCatch(
    {
      # tiny test: 2x2 identity precision, zero rhs → should return a numeric(2)
      chol_sample_precision_cpp(diag(2), rep(0, 2))
      FALSE # call succeeded → no recompile needed
    },
    error = function(e) TRUE
  ) # stale pointer → recompile
}
if (.needs_cpp_compile) {
  cpp_path <- "mnlogit_gibbs_core_sym.cpp"
  if (!file.exists(cpp_path)) cpp_path <- "codes/mnlogit_gibbs_core_sym.cpp"

  # Parallel-safe compilation: Use a lock file to ensure only one process compiles at a time
  lock_file <- paste0(cpp_path, ".lock")

  # Simple retry loop with backoff
  max_retries <- 120 # Wait up to 2 minutes
  retry_count <- 0

  while (file.exists(lock_file) && retry_count < max_retries) {
    Sys.sleep(1 + runif(1, 0, 1)) # Random sleep 1-2s
    retry_count <- retry_count + 1
  }

  if (retry_count >= max_retries) {
    warning("Compilation lock timeout for ", cpp_path, ". Proceeding anyway...")
  }

  # Create lock
  cat(Sys.getpid(), file = lock_file)

  # Execute compilation
  tryCatch({
    Rcpp::sourceCpp(cpp_path)
  }, finally = {
    # Always remove lock
    if (file.exists(lock_file)) {
      pid_in_lock <- tryCatch(as.integer(readLines(lock_file, warn = FALSE, n = 1)), error = function(e) 0)
      if (pid_in_lock == Sys.getpid()) {
        unlink(lock_file)
      }
    }
  })
}
rm(.needs_cpp_compile)



# =============================================================================
# Helper: Draw from N(P^{-1} Pb, P^{-1}) via Cholesky decomposition
# Pure-R fallback for the CAR joint z_c draw (Phase 4, Section C.4.5)
# =============================================================================
chol_sample_precision <- function(P, Pb, d) {
  L <- chol(P) # upper Cholesky: P = L' L
  mu <- backsolve(L, forwardsolve(t(L), Pb))
  z <- rnorm(d)
  mu + backsolve(L, z)
}

# =============================================================================
# Horseshoe helpers (shared by setup / update / storage to avoid duplication)
# =============================================================================
# Regularized (Finnish) horseshoe effective variance: c2 * v / (c2 + v).
.reg_hs_var <- function(v, c2) c2 * v / (c2 + v)
# lambda^2 * tau^2 in the right shape: [k_hs x p] when equation-specific, else a
# length-k_hs vector. (lambda is per-predictor; tau is per-equation or scalar.)
.hs_lam2_tau2 <- function(equation_specific, hs_tau2, hs_lambda2, hs_idx, k_hs, p) {
  tau2_mat <- if (equation_specific) matrix(hs_tau2, k_hs, p, byrow = TRUE) else hs_tau2
  lam2_sub <- if (equation_specific) hs_lambda2[hs_idx, , drop = FALSE] else hs_lambda2[hs_idx]
  lam2_sub * tau2_mat
}

# =============================================================================
# Symmetric HS Pooled Joint Draw
# =============================================================================
draw_beta_symhs_pooled <- function(X, Xt, kappa_w, omega, c_j_mat, prior_P,
                                   prior_Pb, pp, f_bart, use_bart, c_v, Msym, block_sym = NULL, c_block = NULL) {
  k <- nrow(Xt); p <- ncol(omega); K <- k * p
  P_joint  <- matrix(0, K, K)
  Pb_joint <- numeric(K)
  for (ip in seq_len(p)) {                          # data + base prior: block-diagonal by equation
    j <- pp[ip]; om <- omega[, ip]
    A <- prior_P + weighted_crossprod(Xt, X, om)
    tgt <- kappa_w[, j] + om * c_j_mat[, ip]
    if (use_bart) tgt <- tgt - om * f_bart[, ip]
    rows <- ((ip - 1) * k + 1):(ip * k)
    P_joint[rows, rows] <- A
    Pb_joint[rows] <- prior_Pb[, ip] + Xt %*% tgt
  }
  for (v in which(c_v > 0)) {                       # symmetric HS: couples equations at predictor v
    idx <- v + (seq_len(p) - 1) * k
    P_joint[idx, idx] <- P_joint[idx, idx] + c_v[v] * Msym
  }
  if (!is.null(block_sym)) {
      # V3: the block term must be category-symmetric TOO. `Mb` alone is applied within each
      # equation (implicitly Mb %x% I_p), whose quadratic form swings 2.7x with the arbitrary
      # baseline choice; kron(Msym, Mb) is invariant to 7 s.f. and equals the zero-sum penalty
      # exactly (scratchpad/kron_check.R). Index order is column-fastest within equation, which
      # is what kronecker(Msym, M) expects.
    for (bs_id in seq_along(block_sym)) {
      bs <- block_sym[[bs_id]]
      if (!bs$hs_on || c_block[bs_id] <= 0) next
      ix <- as.vector(vapply(seq_len(p), function(ip) (ip - 1) * k + bs$ret, numeric(length(bs$ret))))
      P_joint[ix, ix] <- P_joint[ix, ix] + c_block[bs_id] * kronecker(Msym, bs$M)
    }
  }
  matrix(chol_sample_precision_cpp(P_joint, Pb_joint), k, p)
}



# =============================================================================
# Main sampler — Rcpp-accelerated (Symmetric Adaptive)
# =============================================================================
mnlogit_rcpp_sym <- function(X, Y, intercept = FALSE, baseline = ncol(Y),
                                  niter = 1000, nburn = 500, thin = 1L, A0 = 2.0,
                                  empirical_intercept_prior = FALSE,
                                  y_weight = NULL, use_re = FALSE, group_idx = NULL,
                                  use_bart = FALSE, bart_idx = 1:ncol(X), n_trees_bart = 50, n_threads_bart = 1,
                                  store_bart_trees = TRUE, do_slim_trees = TRUE,
                                  save_bart_to_disk = FALSE,
                                  save_posterior_to_disk = FALSE,
                                  disk_path = tempdir(),
                                  bart_batch_size = 50,
                                  posterior_batch_size = 50,
                                  bart_warmup = 0,
                                  bart_tempering = FALSE,
                                  use_tempering = FALSE, tempering_T0 = 0.05,
                                  bart_base = 0.95,
                                  bart_power = 2.0,
                                  bart_k = 2.0,                # leaf-prior shrinkage (dbarts default 2; larger = tighter f, tames near-separation tails)
                                  bart_symmetric = FALSE,
                                  re_cor_threshold = 0.85,
                                  re_cor_min_groups = 0.5,
                                  store_f = FALSE,
                                  symmetric = FALSE,
                                  # ---- ALTERNATIVE-SPECIFIC term (conditional-logit style) ----
                                  # alt_spec_Z: n x p_all matrix whose column j is an attribute OF
                                  # ALTERNATIVE j (e.g. the spatial Y lag = neighbourhood share of
                                  # class j). Enters every utility with ONE SHARED coefficient delta:
                                  #     V_ij = x_i'beta_j + delta * z_ij
                                  # WHY IT EXISTS. With only CASE-specific covariates (one x_i per
                                  # pixel, identical across alternatives) a nested logit's lambda is
                                  # identified solely by the curvature of the log-sum-exp -- the
                                  # textbook fragile case, and what we measure on the real GLOBIOM
                                  # root (lambda_Cropland -0.131 [-0.245,-0.056] at 25k pixels).
                                  # An alternative-specific regressor is the classical fix: it makes
                                  # IV_c depend on the PATTERN of z across alternatives, which the
                                  # root design cannot reproduce. Validated by maximum likelihood in
                                  # experiments/nested/altspec_identification_proof.R: SE(lambda)
                                  # falls 1.4x-5.5x, the gain GROWING as identification degrades.
                                  # delta's full conditional is a conjugate 1-D Gaussian given the
                                  # PG weights, so this costs one scalar draw per sweep.
                                  alt_spec_Z = NULL, alt_spec_prior_sd = 10, alt_spec_allow_unvalidated = FALSE,
                                  positive_constraints = NULL, negative_constraints = NULL,
                                  use_ncp = TRUE, calc_loo = FALSE,
                                  support_prior_strength = 0, re_asis = FALSE, init_jitter = 0,
                                  collapse_re_var = FALSE, collapse_re_var_validated = FALSE,
                                  re_regularize = FALSE,   # apply the Finnish-HS variance cap to the STANDARD half-Cauchy RE prior too (update_re_precision_hc_sym), so the standard and collapsed paths target the SAME regularised-horseshoe posterior. The principled cure for the heavy-tailed per-(cov,cat) RE scales running off. Off by default (production-neutral); turn on together with collapse for a fair gate.
                                  collapse_slab_c2 = 100,  # shared regularised (Finnish-HS) cap scale for BOTH the collapsed and (when re_regularize=TRUE) the standard RE variance: effective var = c2*s2/(c2+s2) -> sigma <= sqrt(c2) (=10). Stops the unregularised half-Cauchy heavy tail running off (sigma=138). Tune via the gate.
                                  loo_scaling_factor = 1, chain_id = NULL, progress_cb = NULL,
                                  re_idx = 1:ncol(X), linear_idx = 1:ncol(X), prior_mu = NULL, prior_V_inv = NULL,
                                  use_horseshoe = FALSE, horseshoe_idx = NULL,
                                  equation_specific_hs = FALSE,
                                  estimate_c2 = FALSE, slab_df = 4, slab_s2 = 4, slab_s2_re = 4,
                                  # SEPARATE df for the RE slab. `slab_df` is the FE slab's and must
                                  # stay at the permissive P&V default (nu=4), where the slab exists to
                                  # let genuinely large coefficients ESCAPE shrinkage. The RE slab has
                                  # the opposite job -- impose a bound the data do not ask for -- so a
                                  # heavy tail there is self-defeating: at nu=4 the prior on c2 has
                                  # INFINITE variance (a = nu/2 = 2) and a 90% range of [1.7, 22.5],
                                  # i.e. a sigma cap anywhere in [1.3, 4.7]. nu=20 gives c2 90%
                                  # [2.6, 7.4] -> sigma cap [1.60, 2.72]: bounded, still free to move
                                  # if a node needs it. This also DISSOLVES the identifiability worry --
                                  # if the likelihood cannot inform c2 the posterior rests on the prior,
                                  # which IS the intended bound, instead of drifting as the learned tau did.
                                  slab_df_re = 10,   # MEASURED 2026-08-13 (Forests, both c2 starts):
                                  # nu=10 is where identification is cleanest (starts at c2=4 and
                                  # c2=100 both land at 4.12, agreement to 3 s.f.) AND held-out is best
                                  # (-725.5 vs -732.2 at nu=20 and -728.3 at FIXED c2=4). nu=6 gains
                                  # nothing and its starts diverge more (3.76 vs 3.40) with a much
                                  # heavier tail (90% upper 14.1). So the data DO locate c2 ~ 3.5-4.1 --
                                  # this is not the prior answering itself.
                                  p0_mu = NULL, tau0_mu = NULL,
                                  standardize = TRUE, method = c("standardize", "center", "scale", "QR", "none"),
                                  gamma_matched_pg = FALSE,
                                  init_state = NULL, prior_a_re = 0.01, prior_b_re = 0.01,
                                  use_half_cauchy_re = TRUE, re_scale_A = 1.0,
                                  # ── HIERARCHICAL HORSESHOE ON THE RE SCALES (2026-08-13) ──
                                  # OFF: each (cov,cat) cell gets an INDEPENDENT half-Cauchy(0, re_scale_A)
                                  # on sigma -- a local lambda with a FIXED scale and NO global tau, i.e.
                                  # not a horseshoe. Nothing supplies sparsity pressure, so at 2-20
                                  # obs/param the Cauchy tail wins: measured sigma median 2.16, q90 24.9,
                                  # max 87.7 on STANDARDISED covariates (61% above 1.0). Consequence:
                                  # random SLOPES scored worse than no REs at all (McFadden 0.188 vs
                                  # 0.331) -- see the re-block-size-tradeoff memory.
                                  # ON: sigma_k ~ C+(0, tau), tau ~ C+(0, re_hs_tau0) SHARED across all RE
                                  # cells, so most sigma_k are pulled toward 0 while a genuinely varying
                                  # covariate escapes on its heavy-tailed local lambda. The intercept is
                                  # deliberately NOT special-cased: the country level really does vary
                                  # (+36 held-out LL) so its own lambda lets it escape -- that is the
                                  # mechanism working, not a leak.
                                  # Conjugate given the existing aux (xi = 1/tau^2):
                                  #   xi | {a_k}, b ~ Gamma((K+1)/2, sum(a_k) + b);  b | xi ~ Exp(xi + 1/tau0^2)
                                  # and passing re_scale_A = tau into the C++ makes ITS aux draw
                                  # a_k ~ Exp(prec_k + 1/tau^2) the correct conditional -> no C++ change.
                                  # tau0 small = "few drivers truly vary" (measured ~0 of 8 on Forests).
                                  re_hs_global = FALSE, re_hs_tau0 = 0.1,
                                  # FULL BAYES ON THE RE SLAB. estimate_slab_c2 = TRUE samples c2
                                  # instead of fixing it, with the standard Finnish-horseshoe prior
                                  # c2 ~ InvGamma(slab_df/2, slab_df*slab_s2_re/2) (Piironen-Vehtari).
                                  # NOT conjugate: c2 enters through v~ = c2*s2/(c2+s2), so this is a
                                  # 1-D slice on log(c2), the same device the tau_raw draw already uses.
                                  # Motivation: the cap IS part of the estimation, so tuning it by an
                                  # external held-out criterion both burns a data split and biases any
                                  # number reported on that split. Sampling it costs no split, adapts
                                  # PER NODE (Forests carries 54 params/country, Pasture 36) and
                                  # propagates c2 uncertainty into the coefficients.
                                  # RISK TO WATCH: the slab works BECAUSE it imposes a bound the data
                                  # do not ask for. If the likelihood cannot constrain sigma (2-20
                                  # obs/param here) c2 may simply drift up and reproduce the
                                  # unregularised funnel -- exactly how the learned global tau failed
                                  # (re_hs_global: tau tracked the average sigma instead of shrinking
                                  # it). Trace `post_slab_c2` and check it settles rather than drifts.
                                  estimate_slab_c2 = FALSE,
                                  # PARAMETER-COUNT-AWARE SLAB. collapse_slab_c2 = "auto" sets the cap
                                  # from the number of RE covariates: the country deviation in the
                                  # linear predictor is eta_g = sum_v x_v * b_{v,g}, so with K
                                  # STANDARDISED covariates each of deviation SD sigma the TOTAL country
                                  # shift has SD ~ sigma*sqrt(K). Fixing a defensible total shift R (in
                                  # logits) gives sigma <= R/sqrt(K), i.e. c2 = R^2 / K. Adding
                                  # covariates then tightens each one's cap so the AGGREGATE country
                                  # effect stays plausible, instead of K of them each free to reach R.
                                  # R = 3 logits is a generous country-level shift in composition
                                  # log-odds. K=9 -> c2=1; K=1 -> c2=9 (a lone random intercept may be
                                  # larger precisely because it is the only one).
                                  re_slab_range = 3,
                                  # ── Phase 1: Empirical Bayes + Precision-weighted HS ──
                                  use_wls_init = TRUE,
                                  use_precision_hs = TRUE,
                                  # ── Phase 2: Spike-and-Slab RE Inclusions ─────────────
                                  use_spike_slab = TRUE,
                                  disable_separation_detection = FALSE,
                                  separation_as_prior = TRUE,
                                  sep_prior_floor = 0.05,
                                  separation_soft = FALSE,
                                  sep_overlap_hard = 0.0,
                                  sep_overlap_neutral = 0.5,
                                  a_pi = 1.0,
                                  b_pi = 3.0,
                                  store_delta = FALSE,
                                  # ── Phase 4: CAR Prior on Group Deviations ────────────
                                  use_car = TRUE,
                                  country_adjacency = NULL,
                                  car_rho = 0.95,
                                  a_spatial = 0.5,
                                  b_spatial = 0.5,
                                  # ── Reference Coding (static drop-one per const-sum block) ──
                                  const_sum_blocks = NULL,
                                  # ── General rank guard: drop constant / duplicate / exactly-collinear
                                  #    columns (data-quality) for a full-rank fit; reconstruct them as 0.
                                  handle_rank_deficiency = TRUE, rank_tol = 1e-7,
                                  symmetric_hs = FALSE,
                                  # CHANNEL-SPLIT HORSESHOE (symmetric_hs only, default OFF = historical).
                                  # Under symmetric_hs the const-sum block columns are shrunk by their own
                                  # within-block kernel (c_block) and are explicitly EXEMPTED from the
                                  # across-category kernel c_v. They nevertheless remained in hs_idx, so
                                  # they contributed k_block * p to the SHAPE of the global tau2 update
                                  # while their (c_block-crushed) coefficients contributed ~0 to its RATE
                                  # -- a one-way ratchet that collapses tau2 and crushes the unconstrained
                                  # drivers. TRUE gives the block channel its own global scale
                                  # (blk_tau2 / blk_xi) and restricts the covariate channel's pool to the
                                  # non-block columns, so neither channel's shape count includes columns
                                  # the other channel governs.
                                  hs_channel_split = FALSE,
                                  # LIVE FE-HORSESHOE KERNEL (diagonal path only, default OFF =
                                  # historical). hs_prec_kernel is bound ONCE before the Gibbs loop,
                                  # to a hs_prec_mat that is still all zeros; hs_prec_mat is only
                                  # populated inside the loop, and R copies on assignment, so the
                                  # binding never rebinds. The three C++ fast paths
                                  # (gibbs_step_pooled / gibbs_step_re / gibbs_step_re_ncp) therefore
                                  # receive a ZERO FE-horseshoe precision for the whole run; only the
                                  # R-level TMVN and ASIS-CP blocks read the live hs_prec_mat. This
                                  # was previously measured as inert, but that measurement predates
                                  # the NULL -> integer(0) fix, when hs_idx was empty and hs_prec_mat
                                  # never became non-zero. TRUE refreshes the binding each sweep.
                                  # It changes every diagonal-horseshoe fit -- gate on held-out
                                  # log-likelihood before adopting.
                                  hs_kernel_live = FALSE,
                                  # RE-SIDE SYMMETRIC VARIANCE, separable from the FE side.
                                  # NULL = follow symmetric_hs (historical). update_re_precision_hc_sym
                                  # forms its sum of squares from CATEGORY-CENTRED deviations
                                  # (d - mean(d) over p_all) but the resulting sigma is applied to an
                                  # UNCENTRED draw in gibbs_step_re_ncp -- the same CLR deflation the FE
                                  # side had. Set FALSE to keep symmetric FE shrinkage while routing the
                                  # RE variance through the plain update_re_precision_hc.
                                  re_prec_sym = NULL,
                                  # NULL (default) = CORRECT BY CONSTRUCTION: resolves to FALSE
                                  # whenever the symmetric RE updater is in use, TRUE otherwise.
                                  # update_re_precision_hc_sym builds its sum of squares from
                                  # category-CENTRED deviations while gibbs_step_re_ncp draws the REs
                                  # UNCENTRED, and nothing constrains them to the zero-sum subspace,
                                  # so the centred statistic is never the right one there -- it zeroes
                                  # the random effects (measured: RE sd 0.0000, -51.6 nats held-out).
                                  # That is a BUG, not a modelling option, so no caller should have to
                                  # remember to opt out: TRUE is reachable only by asking for it
                                  # explicitly, to reproduce a pre-2026-08-21 fit.
                                  re_prec_center = NULL,
                                  # PER-FAMILY GLOBAL SCALES. NULL = one tau over every shrunk column
                                  # (historical). Otherwise a NAMED LIST OF REGEXES matched against
                                  # colnames(X), e.g. list(topo = "^(Slope_rad|Elevation|Aspect)",
                                  # socio = "^(log1p_|GHM_|CISI)"). Each family gets its OWN (tau2, xi)
                                  # with a shape from its OWN column count, and each const-sum block
                                  # gets its own too. Rationale: one global tau over heterogeneous
                                  # columns lets a family of near-null covariates drag down a family of
                                  # strong ones -- lambda is per-column but cannot rescue a covariate
                                  # once tau has collapsed underneath it. Implies hs_channel_split.
                                  # Columns matching nothing fall into a "_rest" family.
                                  hs_groups = NULL) {
  # --- 1. SETUP (identical to original) ---
  # Store-time self-check (default OFF): options(mnlogit.selfcheck=TRUE) prints,
  # for the first N stored draws, the log-lik recomputed from the exact stored
  # zero-sum coefficient array vs the sampler's own ll$total. Localizes any
  # coefficient/likelihood divergence to the storage/expansion step.
  .selfcheck_store <- isTRUE(getOption("mnlogit.selfcheck", FALSE))
  .selfcheck_n     <- as.integer(getOption("mnlogit.selfcheck_n", 3L))
  if (niter <= nburn) stop("niter must be > nburn.")
  if (use_bart && nrow(X) > 50000 && nburn < 2000) {
    warning("For massive N with BART, consider nburn > 5000.")
  }
  if (any(!is.finite(Y))) stop("Y contains NA/Inf.")
  if (any(!is.finite(X))) stop("X contains NA/Inf.")
  if (any(rowSums(Y) <= 0)) stop("Zero/negative row sums in Y.")

  bart_disk_path <- disk_path
  posterior_disk_path <- disk_path

  if (is.null(y_weight)) {
    y_weight <- 1 / min(rowSums(Y))
    message("Auto y_weight: ", round(y_weight, 6))
  }
  if (missing(loo_scaling_factor) || loo_scaling_factor == 1) {
    loo_scaling_factor <- 1 / y_weight
  }

  nn_weighted <- as.vector(rowSums(Y) * y_weight)
  kappa_weighted <- (Y * y_weight) - (nn_weighted / 2)

  X <- as.matrix(X)
  if (intercept && !any(colSums(X) == nrow(X))) {
    X <- cbind(X, intercept = 1)
  }
  X_orig <- X

  # Convert character indices to integers early for robustness
  if (is.character(linear_idx)) linear_idx <- match(linear_idx, colnames(X))
  if (is.character(bart_idx)) bart_idx <- match(bart_idx, colnames(X))
  if (is.character(re_idx)) re_idx <- match(re_idx, colnames(X))
  if (is.character(horseshoe_idx)) horseshoe_idx <- match(horseshoe_idx, colnames(X))
  if (is.character(positive_constraints)) positive_constraints <- match(positive_constraints, colnames(X))
  if (is.character(negative_constraints)) negative_constraints <- match(negative_constraints, colnames(X))

  # Remove NAs from indexing (variables not found)
  linear_idx <- linear_idx[!is.na(linear_idx)]
  bart_idx <- bart_idx[!is.na(bart_idx)]

  n <- nrow(X)
  k <- ncol(X)
  p_all <- ncol(Y)
  p <- p_all - 1
  pp <- (1:p_all)[-baseline]

  # --- ALTERNATIVE-SPECIFIC design: to BASELINE-REMOVED coordinates -------------------------
  # Utilities are carried as differences from the baseline alternative, so an attribute of
  # alternative j contributes (z_ij - z_i,baseline) to utility j -- the standard conditional-logit
  # transform. Only these differences are identified; a common shift across alternatives cancels.
  use_alt_spec <- !is.null(alt_spec_Z)
  Zt <- NULL
  if (use_alt_spec && !isTRUE(alt_spec_allow_unvalidated)) stop(
    "alt_spec_Z is IMPLEMENTED BUT NOT CORRECT -- delta is not recovered. Do not use for results.\n",
    "  Attempt 1 (conjugate PG draw treating c_j as fixed): delta biased by a clean factor ~2\n",
    "     true 1/2/3 -> 0.489/0.999/1.462 (ratios 2.046/2.002/2.053). Cause: this sampler carries a\n",
    "     ONE-VS-REST decomposition, predictor (psi_j - c_j) with c_j = log sum_{k!=j} exp(psi_k);\n",
    "     delta enters psi_k for EVERY alternative so c_j depends on delta, and the correct effective\n",
    "     regressor is Zt_j - d c_j/d delta, not Zt_j.\n",
    "  Attempt 2 (random-walk Metropolis on the exact likelihood): WORSE -- true 1/2/3 -> -1.301/\n",
    "     -1.109/0.623, no consistent relation, i.e. the utility reconstruction inside the MH step\n",
    "     does not match what the sweep actually uses.\n",
    "  A correct implementation needs the conditional derived properly against the c_j decomposition\n",
    "  (or an MH step that reuses the sweep's own utility path rather than rebuilding it).\n",
    "  The STATISTICAL case for the term is sound and independently validated by maximum likelihood in\n",
    "  experiments/nested/altspec_identification_proof.R (SE(lambda) 1.4x-5.5x tighter).\n",
    "  Pass alt_spec_allow_unvalidated=TRUE only to work ON this feature.")
  if (use_alt_spec) {
    alt_spec_Z <- as.matrix(alt_spec_Z)
    if (nrow(alt_spec_Z) != nrow(Y) || ncol(alt_spec_Z) != p_all)
      stop(sprintf("alt_spec_Z must be %d x %d (n x p_all); got %d x %d",
                   nrow(Y), p_all, nrow(alt_spec_Z), ncol(alt_spec_Z)))
    Zt <- alt_spec_Z[, pp, drop = FALSE] - alt_spec_Z[, baseline]
    Zt[!is.finite(Zt)] <- 0
    cat(sprintf("Alternative-specific term ON: 1 shared delta over %d alternatives (sd(Zt)=%.4f)\n",
                p_all, stats::sd(Zt)))
  }
  curr_delta <- 0
  delta_mh_sd <- 0.5; delta_acc_n <- 0L; delta_acc_k <- 0L   # adaptive RW-Metropolis state for delta

  # --- 2. PREPROCESSING ---
  if (missing(method)) {
    methods_chosen <- if (standardize) "standardize" else "none"
  } else {
    methods_chosen <- match.arg(method, several.ok = TRUE)
  }

  do_std <- "standardize" %in% methods_chosen
  do_cen <- "center" %in% methods_chosen || do_std
  do_scl <- "scale" %in% methods_chosen || do_std
  do_qr <- "QR" %in% methods_chosen

  # Identify continuous variables (needed for rescaling AND BART scaling)
  # Continuous variables are only candidates for centering/scaling if they are not naturally bounded in [-1, 1] (e.g. shares, contrasts, bounded indices)
  cont_idx <- apply(X, 2, function(x) {
    is.numeric(x) &&
      length(unique(x)) > 2 &&
      max(abs(x), na.rm = TRUE) > 1.0
  })

  # ── Native constant-sum handling ──
  # const_sum_blocks = "auto": detect compositional (sum-to-constant) blocks HERE on the raw design
  # (row-sums intact, before centering). Downstream this drops one reference per block for a full-rank
  # fit and reconstructs the block MEAN-CENTERED (sum-to-zero) over ALL original columns (see the
  # do_block_rotation reconstruction), so a RAW X goes in and a clean full-column posterior comes out.
  # Continuous indices/values (e.g. yield_index) are not compositional and are never flagged.
  if (is.character(const_sum_blocks) && length(const_sum_blocks) == 1L && const_sum_blocks == "auto") {
    const_sum_blocks <- tryCatch(detect_constant_sum_blocks(X), error = function(e) NULL)
    if (length(const_sum_blocks)) {
      .bn <- names(const_sum_blocks); if (is.null(.bn)) .bn <- as.character(seq_along(const_sum_blocks))
      cat(sprintf("Auto-detected %d constant-sum block(s): %s\n", length(const_sum_blocks),
                  paste(sprintf("%s(%d cols)", .bn, lengths(const_sum_blocks)), collapse = ", ")))
    } else { const_sum_blocks <- NULL; cat("Auto const-sum: no compositional blocks detected.\n") }
  }

  if (!is.null(const_sum_blocks)) {
    block_col_idx <- unique(unlist(const_sum_blocks))
    # Note: block_col_idx indices are currently wrt full X_mat from the wrapper.
    # The wrapper's X_mat corresponds perfectly to this function's X at this point,
    # because we haven't done any linear_idx subsetting yet.
    cont_idx[block_col_idx] <- FALSE
  }
  
  X_rescaling <- matrix(0, 2, sum(cont_idx))
  rownames(X_rescaling) <- c("mean", "sd")
  X_rescaling[2, ] <- 1 # Default SD=1

  if (do_cen || do_scl) {
    if (sum(cont_idx) > 0) {
      if (do_cen) {
        X_rescaling[1, ] <- colMeans(X[, cont_idx, drop = FALSE])
        X[, cont_idx] <- sweep(X[, cont_idx, drop = FALSE], 2, X_rescaling[1, ], "-")
      }
      if (do_scl) {
        X_names <- colnames(X)
        X_rescaling[2, ] <- apply(X[, cont_idx, drop = FALSE], 2, sd)
        X_rescaling[2, X_rescaling[2, ] == 0] <- 1
        X[, cont_idx] <- sweep(X[, cont_idx, drop = FALSE], 2, X_rescaling[2, ], "/")
        colnames(X) <- X_names
      }
    }
  }

  # --- 2.05 Linear Design Matrix Selection & Index Remapping ---
  X_full_processed <- X
  X_full_for_bart <- X

  if (intercept) {
    # Ensure intercept is in linear_idx if it was added/exists
    int_pos <- which(colnames(X) == "intercept")
    if (length(int_pos) > 0 && !(int_pos %in% linear_idx)) {
      linear_idx <- c(linear_idx, int_pos)
    }
  }

  # Map indices from full X to linear subset (to maintain support for re_idx, horseshoe_idx, etc.)
  map_idx <- function(idx, full_idx) {
    if (is.null(idx)) {
      return(NULL)
    }
    m <- match(idx, full_idx)
    m[!is.na(m)]
  }

  re_idx <- map_idx(re_idx, linear_idx)
  horseshoe_idx <- map_idx(horseshoe_idx, linear_idx)
  positive_constraints <- map_idx(positive_constraints, linear_idx)
  negative_constraints <- map_idx(negative_constraints, linear_idx)
  
  if (!is.null(const_sum_blocks)) {
    const_sum_blocks <- lapply(const_sum_blocks, function(b) map_idx(b, linear_idx))
    # Remove any blocks that don't have at least 2 active elements left
    const_sum_blocks <- const_sum_blocks[sapply(const_sum_blocks, length) >= 2]
    if (length(const_sum_blocks) == 0) const_sum_blocks <- NULL
  }

  # Redefine X and k for the linear part
  X <- X_full_processed[, linear_idx, drop = FALSE]
  k <- ncol(X)
  X_orig <- X_orig[, linear_idx, drop = FALSE]

  bart_scaling <- NULL
  if (use_bart) {
    # Robust bart_scaling: map bart_idx to X_rescaling columns
    # Note: X_rescaling only contains continuous variables
    bart_mu <- rep(0, length(bart_idx))
    bart_sd <- rep(1, length(bart_idx))

    # full_to_rescale maps full column index to X_rescaling column index
    full_to_rescale <- rep(NA, ncol(X_full_processed))
    full_to_rescale[cont_idx] <- 1:sum(cont_idx)

    for (i in seq_along(bart_idx)) {
      full_idx <- bart_idx[i]
      rescale_idx <- full_to_rescale[full_idx]
      if (!is.na(rescale_idx)) {
        bart_mu[i] <- X_rescaling[1, rescale_idx]
        bart_sd[i] <- X_rescaling[2, rescale_idx]
      }
    }
    bart_scaling <- list(mu = bart_mu, sd = bart_sd)
  }

  if (do_qr) {
    if (use_horseshoe || !is.null(positive_constraints) || !is.null(negative_constraints)) {
      warning(
        "QR decomposition ('method=QR') is active alongside Horseshoe shrinkage or sign constraints.\n",
        "Back-transformation (beta = R^-1 * theta) will destroy coefficient sparsity and sign consistency.\n",
        "Recommendation: Use 'method=c(\"center\", \"scale\")' instead of 'QR' for these configurations."
      )
    }
    # Apply QR only to the subsetted linear matrix
    qr_X <- qr(X)
    X_R_mat <- qr.R(qr_X)
    X <- qr.Q(qr_X)
  }

  # Update continuous variable tracking for back-transformation
  if (do_cen || do_scl) {
    # 1. Map each column of the FULL X to its position in the ORIGINAL X_rescaling
    full_to_rescale <- rep(NA, length(cont_idx))
    names(full_to_rescale) <- names(cont_idx)
    full_to_rescale[cont_idx] <- seq_len(sum(cont_idx))

    # 2. For the subsetted X (linear_idx), find which rescale columns we need
    needed_rescale_cols <- full_to_rescale[linear_idx]

    # Validation Check
    if (any(cont_idx[linear_idx] & is.na(needed_rescale_cols))) {
      offenders <- linear_idx[which(cont_idx[linear_idx] & is.na(needed_rescale_cols))]
      offender_names <- if (!is.null(colnames(X_full_processed))) colnames(X_full_processed)[offenders] else "unnamed"
      stop(sprintf(
        "Alignment Error: Continuous variables at indices [%s] (%s) have no corresponding rescaling factors in X_rescaling.",
        paste(offenders, collapse = ","), paste(offender_names, collapse = ",")
      ))
    }

    # Keep only those that are actually continuous
    needed_rescale_cols <- needed_rescale_cols[!is.na(needed_rescale_cols)]

    # 3. Subset X_rescaling and cont_idx
    X_rescaling <- X_rescaling[, needed_rescale_cols, drop = FALSE]
    cont_idx <- cont_idx[linear_idx]

    # Final Sanity Check
    if (sum(cont_idx) != ncol(X_rescaling)) {
      stop(sprintf(
        "Dimension Mismatch: %d continuous variables found but X_rescaling has %d columns.",
        sum(cont_idx), ncol(X_rescaling)
      ))
    }
  }

  # --- 2.1 BART Design Matrix Selection ---
  if (use_bart) {
    # Select subset of processed covariates for BART from the full matrix
    X_bart_final <- X_full_for_bart[, bart_idx, drop = FALSE]
    # Robust calibration references: use empirical range to ensure different terminal nodes
    bart_ref1_mat <- matrix(apply(X_bart_final, 2, min), 1, length(bart_idx))
    bart_ref2_mat <- matrix(apply(X_bart_final, 2, max), 1, length(bart_idx))
  }

  # --- 3. INITIALIZATION ---
  # Hoist matrix conversions and unnaming to avoid overhead in Gibbs loop
  cov_names_save <- colnames(X)
  int_idx <- integer(0)
  if (!is.null(cov_names_save)) {
    int_idx <- which(tolower(cov_names_save) == "intercept")
  }
  if (length(int_idx) == 0) {
    int_idx <- which(apply(X_orig, 2, function(col) all(col == 1)))
  }
  X <- as.matrix(unname(X))
  Xt <- as.matrix(unname(t(X)))
  kappa_weighted <- as.matrix(unname(kappa_weighted))

  curr_beta <- matrix(0, k, p)
  
  horseshoe_idx_pre_drop <- horseshoe_idx
  
  # --- Design reduction: const-sum reference coding (mean-centered reconstruction) + rank guard ---
  # (1) const-sum blocks -> drop one reference each, reconstruct MEAN-CENTERED over the block.
  # (2) rank guard -> drop constant / duplicate / exactly-collinear columns, reconstructed as 0.
  # Both feed ONE keep_mask; the do_block_rotation reconstruction returns a FULL-column posterior.
  block_info <- list(); cs_drop <- integer(0)
  if (!is.null(const_sum_blocks) && length(const_sum_blocks) > 0) {
    if (use_re) for (idx in const_sum_blocks) {
      n_re <- sum(idx %in% re_idx)
      if (n_re > 0 && n_re < length(idx))
        stop("Block contains a mix of fixed and random effects. Blocks must be all fixed or all RE.")
    }
    block_info <- lapply(const_sum_blocks, function(idx) {
      K_blk <- length(idx); vars <- apply(X[, idx, drop = FALSE], 2, var)
      drop_col <- as.integer(which.max(vars))
      # `const` = the block's (near-)constant row sum. The reconstruction re-centres the block
      # (b_j -> b_j - mean_k), which changes each class utility by -mean_k * rowSum(block) =
      # -mean_k * const: a PER-CLASS constant. Without compensating the intercept by
      # +mean_k*const the exported full-column coefficients are NOT softmax-equivalent to the
      # fitted model (verified: mean|P_hat - P_true| 0.27 vs 0.003). See the reconstruction below.
      rs <- rowSums(X[, idx, drop = FALSE]); rs <- rs[is.finite(rs) & rs > 1e-4]
      list(all_idx = idx, K = K_blk, current_drop = drop_col, active = setdiff(1:K_blk, drop_col),
           const = if (length(rs)) stats::median(rs) else 0)
    })
    cs_drop <- vapply(block_info, function(b) b$all_idx[b$current_drop], integer(1))
    # intercept column (pre-drop coordinates) that absorbs the re-centring shift. Without one
    # the shift cannot be compensated -> warn rather than silently export biased coefficients.
    .cs_int <- which(apply(X, 2, function(z) all(is.finite(z) & z == 1)))
    cs_int_col <- if (length(.cs_int)) as.integer(.cs_int[1]) else NA_integer_
    if (is.na(cs_int_col))
      warning("const_sum_blocks: no intercept column found; the zero-sum reconstruction shifts ",
              "each class utility by a constant that cannot be absorbed -> the returned ",
              "coefficients will not reproduce the fitted probabilities. Add an intercept column.")
  } else cs_int_col <- NA_integer_
  rd_drop <- integer(0)
  if (isTRUE(handle_rank_deficiency)) {
    int_cols  <- which(apply(X, 2, function(z) all(z == 1)))                 # intercept(s): never dropped
    protected <- unique(c(int_cols, unlist(const_sum_blocks)))               # const-sum blocks handled above
    cand   <- setdiff(seq_len(k), protected)
    consts <- cand[vapply(cand, function(j) { v <- var(X[, j]); is.na(v) || v < rank_tol }, logical(1))]
    keep0  <- setdiff(cand, consts)
    extra  <- integer(0)
    if (length(keep0) >= 1) {
      allc <- sort(unique(c(int_cols, keep0)))
      qq <- qr(X[, allc, drop = FALSE], tol = rank_tol)
      if (qq$rank < length(allc)) extra <- setdiff(allc[qq$pivot[(qq$rank + 1):length(allc)]], int_cols)
    }
    rd_drop <- setdiff(unique(c(consts, extra)), protected)
    if (length(rd_drop)) {
      nm <- if (!is.null(colnames(X))) colnames(X)[rd_drop] else as.character(rd_drop)
      cat(sprintf("Rank guard: dropped %d degenerate column(s) (constant/duplicate/collinear) -> coef 0: %s\n",
                  length(rd_drop), paste(nm, collapse = ", ")))
    }
  }
  drop_all <- unique(c(cs_drop, rd_drop))
  if (length(drop_all) > 0) {
    do_block_rotation <- TRUE
    keep_mask <- rep(TRUE, k); keep_mask[drop_all] <- FALSE
    X_full <- X; X <- X[, keep_mask, drop = FALSE]; Xt <- as.matrix(unname(t(X))); k_active <- ncol(X)
    full_to_active <- rep(NA, k); full_to_active[keep_mask] <- 1:k_active
    # PRESERVE NULL. `full_to_active[NULL]` yields integer(0), NOT NULL, which silently broke two
    # downstream is.null() guards whenever ANY column was dropped (const-sum ref or rank guard):
    #  (1) hs_idx <- if (is.null(horseshoe_idx)) 1:k else horseshoe_idx -> EMPTY -> k_hs = 0, so the
    #      FE horseshoe applied to NO covariate (log signature: "tau0_pooled=-0.0331 (p0=-0.5)";
    #      p0 = min(p0, k_hs-0.5) = -0.5 and tau0 goes NEGATIVE).
    #  (2) has_constraints <- !is.null(positive_constraints) || ... -> TRUE with none requested, so
    #      the sampler ran the constrained TMVN branch: a few Gibbs sweeps instead of an exact
    #      Cholesky draw -> slower mixing. Measured (Forests n=4000, 600 iter, horseshoe OFF):
    #      held-out -744.7 -> -714.4 (+30.3 LL), 1.11x faster once fixed.
    .remap <- function(idx) if (is.null(idx)) NULL else
                            as.integer(unname(na.omit(full_to_active[idx])))
    re_idx <- .remap(re_idx)
    horseshoe_idx <- .remap(horseshoe_idx)
    positive_constraints <- .remap(positive_constraints)
    negative_constraints <- .remap(negative_constraints)
    int_idx <- .remap(int_idx)
    k <- k_active
  } else {
    do_block_rotation <- FALSE; X_full <- X; full_to_active <- 1:k
  }
  n_blocks <- length(block_info)

  curr_beta <- matrix(0, k, p)
  bart_shifts <- matrix(0, if (bart_symmetric) p_all else p, 1)  # one row per ensemble (p_all if CLR)
  prior_V_inv_originally_null <- is.null(prior_V_inv)
  if (is.null(prior_mu)) prior_mu <- matrix(0, k, p)
  if (!is.matrix(prior_mu)) prior_mu <- matrix(as.vector(prior_mu), k, p)
  bart_alpha <- 1.0 # NEW: Global tempering factor
  if (is.null(prior_V_inv)) {
    # Default: A0 for predictors, but allow wide prior (variance = 100.0, sd = 10) for the intercept (similar to brms/Stan)
    diag_vals <- rep(1 / A0, k)
    if (length(int_idx) > 0) {
      diag_vals[int_idx] <- 1 / 100.0
    } else if (intercept) {
      diag_vals[k] <- 1 / 100.0
    }
    prior_V_inv <- diag(diag_vals, k)
  }

  if (use_re && use_half_cauchy_re) {
    a_re <- matrix(re_scale_A^2 / 2, nrow = k, ncol = p)
  } else {
    a_re <- NULL
  }

  if (empirical_intercept_prior) {
    if (length(int_idx) > 0) {
      Y_smoothed <- colSums(Y) + 1.0
      empirical_log_odds <- log(Y_smoothed / Y_smoothed[baseline])
      empirical_mu <- empirical_log_odds[-baseline]

      for (idx in int_idx) {
        prior_mu[idx, ] <- empirical_mu
        if (prior_V_inv_originally_null) {
          prior_V_inv[idx, idx] <- 0.25
        }
      }
    }
  }

  prior_P <- as.matrix(unname(prior_V_inv))
  prior_Pb <- as.matrix(unname(prior_P %*% prior_mu)) # k x p, static

  # Baseline RE precision: set to effectively infinite (1e12) to lock fixed effects
  # Initialize at the half-Cauchy prior mode for re_scale_A for early iteration exploration
  if (is.character(collapse_slab_c2) && identical(collapse_slab_c2[1], "auto")) {
    .Kre <- max(1L, length(re_idx))
    collapse_slab_c2 <- (re_slab_range^2) / .Kre
    if (isTRUE(use_re)) cat(sprintf("Auto slab: %d RE covariate(s), total country shift R=%.1f -> c2 = R^2/K = %.3f (sigma <= %.3f)\n",
        .Kre, re_slab_range, collapse_slab_c2, sqrt(collapse_slab_c2)))
  }
  prec_init <- 1 / (re_scale_A^2)
  prec_beta_pooled <- matrix(prec_init, k, p)
  sigma_beta_pooled <- matrix(re_scale_A, k, p)
  # global horseshoe scale state: tau_cur is what gets passed as the half-Cauchy scale each sweep.
  # Starts at re_scale_A so re_hs_global = FALSE reproduces the previous behaviour EXACTLY.
  re_tau_cur <- re_scale_A
  re_hs_xi   <- 1 / (re_scale_A^2)     # xi = 1/tau^2
  re_hs_b    <- 1 / (re_hs_tau0^2)     # aux for tau's own half-Cauchy

  if (length(re_idx) < k) {
    fixed_idx <- setdiff(seq_len(k), re_idx)
    prec_beta_pooled[fixed_idx, ] <- 1e12
    sigma_beta_pooled[fixed_idx, ] <- 1e-6
  }

  # Precompute complement of re_idx to avoid recomputing in the Gibbs loops
  re_complement <- setdiff(seq_len(k), re_idx)

  mu_pooled <- curr_beta # Initial reference

  if (use_re) {
    groups <- unique(group_idx)
    n_groups <- length(groups)
    # Pre-convert to unnamed integer lists and matrices for Rcpp bridge speed
    idx_list <- lapply(groups, function(m) as.integer(unname(which(group_idx == m))))
    Xm_list <- lapply(idx_list, function(idx) as.matrix(unname(X[idx, , drop = FALSE])))
    Xmt_list <- lapply(Xm_list, function(m) as.matrix(unname(t(m))))

    # Precompute structural masking matrix for group-specific REs
    re_mask <- matrix(1.0, nrow = k, ncol = n_groups)
    y_mask <- matrix(1.0, nrow = p, ncol = n_groups)
    sparsity_flagged <- matrix(FALSE, nrow = k, ncol = n_groups)

    # Support-aware RE prior: per (covariate x group) SD multiplier in (0,1] = 1/sqrt((n_g/PR)^strength),
    # PR = within-group participation ratio (effective # obs informing that group's slope). Shrinks
    # information-sparse RE slopes (the graduated generalization of the A2 no-variation pin); column m
    # aligns with Xm_list[[m]] / re_mask column m. All-ones when off. Globally-constant cols (random
    # intercept) exempted. Applied as the SD scaler in the C++ NCP core.
    re_support_mat <- matrix(1.0, nrow = k, ncol = n_groups)
    if (support_prior_strength > 0 && length(re_idx) >= 1) {
      re_glob_const <- apply(X[, re_idx, drop = FALSE], 2, function(cc) { s <- sd(cc); !is.finite(s) || s < 1e-8 })
      for (m in seq_len(n_groups)) {
        Xg <- Xm_list[[m]][, re_idx, drop = FALSE]; ng <- nrow(Xg)
        Xc <- sweep(Xg, 2, colMeans(Xg), "-")             # within-group centering -> slope information (PR is not location-invariant)
        s2 <- colSums(Xc^2); s4 <- colSums(Xc^4); pr <- s2^2 / pmax(s4, 1e-12)
        sd_fac <- 1 / sqrt(pmax((ng / pmax(pr, 1))^support_prior_strength, 1e-12))
        sd_fac[re_glob_const] <- 1.0
        re_support_mat[re_idx, m] <- sd_fac
      }
      message(sprintf("RE support prior (participation ratio, strength %.2f): min SD factor %.3f.", support_prior_strength, min(re_support_mat[re_idx, ])))
    }

    # Sketch A: per-cell separation prior offset (<=0), added to the spike-slab
    # inclusion log-odds instead of hard-masking quasi-separated cells. Default 0
    # (neutral). Calibrated so a fully-disjoint cell (overlap<=0) gets prior
    # inclusion ~ sep_prior_floor, and a well-overlapping cell (overlap>=0.5) is
    # neutral (offset 0). Structural non-identification stays a hard mask below.
    sep_logit <- array(0.0, dim = c(k, p, n_groups))
    SEP_A <- log(sep_prior_floor) - log1p(-sep_prior_floor)   # logit(floor)
    SEP_B <- -2 * SEP_A                                        # offset -> 0 at overlap=0.5
    absent_category <- matrix(FALSE, nrow = p, ncol = n_groups)
    is_global_intercept <- sapply(seq_len(k), function(v) {
      all(X[, v] == 1.0)
    })
    group_sizes <- sapply(idx_list, length)

    for (v in seq_len(k)) {
      if (is_global_intercept[v]) next
      var_global <- var(X[, v])
      if (is.na(var_global) || var_global < 1e-12) {
        re_mask[v, ] <- 0.0
        next
      }
      for (m in 1:n_groups) {
        var_local <- var(Xm_list[[m]][, v])
        std_var <- var_local / var_global
        n_nonzero <- sum(Xm_list[[m]][, v] != 0.0)

        # Scale-dependent non-zero threshold (3 to 10 observations depending on group size)
        x_thresh <- min(10L, max(3L, as.integer(ceiling(0.02 * group_sizes[m]))))

        # Hard mask TRUE rank deficiency; soft flag others
        if (is.na(var_local) || var_local < 1e-12) {
          re_mask[v, m] <- 0.0
        } else if (is.na(std_var) || std_var < 0.005 || n_nonzero < x_thresh) {
          sparsity_flagged[v, m] <- TRUE
        }
      }
    }

    # Y-side structural masking [p x n_groups] based on cell counts and standardized variance
    mean_pixel_total <- mean(rowSums(Y))

    for (ip in seq_len(p)) {
      j <- pp[ip]
      var_global <- var(Y[, j])
      if (is.na(var_global) || var_global < 1e-8) {
        y_mask[ip, ] <- 0.0
        next
      }

      for (m in seq_len(n_groups)) {
        group_obs <- Y[idx_list[[m]], j]

        # Calculate number of effective positive observations (scale/unit invariant)
        n_positive_eff <- sum(group_obs) / mean_pixel_total

        var_local <- var(group_obs)
        std_var <- var_local / var_global

        # Scale-dependent category count threshold (2 to 5 observations depending on group size)
        y_thresh <- min(5L, max(2L, as.integer(ceiling(0.01 * group_sizes[m]))))

        # Hard mask TRUE zero variance; soft flag absent category
        n_raw_present <- sum(group_obs > mean_pixel_total * 0.001)
        if (n_raw_present == 0) {
          absent_category[ip, m] <- TRUE
        }
        if (is.na(var_local) || var_local < 1e-8) {
          y_mask[ip, m] <- 0.0
        }
      }
    }



    # Diagnostic summary
    n_y_masked <- sum(y_mask == 0)
    if (n_y_masked > 0) {
      message(sprintf(
        "Y-mask: %d of %d category-group combinations fully masked (constant share or sparse cell).",
        n_y_masked, p * n_groups
      ))
    }

    # =================================================================
    # AUTOMATIC RE COLLINEARITY SCREEN
    # Computes within-group pairwise correlations for all RE predictors.
    # For any pair exceeding cor_threshold in a majority of groups,
    # the variable with lower global variance is removed from re_mask
    # for ALL groups (structural collinearity → pooled-only effect).
    # Variables removed only in specific groups (local collinearity)
    # are masked just for those groups.
    # =================================================================
    global_collinear_idx <- integer(0)   # default: no collinear variables

    if (length(re_idx) >= 2) {
      re_col_names <- colnames(X)[re_idx] # names for messaging
      n_re_vars <- length(re_idx)

      # --- Per-group correlation matrices ---
      # pairwise_hits[i, j] = number of groups where |cor(i,j)| > threshold
      pairwise_hits <- matrix(0L, n_re_vars, n_re_vars)
      pairwise_maxcor <- matrix(0.0, n_re_vars, n_re_vars)
      group_cor_list <- vector("list", n_groups)

      for (m in seq_len(n_groups)) {
        idx_m <- idx_list[[m]]
        if (length(idx_m) < 10L) {
          group_cor_list[[m]] <- matrix(NA_real_, n_re_vars, n_re_vars)
          next
        }
        Xm_re <- Xm_list[[m]][, re_idx, drop = FALSE]

        # Suppress warnings for zero-variance columns (already masked)
        cm <- suppressWarnings(cor(Xm_re))
        cm[is.na(cm)] <- 0.0
        diag(cm) <- 0.0 # ignore self-correlation
        group_cor_list[[m]] <- cm

        hit <- abs(cm) > re_cor_threshold
        pairwise_hits <- pairwise_hits + hit
        pairwise_maxcor <- pmax(pairwise_maxcor, abs(cm))
      }

      # Valid groups (enough obs) for denominator
      n_valid_groups <- sum(sapply(group_cor_list, function(cm) !all(is.na(cm))))
      hit_fraction <- pairwise_hits / max(n_valid_groups, 1L)

      # --- Identify collinear pairs ---
      # Work upper-triangle only to avoid double-processing
      globally_removed <- integer(0) # re_idx positions removed everywhere
      locally_masked <- list() # re_idx positions masked per group

      for (i in seq_len(n_re_vars - 1L)) {
        for (j in (i + 1L):n_re_vars) {
          if (hit_fraction[i, j] < re_cor_min_groups) next # not widespread enough

          # Decide which variable to drop: lower global variance loses
          var_i <- var(X_orig[, re_idx[i]])
          var_j <- var(X_orig[, re_idx[j]])
          drop_pos <- if (var_i <= var_j) i else j # position in re_idx
          keep_pos <- if (drop_pos == i) j else i

          is_global <- hit_fraction[i, j] >= 0.75 # >75% of groups → global drop

          if (is_global) {
            if (!(drop_pos %in% globally_removed)) {
              globally_removed <- c(globally_removed, drop_pos)
              message(sprintf(
                paste0(
                  "RE collinearity [GLOBAL]: removing '%s' from re_idx ",
                  "(r > %.2f in %.0f%% of groups; keeping '%s')."
                ),
                re_col_names[drop_pos],
                re_cor_threshold,
                hit_fraction[i, j] * 100,
                re_col_names[keep_pos]
              ))
            }
          } else {
            # Local: mask only in the offending groups
            bad_groups <- which(sapply(
              group_cor_list,
              function(cm) !is.na(cm[i, j]) && abs(cm[i, j]) > re_cor_threshold
            ))
            for (m in bad_groups) {
              if (re_mask[re_idx[drop_pos], m] != 0.0) {
                re_mask[re_idx[drop_pos], m] <- 0.0
                locally_masked[[length(locally_masked) + 1L]] <-
                  list(
                    var = re_col_names[drop_pos], group = m,
                    cor = group_cor_list[[m]][i, j]
                  )
              }
            }
            if (length(locally_masked) > 0L) {
              message(sprintf(
                paste0(
                  "RE collinearity [LOCAL]: masking '%s' in %d group(s) ",
                  "(r > %.2f; keeping '%s')."
                ),
                re_col_names[drop_pos],
                length(bad_groups),
                re_cor_threshold,
                re_col_names[keep_pos]
              ))
            }
          }
        }
      }

      # --- A2: within-group no-X-variation (the correlation screen misses constants:
      #     constant column -> NA cor -> 0 -> never flagged). Mask such RE slopes locally;
      #     exempt globally-constant columns (the intercept's random intercept is identified). ---
      re_global_const <- apply(X[, re_idx, drop = FALSE], 2, function(cc) { s <- sd(cc); !is.finite(s) || s < 1e-8 })
      n_a2 <- 0L
      for (m in seq_len(n_groups)) {
        sds_m <- apply(Xm_list[[m]][, re_idx, drop = FALSE], 2, sd)
        for (r in which((!is.finite(sds_m) | sds_m < 1e-8) & !re_global_const)) {
          if (re_mask[re_idx[r], m] != 0.0) { re_mask[re_idx[r], m] <- 0.0; n_a2 <- n_a2 + 1L }
        }
      }
      if (n_a2 > 0L) message(sprintf("RE no-variation screen (A2): masked %d (group x covariate) slopes with no within-group X variation.", n_a2))

      # --- Convert global removals to informative prior ---
      if (length(globally_removed) > 0L) {
        global_collinear_idx <- re_idx[unique(globally_removed)]
        message(sprintf(
          "RE collinearity screen complete: %d variable(s) globally demoted to informative pi_inc prior.",
          length(unique(globally_removed))
        ))
      } else {
        global_collinear_idx <- integer(0)
      }

      # --- Summary table ---
      if (length(globally_removed) > 0L || length(locally_masked) > 0L) {
        cat("\n--- RE Collinearity Screen Summary ---\n")
        if (length(globally_removed) > 0L) {
          cat("Globally removed from RE:\n")
          cat(paste0("  ", re_col_names[unique(globally_removed)], "\n"))
        }
        if (length(locally_masked) > 0L) {
          cat("Locally masked (group-specific):\n")
          for (lm_entry in locally_masked) {
            cat(sprintf(
              "  %s in group %d (r = %.3f)\n",
              lm_entry$var, lm_entry$group, lm_entry$cor
            ))
          }
        }
        cat("--------------------------------------\n\n")
      } else {
        message("RE collinearity screen: no problematic pairs detected.")
      }

      rm(
        pairwise_hits, pairwise_maxcor, group_cor_list,
        hit_fraction, n_valid_groups
      )
    }

    # Precompute masked coordinates by group and category to avoid repeating which() in inner loops
    masked_by_group_cat <- lapply(seq_len(n_groups), function(m) {
      idx_m <- idx_list[[m]]
      lapply(seq_len(p), function(ip) {
        j <- pp[ip]
        base_masked <- which((re_mask[, m] == 0.0 | (y_mask[ip, m] == 0.0 & !is_global_intercept)) & seq_len(k) %in% re_idx)

        if (length(base_masked) == length(re_idx) || length(idx_m) < 10L) {
          return(base_masked)
        }

        # Dynamic outcome presence threshold — use absolute threshold
        # based on mean_pixel_total to avoid degeneracy when Y is stabilized
        # (where min(Y[,j]) ≈ stabilization epsilon, making relative threshold trivially pass)
        presence_thresh <- mean_pixel_total * 0.01
        group_obs <- Y[idx_m, j]
        is_j <- (group_obs > presence_thresh)

        # Find other outcomes presence in this group
        is_base <- (rowSums(Y[idx_m, -j, drop=FALSE]) > presence_thresh)

        additional_masks <- integer(0)

        if (disable_separation_detection) {
          additional_masks <- integer(0)
        } else if (any(is_j) && any(is_base)) {
          for (v in setdiff(re_idx, base_masked)) {
            if (is_global_intercept[v]) next
            vals <- Xm_list[[m]][, v]
            vals_1 <- vals[is_j]
            vals_base <- vals[is_base]

            # 1. Quantile Boundaries (Ignore extreme 5% tails)
            q1_low <- quantile(vals_1, 0.05, names = FALSE)
            q1_high <- quantile(vals_1, 0.95, names = FALSE)
            qb_low <- quantile(vals_base, 0.05, names = FALSE)
            qb_high <- quantile(vals_base, 0.95, names = FALSE)

            # 2. Local Variance Check (Loosened threshold)
            var_1 <- var(vals_1)
            if (is.na(var_1)) var_1 <- 0.0

            # 3. Zero-Bounded Functional Separation Check
            # If 95% of one group is practically zero, but the 95th percentile of the other is substantive
            base_is_zero <- (qb_high < 1e-3)
            out_is_zero <- (q1_high < 1e-3)
            functional_separation <- (base_is_zero && q1_high > 0.05) ||
              (out_is_zero && qb_high > 0.05)

            # Trigger mask if distributions don't meaningfully overlap,
            # if variance is critically low, or if zero-bounded separation exists.
            if (q1_low >= qb_high || qb_low >= q1_high || var_1 < 1e-4 || functional_separation) {
              additional_masks <- c(additional_masks, v)
            }
          }
        } else if (any(is_j) && !any(is_base)) {
          # If there is no baseline at all in this group, every RE for this category is unidentifiable
          additional_masks <- c(additional_masks, setdiff(re_idx, base_masked))
        } else if (!any(is_j)) {
          # If the target category is effectively absent (failed presence_thresh),
          # mask all covariates to prevent MCMC drift. We leave the intercept unmasked
          # so it can cleanly drift to -Inf and correctly predict probability ~ 0.
          cov_idx <- setdiff(re_idx, base_masked)
          cov_idx <- cov_idx[!is_global_intercept[cov_idx]]
          additional_masks <- c(additional_masks, cov_idx)
        }

        sort(unique(c(base_masked, additional_masks)))
      })
    })

    curr_beta_c <- array(0, c(k, p, n_groups))
    mu_pooled <- matrix(0, k, p)
    # group_sizes already computed at line 399 — no recomputation needed

    # ---------------------------------------------------------------
    # PRECOMPUTE separation-aware masks (Issues 1-3 from Doc 10)
    # These are static (masked_by_group_cat never changes), so compute
    # once here instead of rebuilding every iteration in the Gibbs loop.
    # ---------------------------------------------------------------

    # 3D category-specific mask for precision updates
    re_mask_cube <- array(1.0, dim = c(k, p, n_groups))
    for (m_d in seq_len(n_groups)) {
      for (ip_d in seq_len(p)) {
        re_mask_cube[, ip_d, m_d] <- re_mask[, m_d]
        # NOTE: Do NOT apply masked_by_group_cat here — delta handles
        # separation-flagged cells. re_mask_cube encodes only structural
        # zeros (rank deficiency, zero variance) so that the precision
        # update in Section D can include learnable cells when delta=1.
      }
    }

    # Issue 3: sep_mask_by_ip — per-category separation mask matrices
    # for the ASIS interweaving step. Each element is a k × n_groups
    # matrix that can be applied via element-wise multiplication.
    sep_mask_by_ip <- lapply(seq_len(p), function(ip) {
      m_mat <- matrix(1.0, k, n_groups)
      for (m in seq_len(n_groups)) {
        sv <- masked_by_group_cat[[m]][[ip]]
        if (length(sv) > 0L) m_mat[sv, m] <- 0.0
      }
      m_mat
    })

    # =========================================================================
    # PHASE 1A: WLS EMPIRICAL BAYES INITIALIZATION
    # Computes fast multinomial log-odds regression to initialize mu_pooled.
    # Removes the scale-discovery phase from burn-in.
    # =========================================================================
    if (use_wls_init && is.null(init_state)) {
      cat(">>> Phase 1: WLS empirical Bayes initialization...\n")

      # Smoothed shares — avoids log(0), preserves sum-to-one property
      Y_smooth <- Y + 0.5 / p_all
      row_sums <- rowSums(Y_smooth)

      for (ip in seq_len(p)) {
        j <- pp[ip]
        log_odds <- log((Y_smooth[, j] / row_sums) / (Y_smooth[, baseline] / row_sums))
        log_odds[!is.finite(log_odds)] <- 0
        log_odds <- pmax(pmin(log_odds, 8), -8)

        # Ridged least squares (not plain lm.fit): on near-collinear blocks (e.g. focal shares)
        # an unregularized fit produces huge cancelling coefficients (range ~[-38,38]) -> an
        # unstable warm start that destabilizes the first PG draws. A small ridge keeps the init
        # bounded and well-conditioned regardless of collinearity; clamp as a final safety net.
        fit <- tryCatch(
          as.vector(solve(crossprod(X) + 1e-2 * diag(k), crossprod(X, log_odds))),
          error = function(e) rep(0, k)
        )
        mu_pooled[, ip] <- pmax(pmin(replace(fit, !is.finite(fit), 0), 8), -8)
      }

      # Overdisperse the (deterministic) WLS fixed-effect init per chain so Rhat is an HONEST
      # convergence check (not chains starting from the same point). Uniform +/- init_jitter
      # around the WLS estimate, drawn in the per-chain RNG stream (future.seed / chain_id).
      # Small (~0.01-0.1) keeps the warm start; larger gives stronger overdispersion for Rhat.
      if (init_jitter > 0) mu_pooled <- mu_pooled + matrix(runif(k * p, -init_jitter, init_jitter), k, p)

      # Warm-start curr_beta_c and z_c consistently
      z_c <- array(0, c(k, p, n_groups))
      sigma_init <- re_scale_A * 0.3 # 30% of prior scale
      for (m in seq_len(n_groups)) {
        for (ip in seq_len(p)) {
          raw_noise <- rnorm(k, 0, 1) * (re_mask[, m] > 0)
          if (length(re_complement) > 0) raw_noise[re_complement] <- 0
          sv <- masked_by_group_cat[[m]][[ip]]
          if (length(sv) > 0) raw_noise[sv] <- 0
          
          if (use_re && use_ncp) {
            z_c[, ip, m] <- raw_noise * 0.3
            curr_beta_c[, ip, m] <- mu_pooled[, ip] + re_scale_A * z_c[, ip, m]
          } else {
            curr_beta_c[, ip, m] <- mu_pooled[, ip] + sigma_init * raw_noise
          }
        }
      }
      cat(sprintf(
        "    WLS init: mu_pooled range [%.2f, %.2f]\n",
        min(mu_pooled), max(mu_pooled)
      ))
    }

    # =========================================================================
    # PHASE 2: SPIKE-AND-SLAB RE INCLUSION SETUP
    #
    # Three classes of (v, ip, m) cells:
    #   [A] Structural zeros  : re_mask[v,m]=0  -> delta permanently 0, never updated
    #   [B] Learnable cells   : separation-flagged but structurally active
    #                           -> delta starts at 0, updated each iteration
    #   [C] Free cells        : structurally active, not separation-flagged
    #                           -> delta starts at 1, updated each iteration
    # =========================================================================
    if (use_spike_slab) {
      # Identify structural zeros (per group × category)
      structural_zero_cells <- lapply(seq_len(n_groups), function(m) {
        lapply(seq_len(p), function(ip) {
          which((re_mask[, m] == 0.0 | (y_mask[ip, m] == 0.0 & !is_global_intercept)) & seq_len(k) %in% re_idx)
        })
      })

      # Derive separation-only flags by subtracting structural zeros
      # from masked_by_group_cat (avoids fragile <<- scoping)
      separation_flagged <- lapply(seq_len(n_groups), function(m) {
        lapply(seq_len(p), function(ip) {
          all_masked  <- masked_by_group_cat[[m]][[ip]]
          struct_z    <- structural_zero_cells[[m]][[ip]]
          setdiff(all_masked, struct_z)
        })
      })

      # Sketch A: severity-scaled separation prior. For each separation-flagged
      # (and now learnable) cell, recompute the 5-95% support overlap and map it to
      # a one-sided inclusion penalty (<=0) added to the delta log-odds in the Gibbs
      # loop. This refines the flat pi_inc prior so more-separated cells are excluded
      # more confidently. Default 0 (neutral); disabled via separation_as_prior=FALSE.
      if (separation_as_prior) {
        for (m in seq_len(n_groups)) {
          idx_m <- idx_list[[m]]
          presence_thresh <- mean_pixel_total * 0.01
          for (ip in seq_len(p)) {
            sv <- separation_flagged[[m]][[ip]]
            sv <- sv[!is_global_intercept[sv]]
            if (length(sv) == 0L) next
            j <- pp[ip]
            is_j    <- Y[idx_m, j] > presence_thresh
            is_base <- rowSums(Y[idx_m, -j, drop = FALSE]) > presence_thresh
            for (v in sv) {
              xv <- Xm_list[[m]][, v]
              v1 <- xv[is_j]; vb <- xv[is_base]
              if (length(v1) < 2L || length(vb) < 2L) { sep_logit[v, ip, m] <- SEP_A; next }
              q1l <- quantile(v1, 0.05, names = FALSE); q1h <- quantile(v1, 0.95, names = FALSE)
              qbl <- quantile(vb, 0.05, names = FALSE); qbh <- quantile(vb, 0.95, names = FALSE)
              inter <- min(q1h, qbh) - max(q1l, qbl)
              uni   <- max(q1h, qbh) - min(q1l, qbl) + 1e-9
              overlap <- inter / uni
              sep_logit[v, ip, m] <- min(SEP_A + SEP_B * max(overlap, 0), 0)
            }
          }
        }
      }

      # Sketch A (soft tier): classify PARTIAL separation — cells whose present- vs
      # baseline-support overlap falls in (sep_overlap_hard, sep_overlap_neutral).
      # These are NOT hard-flagged (never ASIS-masked). When separation_soft = TRUE
      # they are made learnable (delta-sampled, like sparse/collinear cells) and carry
      # the sep_logit prior, so a finite Bayes factor + moderate prior jointly decide
      # inclusion. When FALSE they are only COUNTED (observe mode) — inference unchanged.
      soft_flagged <- lapply(seq_len(n_groups), function(m) vector("list", p))
      for (m in seq_len(n_groups)) {
        idx_m <- idx_list[[m]]
        presence_thresh <- mean_pixel_total * 0.01
        for (ip in seq_len(p)) {
          soft_flagged[[m]][[ip]] <- integer(0)
          j <- pp[ip]
          is_j    <- Y[idx_m, j] > presence_thresh
          is_base <- rowSums(Y[idx_m, -j, drop = FALSE]) > presence_thresh
          if (!any(is_j) || !any(is_base)) next
          cand <- setdiff(re_idx, c(structural_zero_cells[[m]][[ip]],
                                    masked_by_group_cat[[m]][[ip]],
                                    which(is_global_intercept)))
          for (v in cand) {
            xv <- Xm_list[[m]][, v]; v1 <- xv[is_j]; vb <- xv[is_base]
            if (length(v1) < 2L || length(vb) < 2L) next
            q1l <- quantile(v1, 0.05, names = FALSE); q1h <- quantile(v1, 0.95, names = FALSE)
            qbl <- quantile(vb, 0.05, names = FALSE); qbh <- quantile(vb, 0.95, names = FALSE)
            overlap <- (min(q1h, qbh) - max(q1l, qbl)) / (max(q1h, qbh) - min(q1l, qbl) + 1e-9)
            if (overlap > sep_overlap_hard && overlap < sep_overlap_neutral) {
              soft_flagged[[m]][[ip]] <- c(soft_flagged[[m]][[ip]], v)
              if (separation_soft) sep_logit[v, ip, m] <- min(SEP_A + SEP_B * max(overlap, 0), 0)
            }
          }
        }
      }
      n_soft <- sum(vapply(seq_len(n_groups), function(m)
        sum(vapply(seq_len(p), function(ip) length(soft_flagged[[m]][[ip]]), integer(1))), integer(1)))
      cat(sprintf("Phase 2: soft partial-separation band (%.2f, %.2f) flags %d (v,ip,m) cells [%s]\n",
                  sep_overlap_hard, sep_overlap_neutral, n_soft,
                  if (separation_soft) "ACTIVE: learnable + sep_logit" else "observe-only"))

      # Learnable cells: separation-flagged + sparsity-flagged + absent category + collinear, excluding structural zeros
      learnable_cells <- lapply(seq_len(n_groups), function(m) {
        lapply(seq_len(p), function(ip) {
          sep_learnable <- separation_flagged[[m]][[ip]]
          struct_z <- structural_zero_cells[[m]][[ip]]

          # Sparsity-flagged cells (not structurally zero, not already separation-flagged)
          sparse_v <- which(sparsity_flagged[, m] & seq_len(k) %in% re_idx & !is_global_intercept)
          sparse_learnable <- setdiff(sparse_v, union(struct_z, sep_learnable))

          # Absent-category groups: non-intercept covariates have no identification
          absent_l <- integer(0)
          if (absent_category[ip, m]) {
            absent_l <- setdiff(re_idx, c(struct_z, which(is_global_intercept)))
          }

          # Collinear variables
          coll_learnable <- intersect(global_collinear_idx, re_idx)
          coll_learnable <- setdiff(coll_learnable, struct_z)

          # Soft partial-separation cells (opt-in; never structural zeros by construction)
          soft_l <- if (separation_soft) setdiff(soft_flagged[[m]][[ip]], struct_z) else integer(0)

          sort(unique(c(sep_learnable, sparse_learnable, absent_l, coll_learnable, soft_l)))
        })
      })

      n_learnable <- sum(sapply(seq_len(n_groups), function(m) {
        sum(sapply(seq_len(p), function(ip) {
          length(learnable_cells[[m]][[ip]])
        }))
      }))
      cat(sprintf(
        "Phase 2: Spike-and-slab over %d learnable (v,ip,m) cells\n",
        n_learnable
      ))

      # Initialize delta — 1 for free cells, 0 for learnable + structural zeros
      delta <- array(1L, c(k, p, n_groups))
      for (m in seq_len(n_groups)) {
        for (ip in seq_len(p)) {
          sz <- structural_zero_cells[[m]][[ip]]
          if (length(sz) > 0) delta[sz, ip, m] <- 0L
          lc <- learnable_cells[[m]][[ip]]
          if (length(lc) > 0) delta[lc, ip, m] <- 0L
          if (length(re_complement) > 0) delta[re_complement, ip, m] <- 0L
        }
      }

      # Build per-cell hyperparameters for pi_inc
      # Free cells: Beta(2,2) — symmetric, mode=0.5, mildly informative
      # Flagged cells: Beta(1,b) — mode at 0, skeptical toward inclusion
      a_pi_mat <- matrix(2.0, k, p)
      b_pi_mat <- matrix(2.0, k, p)

      if (length(global_collinear_idx) > 0) {
        a_pi_mat[global_collinear_idx, ] <- 1.0
        b_pi_mat[global_collinear_idx, ] <- 19.0  # Beta(1,19) -> 5% prior inclusion
      }

      # Separation-flagged cells
      for (m in seq_len(n_groups)) {
        for (ip in seq_len(p)) {
          sep_vars <- separation_flagged[[m]][[ip]]
          if (length(sep_vars) > 0) {
            a_pi_mat[sep_vars, ip] <- 1.0
            b_pi_mat[sep_vars, ip] <- pmax(b_pi_mat[sep_vars, ip], 19.0)  # Beta(1,19) -> 5%
          }
        }
      }

      # Sparsity-flagged cells
      for (v in re_idx) {
        if (is_global_intercept[v]) next
        for (m in seq_len(n_groups)) {
          if (sparsity_flagged[v, m]) {
            a_pi_mat[v, ] <- pmin(a_pi_mat[v, ], 1.0)
            b_pi_mat[v, ] <- pmax(b_pi_mat[v, ], 9.0)  # Beta(1,9) -> 10%
          }
        }
      }

      pi_inc <- matrix(0.0, k, p)
      for (ip in seq_len(p)) {
        for (v in seq_len(k)) {
          pi_inc[v, ip] <- a_pi_mat[v, ip] / (a_pi_mat[v, ip] + b_pi_mat[v, ip])
        }
      }
      if (length(re_complement) > 0) pi_inc[re_complement, ] <- 0

      # Storage — posterior count of delta accumulated during sampling
      if (!save_posterior_to_disk) {
        post_delta_count <- array(0L, c(k, p, n_groups))
      }
    }

    # =========================================================================
    # PHASE 4: CONDITIONAL AUTOREGRESSIVE (CAR) PRIOR SETUP
    #
    # For each (v, ip), the vector of included group deviations follows
    #   beta_dev[incl] ~ N(0, (tau_spatial * Q_car[incl,incl])^{-1})
    # where Q_car = D - car_rho * W  (proper CAR, always PD for rho in (0,1)).
    #
    # tau_spatial[v, ip] updated via Gamma conjugate from quadratic form.
    # =========================================================================
    car_active <- use_car && !is.null(country_adjacency)

    if (car_active) {
      stopifnot(
        is.matrix(country_adjacency),
        nrow(country_adjacency) == n_groups,
        ncol(country_adjacency) == n_groups
      )

      W_car <- 0.5 * (country_adjacency + t(country_adjacency)) # symmetrize
      diag(W_car) <- 0
      D_car <- diag(rowSums(W_car))
      
      degrees <- rowSums(W_car)
      isolated <- which(degrees == 0)
      D_car_safe <- D_car
      if (length(isolated) > 0) {
        diag(D_car_safe)[isolated] <- 1e-4
      }

      lambda_max_normalized <- max(eigen(
        solve(sqrt(D_car_safe)) %*% W_car %*% solve(sqrt(D_car_safe)),
        symmetric = TRUE, only.values = TRUE
      )$values)
      rho_max_safe <- 0.99 / lambda_max_normalized
      if (car_rho >= rho_max_safe) {
        warning(sprintf(
          "car_rho=%.3f is near or above stability limit %.3f. Clamping to %.3f.",
          car_rho, rho_max_safe, rho_max_safe * 0.99
        ))
        car_rho <- rho_max_safe * 0.99
      }
      Q_car <- D_car_safe - car_rho * W_car # proper CAR precision

      # Verify positive definiteness
      min_ev_Q <- min(eigen(Q_car, symmetric = TRUE, only.values = TRUE)$values)
      if (min_ev_Q <= 0) {
        stop(sprintf("CAR: Q_car not PD (min eigenvalue = %.4g). Reduce car_rho.", min_ev_Q))
      }

      # Initialize spatial precision: one scalar per (v, ip)
      tau_spatial <- matrix(1.0, k, p)
      tau_spatial[re_complement, ] <- 0 # Fixed effects have no spatial prior

      # Average connectivity for diagnostics
      n_isolated <- sum(rowSums(W_car) == 0)
      cat(sprintf(
        "Phase 4: CAR rho=%.2f | avg neighbors=%.1f | isolated groups=%d\n",
        car_rho, mean(rowSums(W_car)), n_isolated
      ))
      if (n_isolated > 0) {
        message(sprintf("CAR: %d groups have no neighbors - spatial prior inactive for them", n_isolated))
      }

      # Storage
      if (!save_posterior_to_disk) {
        post_tau_spatial <- array(0, c(k, p, niter - nburn))
      }
    } else {
      car_active <- FALSE
    }
  } else {
    # Non-RE path: disable adaptive features that require RE
    car_active <- FALSE
  }

  # Pre-compute group-size weights for RE BART intercept absorption
  if (use_re && use_bart) {
    weights <- group_sizes / sum(group_sizes)
  }

  # Initialize z_c unconditionally — Phase 1A warm-start overwrites if use_wls_init = TRUE
  if (use_re && use_ncp) {
    if (!exists("z_c", inherits = FALSE)) z_c <- array(0, c(k, p, n_groups))
  }

  # --- Horseshoe Setup ---
  if (use_horseshoe) {
    # Per-Helmert-contrast shrinkage (equation_specific) is order-dependent AND
    # mis-shapes c_v under symmetric_hs: var_eff becomes [k x p] but c_v[hs_idx] is
    # length k -> recycling ("number of items to replace is not a multiple") and only
    # the first contrast column is used. The symmetric horseshoe must shrink each
    # covariate's WHOLE zero-sum effect by one scale (group/CLR horseshoe from the
    # rotation-invariant magnitude ||beta_v||^2) = the pooled (equation_specific=FALSE)
    # path, which is order-invariant. Force it.
    if (isTRUE(symmetric_hs) && isTRUE(equation_specific_hs)) {
      equation_specific_hs <- FALSE
    }
    hs_idx <- if (is.null(horseshoe_idx)) 1:k else horseshoe_idx
    k_hs <- length(hs_idx)

    # NEW: Tau0 Calibration (Piironen-Vehtari recommended)
    # sigma_util ≈ pi/sqrt(3) for Logit link
    sigma_util <- pi / sqrt(3)
    p0_m <- p0_mu %||% max(1, round(k_hs * 0.3))

    # Ensure p0 < k_hs to avoid Inf or negative global shrinkage
    p0_m <- min(p0_m, k_hs - 0.5)

    # Effective N: pooled uses total observations (Piironen-Vehtari recommendation).
    # Using n_groups systematically under-shrinks when n_groups << n.
    N_eff_pooled <- n

    tau0_pooled <- tau0_mu %||% ((p0_m / (k_hs - p0_m)) * (sigma_util / sqrt(N_eff_pooled)))

    cat(sprintf(
      "Horseshoe Calibration: tau0_pooled=%.4f (p0=%.1f)\n",
      tau0_pooled, p0_m
    ))

    if (equation_specific_hs) {
      hs_lambda2 <- matrix(1, k, p)
      hs_nu <- matrix(1, k, p)
    } else {
      hs_lambda2 <- rep(1, k)
      hs_nu <- rep(1, k)
    }
    hs_tau2 <- if (equation_specific_hs) rep(tau0_pooled^2, p) else tau0_pooled^2
    hs_xi <- if (equation_specific_hs) rep(1, p) else 1
    hs_c2 <- slab_s2
    hs_zeta <- 1
    # Pre-allocate precision corrections for C++
    hs_prec_mat <- matrix(0, k, p)
    c_v <- numeric(k)
  } else {
    hs_prec_mat <- matrix(0, k, p)
    c_v <- numeric(k)
  }

  # =========================================================================
  # PHASE 1B: PRECISION-WEIGHTED HORSESHOE CALIBRATION
  # Computes effective information per predictor using WLS-predicted omega
  # weights. Predictors with low within-group variance get tighter tau0.
  # =========================================================================
  if (use_precision_hs && use_horseshoe && use_re && is.null(init_state)) {
    # Empirical PG weights from WLS-initialized linear predictors
    eta_init <- X %*% mu_pooled # n x p
    omega_init <- matrix(0.25, n, p) # fallback: PG(1,0) mean
    # nn_weighted already computed at L167 with guaranteed non-null y_weight — no redefinition needed
    for (ip in seq_len(p)) {
      z_ip <- eta_init[, ip]
      omega_init[, ip] <- nn_weighted * ifelse(
        abs(z_ip) < 1e-6, 0.25,
        0.5 * tanh(z_ip / 2) / pmax(abs(z_ip), 1e-8)
      )
    }

    # Effective Fisher information per predictor (pooled across categories)
    n_eff_per_var <- sapply(seq_len(k), function(v) {
      if (is_global_intercept[v]) {
        return(n)
      }
      vals <- sapply(seq_len(n_groups), function(m) {
        idx_m <- idx_list[[m]]
        if (length(idx_m) <= 1) return(0)
        om_bar <- mean(omega_init[idx_m, ], na.rm = TRUE)
        var_xv <- var(X[idx_m, v], na.rm = TRUE)
        if (is.na(var_xv) || is.na(om_bar)) return(0)
        om_bar * var_xv * length(idx_m)
      })
      sum(vals, na.rm = TRUE)
    })
    n_eff_per_var <- pmax(n_eff_per_var, 1)

    # Scale factor: sqrt(n_eff / n) -> more info = wider allowed range
    tau_v_scale <- sqrt(n_eff_per_var / n)
    tau_v_scale[!is.finite(tau_v_scale)] <- 1
    tau_v_scale <- pmax(pmin(tau_v_scale, 5), 0.02)

    # Modulate initial lambda2: data-dense predictors start less shrunk
    if (equation_specific_hs) {
      for (ip in seq_len(p)) {
        hs_lambda2[hs_idx, ip] <- tau_v_scale[hs_idx]
      }
    } else {
      hs_lambda2[hs_idx] <- tau_v_scale[hs_idx]
    }

    # Store for post-hoc reporting
    attr(hs_lambda2, "tau_v_scale") <- tau_v_scale

    cat(sprintf(
      "    Precision-weighted HS: tau_v_scale range [%.3f, %.3f]\n",
      min(tau_v_scale[hs_idx]), max(tau_v_scale[hs_idx])
    ))
  }

  if (use_bart) {
    curr_f <- matrix(0, n, p)
    bart_samplers <- list()

    # --- Symmetric (CLR / per-category) BART setup -------------------------
    # Model ONE ensemble per category f_1..f_p_all (NO basis), zero-sum across
    # categories. The likelihood uses the baseline-relative f_tilde = E %*% f with
    # the per-category centered map E[ip,j] = 1{j=cat(ip)} - 1{j=baseline}. f is
    # centered across categories on output (CLR). Order-invariant: permuting the
    # category labels just permutes the ensembles. No Helmert contrasts.
    E_bart <- g_contrib <- NULL
    n_ens <- p
    if (bart_symmetric) {
      E_bart <- matrix(0, p, p_all)
      for (ip in seq_len(p)) E_bart[ip, pp[ip]] <- 1
      E_bart[, baseline] <- E_bart[, baseline] - 1                         # p x p_all map
      g_contrib <- matrix(0, n, p_all)                                     # per-category functions
      n_ens <- p_all
    }

    for (ic in 1:n_ens) {
      ctrl <- dbarts::dbartsControl(
        n.samples = 1, n.burn = 0,
        n.trees = n_trees_bart,
        n.threads = n_threads_bart,
        keepTrees = TRUE,
        n.chains = 1, updateState = TRUE
      )
      init_resp <- if (bart_symmetric) {
        ws0 <- sapply(1:p, function(ip) kappa_weighted[, pp[ip]] / 0.25)   # n x p working response
        as.numeric((ws0 %*% E_bart[, ic]) / max(sum(E_bart[, ic]^2), 1e-8))
      } else {
        as.numeric(kappa_weighted[, pp[ic]] / 0.25)
      }
      bart_samplers[[ic]] <- dbarts::dbarts(
        X_bart_final, init_resp,
        control = ctrl, sigma = 1.0,
        tree.prior = dbarts:::cgm(base = bart_base, power = bart_power),
        node.prior = dbarts:::normal(bart_k),
        resid.prior = dbarts:::chisq(df = 1e10, quant = 0.5)
      )
    }
  }

  has_constraints <- !is.null(positive_constraints) || !is.null(negative_constraints)
  if (has_constraints) {
    pv <- diag(MASS::ginv(prior_P))
    psd <- sqrt(pmax(pv, 1e-12))
    start_lo <- rep(-Inf, k)
    start_hi <- rep(Inf, k)
    if (!is.null(positive_constraints)) start_lo[positive_constraints] <- -2 * psd[positive_constraints]
    if (!is.null(negative_constraints)) start_hi[negative_constraints] <- 2 * psd[negative_constraints]
    decay_rate <- 5
  }

  # =========================================================================
  # HOT-START OVERRIDES
  # =========================================================================
  if (!is.null(init_state)) {
    cat(">>> Overriding initial values with hot-start state...\n")

    # --- Dimension Safety Check ---
    if (use_re && !is.null(init_state$beta)) {
      expected_dim <- as.integer(c(k, p, n_groups))
      actual_dim <- dim(init_state$beta)
      if (!identical(expected_dim, actual_dim)) {
        warning(sprintf(
          paste0(
            "Hot-start SKIPPED: beta dimension mismatch.\n",
            "  Expected: [%s]\n  Got:      [%s]\n",
            "  Likely cause: data filtering changed group count or covariate set."
          ),
          paste(expected_dim, collapse = ","),
          paste(actual_dim, collapse = ",")
        ))
        init_state <- NULL # abort hot-start, proceed cold
      }
    }
    if (!is.null(init_state) && !is.null(init_state$mu)) {
      if (!identical(dim(init_state$mu), as.integer(c(k, p)))) {
        warning("Hot-start SKIPPED: mu dimension mismatch.")
        init_state <- NULL
      }
    }

    if (!is.null(init_state$beta)) {
      if (use_re) curr_beta_c <- init_state$beta else curr_beta <- init_state$beta
    }
    if (!is.null(init_state$mu)) mu_pooled <- init_state$mu
    if (!is.null(init_state$prec_beta)) prec_beta_pooled <- init_state$prec_beta
    if (!is.null(init_state$z_c)) z_c <- init_state$z_c
    if (!is.null(init_state$sigma_re)) sigma_beta_pooled <- init_state$sigma_re
    if (!is.null(init_state$a_re)) a_re <- init_state$a_re

    if (use_horseshoe && !is.null(init_state$horseshoe)) {
      hs <- init_state$horseshoe
      if (!is.null(hs$lambda2)) hs_lambda2 <- hs$lambda2
      if (!is.null(hs$tau2)) hs_tau2 <- hs$tau2
      if (!is.null(hs$nu)) hs_nu <- hs$nu
      if (!is.null(hs$xi)) hs_xi <- hs$xi
      if (!is.null(hs$zeta)) hs_zeta <- hs$zeta
      if (!is.null(hs$c2)) hs_c2 <- hs$c2
    }

    if (use_bart && !is.null(init_state$bart_states)) {
      for (ip in 1:p) {
        if (!is.null(init_state$bart_states[[ip]])) {
          bart_samplers[[ip]]$setState(init_state$bart_states[[ip]])
        }
      }
    }

    # --- Adaptive Phase Hot-Start Restoration ---
    if (use_spike_slab && use_re) {
      if (!is.null(init_state$delta)) delta <- init_state$delta
      if (!is.null(init_state$pi_inc)) pi_inc <- init_state$pi_inc
      if (!is.null(init_state$a_pi_mat)) a_pi_mat <- init_state$a_pi_mat
      if (!is.null(init_state$b_pi_mat)) b_pi_mat <- init_state$b_pi_mat
    }
    if (car_active && !is.null(init_state$tau_spatial)) {
      tau_spatial <- init_state$tau_spatial
    }
  }

  # --- 4. STORAGE & METADATA ---
  thin <- max(1L, as.integer(thin))
  nretain <- (niter - nburn) %/% thin      # thinned # of STORED draws (RAM/disk shrink ~thin x; ESS preserved for thin << autocorr-time)

  if (save_bart_to_disk || save_posterior_to_disk) {
    dir.create(disk_path, recursive = TRUE, showWarnings = FALSE)
    meta_file <- file.path(disk_path, "model_metadata.qs")
    # Always OVERWRITE metadata so it matches THIS run's batches. A stale file from a
    # prior run with a different category set caused recovery to mis-detect p_all
    # ("length of dimnames[2] not equal to array extent"). One dir = one config, so
    # rewriting every run is safe and prevents that mismatch.
    {
      meta <- list(
        cov_names = if (!is.null(cov_names_save)) cov_names_save else paste0("V", 1:k),
        cat_names = if (!is.null(colnames(Y))) colnames(Y) else paste0("C", 1:p),
        baseline_name = if (!is.null(colnames(Y))) colnames(Y)[baseline] else "Base",
        k = k,
        p = p,
        use_re = use_re,
        X_rescaling = if (do_cen || do_scl) X_rescaling else NULL,
        do_cen = do_cen,
        do_scl = do_scl,
        # The per-draw disk-save loop ALREADY back-transforms each batch to physical
        # scale (see "BACK-TRANSFORMATION (Single Sample)" below), so recovery must NOT
        # re-apply the standardization or it double-divides slopes by sd. This flag tells
        # recover_mnlogit_posterior the batches are already physical.
        batches_back_transformed = !("none" %in% methods_chosen),
        cont_idx = cont_idx,
        # --- BART metadata: lets recovery/downstream reconstruct f correctly ---
        use_bart = use_bart,
        bart_symmetric = bart_symmetric,
        bart_idx = if (use_bart) bart_idx else NULL,
        bart_clr = isTRUE(use_bart && bart_symmetric),   # symmetric BART = CLR per-category (f = centered g)
        bart_pp = if (use_bart) pp else NULL,
        p_all = p_all,
        # actual group identities in the SAME appearance order as the postb_total 3rd/slice
        # dimension (groups <- unique(group_idx)); lets recovery label the group dim with real
        # ids instead of positional "1..n", so downstream keys by identity, not position.
        group_levels = if (use_re) as.character(groups) else NULL
      )
      qs2::qs_save(meta, meta_file)
    }
  }

  if (save_posterior_to_disk) {
    postb_total <- NULL
    postb_pooled <- NULL
  } else {
    postb_total <- if (use_re) array(0, c(ncol(X_full), p_all, n_groups, nretain)) else array(0, c(ncol(X_full), p_all, nretain))
    postb_pooled <- array(0, c(ncol(X_full), p_all, nretain))
  }
  if (use_bart) {
    post_f_sum <- matrix(0, n, p_all)
    post_f_sum_sq <- matrix(0, n, p_all)
    if (store_f) post_f <- array(0, c(n, p_all, nretain))
  }
  post_log_lik <- if (save_posterior_to_disk) numeric(0) else numeric(nretain)
  post_ll_pw <- if (calc_loo && !save_posterior_to_disk) matrix(0, nretain, n) else NULL
  post_kappa_pooled <- if (use_horseshoe && !save_posterior_to_disk) matrix(0, k_hs * p, nretain) else NULL
  post_c2 <- if (use_horseshoe && !save_posterior_to_disk) numeric(nretain) else NULL
  post_sigma_re <- if (use_re && !save_posterior_to_disk) matrix(0, ncol(X_full) * p_all, nretain) else NULL
  # global horseshoe scale trace (nretain is only defined here, not at the state init above)
  post_re_tau <- if (isTRUE(use_re) && isTRUE(re_hs_global) && !save_posterior_to_disk) numeric(nretain) else NULL
  post_slab_c2 <- if (isTRUE(use_re) && isTRUE(estimate_slab_c2) && !save_posterior_to_disk) numeric(nretain) else NULL
  post_delta <- if (use_alt_spec) numeric(nretain) else NULL   # alternative-specific coefficient draws
  tree_store <- if ((store_bart_trees || save_bart_to_disk) && use_bart && !save_posterior_to_disk) vector("list", nretain) else NULL

  # --- Batched Disk Buffers ---
  if (save_bart_to_disk && use_bart && !save_posterior_to_disk) {
    bart_batch_buffer <- list()
    bart_batch_files <- character()
    c_label <- if (is.null(chain_id)) "0" else as.character(chain_id)
    old_files <- list.files(bart_disk_path, pattern = sprintf("^bart_batch_[0-9]+_chain_%s\\.qs$", c_label), full.names = TRUE)
    if (length(old_files) > 0) unlink(old_files)
  }

  if (save_posterior_to_disk) {
    posterior_batch_buffer <- list()
    posterior_batch_files <- character()
    c_label <- if (is.null(chain_id)) "0" else as.character(chain_id)
    old_files <- list.files(posterior_disk_path, pattern = sprintf("^posterior_batch_[0-9]+_chain_%s\\.qs$", c_label), full.names = TRUE)
    if (length(old_files) > 0) unlink(old_files)
  }

  # Fix #3: Explicit initialization avoids fragile exists() scoping that could
  # pick up stale variables from a prior call in the parent environment.
  postb_total_std <- NULL
  postb_pooled_std <- NULL

  # Fix #2: Precompute back-transformation indices for disk path.
  # These use dedicated names (bt_*) to avoid shadowing the BART intercept
  # int_idx set on line 312, which is still needed inside Section E.
  if (save_posterior_to_disk && (do_cen || do_scl) && sum(cont_idx) > 0) {
    bt_is_cont <- which(cont_idx)
    bt_int_idx <- which(apply(X_orig, 2, var) == 0)
  }

  # Convert pp to integer vector for C++ (already 1-based, which C++ code expects)
  pp_int <- as.integer(pp)
  group_idx_0 <- if (use_re) as.integer(match(group_idx, groups) - 1L) else integer(0)

  if (is.null(chain_id)) {
    cat("Sampling (Rcpp-accelerated)...\n")
    pb <- utils::txtProgressBar(min = 0, max = niter, style = 3)
  } else {
    cat(sprintf("Chain %d: Sampling (Rcpp-accelerated)...\n", chain_id))
  }

  # =====================================================================
  # SYMMETRIC HORSESHOE SETUP
  # =====================================================================
  if (use_horseshoe && symmetric_hs) {
    # Symmetric (CLR) coupling block in baseline coords: M = I - 11'/p_all is the
    # isotropic precision on the zero-sum subspace; the per-covariate horseshoe scale
    # c_v multiplies it. No basis (Helmert/ILR) is used — shrinkage is on the
    # rotation-invariant magnitude, so it is order-invariant by construction.
    Msym <- diag(1, p) - matrix(1 / p_all, p, p)         # p = p_all - 1
    # Kernels must NOT also apply the diagonal HS when symmetric_hs is on.
    
    # equation_specific is forced FALSE under symmetric_hs (see top of horseshoe setup)
    var_eff <- .reg_hs_var(hs_tau2 * hs_lambda2[hs_idx], hs_c2)
    c_v[hs_idx] <- 1 / pmax(var_eff, 1e-12)
    
    # Within-block symmetric shrinkage setup
    block_sym <- NULL
    if (!is.null(const_sum_blocks) && length(const_sum_blocks) > 0) {
      block_sym <- lapply(block_info, function(b) {
        Kb <- b$K
        if (Kb < 2) return(NULL)
        ret <- as.integer(unname(na.omit(full_to_active[b$all_idx[b$active]])))
        if (length(ret) < 1) return(NULL)
        Mb  <- diag(1, Kb - 1) - matrix(1 / Kb, Kb - 1, Kb - 1)
        hs_on <- if (is.null(horseshoe_idx_pre_drop)) TRUE else all(b$all_idx %in% horseshoe_idx_pre_drop)
        list(ret = ret, Kb = Kb, M = Mb, hs_on = hs_on)
      })
      block_sym <- Filter(Negate(is.null), block_sym)
    }
  } else {
    Msym <- NULL
    block_sym <- NULL
  }
  c_block <- if (!is.null(block_sym)) numeric(length(block_sym)) else numeric(0)
  block_lambda2 <- if (!is.null(block_sym)) rep(1, length(block_sym)) else numeric(0)
  block_nu <- if (!is.null(block_sym)) rep(1, length(block_sym)) else numeric(0)

  # --- CHANNEL SPLIT: two disjoint horseshoe hierarchies ------------------------------
  # hs_pool  = the covariate channel governed by c_v * Msym (across-category coupling)
  # hs_blk_ids / blk_tau2 = the compositional channel governed by c_block * kron(Msym, Mb)
  .hs_idx_safe <- if (use_horseshoe) hs_idx else integer(0)
  hs_blk_ids <- if (!is.null(block_sym)) {
    which(vapply(block_sym, function(b) isTRUE(b$hs_on), logical(1)))
  } else integer(0)
  blk_ret_all <- if (length(hs_blk_ids) > 0) {
    sort(unique(unlist(lapply(block_sym[hs_blk_ids], `[[`, "ret"))))
  } else integer(0)
  # Per-family tau is available on BOTH kernels. The compositional channel only exists where
  # block_sym does (symmetric only), so under the diagonal kernel hs_pool is simply all of hs_idx
  # and the families partition that.
  hs_grouped <- isTRUE(use_horseshoe) && !is.null(hs_groups) && length(hs_groups) > 0
  hs_split_on <- (isTRUE(use_horseshoe) && isTRUE(symmetric_hs) && isTRUE(hs_channel_split) &&
    length(blk_ret_all) > 0 && length(setdiff(.hs_idx_safe, blk_ret_all)) > 0) ||
    (hs_grouped && length(blk_ret_all) > 0)
  hs_pool <- if (hs_split_on) setdiff(.hs_idx_safe, blk_ret_all) else .hs_idx_safe
  if (hs_grouped && length(hs_pool) == 0) { hs_grouped <- FALSE }

  # Compositional channel: ONE shared scale when only hs_channel_split is on, one PER BLOCK when
  # families are requested. blk_tau2 starts at the shared tau0^2 so every variant begins identically.
  .blk_n   <- if (hs_grouped && hs_split_on && !is.null(block_sym)) length(block_sym) else 1L
  blk_tau2 <- rep(if (use_horseshoe) hs_tau2[1] else 1, max(.blk_n, 1L))
  blk_xi   <- rep(1, max(.blk_n, 1L))

  # Covariate channel: partition hs_pool into families by regex against the column names.
  hs_grp <- integer(0); hs_grp_names <- character(0)
  hs_tau2_g <- numeric(0); hs_xi_g <- numeric(0)
  if (hs_grouped) {
    # X is unnamed by the pre-loop hoist; recover the ACTIVE column names by mapping the
    # pre-drop names (cov_names_save) through full_to_active.
    .cn <- colnames(X)
    if (is.null(.cn) && !is.null(cov_names_save) && length(cov_names_save) == length(full_to_active)) {
      .cn <- rep(NA_character_, k)
      .ok <- !is.na(full_to_active)
      .cn[full_to_active[.ok]] <- cov_names_save[.ok]
    }
    if (is.null(.cn)) .cn <- paste0("V", seq_len(k))
    .cn[is.na(.cn)] <- paste0("V", which(is.na(.cn)))
    .nm <- .cn[hs_pool]
    hs_grp <- rep(NA_integer_, length(hs_pool))
    for (gi_ in seq_along(hs_groups)) {
      hit <- is.na(hs_grp) & grepl(hs_groups[[gi_]], .nm, perl = TRUE)
      hs_grp[hit] <- gi_
    }
    hs_grp_names <- names(hs_groups)
    if (is.null(hs_grp_names)) hs_grp_names <- paste0("g", seq_along(hs_groups))
    if (any(is.na(hs_grp))) {                       # unmatched columns get their own family
      hs_grp[is.na(hs_grp)] <- length(hs_grp_names) + 1L
      hs_grp_names <- c(hs_grp_names, "_rest")
    }
    .keep <- sort(unique(hs_grp))                   # drop empty families, renumber densely
    hs_grp <- match(hs_grp, .keep)
    hs_grp_names <- hs_grp_names[.keep]
    hs_tau2_g <- rep(hs_tau2[1], length(hs_grp_names))
    hs_xi_g   <- rep(1, length(hs_grp_names))
    cat(sprintf("Horseshoe PER-FAMILY tau: %d covariate families (%s) + %d const-sum block(s)\n",
                length(hs_grp_names),
                paste(sprintf("%s:%d", hs_grp_names, tabulate(hs_grp, length(hs_grp_names))),
                      collapse = ", "),
                length(hs_blk_ids)))
  } else if (hs_split_on) {
    cat(sprintf(
      "Horseshoe CHANNEL SPLIT: covariate channel %d cols, compositional channel %d cols in %d block(s)\n",
      length(hs_pool), length(blk_ret_all), length(hs_blk_ids)
    ))
  }
  
  block_id_vec <- NULL
  block_size_vec <- NULL
  if (!is.null(block_sym)) {
    block_id_vec <- rep(-1L, k)
    block_size_vec <- integer(length(block_sym))
    for (bs_id in seq_along(block_sym)) {
      bs <- block_sym[[bs_id]]
      block_id_vec[bs$ret] <- bs_id - 1L
      block_size_vec[bs_id] <- bs$Kb
      
      if (bs$hs_on && use_horseshoe && symmetric_hs) {
        # Note: Design choice - block members shrink within-block instead of across-categories, not in addition to.
        # Zeroing out c_v here ensures they don't get double-shrunk.
        c_v[bs$ret] <- 0
        var_eff_blk <- (if (hs_split_on) blk_tau2[min(bs_id, length(blk_tau2))] else hs_tau2[1]) *
          block_lambda2[bs_id]
        var_eff_blk <- .reg_hs_var(var_eff_blk, hs_c2)
        c_block[bs_id] <- 1 / max(var_eff_blk, 1e-12)
      }
    }
  }
  
  re_prec_sym <- if (is.null(re_prec_sym)) isTRUE(symmetric_hs) else isTRUE(re_prec_sym)
  # Resolve the centring default AFTER re_prec_sym is known: the centred statistic is only ever used
  # by the symmetric updater, and it is wrong there, so switch it off exactly when that updater runs.
  if (is.null(re_prec_center)) re_prec_center <- !isTRUE(re_prec_sym)
  if (isTRUE(re_prec_center) && isTRUE(re_prec_sym))
    warning("re_prec_center=TRUE with the symmetric RE updater reproduces the pre-2026-08-21 RE ",
            "collapse (variance estimated on centred deviations, REs drawn uncentred). ",
            "Use only to reproduce an old fit.", call. = FALSE)
  hs_prec_kernel <- if (use_horseshoe && symmetric_hs) matrix(0, k, p) else hs_prec_mat
  # Support-aware FIXED-effect prior: global participation-ratio shrinkage of sparse covariates'
  # fixed effects (mirrors CLR); scales the per-predictor precision c_v / hs_prec rows in-loop.
  fe_support <- fe_support_factor(X, support_prior_strength)

  # =====================================================================
  # GIBBS LOOP
  # =====================================================================
  burn_buffer <- max(50, floor(0.05 * nburn))
  # TEMPERING SETUP
  if (use_tempering) {
    cat(sprintf("Likelihood Tempering Enabled: linear from T=%.3f to 1.0 over first half of burnin (%d iters)\n", tempering_T0, floor(nburn/4)))
  }
  nburn_half <- max(1L, floor(nburn / 4))
  for (iter in 1:niter) {
    # Refresh the kernel's HS precision from the LIVE hs_prec_mat (updated at the end of the
    # previous sweep). Without this the binding above is frozen at its pre-loop value.
    if (isTRUE(hs_kernel_live) && use_horseshoe && !symmetric_hs) hs_prec_kernel <- hs_prec_mat
    if (use_tempering && iter <= nburn_half && nburn > 0) {
      temp_iter <- tempering_T0 + (1.0 - tempering_T0) * (iter / nburn_half)
    } else {
      temp_iter <- 1.0
    }
    nn_weighted_iter <- nn_weighted * temp_iter
    kappa_weighted_iter <- kappa_weighted * temp_iter
    # =====================================================================
    # --- BART Tempering ---
    if (use_bart) {
      if (iter <= bart_warmup) {
        bart_alpha <- 0
      } else if (bart_tempering && iter <= nburn) {
        # Linear ramp from 0.1 to 1.0 during remaining burn-in
        # Buffer: ensures at least 10% of burn-in is spent at full alpha
        tempering_buffer <- floor(nburn / 10)
        ramp_end <- nburn - tempering_buffer

        if (iter <= ramp_end && ramp_end > bart_warmup) {
          progress <- (iter - bart_warmup) / (ramp_end - bart_warmup)
          bart_alpha <- 0.1 + 0.9 * max(0, min(1, progress))
        } else {
          bart_alpha <- 1.0
        }
      } else {
        bart_alpha <- 1.0
      }
    }

    # --- Annealing walls ---
    if (has_constraints) {
      if (iter <= nburn) {
        df <- exp(-decay_rate * (iter / nburn))
        cur_lo <- start_lo * df
        cur_hi <- start_hi * df
        if (iter > nburn * 0.95) {
          if (!is.null(positive_constraints)) cur_lo[positive_constraints] <- 0
          if (!is.null(negative_constraints)) cur_hi[negative_constraints] <- 0
        }
      } else {
        cur_lo <- start_lo
        cur_hi <- start_hi
        if (!is.null(positive_constraints)) cur_lo[positive_constraints] <- 0
        if (!is.null(negative_constraints)) cur_hi[negative_constraints] <- 0
      }
    }

    # =================================================================
    # A. UPDATE UTILITIES + C_J  —  delegated to C++
    # =================================================================
    # The additive n x p utility channel carries BART and/or the alternative-specific term.
    f_bart_mat <- if (use_bart) bart_alpha * curr_f else matrix(0, n, p)
    if (use_alt_spec) f_bart_mat <- f_bart_mat + curr_delta * Zt

    if (use_re) {
      uc <- update_utilities_and_cj_re(
        Xm_list, idx_list, n, curr_beta_c, f_bart_mat, pp_int,
        baseline, p_all, use_bart
      )
      U <- uc$U
      c_j_mat <- uc$c_j_mat
    } else {
      # POOLED: full C++ path
      uc <- update_utilities_and_cj(
        X, curr_beta, f_bart_mat, pp_int,
        baseline, p_all, use_bart
      )
      U <- uc$U
      c_j_mat <- uc$c_j_mat
    }

    # =================================================================
    # B. POLYA-GAMMA DRAWS  —  batched in R (already compiled C under the hood)
    # =================================================================
    all_z <- as.vector(U[, pp, drop = FALSE] - c_j_mat)
    all_omega <- fast_rpg(n * p, rep(nn_weighted_iter, p), all_z, gamma_matched_pg = gamma_matched_pg)
    omega <- matrix(pmax(all_omega, 1e-6), n, p)

    # =================================================================
    # C. COEFFICIENT SAMPLING  —  delegated to C++
    # =================================================================
    if (!use_re && !has_constraints) {
      if (use_horseshoe && symmetric_hs) {
        curr_beta <- draw_beta_symhs_pooled(
          X, Xt, kappa_weighted_iter, omega, c_j_mat,
          prior_P, prior_Pb, pp_int, f_bart_mat, use_bart,
          c_v, Msym, block_sym, c_block
        )
      } else {
        # >>> FAST PATH: entire category loop in one C++ call <<<
        curr_beta <- gibbs_step_pooled(
          X, Xt, kappa_weighted_iter, omega, c_j_mat,
          prior_P, hs_prec_kernel, prior_Pb, pp_int,
          f_bart_mat, use_bart
        )
      }
    } else if (use_re && !has_constraints) {
      # >>> RE path in C++ <<<
      sigma_mat <- matrix(
        1 / sqrt(pmax(prec_beta_pooled, 1e-8)),
        nrow = k, ncol = p
      )

      if (use_ncp) {
        # NOTE (Issue 6, Doc 10): gibbs_step_re_ncp uses re_mask/y_mask only.
        # Dynamic separation masks (masked_by_group_cat) are applied one
        # iteration later in ASIS step C.5. This one-step lag is acceptable —
        # the ASIS correction is deterministic and takes effect before the
        # next iteration's utility update.
        re_res <- gibbs_step_re_ncp(
          X, Xt, kappa_weighted_iter,
          omega, c_j_mat,
          prior_P, hs_prec_kernel, prior_Pb,
          mu_pooled, sigma_mat,
          z_c, pp_int, idx_list,
          Xm_list, Xmt_list, n_groups, f_bart_mat,
          use_bart, as.integer(re_idx - 1L),
          re_mask, y_mask, as.integer(is_global_intercept),
          re_support_mat,
          if (use_half_cauchy_re && !is.null(a_re)) a_re else matrix(1, k, p),   # half-Cauchy aux for RE-scale ASIS
          isTRUE(re_asis) && use_half_cauchy_re && !is.null(a_re),                # RE-scale ASIS interweave
          # return X'Omega X + linear term per equation so the NON-RE covariates can get the symmetric
          # coupling below (this kernel draws mu one equation at a time and cannot couple them itself)
          isTRUE(use_horseshoe) && isTRUE(symmetric_hs) && length(re_complement) > 0
        )

        curr_beta_c <- re_res$beta_c
        mu_pooled <- re_res$mu
        z_c <- re_res$z_c
        P_lik_ip <- re_res$P_lik; Pb_lik_ip <- re_res$Pb_lik   # NULL unless requested above
      } else {
        re_res <- gibbs_step_re(
          X, Xt, kappa_weighted_iter, omega, c_j_mat,
          prior_P, hs_prec_kernel, prior_Pb, mu_pooled, prec_beta_pooled,
          pp_int, idx_list, Xm_list, Xmt_list,
          n_groups, f_bart_mat, use_bart,
          as.integer(re_idx - 1L), # 0-based for C++
          re_mask, y_mask, as.integer(is_global_intercept)
        )
        curr_beta_c <- re_res$beta_c
        mu_pooled <- re_res$mu
      }
    } else {
      # =================================================================
      # C. CONSTRAINED PATH
      # Branches are evaluated ONCE per iteration, not per category.
      # sigma_mat and re_complement are hoisted outside the ip loop.
      # =================================================================

      if (!use_re) {
        # -----------------------------------------------------------------
        # C1. POOLED CONSTRAINED (no RE)
        # -----------------------------------------------------------------
        for (ip in 1:p) {
          j <- pp[ip]
          om_p <- omega[, ip]
          c_j <- c_j_mat[, ip]
          target <- kappa_weighted_iter[, j] + om_p * c_j
          if (use_bart) target <- target - om_p * (bart_alpha * curr_f[, ip])

          P <- prior_P + weighted_crossprod(Xt, X, om_p)
          if (use_horseshoe) diag(P) <- diag(P) + hs_prec_mat[, ip]
          Pb <- prior_Pb[, ip] + Xt %*% target

          lo <- rep(-Inf, k)
          hi <- rep(Inf, k)
          if (!is.null(positive_constraints)) lo[positive_constraints] <- cur_lo[positive_constraints]
          if (!is.null(negative_constraints)) hi[negative_constraints] <- cur_hi[negative_constraints]

          curr_beta[, ip] <- sample_tmvn_precision_gibbs_cpp(
            P       = P,
            Pb      = Pb,
            lo      = lo,
            hi      = hi,
            init    = curr_beta[, ip],
            n_steps = if (iter <= nburn) 2L else 1L
          )
        }
      } else if (use_ncp) {
        # -----------------------------------------------------------------
        # C2. CONSTRAINED RE — NCP PATH
        # beta_c[, m] = mu + sigma * z_c[, m],  z_c ~ N(0, I)
        # Constraint beta >= L  =>  z >= (L - mu) / sigma
        # mu bound accumulated as max over groups of (L - sigma * z_g)
        # -----------------------------------------------------------------
        sigma_mat <- matrix(1 / sqrt(pmax(prec_beta_pooled, 1e-8)), nrow = k, ncol = p)


        for (ip in 1:p) {
          j <- pp[ip]
          om_p <- omega[, ip]
          c_j <- c_j_mat[, ip]
          sig <- sigma_mat[, ip]

          P_mu_acc <- matrix(0, k, k)
          Pb_mu_acc <- rep(0, k)

          # Initialise mu walls at the hard constraint bounds
          mu_lo <- rep(-Inf, k)
          mu_hi <- rep(Inf, k)
          if (!is.null(positive_constraints)) mu_lo[positive_constraints] <- cur_lo[positive_constraints]
          if (!is.null(negative_constraints)) mu_hi[negative_constraints] <- cur_hi[negative_constraints]

          for (m in 1:n_groups) {
            idx <- idx_list[[m]]
            if (length(idx) == 0) next
            Xm <- Xm_list[[m]]
            Xmt <- Xmt_list[[m]]
            om_m <- om_p[idx]

            # Localize RE scales with group-specific structural mask
            sig_m <- sig * re_mask[, m]

            # Base precision term for this group (used for both z and mu)
            Xmt_om_Xm <- Xmt %*% (Xm * om_m)

            # Working residual after removing mu contribution
            r_z <- kappa_weighted_iter[idx, j] + om_m * (c_j[idx] - Xm %*% mu_pooled[, ip])
            if (use_bart) r_z <- r_z - om_m * (bart_alpha * curr_f[idx, ip])

            # CORRECT (Issue 5, Doc 10): element-wise Hadamard, not matrix product.
            # NCP precision = I + diag(sig) %*% X'OmX %*% diag(sig)
            #               = I + outer(sig, sig) * X'OmX
            P_z <- diag(1, k) + outer(sig_m, sig_m) * Xmt_om_Xm
            Pb_z <- sig_m * as.vector(Xmt %*% r_z)

            # Correct bound transformation: z >= (L - mu) / sigma_m
            lo_z <- rep(-Inf, k)
            hi_z <- rep(Inf, k)
            if (!is.null(positive_constraints)) {
              raw_lo <- (cur_lo[positive_constraints] - mu_pooled[positive_constraints, ip]) /
                pmax(sig_m[positive_constraints], 1e-8)
              lo_z[positive_constraints] <- pmax(raw_lo, -1e6)
            }
            if (!is.null(negative_constraints)) {
              raw_hi <- (cur_hi[negative_constraints] - mu_pooled[negative_constraints, ip]) /
                pmax(sig_m[negative_constraints], 1e-8)
              hi_z[negative_constraints] <- pmin(raw_hi, 1e6)
            }

            z_c[, ip, m] <- sample_tmvn_precision_gibbs_cpp(
              P       = P_z,
              Pb      = Pb_z,
              lo      = lo_z,
              hi      = hi_z,
              init    = z_c[, ip, m],
              n_steps = if (iter <= nburn) 2L else 1L
            )

            # Recover beta_c and hard clip for one-iteration lag safety
            curr_beta_c[, ip, m] <- mu_pooled[, ip] + sig_m * z_c[, ip, m]
            if (!is.null(positive_constraints)) {
              curr_beta_c[positive_constraints, ip, m] <- pmax(
                curr_beta_c[positive_constraints, ip, m], cur_lo[positive_constraints]
              )
            }
            if (!is.null(negative_constraints)) {
              curr_beta_c[negative_constraints, ip, m] <- pmin(
                curr_beta_c[negative_constraints, ip, m], cur_hi[negative_constraints]
              )
            }

            # Explicitly pin fixed effects and masked random effects
            if (length(re_complement) > 0) {
              curr_beta_c[re_complement, ip, m] <- mu_pooled[re_complement, ip]
              z_c[re_complement, ip, m] <- 0
            }
            if (!use_spike_slab) {
              masked_re <- masked_by_group_cat[[m]][[ip]]
              if (length(masked_re) > 0) {
                curr_beta_c[masked_re, ip, m] <- mu_pooled[masked_re, ip]
                z_c[masked_re, ip, m] <- 0
              }
            }

            # Accumulate tighter mu walls from this group's z draw
            if (!is.null(positive_constraints)) {
              mu_lo[positive_constraints] <- pmax(
                mu_lo[positive_constraints],
                cur_lo[positive_constraints] - sig_m[positive_constraints] * z_c[positive_constraints, ip, m]
              )
            }
            if (!is.null(negative_constraints)) {
              mu_hi[negative_constraints] <- pmin(
                mu_hi[negative_constraints],
                cur_hi[negative_constraints] - sig_m[negative_constraints] * z_c[negative_constraints, ip, m]
              )
            }

            # z-residualised sufficient stats for mu update
            r_mu <- kappa_weighted_iter[idx, j] + om_m * (c_j[idx] - Xm %*% (sig_m * z_c[, ip, m]))
            if (use_bart) r_mu <- r_mu - om_m * (bart_alpha * curr_f[idx, ip])

            P_mu_acc <- P_mu_acc + Xmt_om_Xm
            Pb_mu_acc <- Pb_mu_acc + as.vector(Xmt %*% r_mu)
          }

          # mu draw — uses walls accumulated from ALL groups
          P_mu <- prior_P + P_mu_acc
          diag(P_mu) <- diag(P_mu) + hs_prec_mat[, ip]

          mu_pooled[, ip] <- sample_tmvn_precision_gibbs_cpp(
            P       = P_mu + diag(1e-10, k),
            Pb      = prior_Pb[, ip] + Pb_mu_acc,
            lo      = mu_lo, # accumulated wall, not raw cur_lo
            hi      = mu_hi,
            init    = mu_pooled[, ip],
            n_steps = if (iter <= nburn) 3L else 2L
          )
        }
      } else {
        # -----------------------------------------------------------------
        # C3. CONSTRAINED RE — CP PATH (fallback when use_ncp = FALSE)
        # -----------------------------------------------------------------


        for (ip in 1:p) {
          j <- pp[ip]
          om_p <- omega[, ip]
          c_j <- c_j_mat[, ip]

          resid_Pb_sum <- rep(0, k)
          P_fixed_sum <- matrix(0, k, k)

          for (m in 1:n_groups) {
            idx <- idx_list[[m]]
            if (length(idx) == 0) next
            Xm <- Xm_list[[m]]
            Xmt <- Xmt_list[[m]]
            om_m <- om_p[idx]

            target_m <- kappa_weighted_iter[idx, j] + om_m * c_j[idx]
            if (use_bart) target_m <- target_m - om_m * (bart_alpha * curr_f[idx, ip])

            P_m <- weighted_crossprod(Xmt, Xm, om_m)
            P <- diag(prec_beta_pooled[, ip], k) + P_m
            Pb <- as.vector(prec_beta_pooled[, ip] * mu_pooled[, ip]) + Xmt %*% target_m

            lo <- rep(-Inf, k)
            hi <- rep(Inf, k)
            if (!is.null(positive_constraints)) lo[positive_constraints] <- cur_lo[positive_constraints]
            if (!is.null(negative_constraints)) hi[negative_constraints] <- cur_hi[negative_constraints]

            curr_beta_c[, ip, m] <- sample_tmvn_precision_gibbs_cpp(
              P       = P,
              Pb      = Pb,
              lo      = lo,
              hi      = hi,
              init    = curr_beta_c[, ip, m],
              n_steps = if (iter <= nburn) 2L else 1L
            )
            if (length(re_complement) > 0) curr_beta_c[re_complement, ip, m] <- mu_pooled[re_complement, ip]
            if (!use_spike_slab) {
              masked_re <- masked_by_group_cat[[m]][[ip]]
              if (length(masked_re) > 0) curr_beta_c[masked_re, ip, m] <- mu_pooled[masked_re, ip]
            }

            # Residual-based accumulation for mu update
            P_fixed_sum <- P_fixed_sum + P_m
            delta_g <- curr_beta_c[, ip, m] - mu_pooled[, ip]
            resid_target <- target_m - (Xm %*% delta_g) * om_m
            resid_Pb_sum <- resid_Pb_sum + as.vector(Xmt %*% resid_target)
          }

          cur_prior_P <- prior_P
          if (use_horseshoe) {
            cur_prior_P[hs_idx, hs_idx] <- cur_prior_P[hs_idx, hs_idx] +
              diag(hs_prec_mat[hs_idx, ip], length(hs_idx))
          }

          P_mu <- cur_prior_P + P_fixed_sum + diag(1e-10, k)
          Pb_mu <- as.vector(cur_prior_P %*% prior_mu[, ip]) + resid_Pb_sum

          lo <- rep(-Inf, k)
          hi <- rep(Inf, k)
          if (!is.null(positive_constraints)) lo[positive_constraints] <- cur_lo[positive_constraints]
          if (!is.null(negative_constraints)) hi[negative_constraints] <- cur_hi[negative_constraints]

          mu_pooled[, ip] <- sample_tmvn_precision_gibbs_cpp(
            P       = P_mu,
            Pb      = Pb_mu,
            lo      = lo,
            hi      = hi,
            init    = mu_pooled[, ip],
            n_steps = if (iter <= nburn) 3L else 2L
          )
        }
      }
    }

    # =================================================================
    # C.5 ASIS INTERWEAVING (NCP -> CP -> NCP)  [modified for Phase 2]
    # =================================================================
    if (use_re && use_ncp) {
      # Vectorized integer→double conversion for delta (one R call, no per-ip allocation)
      if (use_spike_slab) {
        delta_dbl <- delta + 0.0
        dim(delta_dbl) <- c(k, p, n_groups)
      }
      if (use_horseshoe && symmetric_hs) {
        kr <- length(re_idx)
        P_blk  <- matrix(0, kr * p, kr * p)
        Pb_blk <- numeric(kr * p)
        
        lo_k <- rep(-Inf, k); hi_k <- rep(Inf, k)
        if (has_constraints) {
          if (!is.null(positive_constraints)) lo_k[positive_constraints] <- cur_lo[positive_constraints]
          if (!is.null(negative_constraints)) hi_k[negative_constraints] <- cur_hi[negative_constraints]
        }
        lo_blk <- rep(lo_k[re_idx], p); hi_blk <- rep(hi_k[re_idx], p)

        for (ip in seq_len(p)) {
          prec_re <- pmax(prec_beta_pooled[, ip], 1e-12)
          ym <- y_mask[ip, ]; re_ym <- re_mask
          if (any(!is_global_intercept))
            re_ym[!is_global_intercept, ] <- t(t(re_ym[!is_global_intercept, , drop = FALSE]) * ym)
          if (use_spike_slab) re_ym <- re_ym * delta_dbl[, ip, ] else re_ym <- re_ym * sep_mask_by_ip[[ip]]
          if (iter <= floor(nburn * 0.1)) for (m in seq_len(n_groups))
            if (absent_category[ip, m]) re_ym[!is_global_intercept, m] <- 0.0

          n_active <- rowSums(re_ym)
          beta_sum <- rowSums(matrix(curr_beta_c[, ip, ], k, n_groups) * re_ym)

          P_cp <- prior_P
          diag(P_cp)[re_idx] <- diag(P_cp)[re_idx] + n_active[re_idx] * prec_re[re_idx]
          Pb_cp <- prior_Pb[, ip]
          Pb_cp[re_idx] <- Pb_cp[re_idx] + prec_re[re_idx] * beta_sum[re_idx]

          if (length(re_complement) > 0) {
            P_cond  <- P_cp[re_idx, re_idx, drop = FALSE]
            Pb_cond <- Pb_cp[re_idx] - P_cp[re_idx, re_complement, drop = FALSE] %*% mu_pooled[re_complement, ip]
          } else { P_cond <- P_cp[re_idx, re_idx, drop = FALSE]; Pb_cond <- Pb_cp[re_idx] }

          rows <- ((ip - 1) * kr + 1):(ip * kr)
          P_blk[rows, rows] <- P_cond
          Pb_blk[rows] <- Pb_cond
        }

        # Symmetric HS coupling across equations, at each HS predictor that is an RE
        hs_re <- intersect(hs_idx, re_idx)
        for (v in hs_re) {
          if (is.na(c_v[v]) || c_v[v] <= 0) next
          r  <- match(v, re_idx)                 # position within re_idx
          ix <- r + (seq_len(p) - 1) * kr
          P_blk[ix, ix] <- P_blk[ix, ix] + c_v[v] * Msym
        }
        
        if (!is.null(block_sym)) {
          for (bs_id in seq_along(block_sym)) {
            bs <- block_sym[[bs_id]]
            if (!bs$hs_on || c_block[bs_id] <= 0) next
            
            ret_re <- intersect(bs$ret, re_idx)
            if (length(ret_re) == 0) next
            
            sub_idx <- match(ret_re, bs$ret)
            M_sub <- bs$M[sub_idx, sub_idx, drop = FALSE]
            # V3: couple columns AND categories (see draw_beta_symhs_pooled)
            ixb <- as.vector(vapply(seq_len(p), function(ip) (ip - 1) * kr + match(ret_re, re_idx),
                                    numeric(length(ret_re))))
            P_blk[ixb, ixb] <- P_blk[ixb, ixb] + c_block[bs_id] * kronecker(Msym, M_sub)
          }
        }

        if (has_constraints) {
          asis_n_steps <- if (nburn == 0L || iter > nburn) 1L else if (iter <= max(1L, floor(nburn * 0.1))) 5L else if (iter <= max(2L, floor(nburn * 0.5))) 3L else 2L
          
          mu_draw <- sample_tmvn_precision_gibbs_cpp(
            P = P_blk, Pb = Pb_blk, lo = lo_blk, hi = hi_blk,
            init = as.vector(mu_pooled[re_idx, , drop = FALSE]),
            n_steps = asis_n_steps
          )
        } else {
          mu_draw <- chol_sample_precision_cpp(P_blk, Pb_blk)
        }
        # V4 = V3 + stale-conditioning fix. Pb_lik was accumulated with the bdev implied by the
        # PRE-draw mu_R (ASIS holds beta fixed and recomputes z in the tail), so the complement
        # redraw must remove X_R %*% mu_R_OLD, not the freshly drawn value.
        mu_re_prev <- mu_pooled[re_idx, , drop = FALSE]
        mu_pooled[re_idx, ] <- matrix(mu_draw, kr, p)

        # ---- SYMMETRIC COUPLING FOR THE NON-RE COVARIATES ------------------------------------
        # The block above couples equations only at intersect(hs_idx, re_idx). Everything else was
        # drawn by gibbs_step_re_ncp, which works one equation at a time (k x k per ip) and takes a
        # DIAGONAL precision, so it cannot apply c_v * Msym. Redraw the complement jointly across
        # equations here, conditioning on the mu[re_idx] just drawn -- the mirror image of the block
        # above. `r_mu` inside the kernel subtracts the RE deviation and BART but NOT mu, so
        # Pb_lik is the correct linear term for the full mu vector and the conditioning is standard.
        # bdev depends on z_c, which the block above does not change, so Pb_lik stays valid.
        if (!is.null(P_lik_ip) && length(re_complement) > 0) {
          nc <- length(re_complement)
          Pc_blk <- matrix(0, nc * p, nc * p); Pbc_blk <- numeric(nc * p)
          for (ip in seq_len(p)) {
            P_full  <- prior_P + P_lik_ip[, , ip]
            Pb_full <- prior_Pb[, ip] + Pb_lik_ip[, ip]
            rows <- ((ip - 1) * nc + 1):(ip * nc)
            Pc_blk[rows, rows] <- P_full[re_complement, re_complement, drop = FALSE]
            Pbc_blk[rows] <- Pb_full[re_complement] -
              P_full[re_complement, re_idx, drop = FALSE] %*% mu_re_prev[, ip]   # pre-draw mu_R
          }
          for (v in intersect(hs_idx, re_complement)) {
            if (is.na(c_v[v]) || c_v[v] <= 0) next
            r <- match(v, re_complement); ix <- r + (seq_len(p) - 1) * nc
            Pc_blk[ix, ix] <- Pc_blk[ix, ix] + c_v[v] * Msym
          }
          if (!is.null(block_sym)) for (bs_id in seq_along(block_sym)) {
            bs <- block_sym[[bs_id]]
            if (!bs$hs_on || c_block[bs_id] <= 0) next
            ret_c <- intersect(bs$ret, re_complement); if (!length(ret_c)) next
            M_sub <- bs$M[match(ret_c, bs$ret), match(ret_c, bs$ret), drop = FALSE]
            # V3: couple columns AND categories (see draw_beta_symhs_pooled)
            ixc <- as.vector(vapply(seq_len(p), function(ip) (ip - 1) * nc + match(ret_c, re_complement),
                                    numeric(length(ret_c))))
            Pc_blk[ixc, ixc] <- Pc_blk[ixc, ixc] + c_block[bs_id] * kronecker(Msym, M_sub)
          }
          mu_pooled[re_complement, ] <- matrix(chol_sample_precision_cpp(Pc_blk, Pbc_blk), nc, p)
          for (ip in seq_len(p)) for (m in seq_len(n_groups))
            curr_beta_c[re_complement, ip, m] <- mu_pooled[re_complement, ip]   # no RE => beta == mu
        }

        # recover z_c / pin spiked & complement cells (unchanged from current ASIS tail)
        for (ip in seq_len(p)) {
          sig <- 1 / sqrt(pmax(prec_beta_pooled[, ip], 1e-12))
          for (m in seq_len(n_groups)) {
            z_c[re_idx, ip, m] <- (curr_beta_c[re_idx, ip, m] - mu_pooled[re_idx, ip]) / sig[re_idx]
            if (!use_spike_slab) {
              active_masked <- masked_by_group_cat[[m]][[ip]]
              if (length(active_masked) > 0) {
                z_c[active_masked, ip, m] <- 0.0
                curr_beta_c[active_masked, ip, m] <- mu_pooled[active_masked, ip]
              }
            } else {
              spiked <- which(delta[, ip, m] == 0L & seq_len(k) %in% re_idx)
              if (length(spiked) > 0) { z_c[spiked, ip, m] <- 0.0; curr_beta_c[spiked, ip, m] <- mu_pooled[spiked, ip] }
            }
            if (length(re_complement) > 0) { z_c[re_complement, ip, m] <- 0; curr_beta_c[re_complement, ip, m] <- mu_pooled[re_complement, ip] }
          }
        }
      } else {
        for (ip in seq_len(p)) {
          prec_re <- pmax(prec_beta_pooled[, ip], 1e-12)
          ym <- y_mask[ip, ]
          re_ym <- re_mask
  
          if (any(!is_global_intercept)) {
            re_ym[!is_global_intercept, ] <-
              t(t(re_ym[!is_global_intercept, , drop = FALSE]) * ym)
          }
  
          # -- MODIFIED: dynamic delta mask replaces static sep_mask_by_ip --
          if (use_spike_slab) {
            re_ym <- re_ym * delta_dbl[, ip, ]   # use precomputed, no allocation
          } else {
            re_ym <- re_ym * sep_mask_by_ip[[ip]]
          }
  
          # Delayed ASIS activation for absent-category groups
          # Allow intercept to participate once it has moved significantly negative;
          # block non-intercept covariates which have no signal when category is absent
          if (iter <= floor(nburn * 0.1)) {
            for (m in seq_len(n_groups)) {
              if (absent_category[ip, m]) {
                re_ym[!is_global_intercept, m] <- 0.0
                # Intercept is NOT zeroed out, allowing it to learn the negative baseline natively
              }
            }
          }
  
          n_active_vec_ip <- rowSums(re_ym)
          beta_ip_mat <- matrix(curr_beta_c[, ip, ], nrow = k, ncol = n_groups)
          beta_sum_active <- rowSums(beta_ip_mat * re_ym)
  
          P_cp <- prior_P
          diag(P_cp)[re_idx] <- diag(P_cp)[re_idx] +
            n_active_vec_ip[re_idx] * prec_re[re_idx]
          if (use_horseshoe) diag(P_cp) <- diag(P_cp) + hs_prec_mat[, ip]
  
          Pb_cp <- prior_Pb[, ip]
          Pb_cp[re_idx] <- Pb_cp[re_idx] +
            prec_re[re_idx] * beta_sum_active[re_idx]
  
          if (length(re_complement) > 0) {
            P_cond <- P_cp[re_idx, re_idx, drop = FALSE]
            Pb_cond <- Pb_cp[re_idx] -
              P_cp[re_idx, re_complement, drop = FALSE] %*% mu_pooled[re_complement, ip]
          } else {
            P_cond <- P_cp
            Pb_cond <- Pb_cp
          }
  
          lo_k <- rep(-Inf, k)
          hi_k <- rep(Inf, k)
          if (has_constraints) {
            if (!is.null(positive_constraints)) lo_k[positive_constraints] <- cur_lo[positive_constraints]
            if (!is.null(negative_constraints)) hi_k[negative_constraints] <- cur_hi[negative_constraints]
          }
  
          asis_n_steps <- if (nburn == 0L || iter > nburn) 1L else if (iter <= max(1L, floor(nburn * 0.1))) 5L else if (iter <= max(2L, floor(nburn * 0.5))) 3L else 2L
  
          mu_pooled[re_idx, ip] <- sample_tmvn_precision_gibbs_cpp(
            P       = P_cond,
            Pb      = Pb_cond,
            lo      = lo_k[re_idx],
            hi      = hi_k[re_idx],
            init    = mu_pooled[re_idx, ip],
            n_steps = asis_n_steps
          )
  
          sig <- 1 / sqrt(prec_re)
          for (m in seq_len(n_groups)) {
            z_c[re_idx, ip, m] <- (curr_beta_c[re_idx, ip, m] -
              mu_pooled[re_idx, ip]) / sig[re_idx]
  
            if (!use_spike_slab) {
              # Static masking path (spike-and-slab already handled this above)
              active_masked <- masked_by_group_cat[[m]][[ip]]
              if (length(active_masked) > 0) {
                z_c[active_masked, ip, m] <- 0.0
                curr_beta_c[active_masked, ip, m] <- mu_pooled[active_masked, ip]
              }
            } else {
              # Ensure spiked cells are pinned after ASIS
              spiked <- which(delta[, ip, m] == 0L & seq_len(k) %in% re_idx)
              if (length(spiked) > 0) {
                z_c[spiked, ip, m] <- 0.0
                curr_beta_c[spiked, ip, m] <- mu_pooled[spiked, ip]
              }
            }
  
            if (length(re_complement) > 0) {
              z_c[re_complement, ip, m] <- 0
              curr_beta_c[re_complement, ip, m] <- mu_pooled[re_complement, ip]
            }
          }
        }
      }
    }

    # =================================================================
    # C.4  SPIKE-AND-SLAB INCLUSION UPDATE
    # Updates delta[v, ip, m] for each learnable cell via a Laplace
    # Bayes factor, then updates pi_inc[v, ip] via Beta conjugate.
    # =================================================================
    if (use_spike_slab && use_re && use_ncp) {
      for (ip in seq_len(p)) {
        j <- pp[ip]
        om_p <- omega[, ip]

        # -- A. Update delta per learnable cell --
        for (m in seq_len(n_groups)) {
          lc <- learnable_cells[[m]][[ip]]
          if (length(lc) == 0) next

          idx_m <- idx_list[[m]]
          if (length(idx_m) < 5L) next # too few obs - skip BF

          om_m <- om_p[idx_m]

          for (v in lc) {
            if (is_global_intercept[v]) next

            # Recompute utilities to avoid BF hysteresis and W_m bias
            U_m <- as.vector(Xm_list[[m]] %*% curr_beta_c[, ip, m])
            if (use_bart) U_m <- U_m + bart_alpha * curr_f[idx_m, ip]
            
            xv_m <- Xm_list[[m]][, v]
            lp_excl_v <- U_m - xv_m * curr_beta_c[v, ip, m]
            
            # Correct Fisher info using multinomial probabilities from clean utilities
            psi_m_clean <- lp_excl_v - c_j_mat[idx_m, ip]
            p_excl <- plogis(psi_m_clean)
            W_m_v_clean <- as.vector(p_excl * (1 - p_excl))
            
            info_vm <- sum(W_m_v_clean * xv_m^2) # Fisher info for variable v in group m
            if (info_vm < 1e-8) next # uninformative - keep current delta

            # Would-be MLE for residual deviation
            resid_v <- kappa_weighted_iter[idx_m, pp[ip]] + W_m_v_clean * c_j_mat[idx_m, ip] - W_m_v_clean * lp_excl_v
            
            dev_vm <- sum(W_m_v_clean * xv_m * resid_v) / info_vm

            # Use sigma consistent with current-iteration precision (same as NCP draw)
            sigma_v <- 1 / sqrt(max(prec_beta_pooled[v, ip], 1e-8))

            # Laplace log Bayes factor: slab vs spike
            denom <- 1 + info_vm * sigma_v^2
            shrink <- info_vm * sigma_v^2 / denom
            log_BF <- 0.5 * shrink * info_vm * dev_vm^2 - 0.5 * log(denom)

            # Posterior inclusion probability
            log_odds_prior <- log(pi_inc[v, ip] + 1e-10) -
              log(1 - pi_inc[v, ip] + 1e-10) +
              (if (separation_as_prior) sep_logit[v, ip, m] else 0)
            
            p_inc <- plogis(log_odds_prior + log_BF)
            p_inc <- pmax(pmin(p_inc, 1 - 1e-6), 1e-6)

            delta[v, ip, m] <- rbinom(1L, 1L, p_inc)

            # Immediately apply: pin if spike
            if (delta[v, ip, m] == 0L) {
              curr_beta_c[v, ip, m] <- mu_pooled[v, ip]
              z_c[v, ip, m] <- 0
            }
          }
        }

        # -- B. Update pi_inc[v, ip] via Beta conjugate --
        for (v in re_idx) {
          if (is_global_intercept[v]) next

          # Only count groups where cell is learnable (not permanently masked)
          elig <- which(vapply(
            seq_len(n_groups),
            function(m) v %in% learnable_cells[[m]][[ip]],
            logical(1)
          ))
          if (length(elig) == 0) next

          n_inc <- sum(delta[v, ip, elig])
          n_elig <- length(elig)

          pi_inc[v, ip] <- rbeta(1L,
            shape1 = a_pi_mat[v, ip] + n_inc,
            shape2 = b_pi_mat[v, ip] + n_elig - n_inc
          )
        }
      }
    }

    # =================================================================
    # C.4.5  CAR JOINT Z_C CORRECTION
    # For each (v, ip) with >=2 included groups, replaces the
    # independent NCP draws with a joint Gaussian draw that
    # incorporates the spatial CAR precision matrix.
    # =================================================================
    if (car_active && use_spike_slab && use_ncp) {
      for (ip in seq_len(p)) {
        j <- pp[ip]
        om_p <- omega[, ip]

        # Precompute base residuals per group using CURRENT curr_beta_c
        # (not stale U from Section A) so the leave-v-out residual is consistent
        group_resid <- lapply(seq_len(n_groups), function(m) {
          idx_m <- idx_list[[m]]
          om_m <- om_p[idx_m]
          lp_m_curr <- as.vector(Xm_list[[m]] %*% curr_beta_c[, ip, m])
          if (use_bart && iter > bart_warmup) lp_m_curr <- lp_m_curr + bart_alpha * curr_f[idx_m, ip]
          kappa_weighted_iter[idx_m, j] + om_m * c_j_mat[idx_m, ip] - om_m * lp_m_curr
        })

        for (v in re_idx) {
          if (is_global_intercept[v]) next

          # Use sigma consistent with current-iteration precision
          sigma_v <- 1 / sqrt(max(prec_beta_pooled[v, ip], 1e-8))
          if (sigma_v < 1e-8) next

          tau_v <- tau_spatial[v, ip]
          if (tau_v < 1e-8) next

          # Groups where v is included AND structurally active
          incl <- which(delta[v, ip, ] == 1L & re_mask[v, ] > 0)
          if (length(incl) < 2) next

          Q_sub <- Q_car[incl, incl, drop = FALSE]
          r_inc <- length(incl)

          # Sufficient stats for joint z_c[v, ip, incl] draw
          Pb_joint <- numeric(r_inc)
          P_diag_data <- numeric(r_inc)

          for (i_inc in seq_along(incl)) {
            m <- incl[i_inc]
            idx_m <- idx_list[[m]]
            xv_m <- Xm_list[[m]][, v]
            om_m <- om_p[idx_m]

            # Leave-DEVIATION-out residual: add back only v's deviation from the
            # pooled mean, NOT its full coefficient. The CAR draw produces beta_dev
            # (z is centered at 0 by the prior; beta_c = mu + sigma*z), so mu_v must
            # remain a fixed offset. Adding the full curr_beta_c[v] would let the
            # drawn deviation absorb mu_v and double-count it in beta_c, biasing
            # every CAR cell by ~mu_v and destabilizing the spatial estimates.
            r_minus_v <- group_resid[[m]] +
              om_m * xv_m * (curr_beta_c[v, ip, m] - mu_pooled[v, ip])

            # In NCP space: z = beta_dev / sigma
            P_diag_data[i_inc] <- sum(om_m * xv_m^2) * sigma_v^2
            Pb_joint[i_inc] <- sigma_v * sum(xv_m * r_minus_v)
          }

          # Joint precision: data information (diagonal) + CAR prior (off-diagonal).
          # The CAR prior beta_dev ~ N(0, (tau*Q)^-1) is on beta_dev; the draw is in
          # NCP z-space (beta_dev = sigma*z), where that prior precision becomes
          # sigma^2 * tau * Q. Without the sigma^2 factor the draw is inconsistent
          # with the tau_spatial update (Section D.3, which assumes precision tau*Q
          # on beta_dev) by a factor of sigma^2, under-regularizing CAR cells. With
          # it, the recovered beta_dev draw is sigma-independent and consistent.
          P_joint <- diag(P_diag_data) + sigma_v^2 * tau_v * Q_sub

          # Joint draw from N(P_joint^{-1} Pb_joint, P_joint^{-1})
          z_new <- tryCatch(
            chol_sample_precision_cpp(P_joint, Pb_joint),
            error = function(e) z_c[v, ip, incl] # fallback: keep current
          )

          z_c[v, ip, incl] <- z_new
          for (i_inc in seq_along(incl)) {
            m <- incl[i_inc]
            curr_beta_c[v, ip, m] <- mu_pooled[v, ip] + sigma_v * z_new[i_inc]
          }
        }
      }
    }

    # =================================================================
    # NOTE: Placed BEFORE Sections D/D2 so that BART intercept absorption
    #       is reflected in the subsequent sigma/horseshoe updates.
    # =================================================================
    # =================================================================
    # ALTERNATIVE-SPECIFIC delta  —  one shared coefficient on Zt (n x p)
    # =================================================================
    # RANDOM-WALK METROPOLIS on the EXACT multinomial likelihood -- NOT a conjugate PG draw.
    # WHY NOT CONJUGATE: this sampler carries a one-vs-rest decomposition, so category j's predictor
    # is (psi_j - c_j) with c_j = log sum_{k != j} exp(psi_k). delta enters psi_k for EVERY
    # alternative, so c_j depends on delta as well; a Gibbs step that treats c_j as fixed uses the
    # effective regressor Zt_j instead of Zt_j - d c_j/d delta. Measured: that biases delta by a
    # clean factor of ~2 (true 1/2/3 -> 0.489/0.999/1.462, ratio 2.046/2.002/2.053). MH sidesteps
    # the whole issue: delta is one scalar, so two likelihood evaluations per sweep is cheap and
    # correct by construction. Step size adapts during burn-in toward a ~0.30 acceptance rate.
    if (use_alt_spec) {
      .lp_all <- function(dl) {                       # utilities (n x p) at a candidate delta
        f_add <- (if (use_bart) bart_alpha * curr_f else matrix(0, n, p)) + dl * Zt
        uu <- if (use_re) update_utilities_and_cj_re(Xm_list, idx_list, n, curr_beta_c, f_add,
                                                     pp_int, baseline, p_all, use_bart)
              else        update_utilities_and_cj(X, curr_beta, f_add, pp_int, baseline, p_all, use_bart)
        compute_loglik(uu$U, Y, as.numeric(y_weight))$total
      }
      d_prop <- curr_delta + rnorm(1, 0, delta_mh_sd)
      lp_cur <- .lp_all(curr_delta) + dnorm(curr_delta, 0, alt_spec_prior_sd, log = TRUE)
      lp_new <- .lp_all(d_prop)     + dnorm(d_prop,     0, alt_spec_prior_sd, log = TRUE)
      acc <- is.finite(lp_new) && (log(runif(1)) < (lp_new - lp_cur))
      if (acc) curr_delta <- d_prop
      delta_acc_n <- delta_acc_n + 1L; delta_acc_k <- delta_acc_k + as.integer(acc)
      if (iter <= nburn && delta_acc_n >= 50L) {      # adapt only during burn-in
        rate <- delta_acc_k / delta_acc_n
        delta_mh_sd <- delta_mh_sd * exp((rate - 0.30) * 1.5)
        delta_mh_sd <- min(max(delta_mh_sd, 1e-4), 10)
        delta_acc_n <- 0L; delta_acc_k <- 0L
      }
    }

    if (use_bart && iter > bart_warmup && !bart_symmetric) {
      for (ip in 1:p) {
        j <- pp[ip]
        lp <- if (use_re) rowSums(X * t(curr_beta_c[, ip, group_idx_0 + 1L])) else X %*% curr_beta[, ip]
        # net out the alternative-specific term so BART does not re-absorb it
        if (use_alt_spec) lp <- lp + curr_delta * Zt[, ip]

        y_star <- (kappa_weighted_iter[, j] / omega[, ip]) + c_j_mat[, ip] - lp
        # Robustness: Cap y_star to prevent exploding intercepts during burn-in
        y_star <- pmax(pmin(y_star, 100), -100)

        # NOTE: BART sampler sees 'raw' residuals, its own scale is internal.
        # We temper BARTs CONTRIBUTION to the utility, not the sampler itself.

        bart_samplers[[ip]]$setResponse(as.numeric(y_star))
        bart_samplers[[ip]]$setWeights(as.numeric(omega[, ip]))
        curr_f[, ip] <- bart_samplers[[ip]]$run(0L, 1L)$train

        # Robustness: Only perform this after some burn-in to avoid exploding intercepts from initial y_star spikes
        if (length(int_idx) > 0 && iter > burn_buffer) {
          shift_val <- mean(curr_f[, ip]) # Capture the mean BEFORE centering
          curr_f[, ip] <- curr_f[, ip] - shift_val

          # Bug 1 Fix: Only absorb the portion that actually affected the utility
          actual_utility_shift <- bart_alpha * shift_val

          if (use_re) {
            curr_beta_c[int_idx, ip, ] <- curr_beta_c[int_idx, ip, ] + actual_utility_shift
            # Robust weighted average for intercept(s)
            if (length(int_idx) == 1) {
              mu_pooled[int_idx, ip] <- sum(weights * curr_beta_c[int_idx, ip, ])
            } else {
              # Handle multiple constant columns if they exist
              for (ii in seq_along(int_idx)) {
                mu_pooled[int_idx[ii], ip] <- sum(weights * curr_beta_c[int_idx[ii], ip, ])
              }
            }
          } else {
            curr_beta[int_idx, ip] <- curr_beta[int_idx, ip] + actual_utility_shift
          }
          # Bug B Fix: Trees are calibrated to shift_val (what left curr_f),
          # not actual_utility_shift (what entered beta). Do NOT use +=.
          bart_shifts[ip, 1] <- shift_val
        } else {
          shift_val <- 0
          bart_shifts[ip, 1] <- 0
        }
      }
    }

    # --- Symmetric (CLR / per-category) BART update -----------------------
    # Fits p_all per-category ensembles f_j; the likelihood uses the baseline-
    # relative f_tilde = E %*% f. f is centered across categories (zero-sum) on
    # output/recovery. No basis -> order-invariant.
    if (use_bart && iter > bart_warmup && bart_symmetric) {
      # 1. Per-category PG working residual after the linear predictor
      Ystar <- matrix(0, n, p)
      for (ip in 1:p) {
        j <- pp[ip]
        lp <- if (use_re) rowSums(X * t(curr_beta_c[, ip, group_idx_0 + 1L])) else X %*% curr_beta[, ip]
        Ystar[, ip] <- pmax(pmin((kappa_weighted_iter[, j] / omega[, ip]) + c_j_mat[, ip] - lp, 100), -100)
      }
      # 2. Backfit each per-category ensemble f_j on its pooled pseudo-response
      #    (sufficient statistic of the conditional Gaussian => valid Gibbs update)
      for (ic in 1:n_ens) {
        r_excl <- Ystar - g_contrib %*% t(E_bart) + outer(g_contrib[, ic], E_bart[, ic])
        wd     <- sweep(omega, 2, E_bart[, ic], "*")
        den    <- pmax(as.numeric(rowSums(sweep(omega, 2, E_bart[, ic]^2, "*"))), 1e-8)
        ytilde <- as.numeric(rowSums(wd * r_excl)) / den
        bart_samplers[[ic]]$setResponse(as.numeric(ytilde))
        bart_samplers[[ic]]$setWeights(as.numeric(den))
        g_contrib[, ic] <- bart_samplers[[ic]]$run(0L, 1L)$train
      }
      # 3. Across-observation intercept absorption (per-category level vs intercept).
      #    Intercept absorbs the category-space shift E %*% m_b; stored trees record
      #    the per-ensemble mean m_b. (Across-CATEGORY zero-sum centering is a gauge
      #    applied at storage/recovery; it leaves f_tilde = E %*% f unchanged.)
      if (length(int_idx) > 0 && iter > burn_buffer) {
        m_b <- colMeans(g_contrib)                                   # length p_all
        g_contrib <- sweep(g_contrib, 2, m_b, "-")
        shift_ip <- as.numeric(E_bart %*% m_b)                       # length p
        for (ip in 1:p) {
          actual_shift <- bart_alpha * shift_ip[ip]
          if (use_re) {
            curr_beta_c[int_idx, ip, ] <- curr_beta_c[int_idx, ip, ] + actual_shift
            for (ii in seq_along(int_idx)) {
              mu_pooled[int_idx[ii], ip] <- sum(weights * curr_beta_c[int_idx[ii], ip, ])
            }
          } else {
            curr_beta[int_idx, ip] <- curr_beta[int_idx, ip] + actual_shift
          }
        }
        bart_shifts[, 1] <- m_b          # per-ensemble (per-category) mean
      } else {
        bart_shifts[, 1] <- 0
      }
      # 4. Reconstruct baseline-relative curr_f used by the rest of the sampler
      curr_f <- g_contrib %*% t(E_bart)
    }

    # =================================================================
    # D. SIGMA UPDATE  [modified for Phase 2]
    # NOTE: Placed AFTER Section E so that BART intercept absorption is
    #       reflected in curr_beta_c / mu_pooled before computing RE precision.
    # =================================================================
    if (use_re) {
      target_for_prec <- curr_beta_c

      # Build the activity mask: delta incorporates all masking when
      # spike-and-slab is active; otherwise use precomputed re_mask_cube
      if (use_spike_slab && use_re) {
        re_mask_for_prec <- array(as.numeric(delta), dim = c(k, p, n_groups))
        # Enforce hard structural zeros regardless of delta state
        # (defends against hot-start where delta could be 1 for rank-deficient cells)
        re_mask_for_prec <- re_mask_for_prec * re_mask_cube
      } else {
        re_mask_for_prec <- re_mask_cube
      }

      if (collapse_re_var) {
        # ============ collapsed RE-variance draw (per-covariate, zero-sum POOLED) ============
        # Partially-collapsed Gibbs validated against the REAL update_re_precision_hc_sym oracle
        # (experiments/mixing/exp7 mechanism + exp8 coupling/masking/support): ONE tau per covariate
        # shared across categories, with the SYMMETRIC zero-sum coupling. Per covariate v: gather the
        # per-(active category, group) PG sufficient stats (M, c) conditioning on the other covariates'
        # current REs, marginalise the zero-sum-coupled group RE block out of the variance (1-D slice
        # on log tau via the per-group eigendecomposition of the prior metric Q_m), then REFILL the
        # block immediately (van Dyk-Park). The zero-sum prior metric Q_m = I_{s} - kappa_m 11' over the
        # s active categories (kappa_m=(2 p_all - s - 1)/p_all^2) reproduces the C++ ss EXACTLY; the
        # support prior whitens it (Q_m/rs^2), matching the support-scaled RE draw + the fixed
        # support-weighted variance step. NO cap (the per-covariate pooling is self-regularising; the
        # Finnish-HS cap biased the marginal). A FRESH PG augmentation is drawn once per (category,
        # group) and reused across covariates; covariates are swept one-at-a-time. Mutually exclusive
        # with re_asis. The guard blocks it until the production collapsed==standard gate passes.
        if (!collapse_re_var_validated)
          stop("collapse_re_var is guarded. Run experiments/mixing/gate_collapse_production.R (both arms re_regularize=FALSE) on a SUBSET of THIS model, confirm collapsed==standard within the standard's own MC error for RE-SD and FE, THEN set collapse_re_var_validated = TRUE.")
        prec_new <- prec_beta_pooled; a_new_mat <- a_re
        # regularised horseshoe (consistent slice, exp9): effective precision tau_eff = tau_raw + 1/c2
        # used EVERYWHERE the RE precision appears (marginal + refill); we slice tau_raw, store tau_eff.
        # inv_c2 = 0 -> exact unregularised collapse (matches the C++ standard's Gamma branch).
        inv_c2 <- if (isTRUE(re_regularize) && collapse_slab_c2 > 0 && is.finite(collapse_slab_c2)) 1 / collapse_slab_c2 else 0
        # Active-cell predicate MUST match what the STANDARD pins, or the collapse inflates the
        # variance on cells the standard holds at mu. Beyond the re_mask_cube screen + y_mask
        # separation, the standard also pins the HARD separation mask masked_by_group_cat (via the
        # C.5 ASIS step, curr_beta_c[masked]<-mu) -- which re_mask_cube DELIBERATELY excludes
        # (line ~909). Without this term the collapse drew REs for quasi-separated cells the standard
        # holds at 0 -> runaway RE-variance (the production-gate drift). Mirror it here.
        sep_ok <- function(v, ip, m) { mm <- masked_by_group_cat[[m]][[ip]]; is.null(mm) || !(v %in% mm) }
        re_active_cell <- function(v, ip, m)
          re_mask_cube[v, ip, m] > 0.5 && (is_global_intercept[v] == 1 || y_mask[ip, m] > 0.5) && sep_ok(v, ip, m)
        # zero-sum prior metric per active-count s: Q_s = I - kappa_s 11', Ri_s = chol(Q_s)^-1
        kap <- function(s) (2 * p_all - s - 1) / p_all^2
        Ris <- vector("list", p); Qs_l <- vector("list", p)
        for (s in seq_len(p)) { Qs <- diag(s) - matrix(kap(s), s, s); Qs_l[[s]] <- Qs; Ris[[s]] <- backsolve(chol(Qs), diag(s)) }
        # 1) fresh PG augmentation + base residual ONCE per (category, group)
        om_l <- vector("list", p); base_l <- vector("list", p)
        for (ip in seq_len(p)) {
          jc <- pp[ip]; mu_ip <- mu_pooled[, ip]
          om_l[[ip]] <- vector("list", n_groups); base_l[[ip]] <- vector("list", n_groups)
          for (m in seq_len(n_groups)) { ii <- idx_list[[m]]; Xm <- Xm_list[[m]]
            eta <- as.vector(Xm %*% curr_beta_c[, ip, m])
            w <- fast_rpg(length(ii), nn_weighted_iter[ii], eta - c_j_mat[ii, ip], gamma_matched_pg = gamma_matched_pg)
            om_l[[ip]][[m]] <- w
            base_l[[ip]][[m]] <- kappa_weighted_iter[ii, jc] + w * (c_j_mat[ii, ip] - as.vector(Xm %*% mu_ip)) }
        }
        # 2) per covariate: pool one tau across categories x groups (zero-sum coupled), refill
        for (ri in seq_along(re_idx)) {
          v <- re_idx[ri]
          Ml <- vector("list", n_groups); cl <- vector("list", n_groups); Sl <- vector("list", n_groups)
          g_l <- vector("list", n_groups); w_l <- vector("list", n_groups); ntot <- 0L
          for (m in seq_len(n_groups)) {
            Xm <- Xm_list[[m]]; S <- integer(0); Mv <- numeric(0); cv <- numeric(0)
            for (ip in seq_len(p)) {
              if (!re_active_cell(v, ip, m)) next
              w <- om_l[[ip]][[m]]; mu_ip <- mu_pooled[, ip]; dev <- curr_beta_c[, ip, m] - mu_ip
              rz <- base_l[[ip]][[m]] - w * (as.vector(Xm %*% dev) - Xm[, v] * dev[v])   # residual excl cov v's RE
              S <- c(S, ip); Mv <- c(Mv, sum(w * Xm[, v]^2)); cv <- c(cv, sum(Xm[, v] * rz))
            }
            s <- length(S); Sl[[m]] <- S
            if (s == 0) next
            Ml[[m]] <- Mv; cl[[m]] <- cv
            rs2 <- max(re_support_mat[v, m]^2, 1e-12)
            Ri <- Ris[[s]] * sqrt(rs2)                                # chol(Q_s/rs2)^-1 = chol(Q_s)^-1 * rs
            B <- crossprod(Ri, Mv * Ri); eg <- eigen(B, symmetric = TRUE)  # Ri' diag(M) Ri
            g_l[[m]] <- eg$values; w_l[[m]] <- as.vector(crossprod(eg$vectors, crossprod(Ri, cv)))
            ntot <- ntot + s
          }
          if (ntot == 0L) next
          tau_prev_eff <- max(mean(prec_beta_pooled[v, ]), 1e-12)
          tau_prev_raw <- max(tau_prev_eff - inv_c2, 1e-12)
          a_v <- rexp(1, tau_prev_raw + 1 / re_scale_A^2)             # half-Cauchy aux | tau_raw
          lp <- function(lt) { te <- exp(lt) + inv_c2; ssum <- 0     # marginal in tau_eff; prior in tau_raw
            for (m in seq_len(n_groups)) { s <- length(Sl[[m]]); if (s == 0) next
              ssum <- ssum + 0.5 * s * log(te) - 0.5 * sum(log(g_l[[m]] + te)) + 0.5 * sum(w_l[[m]]^2 / (g_l[[m]] + te)) }
            ssum - 0.5 * lt - a_v * exp(lt) + lt }                    # Gamma(1/2,a_v) prior on tau_raw + log-tau jacobian
          lt <- log(tau_prev_raw); y0 <- lp(lt) - rexp(1); L <- lt - runif(1); R <- L + 1; gi <- 0L
          while (lp(L) > y0 && gi < 60L) { L <- L - 1; gi <- gi + 1L }; gi <- 0L
          while (lp(R) > y0 && gi < 60L) { R <- R + 1; gi <- gi + 1L }
          repeat { q2 <- runif(1, L, R); if (lp(q2) > y0) { lt <- q2; break }; if (q2 < log(tau_prev_raw)) L <- q2 else R <- q2 }
          tau_v <- min(max(exp(lt) + inv_c2, 1e-4), 1e6)   # tau_eff = tau_raw + 1/c2 (parity clamp; bounded by c2 when regularised)
          # refill the zero-sum RE block per group with tau_eff; inactive cells -> pinned to mu
          for (m in seq_len(n_groups)) {
            S <- Sl[[m]]; s <- length(S); if (s == 0) next
            rs2 <- max(re_support_mat[v, m]^2, 1e-12)
            Pmat <- diag(Ml[[m]], s) + (tau_v / rs2) * Qs_l[[s]]
            Lc <- tryCatch(chol(Pmat), error = function(e) chol(Pmat + diag(1e-6 * max(diag(Pmat)) + 1e-12, s)))
            d <- as.vector(backsolve(Lc, forwardsolve(t(Lc), cl[[m]])) + backsolve(Lc, rnorm(s)))
            for (a in seq_len(s)) curr_beta_c[v, S[a], m] <- mu_pooled[v, S[a]] + d[a]
          }
          for (ip in seq_len(p)) for (m in seq_len(n_groups))
            if (!re_active_cell(v, ip, m))
              curr_beta_c[v, ip, m] <- mu_pooled[v, ip]
          prec_new[v, ] <- tau_v; a_new_mat[v, ] <- a_v
        }
        a_re <- a_new_mat
        re_prec <- list(prec = prec_new, sigma = 1 / sqrt(prec_new + 1e-12), a_aux = a_re)
      } else if (use_half_cauchy_re) {
        if (re_prec_sym) {
          re_prec <- update_re_precision_hc_sym(
            beta_c = target_for_prec,
            mu_pooled = mu_pooled,
            re_idx = as.integer(re_idx - 1L),
            n_groups = n_groups,
            re_mask = re_mask_for_prec,
            y_mask = y_mask,
            is_intercept = as.integer(is_global_intercept),
            prec_prev = prec_beta_pooled,
            a_aux_prev = a_re,
            p_all = p_all,
            re_scale_A = re_tau_cur,
            block_id_opt = block_id_vec,
            block_size_opt = block_size_vec,
            re_regularize = isTRUE(re_regularize),
            slab_c2 = collapse_slab_c2,
            re_support_opt = re_support_mat,  # support-consistent variance: whiten ss by 1/rs^2 (matches the support-scaled RE draw)
            center_ss = isTRUE(re_prec_center)
          )
        } else {
          re_prec <- update_re_precision_hc(
            beta_c = target_for_prec,
            mu_pooled = mu_pooled,
            re_idx = as.integer(re_idx - 1L),
            n_groups = n_groups,
            re_mask = re_mask_for_prec,
            y_mask = y_mask,
            is_intercept = as.integer(is_global_intercept),
            prec_prev = prec_beta_pooled,
            a_aux_prev = a_re,
            re_scale_A = re_tau_cur,
            re_regularize = isTRUE(re_regularize), slab_c2 = collapse_slab_c2
          )
        }
        a_re <- re_prec$a_aux
        # ---- FULL BAYES: sample the slab c2 (1-D slice on log c2) --------------------------------
        if (isTRUE(estimate_slab_c2) && isTRUE(re_regularize) && collapse_slab_c2 > 0) {
          .inv_c2 <- 1 / collapse_slab_c2
          # prec_beta_pooled stores tau_EFF = tau_raw + 1/c2; recover the raw half-Cauchy precision
          .tau_raw <- pmax(prec_beta_pooled[re_idx, , drop = FALSE] - .inv_c2, 1e-10)
          # per-cell sufficient statistics of the group deviations, over ACTIVE cells only
          .SS <- .NN <- matrix(0, length(re_idx), p)
          for (.ip in seq_len(p)) {
            .act_ip <- if (length(y_mask)) y_mask[.ip, ] > 0.5 else rep(TRUE, n_groups)
            for (.i in seq_along(re_idx)) { .v <- re_idx[.i]
              .act <- (re_mask[.v, ] > 0.5) & (is_global_intercept[.v] == 1 | .act_ip)
              if (!any(.act)) next
              .d <- curr_beta_c[.v, .ip, ] - mu_pooled[.v, .ip]
              .SS[.i, .ip] <- sum(.d[.act]^2); .NN[.i, .ip] <- sum(.act) }
          }
          .keep <- .NN > 0
          if (any(.keep)) {
            .tr <- .tau_raw[.keep]; .ss <- .SS[.keep]; .nn <- .NN[.keep]
            .nu <- slab_df_re; .s0 <- slab_s2_re
            .lp <- function(lc) { c2 <- exp(lc); ve <- 1 / (.tr + 1 / c2)
              -0.5 * sum(.nn * log(ve) + .ss / ve) -            # N(0, v~) over active cells
                (.nu / 2 + 1) * lc - (.nu * .s0 / 2) / c2 + lc } # InvGamma prior + log-scale Jacobian
            .lt <- log(collapse_slab_c2); .y0 <- .lp(.lt) - stats::rexp(1)
            .L <- .lt - stats::runif(1); .R <- .L + 1; .g <- 0L
            while (.lp(.L) > .y0 && .g < 60L) { .L <- .L - 1; .g <- .g + 1L }; .g <- 0L
            while (.lp(.R) > .y0 && .g < 60L) { .R <- .R + 1; .g <- .g + 1L }
            .ln <- .lt
            for (.s2i in 1:60) { .q <- .L + stats::runif(1) * (.R - .L)
              if (.lp(.q) > .y0) { .ln <- .q; break }
              if (.q < .lt) .L <- .q else .R <- .q }
            collapse_slab_c2 <- min(max(exp(.ln), 1e-4), 1e6)
          }
        }
        # ---- global horseshoe scale tau (shared across ALL RE cells) --------------------
        # xi = 1/tau^2 | {a_k}, b ~ Gamma((K+1)/2, sum(a_k) + b);  b | xi ~ Exp(xi + 1/tau0^2).
        # a_k are the half-Cauchy auxiliaries the C++ just drew AT the current tau, so this is a
        # valid Gibbs sweep. Only ACTIVE RE cells contribute (finite, positive aux).
        if (isTRUE(re_hs_global)) {
          aa <- a_re[re_idx, , drop = FALSE]
          aa <- aa[is.finite(aa) & aa > 0]
          if (length(aa)) {
            re_hs_xi <- stats::rgamma(1, shape = 0.5 * (length(aa) + 1), rate = sum(aa) + re_hs_b)
            re_hs_xi <- min(max(re_hs_xi, 1e-10), 1e10)
            re_hs_b  <- stats::rexp(1, rate = re_hs_xi + 1 / (re_hs_tau0^2))
            re_tau_cur <- min(max(1 / sqrt(re_hs_xi), 1e-4), 1e3)
          }
        }
      } else {
        re_prec <- update_re_precision(
          target_for_prec, mu_pooled,
          as.integer(re_idx - 1L), n_groups,
          re_mask_for_prec, y_mask,
          as.integer(is_global_intercept),
          prior_a_re, prior_b_re
        )
      }
      prec_beta_pooled <- re_prec$prec
      sigma_beta_pooled <- re_prec$sigma
    }

    # =================================================================
    # D.3  TAU_SPATIAL UPDATE (CAR precision — Phase 4)
    # Conjugate Gamma update for spatial precision per (v, ip).
    # Uses included groups (delta=1) and the proper CAR quadratic form.
    # =================================================================
    if (car_active) {
      for (ip in seq_len(p)) {
        for (v in re_idx) {
          if (is_global_intercept[v]) next

          if (use_spike_slab) {
            incl <- which(delta[v, ip, ] == 1L & re_mask[v, ] > 0)
          } else {
            incl <- which(re_mask_cube[v, ip, ] > 0.5)
          }
          if (length(incl) < 2) next # need >=2 for spatial prior to be meaningful

          beta_dev <- curr_beta_c[v, ip, incl] - mu_pooled[v, ip]
          Q_sub <- Q_car[incl, incl, drop = FALSE]

          # Quadratic form: beta_dev' Q_sub beta_dev
          qf <- as.numeric(crossprod(beta_dev, Q_sub %*% beta_dev))
          if (!is.finite(qf) || qf < 0) qf <- 0

          # Gamma conjugate: tau | data ~ Gamma(a + |incl|/2, b + qf/2)
          shape_sp <- a_spatial + length(incl) / 2
          rate_sp <- b_spatial + qf / 2

          tau_spatial[v, ip] <- rgamma(1L, shape = shape_sp, rate = rate_sp)
          tau_spatial[v, ip] <- pmax(pmin(tau_spatial[v, ip], 1e4), 1e-4)
        }
      }
    }

    # =================================================================
    # (Zero-sum projection removed - using reference coding + block rotation instead)

    # =================================================================
    # D2. HORSESHOE UPDATE
    # NOTE: Placed AFTER Section E for the same reason as D.
    # =================================================================
    if (use_horseshoe) {
      beta_for_hs <- if (use_re) mu_pooled else curr_beta          # k x p (baseline coded)
      if (symmetric_hs) {
        # Basis-free CLR group horseshoe: shrink each covariate by its zero-sum
        # magnitude ||beta_v||^2. The pooled update uses rowSums(beta^2), which is
        # rotation-invariant, so the centered (CLR) coefficients suffice — no Helmert.
        full <- matrix(0, k, p_all); full[, pp] <- beta_for_hs
        beta_for_hs <- sweep(full, 1, rowMeans(full), "-")         # zero-sum (CLR) centering
      }
      # COVARIATE CHANNEL. Under hs_channel_split, hs_pool excludes the const-sum block
      # columns: they are governed by c_block, contribute nothing to this channel's rate,
      # and must therefore not inflate its shape either (see hs_channel_split above).
      if (hs_grouped) {
        # PER-FAMILY GLOBAL SCALES. lambda stays per-column but is driven by its OWN family's tau,
        # and each family's tau draw counts only its own columns in the shape. Same Makalic-Schmidt
        # auxiliaries; the slab stays shared (it is a cap, not a pooling scale).
        tau_col <- hs_tau2_g[hs_grp]
        hs_nu[hs_pool] <- 1 / rgamma(length(hs_pool), shape = 1, rate = 1 + 1 / hs_lambda2[hs_pool])
        ss_k <- rowSums(beta_for_hs[hs_pool, , drop = FALSE]^2)
        hs_lambda2[hs_pool] <- 1 / rgamma(length(hs_pool), shape = (p + 1) / 2,
                                          rate = 1 / hs_nu[hs_pool] + ss_k / (2 * tau_col))
        hs_lambda2[hs_pool] <- pmin(pmax(hs_lambda2[hs_pool], 1e-10), 1e4)

        ss_col <- ss_k / hs_lambda2[hs_pool]
        for (g_ in seq_along(hs_tau2_g)) {
          ii_ <- which(hs_grp == g_)
          if (!length(ii_)) next
          hs_xi_g[g_] <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / hs_tau2_g[g_])
          hs_tau2_g[g_] <- 1 / rgamma(1, shape = (length(ii_) * p + 1) / 2,
                                      rate = 1 / hs_xi_g[g_] + sum(ss_col[ii_]) / 2)
          hs_tau2_g[g_] <- min(max(hs_tau2_g[g_], 1e-10), 1e4)
        }
        if (estimate_c2) {
          hs_zeta <- rig(1, 1 + 1 / hs_c2)
          hs_c2 <- rig((slab_df + length(hs_pool) * p) / 2,
                       (slab_df * slab_s2) / 2 + sum(beta_for_hs[hs_pool, ]^2) / 2)
        } else hs_c2 <- slab_s2
        hs_tau2 <- mean(hs_tau2_g)                      # reported/stored summary only
        var_eff <- .reg_hs_var(hs_tau2_g[hs_grp] * hs_lambda2[hs_pool], hs_c2)
      } else {
      hs <- update_horseshoe(
        beta_for_hs, hs_lambda2, hs_nu, hs_tau2, hs_xi,
        hs_c2, hs_zeta, k, p, hs_pool,
        estimate_c2 = estimate_c2, slab_df = slab_df, slab_s2 = slab_s2,
        equation_specific = equation_specific_hs
      )
      hs_lambda2 <- hs$lambda2; hs_nu <- hs$nu; hs_tau2 <- hs$tau2
      hs_xi <- hs$xi; hs_c2 <- hs$c2; hs_zeta <- hs$zeta

      var_eff <- .reg_hs_var(.hs_lam2_tau2(equation_specific_hs, hs_tau2, hs_lambda2, hs_pool, length(hs_pool), p), hs_c2)
      }
        
      c_v <- numeric(k)
      c_v[hs_pool] <- if (is.matrix(var_eff)) 1 / pmax(var_eff[, 1], 1e-12) else 1 / pmax(var_eff, 1e-12)   # per-predictor precision
      if (support_prior_strength > 0) c_v[hs_pool] <- c_v[hs_pool] * fe_support[hs_pool]   # FE support-aware shrinkage (feeds the symmetric kernel)

      # Under symmetric_hs we use c_v + Msym below; the diagonal hs_prec_mat is unused.
      hs_prec_mat[hs_pool, ] <- 1 / pmax(var_eff, 1e-12)            # kept only for !symmetric_hs paths
      if (support_prior_strength > 0) hs_prec_mat[hs_pool, ] <- hs_prec_mat[hs_pool, ] * fe_support[hs_pool]
      
      if (!is.null(block_sym)) {
        # COMPOSITIONAL CHANNEL. Pass 1 draws the per-block local scales from the CURRENT
        # global scale; pass 2 (split only) draws that global scale from the new locals;
        # pass 3 rebuilds c_block. Splitting the passes leaves the RNG stream unchanged
        # when hs_channel_split = FALSE, so that path stays bit-identical.
        blk_sq <- numeric(length(block_sym))
        blk_df <- numeric(length(block_sym))
        .blk_tau_of <- function(b) if (!hs_split_on) hs_tau2[1] else blk_tau2[min(b, length(blk_tau2))]
        blk_tau_cur <- .blk_tau_of(1L)
        for (bs_id in seq_along(block_sym)) {
          bs <- block_sym[[bs_id]]
          if (!bs$hs_on) next
          
          blk_tau_cur <- .blk_tau_of(bs_id)
          gamma_ret <- beta_for_hs[bs$ret, , drop = FALSE]
          sq_norm <- sum(gamma_ret * (bs$M %*% gamma_ret))
          blk_sq[bs_id] <- sq_norm
          # effective d.o.f.: each block of Kb columns is rank Kb-1 after the const-sum
          # constraint, times p equations.
          blk_df[bs_id] <- p * (bs$Kb - 1)
          
          nu_rate <- 1 + 1 / block_lambda2[bs_id]
          block_nu[bs_id] <- 1 / rgamma(1, shape = 1, rate = nu_rate)
          
          lam_rate <- 1 / block_nu[bs_id] + sq_norm / (2 * blk_tau_cur)
          shape_lam <- (blk_df[bs_id] + 1) / 2
          block_lambda2[bs_id] <- 1 / rgamma(1, shape = shape_lam, rate = lam_rate)
          block_lambda2[bs_id] <- min(max(block_lambda2[bs_id], 1e-10), 1e4)
        }
        
        if (hs_split_on) {
          # Makalic-Schmidt half-Cauchy auxiliary for the compositional global scale(s):
          # one shared scale under hs_channel_split, one PER BLOCK under hs_groups.
          if (length(blk_tau2) > 1L) {
            for (b_ in hs_blk_ids) {
              blk_xi[b_] <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / blk_tau2[b_])
              blk_tau2[b_] <- 1 / rgamma(1, shape = (blk_df[b_] + 1) / 2,
                                         rate = 1 / blk_xi[b_] + (blk_sq[b_] / block_lambda2[b_]) / 2)
              blk_tau2[b_] <- min(max(blk_tau2[b_], 1e-10), 1e4)
            }
          } else {
            blk_xi[1] <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / blk_tau2[1])
            ss_blk <- sum(blk_sq[hs_blk_ids] / block_lambda2[hs_blk_ids])
            blk_tau2[1] <- 1 / rgamma(1, shape = (sum(blk_df[hs_blk_ids]) + 1) / 2,
                                      rate = 1 / blk_xi[1] + ss_blk / 2)
            blk_tau2[1] <- min(max(blk_tau2[1], 1e-10), 1e4)
          }
        }
        
        for (bs_id in seq_along(block_sym)) {
          bs <- block_sym[[bs_id]]
          if (!bs$hs_on) next
          
          blk_tau_cur <- .blk_tau_of(bs_id)
          var_eff_bs <- .reg_hs_var(block_lambda2[bs_id] * blk_tau_cur, hs_c2)
          c_block[bs_id] <- 1 / max(var_eff_bs, 1e-12)
          
          # We zero out the independent HS precision for these rows so they are ONLY shrunk symmetrically 
          # on the covariate axis (matching Tier 3 logic).
          hs_prec_mat[bs$ret, ] <- 0
          if (symmetric_hs) c_v[bs$ret] <- 0
        }
      }
    }

    # =================================================================
    # F. STORE POSTERIOR  —  log-likelihood via C++
    # =================================================================
    if (iter > nburn && (iter - nburn) %% thin == 0L && (iter - nburn) %/% thin <= nretain) {
      s <- (iter - nburn) %/% thin           # thinned store index (skips the utility/LL/reconstruction work for non-stored draws too)

      # Recompute utilities with end-of-iteration state for accurate LL
      f_bart_final <- if (use_bart) bart_alpha * curr_f else matrix(0, n, p)
      uc_ll <- if (use_re) {
        update_utilities_and_cj_re(Xm_list, idx_list, n, curr_beta_c, f_bart_final, pp_int, baseline, p_all, use_bart)
      } else {
        update_utilities_and_cj(X, curr_beta, f_bart_final, pp_int, baseline, p_all, use_bart)
      }

      # Compute log-likelihood (needed for both RAM and Disk)
      ll <- compute_loglik(uc_ll$U, Y, as.numeric(y_weight))

      # --- RANDOMIZED REFERENCE CODING RECONSTRUCTION ---
      if (do_block_rotation) {
        k_full <- ncol(X_full)
        mu_pooled_full <- matrix(0, k_full, p)
        curr_beta_full <- matrix(0, k_full, p)
        if (use_re) {
          curr_beta_c_full <- array(0, dim = c(k_full, p, n_groups))
          sigma_beta_pooled_full <- matrix(0, k_full, p)
        }
        
        # Populate active/dropped components
        for (i in 1:k_full) {
           a_idx <- full_to_active[i]
           if (!is.na(a_idx)) {
             mu_pooled_full[i, ] <- mu_pooled[a_idx, ]
             curr_beta_full[i, ] <- curr_beta[a_idx, ]
             if (use_re) {
               curr_beta_c_full[i, , ] <- curr_beta_c[a_idx, , ]
               sigma_beta_pooled_full[i, ] <- sigma_beta_pooled[a_idx, ]
             }
           }
        }
        
        for (b in seq_along(block_info)) {
           blk <- block_info[[b]]
           drop_col <- blk$all_idx[blk$current_drop]
           state_idx <- full_to_active[blk$all_idx[blk$active]]
           
           mu_mean <- colSums(mu_pooled[state_idx, , drop=FALSE]) / blk$K
           mu_pooled_full[drop_col, ] <- -mu_mean
           for(a_i in seq_along(state_idx)) {
             mu_pooled_full[blk$all_idx[blk$active[a_i]], ] <- mu_pooled[state_idx[a_i], ] - mu_mean
           }
           
           beta_mean <- colSums(curr_beta[state_idx, , drop=FALSE]) / blk$K
           curr_beta_full[drop_col, ] <- -beta_mean
           for(a_i in seq_along(state_idx)) {
             curr_beta_full[blk$all_idx[blk$active[a_i]], ] <- curr_beta[state_idx[a_i], ] - beta_mean
           }
           
           if (use_re) {
             beta_c_mean <- apply(curr_beta_c[state_idx, , , drop=FALSE], c(2,3), sum) / blk$K
             curr_beta_c_full[drop_col, , ] <- -beta_c_mean
             for(a_i in seq_along(state_idx)) {
               # sweep for array assignment
               curr_beta_c_full[blk$all_idx[blk$active[a_i]], , ] <- curr_beta_c[state_idx[a_i], , ] - beta_c_mean
             }
             
             sigma_mean <- colMeans(sigma_beta_pooled[state_idx, , drop=FALSE])
             sigma_beta_pooled_full[drop_col, ] <- sigma_mean
           }

           # --- UTILITY-PRESERVING COMPENSATION -------------------------------------------
           # Re-centring subtracted mean_k from every block coefficient, which subtracts
           # mean_k * rowSum(block) = mean_k * const from each class utility. Put that
           # per-class constant back into the intercept so X %*% beta_full reproduces the
           # fitted probabilities exactly (the block coefficients stay zero-sum, which is the
           # point of the reference coding). Skipped when there is no intercept (warned above).
           if (!is.na(cs_int_col) && blk$const != 0) {
             mu_pooled_full[cs_int_col, ] <- mu_pooled_full[cs_int_col, ] + mu_mean   * blk$const
             curr_beta_full[cs_int_col, ] <- curr_beta_full[cs_int_col, ] + beta_mean * blk$const
             if (use_re)
               curr_beta_c_full[cs_int_col, , ] <- curr_beta_c_full[cs_int_col, , ] + beta_c_mean * blk$const
           }
        }
        
        if (use_horseshoe) {
          if (is.matrix(hs_lambda2)) {
             hs_lambda2_full <- matrix(1, k_full, p)
             for (i in 1:k_full) {
               a_idx <- full_to_active[i]
               if (!is.na(a_idx)) hs_lambda2_full[i, ] <- hs_lambda2[a_idx, ]
             }
          } else {
             hs_lambda2_full <- rep(1, k_full)
             for (i in 1:k_full) {
               a_idx <- full_to_active[i]
               if (!is.na(a_idx)) hs_lambda2_full[i] <- hs_lambda2[a_idx]
             }
          }
          for (b in seq_along(block_info)) {
             blk <- block_info[[b]]
             drop_col <- blk$all_idx[blk$current_drop]
             state_idx <- full_to_active[blk$all_idx[blk$active]]
             if (is.matrix(hs_lambda2)) {
                hs_lambda2_full[drop_col, ] <- colMeans(hs_lambda2[state_idx, , drop=FALSE])
             } else {
                hs_lambda2_full[drop_col] <- mean(hs_lambda2[state_idx])
             }
          }
        }
      } else {
        mu_pooled_full <- mu_pooled
        curr_beta_full <- curr_beta
        if (use_re) {
           curr_beta_c_full <- curr_beta_c
           sigma_beta_pooled_full <- sigma_beta_pooled
        }
        if (use_horseshoe) hs_lambda2_full <- hs_lambda2
      }
      
      # Expand to P_all zero-sum dimensions so columns remain perfectly aligned 
      # regardless of which category is currently acting as the rotated baseline.
      expand_to_zs <- function(mat, pp_idx, p_all_dim) {
         full <- matrix(0, nrow(mat), p_all_dim)
         full[, pp_idx] <- mat
         sweep(full, 1, rowMeans(full), "-")
      }
      
      mu_zs <- expand_to_zs(mu_pooled_full, pp, p_all)
      beta_zs <- expand_to_zs(curr_beta_full, pp, p_all)
      
      if (use_re) {
         beta_c_zs <- array(0, c(ncol(X_full), p_all, n_groups))
         beta_c_zs[, pp, ] <- curr_beta_c_full
         for (m in 1:n_groups) {
            beta_c_zs[, , m] <- sweep(beta_c_zs[, , m, drop=FALSE], 1, rowMeans(beta_c_zs[, , m, drop=FALSE]), "-")
         }
         # Standard deviation computed directly over groups on the zero-sum subspace
         sigma_zs <- apply(beta_c_zs, c(1, 2), sd)
      }

      if (!save_posterior_to_disk) {
        if (use_re) {
          postb_total[, , , s] <- beta_c_zs
          postb_pooled[, , s] <- mu_zs
          post_sigma_re[, s] <- as.vector(sigma_zs)
          if (!is.null(post_re_tau)) post_re_tau[s] <- re_tau_cur   # global HS scale trace
          if (!is.null(post_slab_c2)) post_slab_c2[s] <- collapse_slab_c2
        } else {
          postb_total[, , s] <- beta_zs
          postb_pooled[, , s] <- beta_zs
        }

        # --- STORE-TIME SELF-CHECK (localizes coef/log-lik divergence) ---------
        # Recompute the log-lik from the EXACT stored zero-sum array using the SAME
        # per-group X the likelihood used. If this diverges from ll$total, the
        # zero-sum expansion / RE reconstruction (not the back-transform, not any
        # downstream alignment) is the corrupting step.
        if (isTRUE(.selfcheck_store) && s <= .selfcheck_n) {
          U_chk <- matrix(0, n, p_all)
          if (use_re) {
            for (m in seq_len(n_groups)) {
              ii <- idx_list[[m]]
              U_chk[ii, ] <- Xm_list[[m]] %*% beta_c_zs[, , m]
            }
          } else {
            U_chk <- X %*% beta_zs
          }
          ll_chk <- compute_loglik(U_chk, Y, as.numeric(y_weight))
          # also the pooled-only reconstruction (ignores RE) for reference
          U_pool <- if (use_re) {
            up <- matrix(0, n, p_all)
            for (m in seq_len(n_groups)) up[idx_list[[m]], ] <- Xm_list[[m]] %*% mu_zs
            up
          } else X %*% beta_zs
          ll_pool <- compute_loglik(U_pool, Y, as.numeric(y_weight))
          cat(sprintf("[SELFCHECK s=%d] ll$total=%.2f | from-STORED beta_zs=%.2f (d=%.2f) | pooled-only=%.2f\n",
                      s, ll$total, ll_chk$total, ll_chk$total - ll$total, ll_pool$total))
        }

        if (use_horseshoe) {
          # Level 1: Pooled Mean Kappa
          kappa_pooled <- 1 / (1 + .hs_lam2_tau2(equation_specific_hs, hs_tau2, hs_lambda2, hs_idx, k_hs, p))
          post_kappa_pooled[, s] <- as.vector(kappa_pooled)
          post_c2[s] <- hs_c2
        }
        post_log_lik[s] <- ll$total
        if (use_alt_spec) post_delta[s] <- curr_delta
        if (calc_loo) post_ll_pw[s, ] <- ll$pointwise

        # Adaptive Phase Storage
        if (use_spike_slab && use_re && store_delta && !save_posterior_to_disk) {
          post_delta_count <- post_delta_count + delta
        }

        if (car_active) {
          post_tau_spatial[, , s] <- tau_spatial
        }
      }
      if (use_bart) {
        if (bart_symmetric) {
          # CLR storage: per-category functions centered across categories (zero-sum).
          # Trees store the raw per-category g_j; recovery applies the same centering.
          f_baseline <- sweep(g_contrib, 1, rowMeans(g_contrib), "-")
        } else {
          # Store curr_f in baseline-coded form (baseline = 0) to match likelihood
          f_baseline <- matrix(0, n, p_all)
          f_baseline[, pp] <- curr_f
        }

        post_f_sum <- post_f_sum + bart_alpha * f_baseline
        post_f_sum_sq <- post_f_sum_sq + (bart_alpha * f_baseline)^2
        if (store_f) post_f[, , s] <- bart_alpha * f_baseline
        if (store_bart_trees || save_bart_to_disk || save_posterior_to_disk) {
          state_m <- vector("list", n_ens)
          for (ip in 1:n_ens) {
            if (do_slim_trees) {
              # SLIM TREES: Extract structures and calibrate
              df <- bart_samplers[[ip]]$getTrees()

              # 1. CALIBRATION: Recover hidden dbarts Y-scaling and offset
              # (Uses standardized trees and empirical reference points)
              s1 <- as.numeric(bart_samplers[[ip]]$predict(bart_ref1_mat))
              s2 <- as.numeric(bart_samplers[[ip]]$predict(bart_ref2_mat))
              t1 <- as.numeric(predict_slim_bart_cpp(bart_ref1_mat, df))
              t2 <- as.numeric(predict_slim_bart_cpp(bart_ref2_mat, df))

              t_diff <- t1 - t2
              y_scale <- if (abs(t_diff) > 1e-10) (s1 - s2) / t_diff else 1.0

              y_offset <- s1 - (t1 * y_scale)

              # 2. TRANSFORM SPLIT POINTS: Move to raw scale (if standardized)
              if (!is.null(bart_scaling)) {
                split_rows <- which(df$var > 0)
                if (length(split_rows) > 0) {
                  v_indices <- df$var[split_rows]
                  df$value[split_rows] <- df$value[split_rows] * bart_scaling$sd[v_indices] + bart_scaling$mu[v_indices]
                }
              }

              # 3. APPLY CALIBRATION: Scale leaves and bake in offset
              df$value[df$var == -1] <- df$value[df$var == -1] * y_scale
              total_offset <- y_offset - as.numeric(bart_shifts[ip, 1])
              df$value[df$tree == 1 & df$var == -1] <- df$value[df$tree == 1 & df$var == -1] + total_offset

              # LIVE CHECK: test_pred should match target_pred at ref1
              test_pred <- as.numeric(predict_slim_bart_cpp(bart_ref1_mat, df))
              target_pred <- s1 - as.numeric(bart_shifts[ip, 1])

              if (abs(test_pred - target_pred) > 1e-5) {
                warning(sprintf("Calibration mismatch in sample %d cat %d: Slim=%.4f, Target=%.4f", s, ip, test_pred, target_pred))
              }

              state_m[[ip]] <- df
            } else {
              # FULL STORAGE: Serialized sampler state
              invisible(bart_samplers[[ip]]$state)
              state_m[[ip]] <- serialize(bart_samplers[[ip]], NULL)
            }
          }

          if (save_bart_to_disk && !save_posterior_to_disk) {
            # Add to buffer and write in batches
            bart_batch_buffer[[length(bart_batch_buffer) + 1]] <- state_m
            if (length(bart_batch_buffer) >= bart_batch_size || s == nretain) {
              b_idx <- length(bart_batch_files) + 1
              c_label <- if (is.null(chain_id)) "0" else as.character(chain_id)
              f_name <- file.path(bart_disk_path, sprintf("bart_batch_%d_chain_%s.qs", b_idx, c_label))
              qs2::qs_save(bart_batch_buffer, f_name, compress_level = 1)
              bart_batch_files[b_idx] <- f_name
              # Store the filename in tree_store for the first sample of this batch
              tree_store[[s - length(bart_batch_buffer) + 1]] <- f_name
              bart_batch_buffer <- list() # Clear RAM
            }
          } else if (!save_posterior_to_disk) {
            # RAM Storage
            tree_store[[s]] <- state_m
          }
        }
      }

      if (save_posterior_to_disk) {
        # Construct the whole state sample
        state_sample <- list(
          iter = iter,
          beta = if (use_re) curr_beta_c_full else curr_beta_full,
          mu = if (use_re) mu_pooled_full else curr_beta_full,
          sigma_re = if (use_re) sigma_beta_pooled_full else NULL,
          log_lik = ll$total,
          log_lik_pw = if (calc_loo) ll$pointwise else NULL,
          horseshoe = if (use_horseshoe) {
            list(
              lambda2 = hs_lambda2_full, tau2 = hs_tau2,
              c2 = hs_c2,
              nu = hs_nu, xi = hs_xi, zeta = hs_zeta
            )
          } else {
            NULL
          },
          bart_trees = if (use_bart) state_m else NULL
        )

        # BACK-TRANSFORMATION (Single Sample)
        if (!("none" %in% methods_chosen)) {
          # 1. Reverse QR
          if (do_qr) {
            if (use_re) {
              for (m in 1:n_groups) {
                state_sample$beta[, , m] <- backsolve(X_R_mat, state_sample$beta[, , m, drop = FALSE])
              }
              state_sample$mu <- backsolve(X_R_mat, state_sample$mu, drop = FALSE)
            } else {
              state_sample$beta <- backsolve(X_R_mat, state_sample$beta)
              state_sample$mu <- state_sample$beta
            }
          }
          # 2. Reverse Standardization (uses precomputed bt_is_cont / bt_int_idx
          #    from line ~579 to avoid shadowing the BART int_idx on line 312)
          if ((do_cen || do_scl) && sum(cont_idx) > 0) {
            if (use_re) {
              if (do_scl) {
                for (ci in seq_along(bt_is_cont)) {
                  state_sample$beta[bt_is_cont[ci], , ] <- state_sample$beta[bt_is_cont[ci], , ] / X_rescaling[2, ci]
                  state_sample$mu[bt_is_cont[ci], ] <- state_sample$mu[bt_is_cont[ci], ] / X_rescaling[2, ci]
                }
              }
              if (do_cen && length(bt_int_idx) > 0) {
                for (m in 1:n_groups) {
                  shift <- colSums(state_sample$beta[bt_is_cont, , m, drop = FALSE] * X_rescaling[1, ])
                  for (ii in bt_int_idx) state_sample$beta[ii, , m] <- state_sample$beta[ii, , m] - shift
                }
                shift_p <- colSums(state_sample$mu[bt_is_cont, , drop = FALSE] * X_rescaling[1, ])
                for (ii in bt_int_idx) state_sample$mu[ii, ] <- state_sample$mu[ii, ] - shift_p
              }
            } else {
              if (do_scl) {
                for (ci in seq_along(bt_is_cont)) state_sample$beta[bt_is_cont[ci], ] <- state_sample$beta[bt_is_cont[ci], ] / X_rescaling[2, ci]
              }
              if (do_cen && length(bt_int_idx) > 0) {
                shift <- colSums(state_sample$beta[bt_is_cont, , drop = FALSE] * X_rescaling[1, ])
                for (ii in bt_int_idx) state_sample$beta[ii, ] <- state_sample$beta[ii, ] - shift
              }
              state_sample$mu <- state_sample$beta
            }
          }
        }

        # Add to buffer
        posterior_batch_buffer[[length(posterior_batch_buffer) + 1]] <- state_sample
        if (length(posterior_batch_buffer) >= posterior_batch_size || s == nretain) {
          b_idx <- length(posterior_batch_files) + 1
          c_label <- if (is.null(chain_id)) "0" else as.character(chain_id)
          f_name <- file.path(posterior_disk_path, sprintf("posterior_batch_%d_chain_%s.qs", b_idx, c_label))
          qs2::qs_save(posterior_batch_buffer, f_name, compress_level = 1)
          posterior_batch_files[b_idx] <- f_name
          posterior_batch_buffer <- list() # Clear RAM
        }
      }

      if (!save_posterior_to_disk) {
        # ll was already computed and stored at the top of the loop (F. STORE POSTERIOR)
      }
    }

    if (!is.null(progress_cb)) {
      if (iter %% 10 == 0 || iter == niter) {
        phase <- if (iter <= nburn) "[Burn-in]" else "[Sampling]"
        p_msg <- sprintf("Chain %s: Iteration %d / %d %s", ifelse(is.null(chain_id), "?", chain_id), iter, niter, phase)
        progress_cb(message = p_msg)
      } else {
        progress_cb() # advance without changing message
      }
    } else if (is.null(chain_id)) {
      utils::setTxtProgressBar(pb, iter)
    } else if (iter %% 100 == 0 || iter == niter) {
      phase <- if (iter <= nburn) "[Burn-in]" else "[Sampling]"
      cat(sprintf("Chain %d: Iteration %d / %d %s\n", chain_id, iter, niter, phase))
    }
  }

  # --- Flush remaining BART ensembles ---
  if (save_bart_to_disk && use_bart && !save_posterior_to_disk && length(bart_batch_buffer) > 0) {
    b_idx <- length(bart_batch_files) + 1
    c_label <- if (is.null(chain_id)) "0" else as.character(chain_id)
    f_name <- file.path(bart_disk_path, sprintf("bart_batch_%d_chain_%s.qs", b_idx, c_label))
    qs2::qs_save(bart_batch_buffer, f_name, compress_level = 1)
    bart_batch_files[b_idx] <- f_name
    bart_batch_buffer <- list()
  }

  # --- Flush remaining Posterior samples ---
  if (save_posterior_to_disk && length(posterior_batch_buffer) > 0) {
    b_idx <- length(posterior_batch_files) + 1
    c_label <- if (is.null(chain_id)) "0" else as.character(chain_id)
    f_name <- file.path(posterior_disk_path, sprintf("posterior_batch_%d_chain_%s.qs", b_idx, c_label))
    qs2::qs_save(posterior_batch_buffer, f_name, compress_level = 1)
    posterior_batch_files[b_idx] <- f_name
    posterior_batch_buffer <- list()
  }

  if (is.null(chain_id)) {
    close(pb)
    cat("\nSampling finished.\n")
  } else {
    cat(sprintf("Chain %d: Sampling finished.\n", chain_id))
  }

  # --- Capture and Save Final State for Hot-Start ---
  horseshoe_state <- NULL
  if (use_horseshoe) {
    horseshoe_state <- list(
      lambda2 = hs_lambda2, tau2 = hs_tau2, nu = hs_nu, xi = hs_xi, zeta = hs_zeta,
      c2 = hs_c2,
      blk_tau2 = blk_tau2, blk_xi = blk_xi, block_lambda2 = block_lambda2, block_nu = block_nu
    )
  }

  final_state_raw <- list(
    beta = if (use_re) curr_beta_c else curr_beta,
    mu = if (use_re) mu_pooled else NULL,
    prec_beta = if (use_re) prec_beta_pooled else NULL,
    z_c = if (use_re && use_ncp) z_c else NULL,
    sigma_re = if (use_re) sigma_beta_pooled else NULL,
    a_re = if (use_re && use_half_cauchy_re) a_re else NULL,
    horseshoe = horseshoe_state
  )
  if (use_bart) {
    final_state_raw$bart_states <- lapply(bart_samplers, function(s) s$state)
  }

  # Adaptive Phase final state
  final_state_raw$delta <- if (use_spike_slab && use_re) delta else NULL
  final_state_raw$pi_inc <- if (use_spike_slab && use_re) pi_inc else NULL
  final_state_raw$a_pi_mat <- if (use_spike_slab && use_re) a_pi_mat else NULL
  final_state_raw$b_pi_mat <- if (use_spike_slab && use_re) b_pi_mat else NULL
  final_state_raw$tau_spatial <- if (car_active) tau_spatial else NULL

  if (save_posterior_to_disk) {
    c_label <- if (is.null(chain_id)) "1" else as.character(chain_id)
    f_name_fs <- file.path(posterior_disk_path, sprintf("final_state_chain_%s.qs", c_label))
    qs2::qs_save(final_state_raw, f_name_fs, compress_level = 1)
  }

  # --- 5. BACK-TRANSFORMATION (Only if not already done via Disk path) ---
  if (!save_posterior_to_disk && !("none" %in% methods_chosen)) {
    if (is.null(chain_id)) cat("Back-transforming...\n") else cat(sprintf("Chain %d: Back-transforming...\n", chain_id))

    # 1. Reverse QR FIRST
    if (do_qr) {
      if (use_re) {
        # Vectorized QR backsolve across all groups and samples by flattening extra dimensions
        dim_total <- dim(postb_total)
        dim(postb_total) <- c(dim_total[1], prod(dim_total[-1]))
        postb_total <- backsolve(X_R_mat, postb_total)
        dim(postb_total) <- dim_total

        dim_pooled <- dim(postb_pooled)
        dim(postb_pooled) <- c(dim_pooled[1], prod(dim_pooled[-1]))
        postb_pooled <- backsolve(X_R_mat, postb_pooled)
        dim(postb_pooled) <- dim_pooled
      } else {
        dim_total <- dim(postb_total)
        dim(postb_total) <- c(dim_total[1], prod(dim_total[-1]))
        postb_total <- backsolve(X_R_mat, postb_total)
        dim(postb_total) <- dim_total
        postb_pooled <- postb_total
      }
    }

    # Save standardized parameters before unscaling
    postb_total_std <- postb_total
    postb_pooled_std <- postb_pooled

    # 2. Reverse Standardization LATER
    if ((do_cen || do_scl) && sum(cont_idx) > 0) {
      is_cont <- which(cont_idx)
      int_idx <- which(apply(X_orig, 2, var) == 0) # Find constant columns (usually intercept)

      if (use_re) {
        if (do_scl) {
          for (ci in seq_along(is_cont)) {
            postb_total[is_cont[ci], , , ] <- postb_total[is_cont[ci], , , ] / X_rescaling[2, ci]
            postb_pooled[is_cont[ci], , ] <- postb_pooled[is_cont[ci], , ] / X_rescaling[2, ci]
          }
        }
        if (do_cen && length(int_idx) > 0) {
          # Fully vectorized shift subtraction
          shift <- colSums(postb_total[is_cont, , , , drop = FALSE] * X_rescaling[1, ], dims = 1)
          for (ii in int_idx) {
            target_dim <- dim(postb_total[ii, , , ])
            if (is.null(target_dim)) target_dim <- length(postb_total[ii, , , ])
            postb_total[ii, , , ] <- postb_total[ii, , , ] - array(shift, dim = target_dim)
          }
          shift_p <- colSums(postb_pooled[is_cont, , , drop = FALSE] * X_rescaling[1, ], dims = 1)
          for (ii in int_idx) {
            target_dim <- dim(postb_pooled[ii, , ])
            if (is.null(target_dim)) target_dim <- length(postb_pooled[ii, , ])
            postb_pooled[ii, , ] <- postb_pooled[ii, , ] - array(shift_p, dim = target_dim)
          }
        }
      } else {
        if (do_scl) {
          for (ci in seq_along(is_cont)) {
            postb_total[is_cont[ci], , ] <- postb_total[is_cont[ci], , ] / X_rescaling[2, ci]
          }
        }
        if (do_cen && length(int_idx) > 0) {
          shift <- colSums(postb_total[is_cont, , , drop = FALSE] * X_rescaling[1, ], dims = 1)
          for (ii in int_idx) {
            target_dim <- dim(postb_total[ii, , ])
            if (is.null(target_dim)) target_dim <- length(postb_total[ii, , ])
            postb_total[ii, , ] <- postb_total[ii, , ] - array(shift, dim = target_dim)
          }
        }
        postb_pooled <- postb_total
      }
    }
  }

  # --- 6. ATTACH NAMES ---
  cov_names <- cov_names_save
  cat_names <- colnames(Y)

  if (!save_posterior_to_disk) {
    if (use_re) {
      # Label the group dimension with the ACTUAL group identities in slice order
      # (groups <- unique(group_idx), appearance order), NOT positional "1..n" — the slice
      # index is a POSITION, and mislabeling it as the group value is exactly what breaks
      # any postb_total[,,g] reconstruction. Downstream must key via match(g, group_levels).
      g_names <- as.character(groups)
      dimnames(postb_total) <- list(cov_names, cat_names, g_names, NULL)
      dimnames(postb_pooled) <- list(cov_names, cat_names, NULL)
      if (!is.null(postb_total_std)) {
        dimnames(postb_total_std) <- list(cov_names, cat_names, g_names, NULL)
        dimnames(postb_pooled_std) <- list(cov_names, cat_names, NULL)
      }
    } else {
      dimnames(postb_total) <- list(cov_names, cat_names, NULL)
      dimnames(postb_pooled) <- list(cov_names, cat_names, NULL)
      if (!is.null(postb_total_std)) {
        dimnames(postb_total_std) <- list(cov_names, cat_names, NULL)
        dimnames(postb_pooled_std) <- list(cov_names, cat_names, NULL)
      }
    }
  }

  # --- END-OF-SAMPLER SELF-CHECK: does the RETURNED (back-transformed, raw-space)
  # postb_total reproduce post_log_lik using the ORIGINAL raw X? Isolates the
  # back-transform (scale/QR reversal) as a suspect independent of the store-time check.
  if (.selfcheck_store && !save_posterior_to_disk && exists("X_orig")) {
    tryCatch({
      Xr <- as.matrix(X_orig)
      for (d in seq_len(min(.selfcheck_n, dim(postb_total)[length(dim(postb_total))]))) {
        Ud <- matrix(0, n, p_all)
        if (use_re) {
          for (m in seq_len(n_groups)) {
            ii <- idx_list[[m]]
            Ud[ii, ] <- Xr[ii, , drop = FALSE] %*% postb_total[, , m, d]
          }
        } else {
          Ud <- Xr %*% postb_total[, , d]
        }
        llr <- compute_loglik(Ud, Y, as.numeric(y_weight))
        cat(sprintf("[SELFCHECK-END d=%d] post_log_lik=%.2f | RETURNED postb_total x rawX=%.2f (d=%.2f)\n",
                    d, post_log_lik[d], llr$total, llr$total - post_log_lik[d]))
      }
    }, error = function(e) cat("[SELFCHECK-END] skipped:", conditionMessage(e), "\n"))
  }

  # --- 7. SUMMARY & DIAGNOSTICS ---
  diagnostics <- list(
    loglik = mean(post_log_lik)
  )
  res <- list(
    post_delta = post_delta,      # alternative-specific coefficient (NULL unless alt_spec_Z given)
    postb = postb_total, postb_total = postb_total, postb_pooled = postb_pooled,
    postb_total_std = postb_total_std,
    postb_pooled_std = postb_pooled_std,
    post_sigma_re = post_sigma_re,
    post_re_tau = post_re_tau,
    post_slab_c2 = post_slab_c2,
    sigma_beta_pooled = sigma_beta_pooled,
    post_log_lik = post_log_lik, post_log_lik_pointwise = post_ll_pw,
    diagnostics = diagnostics, baseline = baseline,
    tree_store = if (save_posterior_to_disk) "Integrated in posterior_store" else (if (save_bart_to_disk) bart_batch_files else tree_store),
    posterior_store = if (save_posterior_to_disk) posterior_batch_files else NULL,
    bart_batch_size = if (save_bart_to_disk || save_posterior_to_disk) bart_batch_size else NULL,
    post_f_mean = if (use_bart) post_f_sum / nretain else NULL,
    post_f_sd = if (use_bart) sqrt(pmax(post_f_sum_sq / nretain - (post_f_sum / nretain)^2, 0)) else NULL,
    post_f = if (use_bart && store_f) post_f else NULL,
    bart_scaling = if (use_bart & !do_slim_trees) bart_scaling else NULL,
    bart_idx = if (use_bart & !do_slim_trees) bart_idx else NULL,
    # Adaptive Phase outputs
    post_delta_mean = if (use_spike_slab && use_re && store_delta && !save_posterior_to_disk) post_delta_count / nretain else NULL,
    pi_inc = if (use_spike_slab && use_re) pi_inc else NULL,
    post_tau_spatial = if (car_active && !save_posterior_to_disk) post_tau_spatial else NULL
  )

  # --- Capture Final State for Hot-Start ---
  res$final_state <- final_state_raw

  if (calc_loo && !save_posterior_to_disk) {
    res$waic <- loo::waic(post_ll_pw)
    res$loo <- loo::loo(post_ll_pw)
  }
  if (use_horseshoe) {
    res$horseshoe <- list(
      lambda2 = hs_lambda2, tau2 = hs_tau2,
      c2 = hs_c2,
      post_kappa_pooled = post_kappa_pooled,
      post_c2 = post_c2,
      nu = hs_nu, xi = hs_xi, zeta = hs_zeta,
      hs_pool = hs_pool, channel_split = hs_split_on,
      blk_tau2 = blk_tau2, block_lambda2 = block_lambda2,
      grouped = hs_grouped, tau2_by_family = if (hs_grouped) setNames(hs_tau2_g, hs_grp_names) else NULL,
      family_of_col = if (hs_grouped) setNames(hs_grp_names[hs_grp], .nm) else NULL
    )
  }
  return(res)
}

# =============================================================================
# Recovery Utility: Reconstructs the full RAM-based posterior list from disk batches
# =============================================================================
recover_mnlogit_posterior <- function(path_or_files, chain_id = NULL) {
  # 1. Handle Input: File List or Directory Path
  if (length(path_or_files) == 1 && dir.exists(path_or_files)) {
    cat(sprintf("Scanning directory for posterior batches: %s\n", path_or_files))

    # Define pattern: if chain_id is provided, filter for it
    pattern <- if (is.null(chain_id)) {
      "posterior_batch_.*\\.qs$"
    } else {
      sprintf("posterior_batch_.*_chain_%s\\.qs$", as.character(chain_id))
    }

    posterior_files <- list.files(path_or_files, pattern = pattern, full.names = TRUE)

    # Sort files numerically to ensure iteration order
    # (Extracts the batch number from 'posterior_batch_X_...')
    batch_nums <- as.integer(gsub(".*batch_([0-9]+)_.*", "\\1", basename(posterior_files)))
    posterior_files <- posterior_files[order(batch_nums)]
  } else {
    posterior_files <- path_or_files
  }

  if (is.null(posterior_files) || length(posterior_files) == 0) {
    if (!is.null(chain_id)) {
      message(sprintf("No posterior files found for chain %s.", chain_id))
      return(NULL)
    }
    stop("No posterior files provided or found.")
  }

  cat(sprintf("Recovering posterior from %d batch files...\n", length(posterior_files)))

  # 2. Load all batches into a flat list of samples
  all_samples <- list()
  for (f in posterior_files) {
    if (!file.exists(f)) {
      warning("File not found, skipping: ", f)
      next
    }
    batch_data <- qs2::qs_read(f)
    all_samples <- c(all_samples, batch_data)
  }

  nretain <- length(all_samples)
  if (nretain == 0) {
    return(NULL)
  }

  # 2. Inspect first sample to determine dimensions and features
  # 2. Inspect first sample to determine dimensions and features
  s1 <- all_samples[[1]]
  beta_dim <- dim(s1$beta)
  k <- beta_dim[1]
  p <- beta_dim[2]
  use_re <- length(beta_dim) == 3
  n_groups <- if (use_re) beta_dim[3] else NULL

  has_sigma_re <- !is.null(s1$sigma_re)
  has_loo <- !is.null(s1$log_lik_pw)
  has_hs <- !is.null(s1$horseshoe)
  has_bart <- !is.null(s1$bart_trees)

  # 3. Restore metadata to determine symmetric mode and categories
  base_dir <- if (dir.exists(path_or_files[1])) path_or_files[1] else dirname(path_or_files[1])
  meta_file <- file.path(base_dir, "model_metadata.qs")
  meta <- if (length(meta_file) == 1 && file.exists(meta_file)) qs2::qs_read(meta_file) else NULL

  if (!is.null(meta)) {
    cov_names <- meta$cov_names
    cat_names_meta <- meta$cat_names
    # If cat_names has p+1 elements, symmetric mode was used. If p, it wasn't.
    p_all <- length(cat_names_meta)
    symmetric_mode <- (p_all == p + 1)
    
    cat_names <- cat_names_meta
    baseline_name <- meta$baseline_name %||% "Base"
    baseline <- if (baseline_name %in% cat_names) which(cat_names == baseline_name)[1] else p_all
  } else {
    cov_names <- dimnames(s1$beta)[[1]] %||% paste0("V", 1:k)
    cat_names_meta <- dimnames(s1$beta)[[2]]
    if (!is.null(cat_names_meta) && length(cat_names_meta) == p + 1) {
       symmetric_mode <- TRUE
       p_all <- p + 1
       cat_names <- cat_names_meta
    } else {
       symmetric_mode <- FALSE
       p_all <- p
       cat_names <- cat_names_meta %||% paste0("C", 1:p)
    }
    baseline <- p_all
  }
  
  if (symmetric_mode) {
    pp <- (1:p_all)[-baseline]
  } else {
    pp <- 1:p
  }

  cat(sprintf(
    "Dimensions detected: %d covariates, %d categories, %s\n",
    k, p_all, ifelse(use_re, paste0(n_groups, " groups"), "pooled model")
  ))

  # 4. Pre-allocate RAM storage
  res_list <- list()
  res_list$postb_total <- if (use_re) array(0, c(k, p_all, n_groups, nretain)) else array(0, c(k, p_all, nretain))
  res_list$postb_pooled <- array(0, c(k, p_all, nretain))
  res_list$post_log_lik <- numeric(nretain)

  if (has_sigma_re) res_list$post_sigma_re <- matrix(0, k * p_all, nretain)
  if (has_loo) res_list$post_log_lik_pointwise <- matrix(0, nretain, length(s1$log_lik_pw))
  if (has_bart) res_list$tree_store <- vector("list", nretain)

  if (has_hs) {
    k_hs <- dim(s1$horseshoe$lambda2)[1] %||% length(s1$horseshoe$lambda2)
    p_hs <- if (!is.null(dim(s1$horseshoe$tau2))) dim(s1$horseshoe$tau2)[2] else p
    res_list$horseshoe <- list(
      post_kappa_pooled = matrix(0, k_hs * p_hs, nretain),
      post_kappa_re     = if (use_re && !is.null(s1$horseshoe$phi2)) matrix(0, k_hs * p_hs, nretain) else NULL,
      post_c2           = numeric(nretain),
      post_c2_re        = if (use_re && !is.null(s1$horseshoe$c2_re)) numeric(nretain) else NULL
    )
  }

  expand_to_zs <- function(mat, pp_idx, p_all_dim) {
     full <- matrix(0, nrow(mat), p_all_dim)
     full[, pp_idx] <- mat
     sweep(full, 1, rowMeans(full), "-")
  }

  # 5. Fill storage arrays
  for (s in 1:nretain) {
    samp <- all_samples[[s]]

    if (use_re) {
      if (symmetric_mode) {
         beta_c_zs <- array(0, c(k, p_all, n_groups))
         beta_c_zs[, pp, ] <- samp$beta
         for (m in 1:n_groups) {
            beta_c_zs[, , m] <- sweep(beta_c_zs[, , m, drop=FALSE], 1, rowMeans(beta_c_zs[, , m, drop=FALSE]), "-")
         }
         res_list$postb_total[, , , s] <- beta_c_zs
         res_list$postb_pooled[, , s] <- expand_to_zs(samp$mu, pp, p_all)
         if (has_sigma_re) {
            res_list$post_sigma_re[, s] <- as.vector(apply(beta_c_zs, c(1, 2), sd))
         }
      } else {
         res_list$postb_total[, , , s] <- samp$beta
         res_list$postb_pooled[, , s] <- samp$mu
         if (has_sigma_re) res_list$post_sigma_re[, s] <- as.vector(samp$sigma_re)
      }
    } else {
      if (symmetric_mode) {
         beta_zs <- expand_to_zs(samp$beta, pp, p_all)
         res_list$postb_total[, , s] <- beta_zs
         res_list$postb_pooled[, , s] <- expand_to_zs(samp$mu, pp, p_all)
      } else {
         res_list$postb_total[, , s] <- samp$beta
         res_list$postb_pooled[, , s] <- samp$mu
      }
    }

    res_list$post_log_lik[s] <- samp$log_lik
    if (has_sigma_re && !symmetric_mode) res_list$post_sigma_re[, s] <- as.vector(samp$sigma_re)
    if (has_loo) res_list$post_log_lik_pointwise[s, ] <- samp$log_lik_pw
    if (has_bart) res_list$tree_store[[s]] <- samp$bart_trees

    if (has_hs) {
      k_hs <- dim(samp$mu)[1]
      p_hs <- dim(samp$mu)[2]

      # Level 1: Finnish-regularized pooled kappa (matches RAM-path formula)
      hs_c2_rec <- samp$horseshoe$c2 %||% 15
      hs_tau2_mat_rec <- if (length(samp$horseshoe$tau2) > 1) matrix(samp$horseshoe$tau2, k_hs, p_hs, byrow = TRUE) else samp$horseshoe$tau2
      lam2_tau2 <- samp$horseshoe$lambda2 * hs_tau2_mat_rec
      var_eff_rec <- (hs_c2_rec * lam2_tau2) / (hs_c2_rec + lam2_tau2)
      res_list$horseshoe$post_kappa_pooled[, s] <- as.vector(1 - var_eff_rec / lam2_tau2)
      res_list$horseshoe$post_c2[s] <- hs_c2_rec

      if (!is.null(res_list$horseshoe$post_kappa_re)) {
        # Level 2: Finnish-regularized RE deviation kappa
        hs_c2_re_rec <- samp$horseshoe$c2_re %||% 3
        hs_psi2_mat_rec <- if (length(samp$horseshoe$psi2) > 1) matrix(samp$horseshoe$psi2, k_hs, p_hs, byrow = TRUE) else samp$horseshoe$psi2
        phi2_psi2 <- samp$horseshoe$phi2 * hs_psi2_mat_rec
        var_re_rec <- (hs_c2_re_rec * phi2_psi2) / (hs_c2_re_rec + phi2_psi2)
        res_list$horseshoe$post_kappa_re[, s] <- as.vector(1 - var_re_rec / phi2_psi2)
        res_list$horseshoe$post_c2_re[s] <- hs_c2_re_rec
      }
    }
  }

  # Defensive: ensure name vectors match array extents. Guards against stale/mismatched
  # metadata so recovery degrades gracefully (warns + regenerates) instead of erroring.
  pa <- dim(res_list$postb_total)[2]
  if (length(cat_names) != pa) {
    warning(sprintf("recover: cat_names length (%d) != category extent (%d); regenerating. Check model_metadata.qs is fresh.",
                    length(cat_names), pa))
    cat_names <- paste0("cat", seq_len(pa))
  }
  if (length(cov_names) != dim(res_list$postb_total)[1]) {
    cov_names <- paste0("V", seq_len(dim(res_list$postb_total)[1]))
  }

  if (use_re) {
    # Label the group (3rd) dimension with real identities from metadata (appearance order,
    # matching the sampler's groups <- unique(group_idx)); fall back to NULL if unavailable
    # or length-mismatched. The slice is a POSITION — downstream keys via match(g, group_levels).
    g_lab <- if (!is.null(meta) && !is.null(meta$group_levels) &&
                 length(meta$group_levels) == dim(res_list$postb_total)[3]) meta$group_levels else NULL
    dimnames(res_list$postb_total) <- list(cov_names, cat_names, g_lab, NULL)
    dimnames(res_list$postb_pooled) <- list(cov_names, cat_names, NULL)
  } else {
    dimnames(res_list$postb_total) <- list(cov_names, cat_names, NULL)
    dimnames(res_list$postb_pooled) <- list(cov_names, cat_names, NULL)
  }

  res_list$postb <- res_list$postb_total
  res_list$var_names <- cov_names
  res_list$cat_names <- cat_names

  # 6. Apply back-transformation ONLY if the on-disk batches are still in scaled space.
  # The per-draw disk-save loop already back-transforms every batch to physical scale
  # (state_sample$beta is divided by sd + intercept-shifted before qs_save), so re-applying
  # here double-divides continuous slopes by sd (flattening them ~sd-fold) while post_log_lik
  # stays correct — the exact symptom that made downstream X*beta reconstruction fail. The
  # batches_back_transformed flag (default TRUE for any batch written by this sampler; absent
  # only in stale pre-flag metadata, which was also physical) gates this off. Only a caller
  # that deliberately saves SCALED batches would set it FALSE to re-enable this block.
  if (!is.null(meta) && !is.null(meta$X_rescaling) && isFALSE(meta$batches_back_transformed)) {
    cat("Back-transforming posterior to physical scale based on metadata...\n")
    X_rescaling <- meta$X_rescaling
    cont_idx <- meta$cont_idx
    do_cen <- meta$do_cen
    do_scl <- meta$do_scl

    if ((do_cen || do_scl) && sum(cont_idx) > 0) {
      is_cont <- which(cont_idx)
      int_idx <- which(cov_names == "intercept")

      if (use_re) {
        if (do_scl) {
          for (ci in seq_along(is_cont)) {
            res_list$postb_total[is_cont[ci], , , ] <- res_list$postb_total[is_cont[ci], , , ] / X_rescaling[2, ci]
            res_list$postb_pooled[is_cont[ci], , ] <- res_list$postb_pooled[is_cont[ci], , ] / X_rescaling[2, ci]
          }
        }
        if (do_cen && length(int_idx) > 0) {
          # Vectorized shift: sum across all continuous variables in one shot
          # (matches main sampler approach — clearer and faster than sequential per-variable loop)
          shift <- colSums(res_list$postb_total[is_cont, , , , drop = FALSE] * X_rescaling[1, ], dims = 1)
          for (ii in int_idx) {
            target_dim <- dim(res_list$postb_total[ii, , , ])
            if (is.null(target_dim)) target_dim <- length(res_list$postb_total[ii, , , ])
            res_list$postb_total[ii, , , ] <- res_list$postb_total[ii, , , ] - array(shift, dim = target_dim)
          }
          shift_p <- colSums(res_list$postb_pooled[is_cont, , , drop = FALSE] * X_rescaling[1, ], dims = 1)
          for (ii in int_idx) {
            target_dim <- dim(res_list$postb_pooled[ii, , ])
            if (is.null(target_dim)) target_dim <- length(res_list$postb_pooled[ii, , ])
            res_list$postb_pooled[ii, , ] <- res_list$postb_pooled[ii, , ] - array(shift_p, dim = target_dim)
          }
        }
      } else {
        if (do_scl) {
          for (ci in seq_along(is_cont)) res_list$postb_total[is_cont[ci], , ] <- res_list$postb_total[is_cont[ci], , ] / X_rescaling[2, ci]
        }
        if (do_cen && length(int_idx) > 0) {
          shift <- colSums(res_list$postb_total[is_cont, , , drop = FALSE] * X_rescaling[1, ], dims = 1)
          for (ii in int_idx) {
            target_dim <- dim(res_list$postb_total[ii, , ])
            if (is.null(target_dim)) target_dim <- length(res_list$postb_total[ii, , ])
            res_list$postb_total[ii, , ] <- res_list$postb_total[ii, , ] - array(shift, dim = target_dim)
          }
        }
        res_list$postb_pooled <- res_list$postb_total
      }
    }
  }

  # NEW: Recover final state if it exists (for hot-starting)
  fs_file <- file.path(base_dir, sprintf("final_state_chain_%s.qs", if (is.null(chain_id)) "1" else as.character(chain_id)))
  if (file.exists(fs_file)) {
    res_list$final_state <- qs2::qs_read(fs_file)
  }

  # Surface BART reconstruction metadata so f = C g can be rebuilt on recovery / downstream
  if (!is.null(meta) && isTRUE(meta$use_bart)) {
    res_list$bart <- list(
      symmetric = isTRUE(meta$bart_symmetric),
      bart_idx  = meta$bart_idx,
      pp        = meta$bart_pp %||% pp,
      p_all     = meta$p_all %||% p_all
    )
  }

  cat("Recovery complete.\n")
  return(res_list)
}

# --- Reconstruct BART posterior-mean f from recovered slim trees ---
# tree_store : per-draw list of p slim-tree data.frames (CONTRAST ensembles if symmetric).
# X_bart     : RAW-scale design matrix for the BART covariates (n x length(bart_idx)).
# bart_meta  : res_list$bart from recover_mnlogit_posterior() (symmetric, C, pp, p_all).
# Returns the posterior-mean f [n x p_all]; zero-sum across categories when symmetric
# (f = C g), baseline-coded otherwise. This is the symmetric-aware reconstruction the
# disk path needs (post_f_mean is not stored in batches).
reconstruct_bart_f_mean <- function(tree_store, X_bart, bart_meta) {
  if (is.null(tree_store) || length(tree_store) == 0 || is.null(bart_meta)) return(NULL)
  X_bart <- as.matrix(X_bart)
  n <- nrow(X_bart); p_all <- bart_meta$p_all
  f_sum <- matrix(0, n, p_all); nd <- 0L
  for (draw in tree_store) {
    if (is.null(draw)) next
    pdraw <- length(draw)
    g <- vapply(seq_len(pdraw),
                function(cc) as.numeric(predict_slim_bart_cpp(X_bart, draw[[cc]])),
                numeric(n))                                    # n x p ensemble predictions
    if (isTRUE(bart_meta$symmetric)) {
      f <- sweep(g, 1, rowMeans(g), "-")                       # CLR: center per-category g (zero-sum)
    } else {
      f <- matrix(0, n, p_all); f[, bart_meta$pp] <- g         # baseline-coded
    }
    f_sum <- f_sum + f; nd <- nd + 1L
  }
  if (nd == 0L) return(NULL)
  f_sum / nd
}

# --- Multi-Chain Recovery Helper ---
# Reconstructs a list of chain results from a directory
recover_mnlogit_chains <- function(path) {
  if (!dir.exists(path)) stop("Directory does not exist: ", path)

  # Identify all unique chain IDs in the directory
  files <- list.files(path, pattern = "posterior_batch_.*_chain_.*\\.qs$")
  if (length(files) == 0) stop("No chain batches found in ", path)

  chain_ids <- unique(gsub(".*_chain_([0-9]+)\\.qs$", "\\1", files))
  cat(sprintf("Found %d chains in directory: %s\n", length(chain_ids), paste(chain_ids, collapse = ", ")))

  res_chains <- lapply(chain_ids, function(id) {
    recover_mnlogit_posterior(path, chain_id = id)
  })

  names(res_chains) <- paste0("chain_", chain_ids)
  return(res_chains)
}
