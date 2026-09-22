# =============================================================================
# mncount_rcpp.R  —  Bayesian Negative-Binomial / Poisson regression
#                    via Pólya-Gamma augmentation (Rcpp-accelerated)
#
# RELATIONSHIP TO MNL:
#   The MNL models P(j|i) = softmax(U_ij).  The count model instead models
#   Y_ij ~ NB(r_j, sigmoid(eta_ij))  independently for each outcome column j.
#
#   Three structural changes from mnlogit_rcpp.R:
#     1. No baseline / pp indexing.  All p = ncol(Y) columns are modelled.
#     2. No c_j (competitor utility) computation.  The PG location is just
#        eta_ij = X_i %*% beta_j (+ BART + RE), not eta_ij - c_j.
#     3. kappa_{ij} = (Y_{ij} - r_j)/2 depends on the dispersion r_j, so it
#        is recomputed inside the Gibbs loop after each r update.
#
#   For Poisson: set family = "poisson".  Internally r is fixed at r_fixed
#   (default 1e5); the PG draws become PG(Y_ij + r_fixed, eta_ij) with
#   r_fixed large → PG(Y_ij, eta_ij) in the limit.  kappa = Y/2.
#
# USAGE:
#   Rcpp::sourceCpp("mncount_gibbs_core.cpp")
#   source("mncount_rcpp.R")
#   fit <- mncount_rcpp(X, Y, family = "negbin", ...)
#
# ARGUMENTS (identical to mnlogit_rcpp unless noted):
#   family         : "negbin" (default) or "poisson"
#   r_init         : initial dispersion vector length p  (NB only)
#   r_prior_shape  : Gamma prior shape for r_j            (NB only)
#   r_prior_rate   : Gamma prior rate  for r_j            (NB only)
#   r_mh_sd        : log-scale MH proposal SD for r_j    (NB only)
#   r_fixed        : fixed dispersion when family="poisson" or use_re=FALSE
#   offset         : n x p matrix of offsets (log-scale), e.g. log(exposure)
#
# All RE / BART / Horseshoe / constraint arguments carry over unchanged.
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
library(RcppArmadillo)

aux_path <- "codes/mnl_aux_func.R"
if (!file.exists(aux_path)) aux_path <- "mnl_aux_func.R"
source(aux_path)

# ---------------------------------------------------------------------------
# INERT ARGUMENTS (audited 2026-08-26). Accepted by the signature but reaching NOTHING in the body.
# Kept rather than removed so existing callers do not break, but setting them now WARNS: a knob that
# is documented, passed in production and silently ignored is exactly how the tau0 calibration went
# unnoticed for weeks.
#   tau0_dev        - RE-side horseshoe scale. Never wired. Currently PASSED by run_ls_count_model.R
#                     and postprocess/build_count_html.R, with no effect.
#   p0_re           - RE-side guess at the number of non-null effects. Never wired.
#   bart_batch_size - vestigial from the MNL BART block.
#   symmetric       - zero-sum output coding; the count outcome is not compositional so it has no
#                     meaning here. (The MNL's own `symmetric` IS real -- do not confuse the two.)
# ---------------------------------------------------------------------------
.count_inert_args <- function(called) {
  defaults <- list(tau0_dev = NULL, p0_re = 2, bart_batch_size = 50, symmetric = FALSE)
  set <- names(defaults)[vapply(names(defaults), function(a)
    a %in% names(called) && !identical(called[[a]], defaults[[a]]), logical(1))]
  if (length(set))
    warning("mncount_rcpp: argument(s) ", paste(set, collapse = ", "),
            " are INERT (signature-only; they reach nothing in the sampler), so setting them has ",
            "no effect. See the inert-argument note in codes/count_rcpp.R.", call. = FALSE)
  invisible(set)
}


# ---------------------------------------------------------------------------
# POOLED RE PRECISION (one scale per PREDICTOR, shared across the p outcome columns).
#
# Mirrors update_re_precision_hc but pools the sufficient statistic over outcome columns, which is
# the count analogue of the MNL's symmetric (category-pooled) RE variance. Makalic-Schmidt
# half-Cauchy auxiliaries, with the same optional regularised (Finnish) cap via a 1-D slice on
# log(tau_raw) so that tau_eff = tau_raw + 1/c2 -- identical parameterisation to the C++ path, so
# re_regularize / re_slab_c2 mean the same thing in both.
#
# ss IS UNCENTRED BY CONSTRUCTION. See the re_prec_pooled note in the signature: the count outcome
# columns are independent rates with no baseline and no zero-sum constraint, so the category-centring
# that broke the MNL's RE variance has no meaning here and is not offered as an option.
# ---------------------------------------------------------------------------
.count_re_precision_pooled <- function(beta_c, mu_pooled, re_idx, n_groups, re_mask, y_mask,
                                       is_intercept, prec_prev, re_scale_A, re_regularize,
                                       slab_c2, re_support, p) {
  prec_out <- prec_prev
  a_out    <- matrix(re_scale_A^2 / 2, nrow(prec_prev), ncol(prec_prev))
  inv_A2   <- 1 / (re_scale_A^2)
  inv_c2   <- if (isTRUE(re_regularize) && slab_c2 > 0 && is.finite(slab_c2)) 1 / slab_c2 else 0

  for (v in re_idx) {
    ss <- 0; n_active <- 0
    for (m in seq_len(n_groups)) {
      for (ip in seq_len(p)) {
        if (re_mask[v, ip, m] > 0.5 && (is_intercept[v] == 1L || y_mask[ip, m] > 0.5)) {
          d <- beta_c[v, ip, m] - mu_pooled[v, ip]              # UNCENTRED: the draw's own coords
          w <- if (!is.null(re_support)) { rs <- re_support[v, m]; if (rs > 1e-12) 1 / (rs * rs) else 1 } else 1
          ss <- ss + d * d * w
          n_active <- n_active + 1
        }
      }
    }
    if (n_active < 1) next

    tau_prev <- max(mean(prec_prev[v, ]), 1e-8)
    a_new    <- max(1 / rgamma(1, shape = 1, rate = max(tau_prev - inv_c2, 1e-8) + inv_A2), 1e-12)

    if (inv_c2 <= 0) {
      rate <- a_new + 0.5 * ss
      if (!is.finite(rate) || rate <= 0) rate <- inv_A2
      pr <- rgamma(1, shape = 0.5 + 0.5 * n_active, rate = rate)
    } else {
      lp <- function(lt) { tr <- exp(lt); te <- tr + inv_c2
        0.5 * n_active * log(te) - 0.5 * te * ss - 0.5 * lt - a_new * tr + lt }
      lt <- log(max(tau_prev - inv_c2, 1e-8)); y0 <- lp(lt) - rexp(1)
      Lb <- lt - runif(1); Rb <- Lb + 1
      gi <- 0L; while (lp(Lb) > y0 && gi < 80L) { Lb <- Lb - 1; gi <- gi + 1L }
      gi <- 0L; while (lp(Rb) > y0 && gi < 80L) { Rb <- Rb + 1; gi <- gi + 1L }
      ln <- lt
      for (.s in 1:100) { q <- Lb + runif(1) * (Rb - Lb)
        if (lp(q) > y0) { ln <- q; break }
        if (q < lt) Lb <- q else Rb <- q }
      pr <- exp(ln) + inv_c2
    }
    prec_out[v, ] <- min(max(pr, 1e-4), 1e6)                    # one scale, shared across columns
    a_out[v, ]    <- a_new
  }
  list(prec = prec_out, sigma = 1 / sqrt(prec_out + 1e-12), a_aux = a_out)
}


# --- Compile core C++ (shared with MNL sampler + count-specific extensions) ---
.needs_cpp_compile <- TRUE
if (exists("compute_utilities_count_cpp", mode = "function")) {
  .needs_cpp_compile <- tryCatch({
    # Test if the compiled symbol is actually loaded and callable in this process
    compute_utilities_count_cpp(matrix(0, 1, 1), matrix(0, 1, 1), matrix(0, 1, 1), matrix(0, 1, 1))
    FALSE
  }, error = function(e) TRUE)
}
if (.needs_cpp_compile) {
  # We need both the shared MNL core utilities and the count-specific extensions.
  # Since count_gibbs_core.cpp relies on functions in mnlogit_gibbs_core.cpp, 
  # we combine them into a single translation unit for Rcpp.
  mnl_cpp <- if (file.exists("codes/mnlogit_gibbs_core.cpp")) "codes/mnlogit_gibbs_core.cpp" else "mnlogit_gibbs_core.cpp"
  cnt_cpp <- if (file.exists("codes/count_gibbs_core.cpp")) "codes/count_gibbs_core.cpp" else "count_gibbs_core.cpp"
  if (file.exists(mnl_cpp) && file.exists(cnt_cpp)) {
    tmp_cpp <- tempfile(fileext = ".cpp")
    writeLines(c(readLines(mnl_cpp), readLines(cnt_cpp)), tmp_cpp)
    Rcpp::sourceCpp(tmp_cpp)
  }
}
rm(.needs_cpp_compile)

# =============================================================================
# Internal: log-posterior of r_j (dispersion) given data and linear predictor
# Used for Metropolis-Hastings accept/reject.
# =============================================================================
.lp_r <- function(r, y_j, eta_j, prior_shape, prior_rate) {
  if (r <= 0) return(-Inf)
  ll   <- sum(lgamma(y_j + r) - lgamma(r) - r * log1p(exp(eta_j)))
  lprior <- (prior_shape - 1) * log(r) - prior_rate * r
  ll + lprior
}

# Metropolis-Hastings step for a single dispersion parameter r_j
.update_r_mh <- function(r_curr, y_j, eta_j,
                         prior_shape = 1, prior_rate = 0.1,
                         mh_sd = 0.2) {
  log_r_prop <- log(r_curr) + rnorm(1, 0, mh_sd)
  r_prop     <- exp(log_r_prop)
  log_acc    <- .lp_r(r_prop, y_j, eta_j, prior_shape, prior_rate) +
                log_r_prop -
                .lp_r(r_curr, y_j, eta_j, prior_shape, prior_rate) -
                log(r_curr)   # Jacobian for log-scale proposal
  if (is.finite(log_acc) && log(runif(1)) < log_acc) r_prop else r_curr
}

# =============================================================================
# EXACT Gibbs update of r_j via Chinese-Restaurant-Table (CRT) augmentation
# (Zhou & Carin 2012). For Y_i ~ NB(r, psi_i), psi_i = sigmoid(eta_i):
#   L_i ~ CRT(Y_i, r) = sum_{j=1}^{Y_i} Bernoulli(r/(r+j-1))
#   r | L, eta ~ Gamma(a0 + sum_i L_i,  b0 + sum_i log1p(exp(eta_i)))   [-log(1-psi)=log1p(e^eta)]
# Replaces the single-site log-RW MH: exact, tuning-free, mixes far better.
# Efficiency: we need only sum_i L_i = sum_{j>=1} Binomial(c_j, r/(r+j-1)),
#   c_j = #{i: Y_i >= j}  -> O(max(Y)) binomial draws, not O(sum Y) Bernoullis.
# =============================================================================
.crt_sum <- function(y, r) {
  my <- max(y)
  if (my <= 0) return(0)
  tab <- tabulate(y[y > 0], nbins = my)      # tab[v] = #{Y_i == v}
  cj  <- rev(cumsum(rev(tab)))               # cj[j]  = #{Y_i >= j}
  # j=1 always has prob = r/r = 1 (every customer opens at least one table) -> deterministic = cj[1].
  # rbinom(size, prob=1) intermittently returns NA in R, so NEVER pass prob>=1 to rbinom; only j>=2.
  s <- cj[1]
  if (my >= 2) {
    jj   <- 2:my
    prob <- pmin(pmax(r / (r + jj - 1), 0), 1 - 1e-12)
    s <- s + sum(rbinom(my - 1, cj[jj], prob))
  }
  s
}
.update_r_crt <- function(r_curr, y_j, eta_j, prior_shape = 1, prior_rate = 0.1) {
  Lsum <- .crt_sum(y_j, r_curr)
  rate <- prior_rate + sum(log1p(exp(eta_j)))   # eta is clipped upstream -> exp is safe
  rgamma(1, shape = prior_shape + Lsum, rate = max(rate, 1e-8))
}

# =============================================================================
# Main sampler
# =============================================================================
mncount_rcpp <- function(
    X, Y,
    family         = c("negbin", "poisson"),
    intercept      = FALSE,
    niter          = 1000,   nburn    = 500,   nsample = NULL,
    thin           = 1L,     # store every `thin`-th post-burn draw (RAM: all posterior arrays shrink ~thin x; keeps iterations for ESS but stores fewer). thin << autocorr-time loses ~no ESS.
    A0             = 4.0,
    # --- NB-specific ---
    r_init         = NULL,   # length-p init; NULL → rep(1, p)
    r_prior_shape  = 1,      # Gamma(shape, rate) prior on r_j
    r_prior_rate   = 0.1,
    r_mh_sd        = 0.2,    # log-scale MH proposal SD (only when r_method="mh")
    r_method       = c("crt", "mh"),  # "crt" = exact CRT-Gibbs (default, mixes far better); "mh" = legacy log-RW Metropolis
    r_warmup       = NULL,   # iteration after which the r-update starts. NULL => 0.45*nburn when use_bart (let the MEAN settle Poisson-like first so the CRT sees mu~Y and lands on the good r, not the collapsed r=0.07 bad mode), else burn_buffer.
    bart_absorb_mean = FALSE,# absorb the BART f-mean into the intercept each iter (identification only). Default OFF: the QR cdir move is not POINTWISE mean-preserving under the exp-link -> it LEAKS level into mu (~9x over-pred). Leaving f un-centered is exact for prediction (only U+f matters); the split is a benign label issue.
    mean_param     = TRUE,   # MEAN parameterization: model log(mu)=X.beta+f+offset (r-FREE) by folding -log(r) into the offset that forms eta_nat everywhere. Removes the intercept<->log(r) ridge at its ROOT (vs r_asis's discrete post-hoc shift, which fights BART). When TRUE, r_asis is auto-disabled (subsumed). This is what makes count-BART stable. Off = legacy natural param (mu=r*e^eta).
    r_asis         = TRUE,   # mean-preserving dispersion interweave: after each r draw, shift the intercept by log(r_old/r_new) so mu=r*e^U is held fixed -> r identified by DISPERSION alone (breaks the r<->intercept ridge that made r drift to 0.2/8 instead of the MLE ~2). Off = legacy behaviour. SUPERSEDED by mean_param (kept for the legacy natural-param path).
    re_center      = FALSE,  # IDENTIFY mu AS THE MEAN COUNTRY EFFECT. The likelihood only sees
                             # beta_g = mu + b_g, so the split between mu and the MEAN of the b_g is
                             # unidentified -- only the sum is. The priors break that tie, and they are
                             # asymmetric: mu carries the HORSESHOE (sparsity, pulls to 0) while the
                             # b_g carry only a variance prior (nothing pulls their mean to 0). In the
                             # NCP-ASIS conditional below the precision is (G*prec_re + hs_prec), so
                             # whenever hs_prec dominates 34 countries' evidence the common effect
                             # MIGRATES OUT OF mu INTO THE RE MEAN. Measured on the BOV livestock fit:
                             # GHM_HI mu = 0.05 while mean_g(beta_g) = 4.85 -- 99% of the average effect
                             # sat in the RE mean; median 79% across the 36 RE covariates. The intercept
                             # was the sole exception, and it is the one column excluded from the
                             # horseshoe -- the mechanism's own control.
                             # TRUE = do not sparsity-shrink mu AGAINST the RE mean (drop hs_prec on the
                             # re_idx rows of the mu conditional), so mu is identified by the group means
                             # and equals the average country effect. Sparsity is then judged on the
                             # TOTAL effect, not on an aliased component. Non-RE covariates keep the
                             # horseshoe untouched. Off = legacy behaviour.
    level_asis     = TRUE,   # location ASIS: CP redraw of the RE means (mu_pooled) given the group totals (curr_beta_c FIXED -> mu_g unchanged). Interweaves with the NCP RE draw -> decouples the pooled LEVEL from the RE deviations (breaks the intercept<->RE-mean ridge that leaves the intercept unconverged / the mean 7x too high). Off = legacy.
    r_fixed        = 1e5,    # used when family="poisson"
    # --- Offset (log scale) ---
    offset         = NULL,   # n x p or n x 1 or NULL
    # --- RE ---
    use_re         = FALSE,  group_idx   = NULL,
    re_idx         = NULL,   # defaults to 1:ncol(X)
    prior_a_re     = 0.01,   prior_b_re  = 0.01,
    use_half_cauchy_re = TRUE, re_scale_A = 1.0,
    # POOLED RE VARIANCE ACROSS OUTCOME COLUMNS (ported from the MNL's symmetric RE block,
    # 2026-08-21). FALSE (default) = historical: update_re_precision_hc draws ONE precision per
    # (predictor, outcome column). TRUE = one precision per PREDICTOR, shared across all p outcome
    # columns, so a covariate's group-level spread is estimated from p times as much information.
    #
    # DELIBERATELY UNCENTRED, and that is the whole point of the port. The MNL's
    # update_re_precision_hc_sym built its sum of squares from CATEGORY-CENTRED deviations while the
    # RE draw stayed uncentred; for an RE that is largely a common shift that zeroes the random
    # effects outright (measured: RE sd 0.0000, -51.6 nats). In a COUNT model there is no baseline
    # category and no zero-sum constraint at all -- the p columns are independent rates, not a
    # composition -- so centring here would be not merely mis-scaled but meaningless. Only the
    # uncentred statistic is defined, and it is the one used below.
    re_prec_pooled = FALSE,
    use_country_shrinkage = FALSE, tau0_country = 1.0, slab_c2_country = 4.0, # Country-Gatekeeper shrinkage tau_m ~ Regularized-Half-Cauchy(0, tau0, slab_c2)
    # Does the per-group factor cover the INTERCEPT RE too, or only the slopes?
    # FALSE (default) = slopes only. The country intercept carries most of the model's skill (~45 nats
    # on the pixel model, essentially all of it the intercept), so shrinking it per-country trades away
    # the part that works. TRUE gives the literal "one factor across ALL RE covariates for each group".
    country_shrink_intercept = FALSE,
    # IDENTIFICATION. sigma_v (per-covariate RE scale) and tau_m (per-group factor) enter the RE
    # prior ONLY as the product sigma_v * tau_m, so the pair is not identified: tau_m can shrink and
    # sigma_v inflate to compensate, which enlarges the RE space until it absorbs the fixed effects
    # (measured on a synthetic panel: b[x1] 0.731 -> 0.256 against truth 0.700, mean RE sd 0.46 ->
    # 3.21). TRUE pins geometric-mean(tau_m) = 1 after each draw and moves that overall scale into
    # sigma_v, which leaves every product sigma_v * tau_m EXACTLY unchanged while fixing the ridge.
    # The RELATIVE gating between groups -- the entire point of the feature -- is preserved.
    #
    # DEFAULT FALSE: PARTIALLY WORKING, DO NOT USE WITHOUT READING THIS. Anchoring does fix the FE
    # absorption (b[x1] 0.256 -> 0.823 against truth 0.700; mean RE sd 3.21 -> 0.19), but the rescale
    # happens OUTSIDE the C++ draw, which applies a REGULARIZED half-Cauchy with slab_c2_country --
    # so the anchored tau_m lands outside the cap the prior enforces (observed arithmetic means of
    # 6.0e3 at geometric mean 1, i.e. tau wandering orders of magnitude) and the next sweep's slice
    # starts from an out-of-range value. The identification has to be imposed INSIDE
    # update_country_shrinkage_hc_cpp (draw tau on the anchored scale), not patched afterwards.
    # Until then: leave FALSE and treat use_country_shrinkage itself as experimental.
    country_shrink_anchor = FALSE,
    re_regularize = FALSE, re_slab_c2 = 100,   # regularised-horseshoe RE variance (slice): tau_eff=tau_raw+1/c2, SD<=sqrt(c2). Tames the heavy half-Cauchy tail (ported from the MNL). Off by default.
    support_prior_strength = 0,                # support-aware RE prior: per-(cov,group) participation-ratio SD shrinkage for info-sparse groups (ported from the MNL). 0 = off.
    init_jitter = 0,                           # per-chain overdispersed init: uniform +/- init_jitter on mu_pooled -> honest Rhat. 0 = off.
    re_cor_threshold = 0.85, re_cor_min_groups = 0.5,
    # --- BART ---
    use_bart       = FALSE,  bart_idx    = NULL,
    n_trees_bart   = 50,     n_threads_bart = 1,
    bart_base      = 0.90,   bart_power  = 3.0,  bart_k = 3.0,  # BART regularization (was UNSET -> dbarts loose defaults 0.95/2.0/2 -> f exploded to +-80). Count exp() link amplifies -> tight depth + STRONG leaf shrinkage (high k) to keep f small (+-3-5). Tune up if f still runs.
    bart_f_cap     = 6.0,    # hard cap on |BART f| (log scale) so exp(eta) can't blow up mu under the count exp-link; 6 -> at most exp(6)~400x density shift. Lower to regularise harder.
    store_bart_trees = TRUE, do_slim_trees  = TRUE,
    save_bart_to_disk = FALSE, bart_disk_path = tempdir(),
    bart_batch_size = 50,
    store_f        = FALSE,
    # --- Constraints ---
    positive_constraints = NULL, negative_constraints = NULL,
    # --- Horseshoe ---
    use_horseshoe  = FALSE,  horseshoe_idx = NULL,
    equation_specific_hs = FALSE,
    estimate_c2    = TRUE,   slab_df = 4, slab_s2 = 15,
    p0             = 5,      p0_re    = 2,
    tau0_mu        = NULL,   tau0_dev = NULL,
    # --- Misc ---
    use_ncp        = TRUE,
    calc_loo       = FALSE,
    symmetric      = FALSE,
    chain_id       = NULL,
    progress_cb    = NULL,
    prior_mu       = NULL,   prior_V_inv = NULL,
    linear_idx     = NULL,   # NULL → all columns
    standardize    = TRUE,
    method         = c("standardize", "center", "scale", "QR", "none"),
    gamma_matched_pg = FALSE,
    init_state     = NULL) {

  family <- match.arg(family)
  r_method <- match.arg(r_method)
  poisson_mode <- (family == "poisson")

  # ── 1. VALIDATION ────────────────────────────────────────────────────────
  # --- nsample: the KEPT draws, with niter derived -----------------------------------------------
  # niter is a TOTAL, so asking for N draws means passing N + nburn and remembering to. Getting it
  # wrong is not a small error: a caller that passed the total where it meant the sample ran
  # "0 sweeps (burn 4000)", and the generic symptom is the stop() just below, raised after the
  # design is already loaded. nsample states the quantity people actually have in mind.
  # Passing BOTH is refused rather than silently preferring one.
  if (!is.null(nsample)) {
    if (!missing(niter))
      stop("give nsample (kept draws) or niter (total), not both: nsample implies niter = nsample + nburn.",
           call. = FALSE)
    if (!is.finite(nsample) || nsample < 1) stop("nsample must be >= 1.", call. = FALSE)
    niter <- as.integer(nsample) + as.integer(nburn)
  }
  if (niter <= nburn) stop(sprintf("niter must be > nburn (niter=%d, nburn=%d). Pass nsample = <kept draws> instead of niter to avoid computing the total.", niter, nburn), call. = FALSE)
  if (any(!is.finite(Y))) stop("Y contains NA/Inf.")
  .count_inert_args(as.list(match.call())[-1])
  if (any(!is.finite(X))) stop("X contains NA/Inf.")
  if (any(Y < 0))         stop("Y contains negative values.")

  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(X)
  p <- ncol(Y)                       # model ALL outcome columns (no baseline)

  if (is.null(linear_idx)) linear_idx <- seq_len(ncol(X))
  if (is.null(re_idx))     re_idx     <- seq_len(length(linear_idx))
  if (is.null(bart_idx))   bart_idx   <- seq_len(ncol(X))

  # Convert character indices
  to_int <- function(idx, nm) {
    if (is.character(idx)) idx <- match(idx, nm)
    idx[!is.na(idx)]
  }
  linear_idx           <- to_int(linear_idx,           colnames(X))
  bart_idx             <- to_int(bart_idx,              colnames(X))
  re_idx               <- to_int(re_idx,                NULL)
  horseshoe_idx        <- to_int(horseshoe_idx,         NULL)
  positive_constraints <- to_int(positive_constraints,  colnames(X))
  negative_constraints <- to_int(negative_constraints,  colnames(X))

  # Offset
  if (is.null(offset)) {
    offset_mat <- matrix(0, n, p)
  } else {
    offset_mat <- matrix(offset, n, p)   # recycles n x 1 to n x p
  }

  # ── 2. PREPROCESSING (identical logic to mnlogit_rcpp) ───────────────────
  if (intercept && !any(colSums(X) == n))
    X <- cbind(X, intercept = 1)
  X_orig <- X

  if (missing(method)) {
    methods_chosen <- if (standardize) "standardize" else "none"
  } else {
    methods_chosen <- match.arg(method, several.ok = TRUE)
  }

  do_std <- "standardize" %in% methods_chosen
  do_cen <- "center"      %in% methods_chosen || do_std
  do_scl <- "scale"       %in% methods_chosen || do_std
  do_qr  <- "QR"          %in% methods_chosen

  # Continuous variables are only candidates for centering/scaling if they are not naturally bounded in [-1, 1] (e.g. shares, contrasts, bounded indices)
  cont_idx <- apply(X, 2, function(x) {
    is.numeric(x) && 
    length(unique(x)) > 2 && 
    max(abs(x), na.rm = TRUE) > 1.0
  })
  X_rescaling <- matrix(0, 2, sum(cont_idx)); rownames(X_rescaling) <- c("mean","sd")
  X_rescaling[2, ] <- 1

  if ((do_cen || do_scl) && sum(cont_idx) > 0) {
    if (do_cen) {
      X_rescaling[1, ] <- colMeans(X[, cont_idx, drop = FALSE])
      X[, cont_idx] <- sweep(X[, cont_idx, drop = FALSE], 2, X_rescaling[1, ], "-")
    }
    if (do_scl) {
      X_rescaling[2, ] <- apply(X[, cont_idx, drop = FALSE], 2, sd)
      X_rescaling[2, X_rescaling[2, ] == 0] <- 1
      X[, cont_idx] <- sweep(X[, cont_idx, drop = FALSE], 2, X_rescaling[2, ], "/")
    }
  }

  X_full_processed <- X
  X_full_for_bart  <- X

  if (intercept) {
    int_pos <- which(colnames(X) == "intercept")
    if (length(int_pos) > 0 && !(int_pos %in% linear_idx))
      linear_idx <- c(linear_idx, int_pos)
  }

  map_idx <- function(idx, full_idx) {
    if (is.null(idx)) return(NULL)
    m <- match(idx, full_idx)
    m[!is.na(m)]
  }
  re_idx               <- map_idx(re_idx,               linear_idx)
  horseshoe_idx        <- map_idx(horseshoe_idx,         linear_idx)
  positive_constraints <- map_idx(positive_constraints,  linear_idx)
  negative_constraints <- map_idx(negative_constraints,  linear_idx)

  X      <- X_full_processed[, linear_idx, drop = FALSE]
  k      <- ncol(X)
  X_orig <- X_orig[, linear_idx, drop = FALSE]

  # BART scaling
  bart_scaling <- NULL
  if (use_bart) {
    full_to_rescale <- rep(NA, ncol(X_full_processed))
    full_to_rescale[cont_idx] <- seq_len(sum(cont_idx))
    bart_mu <- rep(0, length(bart_idx)); bart_sd <- rep(1, length(bart_idx))
    for (i in seq_along(bart_idx)) {
      ri <- full_to_rescale[bart_idx[i]]
      if (!is.na(ri)) { bart_mu[i] <- X_rescaling[1, ri]; bart_sd[i] <- X_rescaling[2, ri] }
    }
    bart_scaling <- list(mu = bart_mu, sd = bart_sd)
  }

  if (do_qr) {
    qr_X <- qr(X); X_R_mat <- qr.R(qr_X); X <- qr.Q(qr_X)
  }

  # Subset X_rescaling to linear columns
  if (do_cen || do_scl) {
    full_to_rescale <- rep(NA, length(cont_idx)); names(full_to_rescale) <- names(cont_idx)
    full_to_rescale[cont_idx] <- seq_len(sum(cont_idx))
    needed_rescale_cols <- full_to_rescale[linear_idx]
    needed_rescale_cols <- needed_rescale_cols[!is.na(needed_rescale_cols)]
    X_rescaling <- X_rescaling[, needed_rescale_cols, drop = FALSE]
    cont_idx    <- cont_idx[linear_idx]
    if (sum(cont_idx) != ncol(X_rescaling))
      stop("Dimension mismatch in X_rescaling after subsetting.")
  }

  if (use_bart) {
    X_bart_final <- X_full_for_bart[, bart_idx, drop = FALSE]
    bart_ref1_mat <- matrix(apply(X_bart_final, 2, min), 1, length(bart_idx))
    bart_ref2_mat <- matrix(apply(X_bart_final, 2, max), 1, length(bart_idx))
  }

  # ── 3. INITIALISATION ─────────────────────────────────────────────────────
  cov_names_save <- colnames(X)
  X  <- as.matrix(unname(X))
  Xt <- as.matrix(unname(t(X)))
  Y  <- as.matrix(unname(Y))

  curr_beta <- matrix(0, k, p)
  bart_shifts <- matrix(0, p, 1)

  if (is.null(prior_mu))    prior_mu    <- matrix(0, k, p)
  if (!is.matrix(prior_mu)) prior_mu    <- matrix(as.vector(prior_mu), k, p)
  if (is.null(prior_V_inv)) {
    # Default: A0 for predictors, but allow wide prior (variance = 100.0, sd = 10) for the intercept (similar to brms/Stan)
    diag_vals <- rep(1 / A0, k)
    int_pos <- which(colnames(X) == "intercept")
    if (length(int_pos) > 0) {
      diag_vals[int_pos] <- 1 / 100.0
    } else if (intercept) {
      diag_vals[k] <- 1 / 100.0
    }
    prior_V_inv <- diag(diag_vals, k)
  }
  prior_P  <- as.matrix(unname(prior_V_inv))
  prior_Pb <- as.matrix(unname(prior_P %*% prior_mu))    # k x p, static

  # Dispersion parameters
  if (poisson_mode) {
    r_disp <- rep(r_fixed, p)
    sample_r <- FALSE
  } else {
    r_disp    <- if (!is.null(r_init)) as.numeric(r_init) else rep(1, p)
    sample_r  <- TRUE
  }

  # kappa will be recomputed each iteration (depends on r_disp)
  .kappa <- function(r) (Y - matrix(r, n, p, byrow = TRUE)) / 2   # n x p

  prec_beta_pooled  <- matrix(1e12, k, p)
  sigma_beta_pooled <- matrix(1e-6, k, p)

  if (use_re) {
    groups    <- unique(group_idx)
    n_groups  <- length(groups)
    idx_list  <- lapply(groups, function(m) as.integer(which(group_idx == m)))
    # POSITION of each observation's group in `groups` (APPEARANCE order). curr_beta_c /
    # postb_total slice m belongs to groups[m], NOT to raw group id m -- these differ whenever the
    # rows are not ordered by group (28 of 34 ids on the NUTS3 livestock panel). Anything that
    # indexes the per-group array by observation MUST use gpos, never group_idx.
    gpos      <- match(group_idx, groups)
    Xm_list   <- lapply(idx_list, function(idx) as.matrix(unname(X[idx, , drop = FALSE])))
    Xmt_list  <- lapply(Xm_list, t)
    group_sizes <- sapply(idx_list, length)

    # Support-aware RE prior (ported from the MNL): per (covariate x group) SD multiplier in (0,1]
    # = 1/sqrt((n_g/PR)^strength), PR = within-group participation ratio (effective # obs informing
    # that group's slope). Shrinks information-sparse RE slopes; all-ones when off; globally-constant
    # cols (random intercept) exempt. Column m aligns with Xm_list[[m]] / re_mask column m.
    re_support_mat <- matrix(1.0, nrow = k, ncol = n_groups)
    if (support_prior_strength > 0 && length(re_idx) >= 1) {
      re_glob_const <- apply(X[, re_idx, drop = FALSE], 2, function(cc) { s <- sd(cc); !is.finite(s) || s < 1e-8 })
      for (m in seq_len(n_groups)) {
        Xg <- Xm_list[[m]][, re_idx, drop = FALSE]; ng <- nrow(Xg)
        Xc <- sweep(Xg, 2, colMeans(Xg), "-")
        s2 <- colSums(Xc^2); s4 <- colSums(Xc^4); pr <- s2^2 / pmax(s4, 1e-12)
        sd_fac <- 1 / sqrt(pmax((ng / pmax(pr, 1))^support_prior_strength, 1e-12))
        sd_fac[re_glob_const] <- 1.0
        re_support_mat[re_idx, m] <- sd_fac
      }
      if (is.null(chain_id) || chain_id == 1)
        message(sprintf("RE support prior (participation ratio, strength %.2f): min SD factor %.3f.", support_prior_strength, min(re_support_mat[re_idx, ])))
    }

    # Precompute structural masking matrix for group-specific REs
    re_mask <- matrix(1.0, nrow = k, ncol = n_groups)
    is_global_intercept <- sapply(seq_len(k), function(v) {
      all(X[, v] == 1.0)
    })
    for (m in 1:n_groups) {
      for (v in seq_len(k)) {
        if (!is_global_intercept[v]) {
          # Check local variance in group m
          var_m <- var(Xm_list[[m]][, v])
          if (is.na(var_m) || var_m < 1e-12) {
            re_mask[v, m] <- 0.0
          }
        }
      }
    }

    # =================================================================
    # RE COLLINEARITY SCREEN (Parity with mnlogit_rcpp.R)
    # =================================================================
    if (length(re_idx) >= 2) {
      re_col_names <- colnames(X)[re_idx]
      n_re_vars <- length(re_idx)
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
        cm <- suppressWarnings(cor(Xm_re))
        cm[is.na(cm)] <- 0.0
        diag(cm) <- 0.0
        group_cor_list[[m]] <- cm
        hit <- abs(cm) > re_cor_threshold
        pairwise_hits <- pairwise_hits + hit
        pairwise_maxcor <- pmax(pairwise_maxcor, abs(cm))
      }

      n_valid_groups <- sum(sapply(group_cor_list, function(cm) !all(is.na(cm))))
      hit_fraction <- pairwise_hits / max(n_valid_groups, 1L)

      globally_removed <- integer(0)
      locally_masked <- list()

      for (i in seq_len(n_re_vars - 1L)) {
        for (j in (i + 1L):n_re_vars) {
          if (hit_fraction[i, j] < re_cor_min_groups) next

          var_i <- var(X[, re_idx[i]])
          var_j <- var(X[, re_idx[j]])
          drop_pos <- if (var_i <= var_j) i else j
          keep_pos <- if (drop_pos == i) j else i

          is_global <- hit_fraction[i, j] >= 0.75

          if (is_global) {
            if (!(drop_pos %in% globally_removed)) {
              globally_removed <- c(globally_removed, drop_pos)
            }
          } else {
            bad_groups <- which(sapply(
              group_cor_list,
              function(cm) !is.na(cm[i, j]) && abs(cm[i, j]) > re_cor_threshold
            ))
            for (m in bad_groups) {
              if (re_mask[re_idx[drop_pos], m] != 0.0) {
                re_mask[re_idx[drop_pos], m] <- 0.0
              }
            }
          }
        }
      }

      if (length(globally_removed) > 0L) {
        re_mask[re_idx[globally_removed], ] <- 0.0
      }
    }

    # =================================================================
    # Y-SIDE STRUCTURAL MASKING
    # =================================================================
    y_mask <- matrix(1.0, nrow = p, ncol = n_groups)
    mean_pixel_total <- mean(rowSums(Y))

    for (ip in seq_len(p)) {
      var_global <- var(Y[, ip])
      if (is.na(var_global) || var_global < 1e-8) {
        y_mask[ip, ] <- 0.0
        next
      }

      for (m in seq_len(n_groups)) {
        group_obs <- Y[idx_list[[m]], ip]
        n_positive_eff <- sum(group_obs) / mean_pixel_total
        var_local <- var(group_obs)
        std_var <- var_local / var_global

        y_thresh <- min(5L, max(2L, as.integer(ceiling(0.01 * group_sizes[m]))))
        if (n_positive_eff < y_thresh || is.na(std_var) || std_var < 1e-4) {
          y_mask[ip, m] <- 0.0
        }
      }
    }

    # =================================================================
    # PRECOMPUTE MASK COORDINATES & ADJ MASK (Parity with mnlogit_rcpp.R)
    # =================================================================
    masked_by_group_cat <- lapply(seq_len(n_groups), function(m) {
      lapply(seq_len(p), function(ip) {
        base_masked <- which((re_mask[, m] == 0.0 | (y_mask[ip, m] == 0.0 & !is_global_intercept)) & seq_len(k) %in% re_idx)
        sort(unique(base_masked))
      })
    })

    # 3D category-specific mask for precision updates
    re_mask_cube <- array(1.0, dim = c(k, p, n_groups))
    for (m_d in seq_len(n_groups)) {
      for (ip_d in seq_len(p)) {
        re_mask_cube[, ip_d, m_d] <- re_mask[, m_d]
        masked_vars <- masked_by_group_cat[[m_d]][[ip_d]]
        if (length(masked_vars) > 0L) {
          re_mask_cube[masked_vars, ip_d, m_d] <- 0.0
        }
      }
    }

    curr_beta_c <- array(0, c(k, p, n_groups))
    mu_pooled   <- matrix(0, k, p)
    # Per-chain overdispersed init (honest Rhat): uniform +/- init_jitter around the (cold) 0 start,
    # drawn in the per-chain RNG stream. A hot-start (init_state) overrides below.
    if (init_jitter > 0) mu_pooled <- mu_pooled + matrix(runif(k * p, -init_jitter, init_jitter), k, p)
    if (use_ncp) z_c <- array(0, c(k, p, n_groups))
    
    if (use_half_cauchy_re) {
      a_re <- matrix(re_scale_A^2 / 2, nrow = k, ncol = p)
    } else {
      a_re <- NULL
    }

    curr_tau_country <- rep(1.0, n_groups)
    curr_nu_country  <- rep(1.0, n_groups)

    weights <- group_sizes / sum(group_sizes)
  }

  # Horseshoe setup (identical to MNL)
  if (use_horseshoe) {
    hs_idx <- if (is.null(horseshoe_idx)) seq_len(k) else horseshoe_idx
    k_hs   <- length(hs_idx)
    
    # --- Piironen-Vehtari Calibration ---
    sigma_util <- pi / sqrt(3)
    p0_m <- p0 %||% max(1, round(k_hs * 0.3))
    p0_m <- min(p0_m, k_hs - 0.5)
    
    N_eff_pooled <- n
    tau0_pooled <- tau0_mu %||% ((p0_m / (k_hs - p0_m)) * (sigma_util / sqrt(N_eff_pooled)))
    
    if (chain_id == 1 || is.null(chain_id)) {
      cat(sprintf(
        "Horseshoe Calibration: tau0_pooled=%.4f (p0=%.1f)\n",
        tau0_pooled, p0_m
      ))
    }

    if (equation_specific_hs) {
      hs_lambda2 <- matrix(1, k, p); hs_nu <- matrix(1, k, p)
    } else {
      hs_lambda2 <- rep(1, k);       hs_nu <- rep(1, k)
    }
    hs_tau2 <- if (equation_specific_hs) rep(tau0_pooled^2, p) else tau0_pooled^2
    hs_xi   <- if (equation_specific_hs) rep(1, p) else 1
    hs_c2   <- slab_s2; hs_zeta <- 1
    hs_prec_mat <- matrix(0, k, p)
  } else {
    hs_prec_mat <- matrix(0, k, p)
  }

  # BART samplers
  if (use_bart) {
    curr_f <- matrix(0, n, p)
    bart_samplers <- vector("list", p)
    for (ip in seq_len(p)) {
      ctrl <- dbarts::dbartsControl(
        n.samples = 1, n.burn = 0, n.trees = n_trees_bart,
        n.threads = n_threads_bart, keepTrees = TRUE,
        n.chains = 1, updateState = TRUE)
      init_resp <- as.numeric(log(pmax(Y[, ip], 0.5)) - offset_mat[, ip])  # log-DENSITY init (subtract exposure); log-count init seeds the level ~11 -> inflates eta -> r collapses -> r_asis blows up U
      bart_samplers[[ip]] <- dbarts::dbarts(
        X_bart_final, init_resp, control = ctrl, sigma = 1.0,
        tree.prior = dbarts:::cgm(base = bart_base, power = bart_power),   # was UNSET -> loose defaults -> f exploded
        node.prior = dbarts:::normal(bart_k),                             # leaf shrinkage (bounds f under the exp link)
        resid.prior = chisq(df = 1e10, quant = 0.5))
    }
  }

  # Constraint setup
  has_constraints <- !is.null(positive_constraints) || !is.null(negative_constraints)
  if (has_constraints) {
    pv       <- diag(MASS::ginv(prior_P))
    psd      <- sqrt(pmax(pv, 1e-12))
    start_lo <- rep(-Inf, k); start_hi <- rep(Inf, k)
    if (!is.null(positive_constraints)) start_lo[positive_constraints] <- -2 * psd[positive_constraints]
    if (!is.null(negative_constraints)) start_hi[negative_constraints] <-  2 * psd[negative_constraints]
    decay_rate <- 5
  }

  # ── 4. STORAGE ────────────────────────────────────────────────────────────
  thin <- max(1L, as.integer(thin))
  nretain <- (niter - nburn) %/% thin        # thinned # of STORED draws
  postb_total  <- if (use_re) array(0, c(k, p, n_groups, nretain)) else array(0, c(k, p, nretain))
  postb_pooled <- array(0, c(k, p, nretain))
  post_r       <- if (sample_r) matrix(0, p, nretain) else NULL
  post_sigma_re <- if (use_re) matrix(0, k * p, nretain) else NULL
  post_tau_country <- if (use_re && isTRUE(use_country_shrinkage)) matrix(0, n_groups, nretain) else NULL
  post_log_lik <- numeric(nretain)
  post_ll_pw   <- if (calc_loo) matrix(0, nretain, n) else NULL
  post_kappa_mean <- if (use_horseshoe) matrix(0, k * p, nretain) else NULL

  if (use_bart) {
    post_f_sum    <- matrix(0, n, p)
    post_f_sum_sq <- matrix(0, n, p)
    if (store_f) post_f <- array(0, c(n, p, nretain))
    tree_store <- if (store_bart_trees || save_bart_to_disk) vector("list", nretain) else NULL
    if (save_bart_to_disk) {
      dir.create(bart_disk_path, recursive = TRUE, showWarnings = FALSE)
      bart_batch_buffer <- list(); bart_batch_files <- character()
    }
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

    if (!is.null(init_state)) {
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
        for (ip in seq_len(p)) {
          if (!is.null(init_state$bart_states[[ip]])) {
            bart_samplers[[ip]]$setState(init_state$bart_states[[ip]])
          }
        }
      }
    }
  }

  burn_buffer <- max(50, floor(0.05 * nburn))
  # r-update warmup: with BART + NB, the r-update kicks in while BART/beta are still
  # settling, so mu is TRANSIENTLY inflated -> the CRT rate sum log1p(mu/r) collapses r
  # -> -log(r) in offset_eff inflates U -> locks in a bad mode (validated: Poisson/r-fixed
  # BART is stable, NB-BART diverges). Delay the r-update until the MEAN has settled
  # (Poisson-like) so the CRT sees mu~Y and lands on the good r. Fixed-r (Poisson) needs no warmup.
  if (is.null(r_warmup)) r_warmup <- if (use_bart) max(burn_buffer, floor(0.45 * nburn)) else burn_buffer

  if (is.null(chain_id)) {
    cat(sprintf("mncount_rcpp [%s]: Sampling...\n", family))
    pb <- utils::txtProgressBar(min = 0, max = niter, style = 3)
  } else {
    cat(sprintf("Chain %d [%s]: Sampling...\n", chain_id, family))
  }

  # =====================================================================
  # GIBBS LOOP
  # =====================================================================
  for (iter in seq_len(niter)) {

    # ── Constraint walls (annealing) ──────────────────────────────────
    if (has_constraints) {
      if (iter <= nburn) {
        df     <- exp(-decay_rate * (iter / nburn))
        cur_lo <- start_lo * df;  cur_hi <- start_hi * df
        if (iter > nburn * 0.95) {
          if (!is.null(positive_constraints)) cur_lo[positive_constraints] <- 0
          if (!is.null(negative_constraints)) cur_hi[negative_constraints] <- 0
        }
      } else {
        cur_lo <- start_lo;  cur_hi <- start_hi
        if (!is.null(positive_constraints)) cur_lo[positive_constraints] <- 0
        if (!is.null(negative_constraints)) cur_hi[negative_constraints] <- 0
      }
    }

    # ── Recompute kappa (changes because r_disp may have been updated) ──
    kappa_mat <- .kappa(r_disp)   # n x p,  = (Y - r_j) / 2

    # ── MEAN-PARAMETERIZATION offset (the r-ridge fix) ──────────────────
    # Natural param: mu = r*exp(eta_nat), eta_nat = X.beta + f + offset.  Here the
    # intercept and log(r) BOTH shift log(mu)=log(r)+eta_nat -> a ridge with two
    # self-consistent CRT fixed points (r~2 sensible, r~0.07 with eta inflated).
    # BART flexibility tips the sampler onto the bad one (validated: r=0.07, U=15-21).
    # Mean param: model log(mu)=X.beta+f+offset (r-FREE) by carrying eta_nat =
    # log(mu) - log(r).  Implement by folding -log(r) into the offset used to FORM
    # eta everywhere (U, the PG tilt, the kappa offset-removal, the CRT).  Then
    # mu = r*exp(eta_nat) = exp(X.beta+f+offset) is r-free -> no ridge -> r is
    # identified by dispersion alone.  Subsumes r_asis (that post-hoc intercept shift
    # is exactly this move done discretely; it fights BART, this does not).
    offset_eff <- if (isTRUE(mean_param))
      sweep(offset_mat, 2, log(pmax(r_disp, 1e-8)), "-") else offset_mat

    # ================================================================
    # A. LINEAR PREDICTOR  eta = X*beta (+ BART + RE + offset)
    # KEY DIFFERENCE FROM MNL: no softmax; no c_j; no baseline removal.
    # ================================================================
    if (use_re) {
      U <- matrix(0, n, p)
      for (ip in seq_len(p)) {
        for (m in seq_len(n_groups)) {
          U[idx_list[[m]], ip] <- Xm_list[[m]] %*% curr_beta_c[, ip, m]
        }
        if (use_bart) U[, ip] <- U[, ip] + curr_f[, ip]
        U[, ip] <- U[, ip] + offset_eff[, ip]     # includes -log(r) under mean_param
        U[, ip] <- pmax(pmin(U[, ip], 30), -30)   # wider clip for counts
      }
    } else {
      # C++ utility (reuse from mnlogit; c_j output is ignored here)
      f_bart_mat <- if (use_bart) curr_f else matrix(0, n, p)
      U <- compute_utilities_count_cpp(X, curr_beta, f_bart_mat, offset_eff)
    }

    # ================================================================
    # B. PÓLYA-GAMMA DRAWS
    # h_ij = Y_ij + r_j  (per-observation, per-outcome)
    # z_ij = eta_ij       (linear predictor — NOT utility minus c_j)
    # ================================================================
    omega <- matrix(0, n, p)
    for (ip in seq_len(p)) {
      h_i   <- Y[, ip] + r_disp[ip]    # observation-specific
      z_i   <- U[, ip]
      omega[, ip] <- pmax(
        as.vector(fast_rpg(n, h_i, z_i, gamma_matched_pg = gamma_matched_pg)),
        1e-6)
    }

    # ================================================================
    # B2. REMOVE THE OFFSET FROM THE PG WORKING RESPONSE
    # eta = X*beta + f_bart + OFFSET, so the offset is a KNOWN additive term in the
    # linear predictor. The Gaussian beta update regresses the working response
    # kappa/omega onto X, so the offset must be subtracted from that target exactly
    # as the BART term is (target = kappa - omega*(f_bart + offset)). Doing it once
    # here -- on kappa_mat, after the omega draw -- makes EVERY downstream path
    # (pooled gibbs_step_count_cpp, RE gibbs_step_re[_ncp], the constrained loops,
    # and the BART working response at kappa_mat/omega - lp) inherit it.
    # WITHOUT this, beta (the intercept) absorbs the offset and runs away to the
    # eta clip (+/-30), collapsing the dispersion r toward 0 and nuking all slopes.
    # ================================================================
    if (any(offset_eff != 0)) kappa_mat <- kappa_mat - omega * offset_eff

    # ================================================================
    # C. COEFFICIENT SAMPLING
    # IDENTICAL structure to MNL — pass c_j_mat = 0 (no competitor term).
    # The target for beta_j becomes: kappa_j - omega_j * f_bart_j
    #   (c_j has been removed; everything else is the same).
    # ================================================================
    c_j_zero <- matrix(0, n, p)    # substitutes for c_j_mat

    # EFFECTIVE PER-(cov, group) RE SD MULTIPLIER for THIS sweep. Hoisted out of the RE branch
    # (2026-08-26) so that the COEFFICIENT DRAW and the VARIANCE UPDATE below both see the SAME
    # coordinates. Previously the draw used base x tau_country while update_re_precision_hc was
    # handed the BASE only, so sigma_v was re-inferred from deviations it did not know had been
    # scaled -- sigma_v silently absorbed 1/tau_country and the country shrinkage was counted twice.
    # Same failure class as the centred-ss RE collapse: a variance estimated in different coordinates
    # from the draw it governs.
    re_supp_eff <- re_support_mat
    if (use_re && isTRUE(use_country_shrinkage)) {
      .slopes <- which(!is_global_intercept & (seq_len(k) %in% re_idx))
      if (!isTRUE(country_shrink_intercept)) {
        if (length(.slopes)) re_supp_eff[.slopes, ] <-
          sweep(re_support_mat[.slopes, , drop = FALSE], 2, curr_tau_country, "*")
      } else {
        .all_re <- intersect(seq_len(k), re_idx)
        if (length(.all_re)) re_supp_eff[.all_re, ] <-
          sweep(re_support_mat[.all_re, , drop = FALSE], 2, curr_tau_country, "*")
      }
    }

    if (!use_re && !has_constraints) {
      curr_beta <- gibbs_step_count_cpp(
        X, Xt, kappa_mat, omega,
        prior_P, hs_prec_mat, prior_Pb,
        if (use_bart) curr_f else matrix(0, n, p),
        use_bart)

    } else if (use_re && !has_constraints) {
      sigma_mat <- matrix(1 / sqrt(pmax(prec_beta_pooled, 1e-8)), k, p)
      if (use_ncp) {
        re_res <- gibbs_step_re_ncp(
          X, Xt, kappa_mat, omega, c_j_zero,
          prior_P, hs_prec_mat, prior_Pb,
          mu_pooled, sigma_mat, z_c,
          # pp_int: in count model all columns 1..p are active
          as.integer(seq_len(p)),
          idx_list, Xm_list, Xmt_list, n_groups,
          if (use_bart) curr_f else matrix(0, n, p),
          use_bart, as.integer(re_idx - 1L),
          re_mask, y_mask, as.integer(is_global_intercept),
          # Shared base core's RE-scale-ASIS args: support prior (participation-ratio SD scaling),
          # the half-Cauchy aux, asis=FALSE (ASIS off in the count model).
          re_supp_eff,
          if (use_half_cauchy_re) a_re else matrix(1.0, k, p),
          FALSE)
        curr_beta_c <- re_res$beta_c
        mu_pooled   <- re_res$mu
        z_c         <- re_res$z_c

        # LOCATION ASIS: CP redraw of the RE means given the group TOTALS (curr_beta_c fixed => mu_g and the
        # mean mu=r*exp(U) unchanged). mu_re | {b_g}, sigma ~ N(mean_g b_g, sigma^2/G). Interweaves with the
        # NCP draw above -> breaks the pooled-level <-> RE-deviation ridge (the intercept mixes / stops drifting).
        if (isTRUE(level_asis)) {
          for (ip in seq_len(p)) {
            bbar  <- rowMeans(curr_beta_c[re_idx, ip, , drop = FALSE], dims = 2)   # per-RE-cov mean over groups
            sd_mu <- pmax(sigma_mat[re_idx, ip] / sqrt(n_groups), 1e-8)
            mu_pooled[re_idx, ip] <- rnorm(length(re_idx), bbar, sd_mu)            # curr_beta_c kept fixed
          }
        }

      } else {
        re_res <- gibbs_step_re(
          X, Xt, kappa_mat, omega, c_j_zero,
          prior_P, hs_prec_mat, prior_Pb, mu_pooled, prec_beta_pooled,
          as.integer(seq_len(p)),
          idx_list, Xm_list, Xmt_list, n_groups,
          if (use_bart) curr_f else matrix(0, n, p),
          use_bart, as.integer(re_idx - 1L),
          re_mask, y_mask, as.integer(is_global_intercept))
        curr_beta_c <- re_res$beta_c
        mu_pooled   <- re_res$mu
      }

    } else {
      # Constrained path
      if (!use_re) {
        # C1. POOLED CONSTRAINED
        for (ip in seq_len(p)) {
          om_p   <- omega[, ip]
          target <- kappa_mat[, ip]
          if (use_bart) target <- target - om_p * curr_f[, ip]
          P  <- prior_P + weighted_crossprod(Xt, X, om_p)
          if (use_horseshoe) diag(P) <- diag(P) + hs_prec_mat[, ip]
          Pb <- prior_Pb[, ip] + Xt %*% target
          lo <- rep(-Inf, k); hi <- rep(Inf, k)
          if (!is.null(positive_constraints)) lo[positive_constraints] <- cur_lo[positive_constraints]
          if (!is.null(negative_constraints)) hi[negative_constraints] <- cur_hi[negative_constraints]
          curr_beta[, ip] <- sample_tmvn_precision_gibbs_cpp(
            P, Pb, lo, hi, curr_beta[, ip],
            n_steps = if (iter <= nburn) 2L else 1L)
        }
      } else if (use_ncp) {
        # C2. CONSTRAINED RE — NCP PATH
        sigma_mat <- matrix(1 / sqrt(pmax(prec_beta_pooled, 1e-8)), nrow = k, ncol = p)
        re_complement <- setdiff(seq_len(k), re_idx)

        for (ip in seq_len(p)) {
          om_p <- omega[, ip]
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
            r_z <- kappa_mat[idx, ip] - om_m * (Xm %*% mu_pooled[, ip])
            if (use_bart) r_z <- r_z - om_m * curr_f[idx, ip]

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
            masked_re <- masked_by_group_cat[[m]][[ip]]
            if (length(masked_re) > 0) {
              curr_beta_c[masked_re, ip, m] <- mu_pooled[masked_re, ip]
              z_c[masked_re, ip, m] <- 0
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
            U_z_m <- Xm %*% (sig_m * z_c[, ip, m])
            r_mu <- kappa_mat[idx, ip] - om_m * U_z_m
            if (use_bart) r_mu <- r_mu - om_m * curr_f[idx, ip]

            P_mu_acc <- P_mu_acc + Xmt_om_Xm
            Pb_mu_acc <- Pb_mu_acc + as.vector(Xmt %*% r_mu)
          }

          # mu draw — uses walls accumulated from ALL groups
          P_mu <- prior_P + P_mu_acc
          diag(P_mu) <- diag(P_mu) + hs_prec_mat[, ip]

          mu_pooled[, ip] <- sample_tmvn_precision_gibbs_cpp(
            P       = P_mu + diag(1e-10, k),
            Pb      = prior_Pb[, ip] + Pb_mu_acc,
            lo      = mu_lo,
            hi      = mu_hi,
            init    = mu_pooled[, ip],
            n_steps = if (iter <= nburn) 3L else 2L
          )
        }
      } else {
        # C3. CONSTRAINED RE — CP PATH
        re_complement <- setdiff(seq_len(k), re_idx)

        for (ip in seq_len(p)) {
          om_p <- omega[, ip]

          resid_Pb_sum <- rep(0, k)
          P_fixed_sum <- matrix(0, k, k)

          for (m in 1:n_groups) {
            idx <- idx_list[[m]]
            if (length(idx) == 0) next
            Xm <- Xm_list[[m]]
            Xmt <- Xmt_list[[m]]
            om_m <- om_p[idx]

            target_m <- kappa_mat[idx, ip]
            if (use_bart) target_m <- target_m - om_m * curr_f[idx, ip]

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
            masked_re <- masked_by_group_cat[[m]][[ip]]
            if (length(masked_re) > 0) curr_beta_c[masked_re, ip, m] <- mu_pooled[masked_re, ip]

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

    # ASIS interweaving for NCP (unchanged from MNL)
    if (use_re && use_ncp) {
      re_complement <- setdiff(seq_len(k), re_idx)
      for (ip in seq_len(p)) {
        beta_sum <- rowSums(curr_beta_c[, ip, , drop = FALSE])
        prec_re  <- pmax(prec_beta_pooled[, ip], 1e-12)
        P_cp     <- prior_P
        diag(P_cp)[re_idx] <- diag(P_cp)[re_idx] + n_groups * prec_re[re_idx]
        if (use_horseshoe) {
          hsv <- hs_prec_mat[, ip]
          # re_center: mu of an RE covariate is the MEAN COUNTRY EFFECT, not a sparsity target --
          # shrinking it here just pushes the common effect into the (unpenalised) RE mean.
          if (isTRUE(re_center)) hsv[re_idx] <- 0
          diag(P_cp) <- diag(P_cp) + hsv
        }
        Pb_cp       <- prior_Pb[, ip]
        Pb_cp[re_idx] <- Pb_cp[re_idx] + prec_re[re_idx] * beta_sum[re_idx]
        if (length(re_complement) > 0) {
          P_cond  <- P_cp[re_idx, re_idx, drop = FALSE]
          Pb_cond <- Pb_cp[re_idx] - P_cp[re_idx, re_complement, drop = FALSE] %*% mu_pooled[re_complement, ip]
        } else { P_cond <- P_cp; Pb_cond <- Pb_cp }
        lo_cond <- rep(-Inf, length(re_idx)); hi_cond <- rep(Inf, length(re_idx))
        if (has_constraints) {
          lo_k <- rep(-Inf, k); hi_k <- rep(Inf, k)
          if (!is.null(positive_constraints)) lo_k[positive_constraints] <- cur_lo[positive_constraints]
          if (!is.null(negative_constraints)) hi_k[negative_constraints] <- cur_hi[negative_constraints]
          lo_cond <- lo_k[re_idx]; hi_cond <- hi_k[re_idx]
        }
        mu_pooled[re_idx, ip] <- sample_tmvn_precision_gibbs_cpp(
          P_cond, Pb_cond, lo_cond, hi_cond, mu_pooled[re_idx, ip],
          n_steps = if (iter <= nburn) 2L else 1L)
        sig <- 1 / sqrt(prec_re)
        for (m in seq_len(n_groups)) {
          z_c[re_idx, ip, m] <- (curr_beta_c[re_idx, ip, m] - mu_pooled[re_idx, ip]) / sig[re_idx]

          # Apply structural and category-specific masking to z_c and beta_c
          masked_re <- masked_by_group_cat[[m]][[ip]]
          if (length(masked_re) > 0) {
            z_c[masked_re, ip, m] <- 0.0
            curr_beta_c[masked_re, ip, m] <- mu_pooled[masked_re, ip]
          }

          if (length(re_complement) > 0) {
            z_c[re_complement, ip, m]         <- 0
            curr_beta_c[re_complement, ip, m] <- mu_pooled[re_complement, ip]
          }
        }
      }
    }

    # ================================================================
    # D. SIGMA UPDATE (RE)
    # ================================================================
    if (use_re) {
      target_for_prec <- curr_beta_c

      if (use_half_cauchy_re && isTRUE(re_prec_pooled)) {
        re_prec <- .count_re_precision_pooled(
          beta_c = target_for_prec, mu_pooled = mu_pooled, re_idx = re_idx,
          n_groups = n_groups, re_mask = re_mask_cube, y_mask = y_mask,
          is_intercept = is_global_intercept, prec_prev = prec_beta_pooled,
          re_scale_A = re_scale_A, re_regularize = isTRUE(re_regularize),
          slab_c2 = re_slab_c2, re_support = re_supp_eff, p = p   # SAME coords as the draw
        )
        a_re <- re_prec$a_aux
      } else if (use_half_cauchy_re) {
        re_prec <- update_re_precision_hc(
          beta_c       = target_for_prec,
          mu_pooled    = mu_pooled,
          re_idx       = as.integer(re_idx - 1L),
          n_groups     = n_groups,
          re_mask      = re_mask_cube,
          y_mask       = y_mask,
          is_intercept = as.integer(is_global_intercept),
          prec_prev    = prec_beta_pooled,
          a_aux_prev   = a_re,
          re_scale_A   = re_scale_A,
          re_regularize = isTRUE(re_regularize),
          slab_c2      = re_slab_c2,
          re_support_opt = re_supp_eff        # SAME coords as the draw (incl. tau_country)
        )
        a_re <- re_prec$a_aux
      } else {
        re_prec <- update_re_precision(
          target_for_prec,
          mu_pooled,
          as.integer(re_idx - 1L),
          n_groups,
          re_mask_cube,
          y_mask,
          as.integer(is_global_intercept),
          prior_a_re,
          prior_b_re
        )
      }
      prec_beta_pooled <- re_prec$prec
      sigma_beta_pooled <- re_prec$sigma

      # Country-Gatekeeper hierarchical shrinkage update
      if (isTRUE(use_country_shrinkage)) {
        c_shrink <- update_country_shrinkage_hc_cpp(
          beta_c = target_for_prec,
          mu_pooled = mu_pooled,
          sigma_mat = sigma_beta_pooled,
          re_idx = as.integer(re_idx - 1L),
          re_mask = re_mask_cube,
          y_mask = y_mask,
          is_intercept = as.integer(is_global_intercept),
          tau_prev = curr_tau_country,
          nu_prev = curr_nu_country,
          tau0_country = tau0_country,
          slab_c2_country = slab_c2_country,
          re_support_opt = re_support_mat,   # BASE on purpose: this kernel divides out sigma*rs to
                                             # ISOLATE tau_m; passing the tau-scaled matrix would
                                             # divide tau out twice and pin tau_m at 1.
          include_intercept = isTRUE(country_shrink_intercept)
        )
        curr_tau_country <- c_shrink$tau_country
        curr_nu_country  <- c_shrink$nu_country
        if (isTRUE(country_shrink_anchor)) {
          # geometric mean to 1; compensate sigma on exactly the rows tau multiplies, so the
          # implied prior SD (sigma_v * tau_m) is invariant to machine precision.
          .lt <- log(pmax(curr_tau_country, 1e-12))
          .gm <- exp(mean(.lt[is.finite(.lt)]))
          if (is.finite(.gm) && .gm > 0) {
            curr_tau_country <- curr_tau_country / .gm
            .rows <- if (isTRUE(country_shrink_intercept)) intersect(seq_len(k), re_idx)
                     else which(!is_global_intercept & (seq_len(k) %in% re_idx))
            if (length(.rows)) {
              prec_beta_pooled[.rows, ] <- prec_beta_pooled[.rows, , drop = FALSE] / (.gm^2)
              sigma_beta_pooled[.rows, ] <- sigma_beta_pooled[.rows, , drop = FALSE] * .gm
            }
          }
        }
      }
    }

    # ================================================================
    # D2. HORSESHOE UPDATE
    # ================================================================
    if (use_horseshoe) {
      hs <- update_horseshoe(
        if (use_re) mu_pooled else curr_beta, hs_lambda2, hs_nu, hs_tau2, hs_xi,
        hs_c2, hs_zeta, k, p, hs_idx,
        estimate_c2 = estimate_c2, slab_df = slab_df, slab_s2 = slab_s2,
        equation_specific = equation_specific_hs
      )
      hs_lambda2 <- hs$lambda2; hs_nu <- hs$nu
      hs_tau2    <- hs$tau2;    hs_xi <- hs$xi
      hs_c2 <- hs$c2; hs_zeta <- hs$zeta
      
      hs_tau2_mat <- if (equation_specific_hs) matrix(hs_tau2, k_hs, p, byrow = TRUE) else hs_tau2
      hs_lam2_sub <- if (equation_specific_hs) hs_lambda2[hs_idx, , drop = FALSE] else hs_lambda2[hs_idx]
      var_eff <- (hs_c2 * hs_lam2_sub * hs_tau2_mat) /
                 (hs_c2 + hs_lam2_sub * hs_tau2_mat)
      hs_prec_mat[hs_idx, ] <- 1 / pmax(var_eff, 1e-12)
    }

    # ================================================================
    # E. DISPERSION UPDATE via Metropolis-Hastings  (NB only)
    # The current linear predictor U[, ip] serves as eta for the r update.
    # ================================================================
    if (sample_r && iter > r_warmup) {
      for (ip in seq_len(p)) {
        r_old <- r_disp[ip]
        r_disp[ip] <- if (r_method == "crt")
          .update_r_crt(r_disp[ip], Y[, ip], U[, ip],
                        prior_shape = r_prior_shape, prior_rate = r_prior_rate)
        else
          .update_r_mh(r_disp[ip], Y[, ip], U[, ip],
                       prior_shape = r_prior_shape, prior_rate = r_prior_rate, mh_sd = r_mh_sd)
        # ASIS interweave: hold the mean mu = r * exp(U) fixed by shifting the linear predictor by
        # log(r_old/r_new) for ALL obs. Decouples r from the intercept level so r is identified by
        # dispersion alone (fixes the r<->intercept ridge -> r stops drifting to 0.2/8). Basis-agnostic
        # (X is QR-orthogonalised): find dbeta with X %*% dbeta = 1 so U shifts by exactly `shift`;
        # add it to every group's beta AND the pooled mean, preserving RE deviations.
        if (isTRUE(r_asis) && !isTRUE(mean_param) && is.finite(r_disp[ip]) && r_disp[ip] > 0 && r_old > 0) {
          shift <- log(r_old) - log(r_disp[ip])
          dvec <- shift * as.numeric(solve(crossprod(X) + diag(1e-8, ncol(X)), colSums(X)))   # dbeta: X %*% dbeta ~= 1
          if (use_re) {
            curr_beta_c[, ip, ] <- curr_beta_c[, ip, ] + dvec
            mu_pooled[, ip]     <- mu_pooled[, ip]     + dvec
          } else {
            curr_beta[, ip]     <- curr_beta[, ip]     + dvec
          }
        }
      }
    }

    # ================================================================
    # F. BART UPDATE
    # Working response for BART: y*_i = kappa_i/omega_i - X_i beta
    #   (same formula as MNL but without the c_j offset)
    # ================================================================
    if (use_bart) {
      for (ip in seq_len(p)) {
        lp <- if (use_re)
          rowSums(X * t(curr_beta_c[, ip, gpos]))   # gpos, NOT group_idx (appearance-order slices)
        else
          X %*% curr_beta[, ip]

        y_star <- pmax(pmin(
          kappa_mat[, ip] / omega[, ip] - lp,
          100), -100)

        bart_samplers[[ip]]$setResponse(as.numeric(y_star))
        bart_samplers[[ip]]$setWeights(as.numeric(omega[, ip]))
        curr_f[, ip] <- pmax(pmin(bart_samplers[[ip]]$run(0L, 1L)$train, bart_f_cap), -bart_f_cap)  # CAP RAW output first (exp-link)

        if (isTRUE(bart_absorb_mean) && iter > burn_buffer) {
          # Absorb the BART f MEAN into the intercept, then CAP the f tails. QR-AWARE: X is QR-orthogonalised,
          # so the old `curr_beta_c[int_idx,]` shift hit the WRONG coefficient -> f-mean never absorbed -> f
          # drifted -> exp(eta) exploded -> r collapsed. Use the constant direction (X %*% cdir = 1), like r_asis.
          # WARNING (default OFF): X %*% cdir = 1 only in a LEAST-SQUARES-AVERAGE sense, not POINTWISE, so under
          # the exp-link this move is NOT mean-preserving -> it LEAKS a net positive level into mu every iteration
          # (validated: ~9x over-prediction even with r fixed). For prediction only U+f matters, so leaving f
          # un-centered (bart_absorb_mean=FALSE) is exact; the intercept<->f-mean split is a benign label issue.
          shift_val        <- mean(curr_f[, ip])   # mean of the ALREADY-CAPPED f -> bounded shift into the intercept
          curr_f[, ip]     <- curr_f[, ip] - shift_val
          cdir <- as.numeric(solve(crossprod(X) + diag(1e-8, ncol(X)), colSums(X)))
          if (use_re) {
            curr_beta_c[, ip, ] <- curr_beta_c[, ip, ] + shift_val * cdir
            mu_pooled[, ip]     <- mu_pooled[, ip]     + shift_val * cdir
          } else {
            curr_beta[, ip]     <- curr_beta[, ip]     + shift_val * cdir
          }
          bart_shifts[ip, 1] <- shift_val
        } else {
          bart_shifts[ip, 1] <- 0
        }
      }
    }

    # ================================================================
    # G. STORE POSTERIOR
    # ================================================================
    if (iter > nburn && (iter - nburn) %% thin == 0L && (iter - nburn) %/% thin <= nretain) {
      s <- (iter - nburn) %/% thin           # thinned store index

      if (use_re) {
        postb_total[, , , s] <- curr_beta_c
        postb_pooled[, , s]  <- mu_pooled
        post_sigma_re[, s]   <- as.vector(sigma_beta_pooled)
        if (isTRUE(use_country_shrinkage)) {
          post_tau_country[, s] <- curr_tau_country
        }
      } else {
        postb_total[, , s]  <- curr_beta
        postb_pooled[, , s] <- curr_beta
      }

      if (sample_r) post_r[, s] <- r_disp

      if (use_horseshoe) {
        if (equation_specific_hs) {
          hs_tau2_mat <- matrix(hs_tau2, k, p, byrow = TRUE)
          post_kappa_mean[, s] <- as.vector(1 / (1 + hs_lambda2 * hs_tau2_mat))
        } else {
          post_kappa_mean[, s] <- rep(1 / (1 + hs_lambda2 * hs_tau2), p)
        }
      }

      if (use_bart) {
        post_f_sum    <- post_f_sum    + curr_f
        post_f_sum_sq <- post_f_sum_sq + curr_f^2
        if (store_f) post_f[, , s] <- curr_f
        # ── SLIM-TREE STORAGE (ported from mnlogit_rcpp_sym) ──────────────
        # Calibrate each per-target tree ensemble back to the RAW dbarts output (recover the
        # hidden Y-scale/offset via two ref points, un-standardize split points via bart_scaling,
        # bake the offset into tree-1 leaves). bart_shifts=0 here (bart_absorb_mean off) so trees
        # reproduce the UNCAPPED f -> predict_count re-applies bart_f_cap (stored in the fit).
        if ((store_bart_trees || save_bart_to_disk) && do_slim_trees) {
          state_m <- vector("list", p)
          for (ip in seq_len(p)) {
            df <- bart_samplers[[ip]]$getTrees()
            s1 <- as.numeric(bart_samplers[[ip]]$predict(bart_ref1_mat))
            s2 <- as.numeric(bart_samplers[[ip]]$predict(bart_ref2_mat))
            t1 <- as.numeric(predict_slim_bart_cpp(bart_ref1_mat, df))
            t2 <- as.numeric(predict_slim_bart_cpp(bart_ref2_mat, df))
            t_diff  <- t1 - t2
            y_scale <- if (abs(t_diff) > 1e-10) (s1 - s2) / t_diff else 1.0
            y_offset <- s1 - (t1 * y_scale)
            if (!is.null(bart_scaling)) {
              split_rows <- which(df$var > 0)
              if (length(split_rows) > 0) {
                v_indices <- df$var[split_rows]
                df$value[split_rows] <- df$value[split_rows] * bart_scaling$sd[v_indices] + bart_scaling$mu[v_indices]
              }
            }
            df$value[df$var == -1] <- df$value[df$var == -1] * y_scale
            total_offset <- y_offset - as.numeric(bart_shifts[ip, 1])
            df$value[df$tree == 1 & df$var == -1] <- df$value[df$tree == 1 & df$var == -1] + total_offset
            state_m[[ip]] <- df
          }
          tree_store[[s]] <- state_m
        }
      }

      # NB log-likelihood
      ll <- compute_loglik_nb_cpp(U, Y, r_disp)
      post_log_lik[s] <- ll$total
      if (calc_loo) post_ll_pw[s, ] <- ll$pointwise
    }

    # ── Progress reporting ─────────────────────────────────────────────
    if (!is.null(progress_cb)) {
      if (iter %% 10 == 0 || iter == niter) {
        phase <- if (iter <= nburn) "[Burn-in]" else "[Sampling]"
        progress_cb(message = sprintf(
          "Chain %s: Iteration %d / %d %s",
          ifelse(is.null(chain_id), "?", chain_id), iter, niter, phase))
      } else progress_cb()
    } else if (is.null(chain_id)) {
      utils::setTxtProgressBar(pb, iter)
    } else if (iter %% 100 == 0 || iter == niter) {
      cat(sprintf("Chain %d: %d / %d %s\n", chain_id, iter, niter,
                  if (iter <= nburn) "[Burn-in]" else "[Sampling]"))
    }
  } # end Gibbs loop

  if (is.null(chain_id)) { close(pb); cat("\nSampling finished.\n") } else
    cat(sprintf("Chain %d: Sampling finished.\n", chain_id))

  # ── 5. BACK-TRANSFORMATION (identical to MNL) ─────────────────────────
  if (!("none" %in% methods_chosen)) {
    if (do_qr) {
      for (s in seq_len(nretain)) {
        if (use_re) {
          for (m in seq_len(n_groups))
            postb_total[, , m, s] <- backsolve(X_R_mat, postb_total[, , m, s])
          postb_pooled[, , s] <- backsolve(X_R_mat, postb_pooled[, , s])
        } else {
          postb_total[, , s] <- backsolve(X_R_mat, postb_total[, , s])
        }
      }
      if (!use_re) postb_pooled <- postb_total
    }
    if ((do_cen || do_scl) && sum(cont_idx) > 0) {
      is_cont <- which(cont_idx)
      int_idx <- which(apply(X_orig, 2, var) == 0)
      if (use_re) {
        if (do_scl) for (ci in seq_along(is_cont)) {
          postb_total[is_cont[ci], , , ] <- postb_total[is_cont[ci], , , ] / X_rescaling[2, ci]
          postb_pooled[is_cont[ci], , ]  <- postb_pooled[is_cont[ci], , ]  / X_rescaling[2, ci]
        }
        if (do_cen && length(int_idx) > 0) for (s in seq_len(nretain)) {
          for (m in seq_len(n_groups)) {
            shift <- colSums(postb_total[is_cont, , m, s, drop = FALSE] * X_rescaling[1, ])
            for (ii in int_idx) postb_total[ii, , m, s] <- postb_total[ii, , m, s] - shift
          }
          shift_p <- colSums(postb_pooled[is_cont, , s, drop = FALSE] * X_rescaling[1, ])
          for (ii in int_idx) postb_pooled[ii, , s] <- postb_pooled[ii, , s] - shift_p
        }
      } else {
        if (do_scl) for (ci in seq_along(is_cont))
          postb_total[is_cont[ci], , ] <- postb_total[is_cont[ci], , ] / X_rescaling[2, ci]
        if (do_cen && length(int_idx) > 0) for (s in seq_len(nretain)) {
          shift <- colSums(postb_total[is_cont, , s, drop = FALSE] * X_rescaling[1, ])
          for (ii in int_idx) postb_total[ii, , s] <- postb_total[ii, , s] - shift
        }
        postb_pooled <- postb_total
      }
    }
  }

  # ── 6. NAMES ──────────────────────────────────────────────────────────────
  cov_names <- cov_names_save
  out_names <- colnames(Y)
  if (is.null(out_names)) out_names <- paste0("Y", seq_len(p))

  if (use_re) {
    dimnames(postb_total)  <- list(cov_names, out_names, as.character(seq_len(n_groups)), NULL)
    dimnames(postb_pooled) <- list(cov_names, out_names, NULL)
    if (!is.null(post_tau_country) && exists("group_names") && !is.null(group_names)) {
      rownames(post_tau_country) <- as.character(group_names)
    }
  } else {
    dimnames(postb_total)  <- list(cov_names, out_names, NULL)
    dimnames(postb_pooled) <- list(cov_names, out_names, NULL)
  }

  # --- Capture and Save Final State for Hot-Start ---
  horseshoe_state <- NULL
  if (use_horseshoe) {
    horseshoe_state <- list(
      lambda2 = hs_lambda2, tau2 = hs_tau2, nu = hs_nu, xi = hs_xi, zeta = hs_zeta,
      c2 = hs_c2
    )
  }

  # ── 7. RETURN ─────────────────────────────────────────────────────────────
  res <- list(
    family         = family,
    postb          = postb_total,
    postb_total    = postb_total,
    postb_pooled   = postb_pooled,
    post_sigma_re  = post_sigma_re,
    sigma_beta_pooled = sigma_beta_pooled,
    post_tau_country = post_tau_country,
    post_r         = post_r,           # <── new: posterior draws of dispersion
    r_disp_final   = r_disp,           # <── last MCMC value (useful for warmstarting)
    mean_param     = mean_param,       # <── TRUE => beta models log(mu) (r-free): predict mu=exp(X.beta+f+offset), NO r factor
    offset_used    = !is.null(offset),
    train_density_q = {                                    # <── training density (Y/exp(offset)) quantiles per target -> predict_count density_cap="auto" ceiling
      dens <- Y / exp(offset_mat)
      apply(dens, 2, quantile, probs = c(0.5, 0.9, 0.99, 0.999, 1.0), na.rm = TRUE)
    },
    tree_store     = if (use_bart) tree_store else NULL,   # <── per-draw slim trees (calibrated, raw-scale) for predict_count grid BART
    bart_f_cap     = if (use_bart) bart_f_cap else NULL,   # <── predict_count re-applies this cap per draw (trees are uncapped)
    bart_idx       = if (use_bart) bart_idx else NULL,
    post_log_lik   = post_log_lik,
    post_log_lik_pointwise = post_ll_pw,
    diagnostics    = list(loglik = mean(post_log_lik)),
    post_f_mean    = if (use_bart) post_f_sum / nretain else NULL,
    post_f_sd      = if (use_bart) sqrt(pmax(post_f_sum_sq / nretain -
                                             (post_f_sum / nretain)^2, 0)) else NULL,
    post_f         = if (use_bart && store_f) post_f else NULL,
    final_state    = list(
      beta = if (use_re) curr_beta_c else curr_beta,
      mu = if (use_re) mu_pooled else curr_beta,
      prec_beta = if (use_re) prec_beta_pooled else NULL,
      z_c = if (use_re && use_ncp) z_c else NULL,
      sigma_re = if (use_re) sigma_beta_pooled else NULL,
      a_re = if (use_re && use_half_cauchy_re) a_re else NULL,
      horseshoe = horseshoe_state,
      bart_states = if (use_bart) lapply(bart_samplers, function(s) s$state) else NULL
    )
  )

  if (calc_loo) {
    res$waic <- loo::waic(post_ll_pw)
    res$loo  <- loo::loo(post_ll_pw)
  }
  if (use_horseshoe) {
    res$horseshoe <- list(
      lambda2 = hs_lambda2, tau2 = hs_tau2,
      c2 = hs_c2,
      post_kappa_pooled = post_kappa_mean,
      post_c2 = rep(hs_c2, nretain), # to match expected structure
      nu = hs_nu, xi = hs_xi, zeta = hs_zeta
    )
  }
  return(res)
}
