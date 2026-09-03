# =============================================================================
# predict_gamble() — ONE entry point from a fitted model to predicted shares.
#
#   P <- predict_gamble(fit, X, group_idx = g)                    # design already assembled
#   P <- predict_gamble(fit, X, group_idx = g, Y_prev = Y0, ...)  # + lag/focal features
#
# Dispatches on what `fit` IS, so the caller does not have to know which sampler produced it:
#
#   class "nested_cut"        -> predict_nested_cut()
#   class "nested_iv"         -> predict_nested_iv()
#   count fit (has $r / NB)   -> predict_count()          [returns counts, not shares]
#   otherwise (flat MNL list) -> predict_shares() + the alternative-specific term
#
# Returns an n x J matrix of shares whose rows sum to 1 (except the count path).
#
# WHY THIS EXISTS beyond the three underlying functions:
#   1. predict_shares() ignores the alternative-specific coefficient entirely (zero references to
#      delta). An alt-spec fit scored through it silently loses delta*Z from every utility. This
#      wrapper adds that term back.
#   2. Mundlak group-mean columns are rebuilt from the fit's stored centring/scaling rather than
#      recomputed on the new sample, which would be a silent train/test mismatch.
#   3. The three predictors take different argument names for the same things.
# =============================================================================

.pg_softmax <- function(U) { m <- apply(U, 1, max); E <- exp(U - m); E / rowSums(E) }

.pg_is_count <- function(fit)
  !is.null(fit$post_r) || !is.null(fit$r) || isTRUE(fit$family %in% c("nb", "poisson"))

# Alternative-specific utility: sum_b delta_b * Z_b, matching how the block was CODED at fit time.
# scale by the STORED z_scale (recomputing it on new data rewrites the coefficient's units).
.pg_alt_eta <- function(fit, alt_spec_Z, n, J) {
  if (is.null(alt_spec_Z) || is.null(fit$post_delta)) return(NULL)
  d <- rowMeans(as.matrix(fit$post_delta))
  blocks <- if (is.matrix(alt_spec_Z)) list(alt_spec_Z) else alt_spec_Z
  meta <- fit$alt_spec_meta %||% (fit$horseshoe$alt_scale %||% NULL)
  U <- matrix(0, n, J); pos <- 1L
  for (bi in seq_along(blocks)) {
    b <- blocks[[bi]]
    Z <- if (is.list(b) && !is.null(b$Z)) as.matrix(b$Z) else as.matrix(b)
    coefk <- if (is.list(b) && !is.null(b$coef)) b$coef else
             (if (!is.null(meta[[bi]]$coef)) meta[[bi]]$coef else "shared")
    sc <- if (!is.null(meta[[bi]]$z_scale)) meta[[bi]]$z_scale else
          (if (is.list(b) && !is.null(b$scale)) b$scale else 1)
    if (nrow(Z) != n) stop(sprintf("predict_gamble: alt-spec block %d has %d rows, design has %d.",
                                   bi, nrow(Z), n))
    Zs <- Z / sc
    nd <- switch(coefk, shared = 1L, per_class = J, symmetric = J,
                 stop("predict_gamble: unknown alt-spec coef type '", coefk, "'"))
    if (pos + nd - 1L > length(d))
      stop("predict_gamble: post_delta is shorter than the declared blocks -- the fit and the blocks disagree.")
    dv <- d[pos:(pos + nd - 1L)]; pos <- pos + nd
    U <- U + if (identical(coefk, "shared")) dv * Zs else sweep(Zs, 2, dv, "*")
  }
  U
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

predict_gamble <- function(fit, X, group_idx = NULL, group_levels = NULL,
                           Y_prev = NULL, lag_builder = NULL,
                           alt_spec_Z = NULL, mundlak_def = NULL,
                           bart_cols = NULL, linear_cols = NULL,
                           type = c("mean", "predictive"), thin = 1L,
                           D = 200L, offset = NULL, ...) {
  type <- match.arg(type)
  X <- as.matrix(X)

  # ---- 1. lag / focal features from the previous state --------------------------------------
  # Temporal and spatial lags are DESIGN, not model: they must be built the same way they were at
  # fit time. Supply `lag_builder(X, Y_prev, ...)` (e.g. a closure over apply_recipe / the focal
  # emitter) rather than having this function guess a naming convention.
  if (!is.null(Y_prev)) {
    if (is.null(lag_builder))
      stop("predict_gamble: Y_prev supplied but no lag_builder. Lag/focal columns must be built ",
           "exactly as at fit time -- pass lag_builder = function(X, Y_prev) <returns X with the ",
           "lag columns>, e.g. wrapping apply_recipe().")
    X <- as.matrix(lag_builder(X, Y_prev))
  }

  # ---- 2. Mundlak group-mean columns ---------------------------------------------------------
  md <- mundlak_def %||% fit$mundlak_def
  if (!is.null(md) && !any(md$names %in% colnames(X))) {
    if (is.null(group_idx)) stop("predict_gamble: mundlak_def present but group_idx is NULL.")
    if (!exists("apply_mundlak_design")) source("codes/mundlak.R")
    X <- apply_mundlak_design(md, X, group_idx)      # reuses the FITTED centring/scaling
  }

  # ---- 3. dispatch ---------------------------------------------------------------------------
  if (inherits(fit, "nested_cut")) {
    return(predict_nested_cut(fit, X, D = D,
                              summary = if (type == "mean") "mean" else "draws",
                              group_idx = group_idx))
  }
  if (inherits(fit, "nested_iv")) return(predict_nested_iv(fit, X))
  if (.pg_is_count(fit)) {
    if (is.null(offset)) stop("predict_gamble: count fit needs `offset`.")
    return(predict_count(fit, X, offset = offset, group_idx = group_idx,
                         group_levels = group_levels, type = type, thin = thin, ...))
  }

  # ---- 4. flat MNL ---------------------------------------------------------------------------
  if (is.null(fit$postb_total)) stop("predict_gamble: fit has no postb_total -- not a fitted MNL.")
  k <- dim(fit$postb_total)[1]
  if (ncol(X) != k)
    stop(sprintf("predict_gamble: X has %d columns, the fit has %d coefficients (%s). The design ",
                 ncol(X), k, paste(head(fit$var_names %||% dimnames(fit$postb_total)[[1]], 3), collapse=", ")),
         "must be the FITTED design -- see the reconstruction gotcha in docs/.")
  if (is.null(linear_cols)) linear_cols <- seq_len(k)
  P <- predict_shares(fit, X, bart_cols = bart_cols, linear_cols = linear_cols,
                      group_idx = group_idx, group_levels = group_levels,
                      type = type, thin = thin, ...)
  # add the alternative-specific term, which predict_shares does not know about
  A <- .pg_alt_eta(fit, alt_spec_Z, nrow(X), ncol(P))
  if (!is.null(A)) P <- .pg_softmax(log(pmax(P, 1e-300)) + A)
  P
}
