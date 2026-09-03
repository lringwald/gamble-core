# =============================================================================
# predict_gamble() — ONE entry point from a fitted model to predicted shares.
#
#   P <- predict_gamble(fit, X, group_idx = g)                    # design already assembled
#   P <- predict_gamble(fit, raw = d, recipe = rc, Y_prev = Y0)   # focal/lags rebuilt from Y_prev
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

predict_gamble <- function(fit, X = NULL, group_idx = NULL, group_levels = NULL,
                           raw = NULL, recipe = NULL,
                           Y_prev = NULL, lag_builder = NULL,
                           alt_spec_Z = NULL, mundlak_def = NULL,
                           bart_cols = NULL, linear_cols = NULL,
                           type = c("mean", "predictive"), thin = 1L,
                           D = 200L, offset = NULL, ...) {
  type <- match.arg(type)

  # A `prior_model` bundle carries its own recipe -- unwrap so callers can pass either.
  if (is.null(recipe) && !is.null(fit$recipe) && !is.null(fit$fit)) { recipe <- fit$recipe; fit <- fit$fit }

  # ---- 1. build the design, recomputing lag/focal from the previous state --------------------
  # The RECIPE is the carrier for how lags were built: lu_classes + coord + res + slice define the
  # spatial neighbourhood, transforms the derived columns, col_order the exact training order. So
  # with a recipe nothing has to be inferred -- inject Y_prev as the LU state and re-run it, which
  # is what project_prior() does each step.
  if (!is.null(recipe)) {
    if (is.null(raw)) stop("predict_gamble: `recipe` given but `raw` is NULL -- the recipe rebuilds ",
                           "the design from raw data (it needs coords and the LU state).")
    d <- data.table::as.data.table(data.table::copy(raw))
    if (!is.null(Y_prev)) {
      lu <- recipe$lu_classes
      Yp <- as.matrix(Y_prev); Yp <- Yp / pmax(rowSums(Yp), 1e-12)      # shares, as at fit time
      if (nrow(Yp) != nrow(d)) stop(sprintf("predict_gamble: Y_prev has %d rows, raw has %d.", nrow(Yp), nrow(d)))
      if (ncol(Yp) != length(lu)) stop(sprintf("predict_gamble: Y_prev has %d columns, the recipe has %d LU classes (%s).",
                                               ncol(Yp), length(lu), paste(head(lu, 3), collapse=", ")))
      d[, (lu) := data.table::as.data.table(Yp)]                        # t-1 state -> focal source
    }
    a <- apply_recipe(recipe, d)
    X <- a$X
    if (is.null(group_idx))    group_idx    <- a$group_idx
    if (is.null(group_levels)) group_levels <- recipe$group_levels_appear
  } else {
    if (is.null(X)) stop("predict_gamble: supply either `X`, or `raw` + `recipe`.")
    X <- as.matrix(X)
    # No recipe: a bare fit does NOT record the focal radius/kernel/coords or the transform
    # pipeline, so those columns cannot be reconstructed from the fit alone. Take an explicit
    # builder rather than guess a naming convention.
    if (!is.null(Y_prev)) {
      if (is.null(lag_builder))
        stop("predict_gamble: Y_prev supplied with neither a `recipe` nor a `lag_builder`. Pass ",
             "recipe = <the fit's recipe> together with raw = <data>, which rebuilds focal/lags ",
             "exactly as at fit time; or lag_builder = function(X, Y_prev) for a design assembled ",
             "outside build_recipe().")
      X <- as.matrix(lag_builder(X, Y_prev))
    }
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
  # With BART, postb_total holds ONLY the linear block, so X legitimately has more columns than
  # there are linear coefficients (the BART covariates carry no linear coefficient). Check the
  # LINEAR selection against k, not the full width.
  if (is.null(linear_cols)) linear_cols <- seq_len(k)
  if (length(linear_cols) != k)
    stop(sprintf("predict_gamble: linear_cols selects %d column(s) but the fit has %d linear ",
                 length(linear_cols), k),
         sprintf("coefficient(s) (%s). ", paste(head(fit$var_names %||% dimnames(fit$postb_total)[[1]], 3), collapse=", ")),
         "The design must be the FITTED design -- see the reconstruction gotcha in docs/.")
  if (max(c(linear_cols, bart_cols)) > ncol(X))
    stop(sprintf("predict_gamble: column index %d exceeds the %d columns supplied.",
                 max(c(linear_cols, bart_cols)), ncol(X)))
  if (!is.null(fit$tree_store) && length(fit$tree_store) && is.null(bart_cols))
    stop("predict_gamble: this fit HAS BART trees but bart_cols is NULL -- the BART utility would ",
         "be silently dropped. Pass bart_cols (the fitted bart_idx).")
  P <- predict_shares(fit, X, bart_cols = bart_cols, linear_cols = linear_cols,
                      group_idx = group_idx, group_levels = group_levels,
                      type = type, thin = thin, ...)
  # add the alternative-specific term, which predict_shares does not know about
  A <- .pg_alt_eta(fit, alt_spec_Z, nrow(X), ncol(P))
  if (!is.null(A)) P <- .pg_softmax(log(pmax(P, 1e-300)) + A)
  P
}
