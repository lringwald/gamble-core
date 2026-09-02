# =============================================================================
# Mundlak device: group-mean predictors for the random-effect prior mean.
#
#   b_c ~ N(mu + gamma' xbar_c, sigma^2)                     (prior-mean form)
#
# Marginalising b_c gives the equivalent EXPLICIT form, which is what this file
# builds because it needs no sampler change at all:
#
#   eta_i = x_i' mu  +  x_iv * (gamma_v' xbar_g(i))  +  x_i' u_g(i)
#
# For the random INTERCEPT (x_iv = 1) the added column is just xbar_g(i) — the
# country mean entering as a pixel-level covariate, i.e. the classic device.
# For a random SLOPE on v it is the interaction x_iv * xbar_g(i,j).
#
# Identification (measured on the GLOBIOM pixel design, 2026-09-02):
#   beta is identified from WITHIN-group variation, (beta+gamma) from BETWEEN.
#   With zero within-variation the two are perfectly aliased; here the RE
#   covariates are 60-93% within, so both are identified. Country means of a
#   globally-constant column (e.g. the intercept) carry NO between information
#   and are refused below.
# =============================================================================

build_mundlak_design <- function(X, group_idx, mean_cols,
                                 re_cols = "intercept",
                                 center = TRUE, scale = TRUE,
                                 prefix = "MDL", verbose = TRUE) {
  X <- as.matrix(X)
  stopifnot(!is.null(colnames(X)), length(group_idx) == nrow(X))
  miss <- setdiff(c(mean_cols, setdiff(re_cols, "intercept")), colnames(X))
  if (length(miss)) stop("build_mundlak_design: column(s) not in X: ", paste(miss, collapse = ", "))

  gl  <- unique(group_idx)                       # appearance order — matches the sampler's slices
  gi  <- match(group_idx, gl)
  xb  <- vapply(mean_cols, function(nm) tapply(X[, nm], group_idx, mean)[as.character(gl)],
                numeric(length(gl)))
  xb  <- matrix(xb, nrow = length(gl), dimnames = list(NULL, mean_cols))

  # a group mean with no between-group variation carries nothing and would be a constant column
  keep <- apply(xb, 2, function(z) { s <- stats::sd(z); is.finite(s) && s > 1e-10 })
  if (any(!keep)) {
    if (verbose) message("Mundlak: dropping group-mean column(s) with no between-group variation: ",
                         paste(mean_cols[!keep], collapse = ", "))
    xb <- xb[, keep, drop = FALSE]; mean_cols <- mean_cols[keep]
  }
  if (!ncol(xb)) stop("build_mundlak_design: no usable group-mean columns.")
  ctr <- if (center) colMeans(xb) else rep(0, ncol(xb))
  scl <- if (scale)  apply(xb, 2, stats::sd) else rep(1, ncol(xb))
  xb_s <- sweep(sweep(xb, 2, ctr, "-"), 2, pmax(scl, 1e-12), "/")

  add <- list(); nms <- character(0)
  for (rc in re_cols) {
    base <- if (identical(rc, "intercept")) rep(1, nrow(X)) else X[, rc]
    for (j in seq_len(ncol(xb_s))) {
      add[[length(add) + 1L]] <- base * xb_s[gi, j]
      nms <- c(nms, if (identical(rc, "intercept"))
                      sprintf("%s_%s", prefix, mean_cols[j])
                    else sprintf("%s_%s_x_%s", prefix, rc, mean_cols[j]))
    }
  }
  A <- matrix(unlist(add), nrow = nrow(X)); colnames(A) <- nms
  dup <- nms %in% colnames(X)
  if (any(dup)) stop("build_mundlak_design: name collision with X: ", paste(nms[dup], collapse=", "))
  Xa <- cbind(X, A)
  if (verbose) message(sprintf("Mundlak: added %d column(s) over %d group(s): %s",
                               ncol(A), length(gl), paste(nms, collapse = ", ")))
  list(X = Xa, mundlak_cols = (ncol(X) + 1L):ncol(Xa), names = nms,
       group_levels = gl, xbar = xb, center = ctr, scale = scl,
       mean_cols = mean_cols, re_cols = re_cols)
}

# Rebuild the SAME columns on new data. Reuses the stored centring/scaling so a
# prediction design matches the fitted one (recomputing means on a new sample is
# a silent train/test mismatch).
apply_mundlak_design <- function(def, X_new, group_idx_new) {
  X_new <- as.matrix(X_new)
  gi <- match(group_idx_new, def$group_levels)
  if (anyNA(gi)) stop("apply_mundlak_design: unseen group(s): ",
                      paste(unique(group_idx_new[is.na(gi)]), collapse = ", "))
  xb_s <- sweep(sweep(def$xbar, 2, def$center, "-"), 2, pmax(def$scale, 1e-12), "/")
  add <- list()
  for (rc in def$re_cols) {
    base <- if (identical(rc, "intercept")) rep(1, nrow(X_new)) else X_new[, rc]
    for (j in seq_len(ncol(xb_s))) add[[length(add)+1L]] <- base * xb_s[gi, j]
  }
  A <- matrix(unlist(add), nrow = nrow(X_new)); colnames(A) <- def$names
  cbind(X_new, A)
}
