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

# =============================================================================
# REPORTING: the total country coefficient, decomposed
# =============================================================================
# The sampler returns b_g = mu + u_c for each covariate. Under the Mundlak expansion that is NOT
# the country's coefficient -- the contextual part sits in SEPARATE columns (MDL_*), whose value is
# the country's own mean. So reporting b_g as "the country slope" reports the country as if its
# national context were at the sample mean:
#
#   beta_c = mu           pooled within-group slope
#          + gamma' xbar_c   contextual, deterministic in the country's macro means
#          + u_c          idiosyncratic residual
#
# A country can have a NEGATIVE u_c and a well above-average total: reporting only mu + u_c inverts
# the conclusion. This returns all three pieces per country and class, with intervals taken across
# draws (not from combining separate summaries, which would understate the correlation between
# gamma and u).
#
#   post_pooled : [cov x class x draw]        (fit$postb_pooled)
#   post_total  : [cov x class x group x draw] (fit$postb_total)
#   cov_names   : rownames for the covariate axis
#   group_levels: the APPEARANCE-ORDER ids matching post_total's 3rd axis
mundlak_decompose <- function(post_pooled, post_total, def, cov_names,
                              re_col = "intercept", group_labels = NULL, probs = c(.025, .975)) {
  stopifnot(re_col %in% def$re_cols)
  i_rc <- if (identical(re_col, "intercept")) match("intercept", cov_names) else match(re_col, cov_names)
  if (is.na(i_rc)) stop("mundlak_decompose: '", re_col, "' is not a covariate in this fit.")
  # the gamma columns belonging to THIS random effect
  nm <- if (identical(re_col, "intercept")) sprintf("MDL_%s", def$mean_cols)
        else sprintf("MDL_%s_x_%s", re_col, def$mean_cols)
  i_g <- match(nm, cov_names)
  if (anyNA(i_g)) stop("mundlak_decompose: gamma column(s) missing from the fit: ",
                       paste(nm[is.na(i_g)], collapse = ", "))
  xb_s <- sweep(sweep(def$xbar, 2, def$center, "-"), 2, pmax(def$scale, 1e-12), "/")  # [group x mean_col]
  G <- dim(post_total)[3]; J <- dim(post_total)[2]; D <- dim(post_total)[4]
  if (nrow(xb_s) != G) stop("mundlak_decompose: ", nrow(xb_s), " group means vs ", G, " fitted groups.")
  labs <- group_labels %||% paste0("g", seq_len(G))
  out <- list()
  for (g in seq_len(G)) for (j in seq_len(J)) {
    mu_d  <- post_pooled[i_rc, j, ]                      # draws of mu
    bg_d  <- post_total[i_rc, j, g, ]                    # draws of mu + u_c
    u_d   <- bg_d - mu_d
    ctx_d <- as.numeric(xb_s[g, ] %*% matrix(post_pooled[i_g, j, ], nrow = length(i_g)))
    tot_d <- bg_d + ctx_d                                # = mu + u_c + gamma'xbar_c
    q <- stats::quantile(tot_d, probs, names = FALSE)
    out[[length(out) + 1L]] <- data.frame(
      group = labs[g], class = j, re_col = re_col,
      mu = mean(mu_d), contextual = mean(ctx_d), u = mean(u_d), total = mean(tot_d),
      total_lo = q[1], total_hi = q[2],
      # how much of this country's departure from the pooled slope is EXPLAINED by its context
      share_contextual = if (abs(mean(ctx_d) + mean(u_d)) > 1e-12)
                           abs(mean(ctx_d)) / (abs(mean(ctx_d)) + abs(mean(u_d))) else NA_real_,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

# Variance components: what the contextual columns EXPLAIN of the between-country spread.
# Reported against the same quantity before and after, so "the RE is smaller now" is demonstrated
# rather than asserted. var_total is the spread of the full country coefficients; var_resid the
# spread of u alone.
mundlak_variance_explained <- function(dec) {
  s <- lapply(split(dec, dec$class), function(d) data.frame(
    class = d$class[1],
    var_total = stats::var(d$total), var_contextual = stats::var(d$contextual),
    var_resid = stats::var(d$u),
    explained = if (stats::var(d$total) > 0) 1 - stats::var(d$u)/stats::var(d$total) else NA_real_))
  do.call(rbind, s)
}
