# =============================================================================
# mnlogit_nested_iv.R  —  RECURSIVE nested (hierarchical) MNL with inclusive values
# =============================================================================
# Sequential nested logit with McFadden inclusive values (IV / logsum), over an
# arbitrary-depth tree, reusing mnlogit_rcpp_sym unchanged. Restores the
# "best sub-alternative feeds the parent choice" coupling (e.g. wheat's attractiveness
# raises P(Cropland)) that plain factorized nesting severs.
#
#   P(fine) = prod over the path of  P(child | parent, X)
#   P(child c | parent) = softmax_c( X'gamma_c + lambda_c * IV_c )
#   IV_c(X)             = log sum over c's own alternatives exp( their utilities )
#
# A TREE node is either:
#   * a character vector  -> a terminal nest of FINE class names (len 1 = singleton),
#   * a named list        -> an internal nest; recurse on its children.
# e.g. AgMIP:
#   list(Cropland = list(Arable = c("Wheat",...,"Cropland_arable_other"),
#                        Permanent = c("Grapes",...,"Cropland_permanent_other")),
#        Forest   = c("Forests_primary","Forests_managed"),
#        Natural  = c("Natural_unmanaged","Natural_other"),
#        Grassland = "Grassland", Built_up_area = "Built_up_area", no_choice = "no_choice")
#
# ESTIMATION (sequential/limited-information): fit bottom-up. Each terminal multi-class
# nest and each internal node is a symmetric mnlogit_rcpp_sym fit on the pixels where that
# nest is present; its IV (logsum of its own choice utilities, defined on ALL pixels) is
# appended as a covariate to its PARENT's design. lambda = the parent's coefficient on that
# IV, rescaled by K/(K-1) out of the zero-sum symmetric parameterization (K = #siblings).
# Consistent, not fully efficient (IV is a generated regressor). Reuses the core unchanged.
# =============================================================================

.nested_iv_sampler <- function(...) mnlogit_rcpp_sym(...)
.softmax_rows <- function(U) { m <- apply(U, 1, max); E <- exp(U - m); E / rowSums(E) }

# all fine-class names under a node (recursive)
.fine_of <- function(node) if (is.character(node)) node else unlist(lapply(node, .fine_of), use.names = FALSE)

# posterior-mean utilities U (n x p_all) = X %*% mean(postb_pooled); postb_pooled is the
# zero-sum p_all coefficient array in symmetric mode -> no baseline handling needed.
.fit_utilities <- function(fit, X) {
  beta_mean <- apply(fit$postb_pooled, c(1, 2), mean)
  if (ncol(X) != nrow(beta_mean)) stop(sprintf(".fit_utilities: X %d cols vs beta %d rows", ncol(X), nrow(beta_mean)))
  X %*% beta_mean
}
.inclusive_value <- function(U) { m <- apply(U, 1, max); m + log(rowSums(exp(U - m))) }

.fit_block <- function(X, Y, group_idx = NULL, use_re = FALSE, re_idx = NULL, niter = 1500, nburn = 500) {
  .nested_iv_sampler(
    X = X, Y = Y, intercept = FALSE, symmetric = TRUE, baseline = which.max(colSums(Y)),
    niter = niter, nburn = nburn, thin = 1L,
    use_re = use_re, group_idx = if (use_re) group_idx else NULL,
    re_idx = if (use_re && !is.null(re_idx)) re_idx else 1:ncol(X),
    save_posterior_to_disk = FALSE, calc_loo = FALSE
  )
}

# ---- recursive fit of one node --------------------------------------------
# in_pixels: logical vector, pixels this node is estimated on (nest present).
# Returns: list(is_internal, singleton, fine, fit, iv[all pixels], children, child_names,
#               iv_children, lambda, macro_classes)
.fit_node <- function(node, X, Y, in_pixels, use_iv, use_re, group_idx, re_idx, niter, nburn, min_pixels) {
  if (is.character(node)) {                                   # terminal nest of fine classes
    fine <- node
    if (length(fine) == 1L) return(list(is_internal = FALSE, singleton = TRUE, fine = fine, fit = NULL, iv = NULL))
    keep <- in_pixels & (rowSums(Y[, fine, drop = FALSE]) > 0)
    if (sum(keep) < min_pixels) return(list(is_internal = FALSE, singleton = TRUE, fine = fine, fit = NULL, iv = NULL))
    Ysub <- Y[keep, fine, drop = FALSE]; Ysub <- Ysub / rowSums(Ysub)
    fit  <- .fit_block(X[keep, , drop = FALSE], Ysub, group_idx = group_idx[keep], use_re = use_re, re_idx = re_idx, niter = niter, nburn = nburn)
    fit$classes <- fine
    return(list(is_internal = FALSE, singleton = FALSE, fine = fine, fit = fit,
                iv = .inclusive_value(.fit_utilities(fit, X))))
  }
  # internal node: recurse on children, then fit the node's choice model over children
  child_names <- names(node)
  children <- setNames(vector("list", length(child_names)), child_names)
  for (cn in child_names) {
    cf <- .fine_of(node[[cn]])
    child_in <- in_pixels & (rowSums(Y[, cf, drop = FALSE]) > 0)
    children[[cn]] <- .fit_node(node[[cn]], X, Y, child_in, use_iv, use_re, group_idx, re_idx, niter, nburn, min_pixels)
  }
  Ynode <- vapply(child_names, function(cn) rowSums(Y[, .fine_of(node[[cn]]), drop = FALSE]), numeric(nrow(Y)))
  colnames(Ynode) <- child_names
  iv_children <- child_names[vapply(children, function(c) !is.null(c$iv), logical(1))]
  # Gate on use_iv HERE, not only at the cbind below. The node STORES iv_children, and both
  # .node_iv() and .predict_node() rebuild the node design from it gating on length() alone.
  # Gating just the fit-time cbind left use_iv=FALSE fits carrying a non-empty iv_children, so
  # predict appended IV columns the fit never saw -> ".fit_utilities: X 5 cols vs beta 4 rows".
  # One source of truth: if no IV column entered the fit, the node must not claim IV children.
  if (!isTRUE(use_iv)) iv_children <- character(0)
  Xnode <- X
  if (length(iv_children)) {
    ivm <- vapply(iv_children, function(cn) children[[cn]]$iv, numeric(nrow(X)))
    colnames(ivm) <- paste0("IV_", iv_children); Xnode <- cbind(X, ivm)
  }
  keep <- in_pixels & (rowSums(Ynode) > 0)
  Yk <- Ynode[keep, , drop = FALSE]; Yk <- Yk / rowSums(Yk)
  fit <- .fit_block(Xnode[keep, , drop = FALSE], Yk, group_idx = group_idx[keep], use_re = use_re, re_idx = re_idx, niter = niter, nburn = nburn)
  # lambda_c for each deep child = zero-sum coef on IV_c for child c, rescaled K/(K-1)
  lambda <- setNames(rep(NA_real_, length(iv_children)), iv_children)
  if (length(iv_children)) {
    bm <- apply(fit$postb_pooled, c(1, 2), mean); K <- length(child_names)
    for (cn in iv_children) lambda[cn] <- bm[paste0("IV_", cn) == colnames(Xnode), which(child_names == cn)] * K / (K - 1)
  }
  list(is_internal = TRUE, singleton = FALSE, fine = .fine_of(node), fit = fit,
       iv = .inclusive_value(.fit_utilities(fit, Xnode)),
       children = children, child_names = child_names, iv_children = iv_children,
       lambda = lambda, node = node)
}

# =============================================================================
# nested_iv_fit(X, Y, tree, ...)  — X must include an "intercept" column (intercept=FALSE)
# =============================================================================
nested_iv_fit <- function(X, Y, tree, use_iv = TRUE, group_idx = NULL, use_re = FALSE,
                          re_idx = NULL, niter = 1500, nburn = 500, min_pixels = 50) {
  X <- as.matrix(X); Y <- as.matrix(Y)
  stopifnot(!is.null(colnames(Y)), all(.fine_of(tree) %in% colnames(Y)))
  if (is.null(group_idx)) group_idx <- rep(1L, nrow(X))
  root <- .fit_node(tree, X, Y, rep(TRUE, nrow(X)), use_iv, use_re, group_idx, re_idx, niter, nburn, min_pixels)
  structure(list(root = root, tree = tree, fine_classes = colnames(Y), use_iv = use_iv, X_cols = colnames(X)),
            class = "nested_iv")
}

# ---- recursive inclusive value of a node on new X (rebuilds its own design) --
.node_iv <- function(nd, X_new) {
  if (!nd$is_internal) {
    if (is.null(nd$fit)) return(rep(0, nrow(X_new)))          # singleton -> degenerate IV
    return(.inclusive_value(.fit_utilities(nd$fit, X_new)))
  }
  X_node <- X_new
  if (length(nd$iv_children)) {
    ivm <- vapply(nd$iv_children, function(cn) .node_iv(nd$children[[cn]], X_new), numeric(nrow(X_new)))
    colnames(ivm) <- paste0("IV_", nd$iv_children); X_node <- cbind(X_new, ivm)
  }
  .inclusive_value(.fit_utilities(nd$fit, X_node))
}

# ---- recursive predict -----------------------------------------------------
.predict_node <- function(nd, X_new, P, share) {
  # share: n-vector of the probability mass flowing into this node (from ancestors)
  if (!nd$is_internal) {
    if (nd$singleton || is.null(nd$fit)) {                    # split evenly if a skipped multi-leaf
      for (fc in nd$fine) P[, fc] <- P[, fc] + share / length(nd$fine)
    } else {
      sub <- .softmax_rows(.fit_utilities(nd$fit, X_new)); colnames(sub) <- nd$fine
      for (fc in nd$fine) P[, fc] <- P[, fc] + share * sub[, fc]
    }
    return(P)
  }
  X_node <- X_new
  if (length(nd$iv_children)) {                               # rebuild the node's design with child IVs
    ivm <- vapply(nd$iv_children, function(cn) .node_iv(nd$children[[cn]], X_new), numeric(nrow(X_new)))
    colnames(ivm) <- paste0("IV_", nd$iv_children); X_node <- cbind(X_new, ivm)
  }
  cp <- .softmax_rows(.fit_utilities(nd$fit, X_node)); colnames(cp) <- nd$child_names
  for (cn in nd$child_names) P <- .predict_node(nd$children[[cn]], X_new, P, share * cp[, cn])
  P
}

predict_nested_iv <- function(fit, X_new) {
  X_new <- as.matrix(X_new)
  P <- matrix(0, nrow(X_new), length(fit$fine_classes), dimnames = list(NULL, fit$fine_classes))
  .predict_node(fit$root, X_new, P, rep(1, nrow(X_new)))
}

# ---- extract per-level symmetric (zero-sum) coefficient tables + lambdas ----
# Returns a named list: one k x K coefficient matrix per fitted node (rows=covariates incl IV,
# cols=that node's choices), plus $lambda per internal node. Feed a matrix to the existing
# heatplot; rotation (if desired) applies per matrix exactly as for the flat model.
nested_iv_params <- function(fit) {
  out <- list()
  walk <- function(nd, path) {
    if (!is.null(nd$fit)) {
      bm <- apply(nd$fit$postb_pooled, c(1, 2), mean)
      cn <- if (nd$is_internal) nd$child_names else nd$fine
      dimnames(bm) <- list(if (!is.null(nd$fit$cov_names)) NULL else NULL, cn)
      out[[path]] <<- list(coef = bm, classes = cn, lambda = if (nd$is_internal) nd$lambda else NULL)
    }
    if (isTRUE(nd$is_internal)) for (c in nd$child_names) walk(nd$children[[c]], paste0(path, "/", c))
  }
  walk(fit$root, "root")
  out
}
