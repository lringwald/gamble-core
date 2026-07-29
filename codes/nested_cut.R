# =============================================================================
# nested_cut.R — CUT / multiple-imputation nested MNL with inclusive values
# =============================================================================
# One clean framework for the nested model VARIANTS, all via (tree, use_iv):
#   * flat MNL        : tree = list(all = c(<every fine class>))         -> one softmax
#   * factorized nest : tree = <hierarchy>, use_iv = FALSE                -> levels independent
#   * IV-nested       : tree = <hierarchy>, use_iv = TRUE                 -> logsum coupling + lambda
#
# Improves on the sequential mean-plug-in (mnlogit_nested_iv.R): instead of feeding
# ONE inclusive value (posterior mean) up the tree, we carry M posterior DRAWS of each
# sub-model through the inclusive value. The parent (gamma, lambda) is then fit under
# each imputed IV and the draws are POOLED (Rubin's multiple imputation) — a "cut"
# posterior: information flows leaf -> parent, NOT back (no feedback), but leaf
# uncertainty is honestly propagated into lambda, gamma and the predictive.
# Reuses mnlogit_rcpp_sym UNCHANGED (no C++/offset changes). PG-fast.
#
# TREE node: character vector = terminal nest of fine-class names (len 1 = singleton);
#            named list       = internal nest, recurse on children.
#
#   P(fine) = prod over path P(child|parent);  P(child c|parent)=softmax(X'gamma_c + lambda_c IV_c)
#   IV_c(X) = logsumexp of c's own choice utilities;   lambda = zero-sum IV coef * K/(K-1)
# =============================================================================

.ncut_softmax <- function(U) { m <- apply(U, 1, max); E <- exp(U - m); E / rowSums(E) }
.ncut_lse     <- function(U) { m <- apply(U, 1, max); m + log(rowSums(exp(U - m))) }
.ncut_fine    <- function(node) if (is.character(node)) node else unlist(lapply(node, .ncut_fine), use.names = FALSE)

# --- subsample M posterior beta draws [k x p_all] from a symmetric fit -------
# postb_pooled is [k, p_all, ndraw(, nchain)] zero-sum coeffs (no baseline handling).
.ncut_beta_draws <- function(fit, M) {
  B <- fit$postb_pooled
  d <- dim(B)
  if (length(d) == 4L) { B <- array(B, c(d[1], d[2], d[3] * d[4])); d <- dim(B) }
  k <- d[1]; K <- d[2]; nd <- d[3]
  idx <- if (nd <= M) rep_len(seq_len(nd), M) else round(seq(1, nd, length.out = M))
  lapply(idx, function(s) matrix(B[, , s], nrow = k, ncol = K))       # each [k x p_all], guaranteed matrix
}

.ncut_fit_block <- function(X, Y, group_idx, use_re, re_idx, niter, nburn, thin, init_state = NULL, disk_path = NULL) {
  mnlogit_rcpp_sym(
    X = X, Y = Y, intercept = FALSE, symmetric = TRUE, baseline = which.max(colSums(Y)),
    niter = niter, nburn = nburn, thin = thin,
    use_re = use_re, group_idx = if (use_re) group_idx else NULL,
    re_idx = if (use_re && !is.null(re_idx)) re_idx else 1:ncol(X),
    save_posterior_to_disk = !is.null(disk_path), disk_path = disk_path %||% tempdir(), calc_loo = FALSE,
    support_prior_strength = 2, re_regularize = TRUE, init_jitter = 0.1,
    init_state = init_state,                              # HOT-START across imputations (see below)
    progress_cb = function(...) invisible(NULL)          # silence per-fit progress bar (many fits)
  )
}

# read `ndraws` zero-sum POOLED coefficient draws from streamed batches. Never materialises the big
# per-group array (peak = one batch, independent of niter). Expands each baseline-removed mu to the
# symmetric p_all space (same as the in-RAM postb_pooled), so downstream is identical.
.ncut_pooled_from_disk <- function(disk_path, ndraws) {
  meta <- qs2::qs_read(file.path(disk_path, "model_metadata.qs"))
  p_all <- length(meta$cat_names); bl <- which(meta$cat_names == (meta$baseline_name %||% ""))[1]
  if (is.na(bl)) bl <- p_all; pp <- (seq_len(p_all))[-bl]
  files <- list.files(disk_path, "posterior_batch_.*\\.qs$", full.names = TRUE)  # one fit per tempdir -> one chain
  files <- files[order(as.integer(gsub(".*batch_([0-9]+)_.*", "\\1", basename(files))))]
  mus <- list(); for (f in files) { b <- qs2::qs_read(f); mus <- c(mus, lapply(b, `[[`, "mu")); rm(b) }  # mu is small; betas discarded per batch
  idx <- if (length(mus) <= ndraws) rep_len(seq_along(mus), ndraws) else round(seq(1, length(mus), length.out = ndraws))
  lapply(idx, function(i) { full <- matrix(0, nrow(mus[[i]]), p_all); full[, pp] <- mus[[i]]; sweep(full, 1, rowMeans(full), "-") })
}

# fit a node's sub-MNL and return `ndraws` zero-sum pooled draws [k x p_all] + final_state (hot-start).
# stream_disk=TRUE keeps the sampler's big per-group array off RAM (safer at scale / high niter).
.ncut_fit_draws <- function(X, Y, group_idx, use_re, re_idx, niter, nburn, thin, ndraws,
                            init_state = NULL, stream_disk = FALSE) {
  if (isTRUE(stream_disk)) {
    td <- tempfile("ncut_"); dir.create(td)
    fit <- .ncut_fit_block(X, Y, group_idx, use_re, re_idx, niter, nburn, thin, init_state, disk_path = td)
    out <- list(draws = .ncut_pooled_from_disk(td, ndraws), final_state = fit$final_state)
    unlink(td, recursive = TRUE); out
  } else {
    fit <- .ncut_fit_block(X, Y, group_idx, use_re, re_idx, niter, nburn, thin, init_state)
    list(draws = .ncut_beta_draws(fit, ndraws), final_state = fit$final_state)
  }
}

# split-Rhat / bulk-ESS across chains for a set of per-chain pooled-coefficient draws.
# chain_draws: list length n_chains; each a list of [k x p_all] matrices (the chain's draws).
.ncut_rhat <- function(chain_draws) {
  nc <- length(chain_draws); nd <- min(vapply(chain_draws, length, 0L)); if (nc < 2L || nd < 2L) return(NULL)
  dm <- dim(chain_draws[[1]][[1]]); k <- dm[1]; p <- dm[2]
  arr <- array(0, c(k, p, nd, nc))
  for (c in seq_len(nc)) for (s in seq_len(nd)) arr[, , s, c] <- chain_draws[[c]][[s]]
  rh <- ess <- numeric(k * p); i <- 0L
  for (a in seq_len(k)) for (b in seq_len(p)) { i <- i + 1L; m <- matrix(arr[a, b, , ], nd, nc)
    rh[i]  <- tryCatch(posterior::rhat(m),     error = function(e) NA_real_)
    ess[i] <- tryCatch(posterior::ess_bulk(m), error = function(e) NA_real_) }
  list(max_rhat = max(rh, na.rm = TRUE), frac_rhat_lt_1.01 = mean(rh < 1.01, na.rm = TRUE),
       frac_rhat_lt_1.1 = mean(rh < 1.1, na.rm = TRUE), median_ess = median(ess, na.rm = TRUE), n_params = i)
}

# run n_chains sub-fits (cold, different seeds) for a CONVERGENCE CHECK: pool their draws for the model
# and compute Rhat/ESS across chains. Used on leaves and on each nest evaluated at the MEAN inclusive value.
.ncut_fit_chains <- function(X, Y, group_idx, use_re, re_idx, niter, nburn, thin, ndraws, n_chains,
                             stream_disk = FALSE, conv_draws = 150L, n_cores = 1L) {
  one <- function(c) .ncut_fit_draws(X, Y, group_idx, use_re, re_idx, niter, nburn, thin, conv_draws,
                                     stream_disk = stream_disk)               # cold, independent chain
  nco <- max(1L, min(as.integer(n_cores), n_chains))
  # fork the independent chains when n_cores>1 (mclapply: COW-shared X/Y, no Windows fork -> serial fallback)
  reslist <- if (nco > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(seq_len(n_chains), one, mc.cores = nco, mc.preschedule = FALSE)
  else lapply(seq_len(n_chains), one)
  bad <- vapply(reslist, function(r) inherits(r, "try-error") || is.null(r$draws), logical(1))
  if (any(bad)) stop(".ncut_fit_chains: ", sum(bad), "/", n_chains, " chain(s) failed (see mclapply warnings)")
  chains <- lapply(reslist, `[[`, "draws"); fs <- reslist[[1]]$final_state
  allp <- unlist(chains, recursive = FALSE)                       # pool all chains' draws
  model <- allp[round(seq(1, length(allp), length.out = ndraws))] # subsample ndraws for the model
  list(draws = model, final_state = fs, conv = .ncut_rhat(chains))
}

# --- moment carrier: sigma-point IV fields (mean + leading principal directions) ---
# The M coherent IV fields are low-rank (driven by the leaf's parameter uncertainty), so
# their variation is captured by a few principal directions. Represent that distribution by
# unscented sigma points: the mean, plus mu +/- sqrt(3)*d_q along the top `Q` 1-SD perturbation
# fields d_q. Pool weights {1-Q/3, 1/6,...} reproduce the mean and leading-direction variances
# exactly (3-pt Gauss-Hermite per direction) => 1 + 2Q fits instead of M, same first two IV moments.
.ncut_sigma_iv <- function(iv_by_m, Q) {
  M <- length(iv_by_m); nr <- nrow(iv_by_m[[1]]); nc <- ncol(iv_by_m[[1]]); cn <- colnames(iv_by_m[[1]])
  Fm <- vapply(iv_by_m, as.vector, numeric(nr * nc))          # [(nr*nc) x M]
  mu <- rowMeans(Fm); Fc <- Fm - mu; Q <- max(1L, min(Q, M - 1L))
  ev <- eigen(crossprod(Fc), symmetric = TRUE)                # M x M gram (cheap)
  fields <- list(matrix(mu, nr, nc)); weights <- (1 - Q / 3)
  for (q in seq_len(Q)) {
    dq <- (Fc %*% ev$vectors[, q]) / sqrt(M - 1)              # [(nr*nc)] 1-SD perturbation field
    fields <- c(fields, list(matrix(mu + sqrt(3) * dq, nr, nc), matrix(mu - sqrt(3) * dq, nr, nc)))
    weights <- c(weights, 1/6, 1/6)
  }
  fields <- lapply(fields, function(f) { colnames(f) <- cn; f })
  list(fields = fields, weights = weights / sum(weights))
}

# --- Gaussianity diagnostic on the IV field (gates "auto" mode) ---------------
# The moment carrier assumes the IV field's LEADING-DIRECTION scores are ~Gaussian.
# Quasi-separation (sparse crop x country cells) breeds skew, so CHECK it (free — the
# M IV fields are already computed): the M scores on PC q are sqrt(lambda_q)*v_q; test
# their skew / excess kurtosis. Gaussian if both within tolerance -> moments is safe;
# else fall back to the exact draws cut for that node. (Sampling SD of skew ~ sqrt(6/M),
# of kurtosis ~ sqrt(24/M); default tols are ~2 SD at M~30 to avoid false fallbacks.)
.ncut_iv_gaussianity <- function(iv_by_m, Q, skew_tol = 1.0, kurt_tol = 2.0) {
  M <- length(iv_by_m); nr <- nrow(iv_by_m[[1]]); nc <- ncol(iv_by_m[[1]])
  Fm <- vapply(iv_by_m, as.vector, numeric(nr * nc)); Fc <- Fm - rowMeans(Fm)
  Q <- max(1L, min(Q, M - 1L)); ev <- eigen(crossprod(Fc), symmetric = TRUE)
  sk <- ku <- numeric(Q)
  for (q in seq_len(Q)) {
    s <- sqrt(pmax(ev$values[q], 0)) * ev$vectors[, q]; s <- s - mean(s); v <- mean(s^2)
    if (v < 1e-12) next
    sk[q] <- mean(s^3) / v^1.5; ku[q] <- mean(s^4) / v^2 - 3
  }
  list(gaussian = max(abs(sk)) <= skew_tol && max(abs(ku)) <= kurt_tol,
       max_skew = max(abs(sk)), max_kurt = max(abs(ku)))
}

# --- nested-level progress (the sub-fits' own bars are silenced) ------------
.ncut_prog <- new.env(parent = emptyenv())
.ncut_prog_init <- function(total, on) { .ncut_prog$on <- isTRUE(on); .ncut_prog$total <- total
  .ncut_prog$done <- 0L; .ncut_prog$t0 <- Sys.time(); .ncut_prog$tlast <- .ncut_prog$t0 }
.ncut_prog_tick <- function(label) {
  if (!isTRUE(.ncut_prog$on)) return(invisible())
  now <- Sys.time(); .ncut_prog$done <- .ncut_prog$done + 1L
  dt <- as.numeric(difftime(now, .ncut_prog$tlast, units = "secs")); .ncut_prog$tlast <- now
  el <- as.numeric(difftime(now, .ncut_prog$t0, units = "mins")); rate <- el / .ncut_prog$done
  eta <- max(0, rate * (.ncut_prog$total - .ncut_prog$done))
  message(sprintf("[nested_cut] fit %d/~%d  %-30s %5.1fs | elapsed %.1fm | ~ETA %.1fm",
                  .ncut_prog$done, .ncut_prog$total, label, dt, el, eta))
}
# estimate total sub-fits for the ETA denominator (upper bound: draws/auto -> M per IV nest, moments -> 1+2Q)
.ncut_count_fits <- function(node, M, iv_mode, rank, n_chains = 1L) {   # counts progress TICKS, not raw fits
  if (is.character(node)) return(if (length(node) > 1L) 1L else 0L)                 # 1 tick per multi-class leaf
  child <- sum(vapply(node, .ncut_count_fits, 0L, M = M, iv_mode = iv_mode, rank = rank, n_chains = n_chains))
  iv_kids <- sum(vapply(node, function(ch) if (is.character(ch)) length(ch) > 1L else TRUE, logical(1)))
  node_ticks <- if (iv_kids == 0L) 1L else if (identical(iv_mode, "moments")) 1L + 2L * rank else as.integer(M)
  child + node_ticks + (if (n_chains > 1L) 1L else 0L)                              # +1 convergence tick per nest
}
# ONE node's own ticks (excludes children) — used to advance the counter when a node loads from cache.
.ncut_own_ticks <- function(node, M, iv_mode, rank, n_chains = 1L) {
  if (is.character(node)) return(if (length(node) > 1L) 1L else 0L)
  iv_kids <- sum(vapply(node, function(ch) if (is.character(ch)) length(ch) > 1L else TRUE, logical(1)))
  nt <- if (iv_kids == 0L) 1L else if (identical(iv_mode, "moments")) 1L + 2L * rank else as.integer(M)
  nt + (if (n_chains > 1L) 1L else 0L)
}
.ncut_prog_skip <- function(k, label) {                       # advance the counter for cache-loaded ticks
  if (!isTRUE(.ncut_prog$on)) return(invisible())
  .ncut_prog$done <- .ncut_prog$done + as.integer(k); .ncut_prog$tlast <- Sys.time()
  message(sprintf("[nested_cut] %-30s  <- resumed from disk (+%d)", label, as.integer(k)))
}

# --- disk store (persist node draws + auto-resume) --------------------------
# store_dir holds nodes/<path>.rds (one per completed leaf/nest, children stripped) + manifest.rds
# (a config fingerprint). On relaunch with a MATCHING config, completed nodes load instead of refitting.
.ncut_path_key <- function(path) gsub("[^A-Za-z0-9._-]+", "__", path)      # "root/Forests" -> "root__Forests"
.ncut_config_hash <- function(X, Y, tree, niter, nburn, thin, M, use_re, iv_mode, n_chains,
                              draws_per_impute, moment_rank) {
  Xf <- X[is.finite(X)]
  key <- list(cols = colnames(X), n = nrow(X), p = ncol(X), ycols = colnames(Y),
              sx = sum(Xf), sx2 = sum(Xf^2), sy = sum(Y), niter = niter, nburn = nburn, thin = thin,
              M = M, dpi = draws_per_impute, use_re = isTRUE(use_re), iv_mode = iv_mode,
              n_chains = n_chains, mrank = moment_rank, tree = paste(sort(.ncut_fine(tree)), collapse = "|"))
  if (requireNamespace("digest", quietly = TRUE)) digest::digest(key)
  else paste0("h", format(sum(utf8ToInt(paste(unlist(lapply(key, paste, collapse = ",")), collapse = "||"))),
                          scientific = FALSE), "_", key$n, "x", key$p, "_",
              format(round(key$sx, 2), scientific = FALSE), "_", format(round(key$sy, 2), scientific = FALSE))
}
# inspect a store: which nodes are already cached (skipped on the next resume)?
nested_cut_store_status <- function(store_dir) {
  nd <- file.path(store_dir, "nodes"); mf <- file.path(store_dir, "manifest.rds")
  if (!dir.exists(nd)) { message("no store at ", store_dir); return(invisible(NULL)) }
  man <- if (file.exists(mf)) readRDS(mf) else NULL
  done <- sub("\\.rds$", "", list.files(nd, "\\.rds$"))
  cat(sprintf("store: %s\n  created: %s\n  cached nodes (%d): %s\n", store_dir,
              man$created %||% "?", length(done), paste(gsub("__", "/", done), collapse = ", ")))
  invisible(done)
}

# --- recursive CUT fit of one node ------------------------------------------
# in_pixels: logical, pixels this node is estimated on. M: # imputations.
# Returns a node object; internal nodes hold, per imputation m, that m's parent fit draws.
.ncut_node <- function(node, X, Y, in_pixels, use_iv, use_re, group_idx, re_idx,
                       M, niter, nburn, thin, min_pixels, draws_per_impute, nburn_warm,
                       iv_mode = "draws", moment_rank = 1L, moment_skew_tol = 1.0, moment_kurt_tol = 2.0,
                       stream_disk = FALSE, n_chains = 1L, n_cores = 1L, path = "root", store = NULL) {
  cache_file <- if (!is.null(store)) file.path(store$dir, paste0(.ncut_path_key(path), ".rds")) else NULL
  if (!is.null(cache_file) && file.exists(cache_file)) {                    # ---- resume: load this node ----
    cached <- readRDS(cache_file)
    if (isTRUE(cached$type == "nest")) {                                     # rebuild children (each cached too)
      kids <- setNames(vector("list", length(cached$child_names)), cached$child_names)
      for (cn in cached$child_names) {
        cf <- .ncut_fine(node[[cn]]); child_in <- in_pixels & (rowSums(Y[, cf, drop = FALSE]) > 0)
        kids[[cn]] <- .ncut_node(node[[cn]], X, Y, child_in, use_iv, use_re, group_idx, re_idx,
                                 M, niter, nburn, thin, min_pixels, draws_per_impute, nburn_warm,
                                 iv_mode, moment_rank, moment_skew_tol, moment_kurt_tol,
                                 stream_disk, n_chains, n_cores, path = paste0(path, "/", cn), store = store)
      }
      cached$children <- kids
    }
    .ncut_prog_skip(.ncut_own_ticks(node, M, iv_mode, moment_rank, n_chains), sprintf("resume {%s}", path))
    return(cached)
  }
  if (is.character(node)) {                                     # ---- terminal nest ----
    fine <- node
    if (length(fine) == 1L)
      return(list(type = "singleton", fine = fine))
    keep <- in_pixels & (rowSums(Y[, fine, drop = FALSE]) > 0)
    if (sum(keep) < min_pixels)                                  # too sparse -> degrade to even split
      return(list(type = "even", fine = fine))
    Ysub <- Y[keep, fine, drop = FALSE]; Ysub <- Ysub / rowSums(Ysub)
    res <- if (n_chains > 1L)
      .ncut_fit_chains(X[keep, , drop = FALSE], Ysub, group_idx[keep], use_re, re_idx, niter, nburn, thin, M, n_chains, stream_disk, n_cores = n_cores)
    else .ncut_fit_draws(X[keep, , drop = FALSE], Ysub, group_idx[keep], use_re, re_idx, niter, nburn, thin, M, stream_disk = stream_disk)
    .ncut_prog_tick(sprintf("leaf {%s}%s", paste(fine, collapse = ","), if (n_chains > 1L) sprintf(" x%d", n_chains) else ""))
    res_leaf <- list(type = "leaf", fine = fine, classes = fine, beta_draws = res$draws, conv = res$conv)
    if (!is.null(cache_file)) saveRDS(res_leaf, cache_file)                  # persist -> resumable
    return(res_leaf)
  }
  # ---- internal node ----
  child_names <- names(node)
  children <- setNames(vector("list", length(child_names)), child_names)
  for (cn in child_names) {
    cf <- .ncut_fine(node[[cn]])
    child_in <- in_pixels & (rowSums(Y[, cf, drop = FALSE]) > 0)
    children[[cn]] <- .ncut_node(node[[cn]], X, Y, child_in, use_iv, use_re, group_idx, re_idx,
                                 M, niter, nburn, thin, min_pixels, draws_per_impute, nburn_warm,
                                 iv_mode, moment_rank, moment_skew_tol, moment_kurt_tol, stream_disk, n_chains, n_cores,
                                 path = paste0(path, "/", cn), store = store)
  }
  # aggregated child membership Y at this node
  Ynode <- vapply(child_names, function(cn) rowSums(Y[, .ncut_fine(node[[cn]]), drop = FALSE]), numeric(nrow(Y)))
  colnames(Ynode) <- child_names
  keep <- in_pixels & (rowSums(Ynode) > 0)
  # which children can supply an IV (multi-class leaf or internal, i.e. non-degenerate)
  iv_children <- child_names[vapply(children, function(c) c$type %in% c("leaf", "nest"), logical(1))]
  K <- length(child_names)

  # per-imputation IV of each iv-child on the train grid (M coherent fields)
  iv_by_m <- if (use_iv && length(iv_children))
    lapply(seq_len(M), function(m) {
      ivm <- vapply(iv_children, function(cn) .ncut_node_iv(children[[cn]], X, m), numeric(nrow(X)))
      if (is.null(dim(ivm))) ivm <- matrix(ivm, ncol = length(iv_children))
      colnames(ivm) <- paste0("IV_", iv_children); ivm }) else NULL

  # ---- resolve this node's mode; "auto" = moments only where the IV field is ~Gaussian ----
  # draws  : one fit per imputation (M fits) -> exact cut.
  # moments: mean + leading principal directions (sigma points) -> 1 + 2*rank fits, 2nd-order.
  # auto   : DIAGNOSE the IV field (skew/kurtosis of leading-direction scores); moments if it
  #          passes, else fall back to draws (e.g. sparse/quasi-separated leaves -> non-Gaussian IV).
  node_mode <- iv_mode
  if (iv_mode == "auto") {
    if (is.null(iv_by_m)) { node_mode <- "moments" } else {
      g <- .ncut_iv_gaussianity(iv_by_m, moment_rank, moment_skew_tol, moment_kurt_tol)
      node_mode <- if (g$gaussian) "moments" else "draws"
      message(sprintf("[nested_cut] node {%s}: IV |skew|=%.2f |exkurt|=%.2f -> %s%s",
        paste(child_names, collapse = ","), g$max_skew, g$max_kurt, node_mode,
        if (!g$gaussian) "  (fallback: non-Gaussian IV)" else ""))
    }
  }
  if (is.null(iv_by_m))            { fit_fields <- list(NULL); fit_w <- 1 }          # no IV -> single fit
  else if (node_mode == "moments") { sg <- .ncut_sigma_iv(iv_by_m, moment_rank); fit_fields <- sg$fields; fit_w <- sg$weights }
  else                             { fit_fields <- iv_by_m; fit_w <- rep(1 / M, M) } # exact draws

  # ---- fit at each field, hot-started across fields (m=1 cold, rest warm) ----
  fit_draws <- vector("list", length(fit_fields)); warm <- NULL
  for (j in seq_along(fit_fields)) {
    Xnode <- if (is.null(fit_fields[[j]])) X else cbind(X, fit_fields[[j]])
    Yk <- Ynode[keep, , drop = FALSE]; Yk <- Yk / rowSums(Yk)
    nb_j <- if (j == 1L || is.null(warm)) nburn else nburn_warm
    res_j <- .ncut_fit_draws(Xnode[keep, , drop = FALSE], Yk, group_idx[keep], use_re, re_idx,
                             nb_j + (niter - nburn), nb_j, thin, draws_per_impute, init_state = warm, stream_disk = stream_disk)
    warm <- res_j$final_state
    Bd <- res_j$draws; kk <- nrow(Bd[[1]]); Kfit <- ncol(Bd[[1]])
    pa <- array(unlist(Bd), c(kk, Kfit, length(Bd)))
    if (Kfit == length(child_names)) dimnames(pa) <- list(colnames(Xnode), child_names, NULL)
    fit_draws[[j]] <- pa
    .ncut_prog_tick(sprintf("nest {%s} %d/%d", paste(child_names, collapse = ","), j, length(fit_fields)))
  }

  # ---- resolve to an M-length parent-draw list (resample sigma fits by weight) ----
  parent_draws <- if (length(fit_draws) == M && node_mode != "moments") fit_draws
    else fit_draws[sample(seq_along(fit_draws), M, replace = TRUE, prob = fit_w)]

  # ---- lambda per iv-child (IV rows follow the X columns, in iv_children order) ----
  lambda_draws <- setNames(vector("list", length(iv_children)), iv_children)
  if (use_iv && length(iv_children)) { Pn <- ncol(X)
    for (cn in iv_children) { r <- Pn + match(cn, iv_children); cc <- which(child_names == cn)
      lambda_draws[[cn]] <- unlist(lapply(parent_draws, function(A) A[r, cc, ])) * K / (K - 1) }
  }
  # ---- convergence check (n_chains>1): fit this nest at the MEAN inclusive value across chains -> Rhat/ESS.
  # (The M imputation fits above stay single-chain; this separate multi-chain fit checks the sampler mixes.)
  node_conv <- NULL
  if (n_chains > 1L) {
    Xc <- if (!is.null(iv_by_m)) { ivm <- Reduce(`+`, iv_by_m) / length(iv_by_m)
                                   colnames(ivm) <- paste0("IV_", iv_children); cbind(X, ivm) } else X
    Yc <- Ynode[keep, , drop = FALSE]; Yc <- Yc / rowSums(Yc)
    rc <- .ncut_fit_chains(Xc[keep, , drop = FALSE], Yc, group_idx[keep], use_re, re_idx,
                           niter, nburn, thin, draws_per_impute, n_chains, stream_disk, n_cores = n_cores)
    node_conv <- rc$conv
    .ncut_prog_tick(sprintf("nest {%s} conv x%d", paste(child_names, collapse = ","), n_chains))
  }
  res_nest <- list(type = "nest", node = node, child_names = child_names, children = children,
       iv_children = iv_children, use_iv = use_iv, K = K, iv_mode_used = node_mode, conv = node_conv,
       parent_draws = parent_draws, x_cols = colnames(X), lambda_draws = lambda_draws)
  if (!is.null(cache_file)) { to_cache <- res_nest; to_cache$children <- NULL; saveRDS(to_cache, cache_file) }
  res_nest                                                                   # children rebuilt from their caches on resume
}

# inclusive value of a node evaluated on X_new for imputation m (recursive) --
.ncut_node_iv <- function(nd, X_new, m) {
  if (nd$type %in% c("singleton", "even")) return(rep(0, nrow(X_new)))     # degenerate -> no IV signal
  if (nd$type == "leaf") return(.ncut_lse(X_new %*% nd$beta_draws[[m]]))
  # internal: rebuild its design with its children's IV (imputation m), pick a parent draw
  Xn <- X_new
  if (nd$use_iv && length(nd$iv_children)) {
    ivm <- vapply(nd$iv_children, function(cn) .ncut_node_iv(nd$children[[cn]], X_new, m), numeric(nrow(X_new)))
    colnames(ivm) <- paste0("IV_", nd$iv_children); Xn <- cbind(X_new, ivm)
  }
  B <- nd$parent_draws[[m]]; b <- B[, , sample.int(dim(B)[3], 1)]
  .ncut_lse(Xn %*% b)
}

# =============================================================================
# nested_cut_fit(X, Y, tree, ...)  — X MUST include an "intercept" column
# =============================================================================
nested_cut_fit <- function(X, Y, tree, use_iv = TRUE, group_idx = NULL, use_re = FALSE, re_idx = NULL,
                           M = 40, draws_per_impute = 20, niter = 1200, nburn = 400, thin = 1L,
                           min_pixels = 50, nburn_warm = NULL,
                           iv_mode = c("draws", "moments", "auto"), moment_rank = 1L,
                           moment_skew_tol = 1.0, moment_kurt_tol = 2.0, stream_disk = FALSE, progress = TRUE,
                           n_chains = 1L, n_cores = 1L, store_dir = NULL) {
  iv_mode <- match.arg(iv_mode)
  n_cores <- max(1L, min(as.integer(n_cores), n_chains))
  if (n_cores > 1L && n_chains > 1L)
    message(sprintf("[nested_cut] parallel chains: %d core(s) across %d chains%s", n_cores, n_chains,
                    if (.Platform$OS.type == "windows") " (Windows: forking unavailable -> serial)" else ""))
  X <- as.matrix(X); Y <- as.matrix(Y)
  stopifnot(!is.null(colnames(Y)), all(.ncut_fine(tree) %in% colnames(Y)))
  if (is.null(group_idx)) group_idx <- rep(1L, nrow(X))
  if (is.null(nburn_warm)) nburn_warm <- max(50L, round(nburn / 4))   # short warm burn-in for m>1
  # ---- disk store: persist node draws + auto-resume (config-hash guarded) ----
  store <- NULL
  if (!is.null(store_dir)) {
    nodes_dir <- file.path(store_dir, "nodes"); dir.create(nodes_dir, recursive = TRUE, showWarnings = FALSE)
    h <- .ncut_config_hash(X, Y, tree, niter, nburn, thin, M, use_re, iv_mode, n_chains, draws_per_impute, moment_rank)
    mf <- file.path(store_dir, "manifest.rds")
    if (file.exists(mf)) {
      old <- readRDS(mf)
      if (!identical(old$hash, h))
        stop(sprintf("store_dir '%s' was built with a DIFFERENT config (data/niter/M/... changed).\n  Point store_dir at a fresh path, or delete that folder to refit from scratch.", store_dir))
      ncached <- length(list.files(nodes_dir, "\\.rds$"))
      message(sprintf("[nested_cut] resume: store '%s' matches config; %d node(s) cached -> skipping those.", store_dir, ncached))
    } else saveRDS(list(hash = h, tree = tree, created = format(Sys.time(), "%Y-%m-%d_%H%M"),
                        M = M, niter = niter, nburn = nburn, thin = thin, use_re = use_re,
                        iv_mode = iv_mode, n_chains = n_chains), mf)
    message(sprintf("[nested_cut] store: node draws persist to '%s' (kill-safe / resumable).", store_dir))
    store <- list(dir = nodes_dir)
  }
  .ncut_prog_init(.ncut_count_fits(tree, M, iv_mode, moment_rank, n_chains), progress)  # nested-level progress + ~ETA
  root <- .ncut_node(tree, X, Y, rep(TRUE, nrow(X)), use_iv, use_re, group_idx, re_idx,
                     M, niter, nburn, thin, min_pixels, draws_per_impute, nburn_warm,
                     iv_mode, moment_rank, moment_skew_tol, moment_kurt_tol, stream_disk, n_chains, n_cores,
                     path = "root", store = store)
  fit <- structure(list(root = root, tree = tree, fine_classes = colnames(Y), use_iv = use_iv,
                 X_cols = colnames(X), M = M, draws_per_impute = draws_per_impute,
                 iv_mode = iv_mode, moment_rank = moment_rank, stream_disk = stream_disk, n_chains = n_chains,
                 store_dir = store_dir),
            class = "nested_cut")
  if (n_chains > 1L) { cv <- nested_cut_convergence(fit)
    if (!is.null(cv)) message(sprintf("[nested_cut] convergence (%d chains): worst Rhat %.3f | min %%(Rhat<1.01) %.0f%% | median ESS %.0f  (per-node -> fit$convergence)",
                    n_chains, max(cv$max_rhat, na.rm = TRUE), min(cv$pct_rhat_lt_1.01, na.rm = TRUE), median(cv$median_ess, na.rm = TRUE)))
    fit$convergence <- cv }
  fit
}

# report the IV mode actually used at each internal node (esp. auto fallbacks)
nested_cut_modes <- function(fit) {
  out <- list(); walk <- function(nd, path) {
    if (isTRUE(nd$type == "nest")) {
      out[[path]] <<- if (length(nd$iv_children)) nd$iv_mode_used %||% "n/a" else "no-IV"
      for (cn in nd$child_names) walk(nd$children[[cn]], paste0(path, "/", cn))
    }
  }
  walk(fit$root, "root"); unlist(out)
}

# per-node convergence table (only populated when fit with n_chains > 1): leaves + each nest
# at the mean IV. Split-Rhat / bulk-ESS across chains on the pooled coefficients.
nested_cut_convergence <- function(fit) {
  rows <- list(); walk <- function(nd, path) {
    if (!is.null(nd$conv)) rows[[path]] <<- data.frame(node = path,
        max_rhat = nd$conv$max_rhat, pct_rhat_lt_1.01 = 100 * nd$conv$frac_rhat_lt_1.01,
        pct_rhat_lt_1.1 = 100 * nd$conv$frac_rhat_lt_1.1, median_ess = nd$conv$median_ess,
        n_params = nd$conv$n_params, stringsAsFactors = FALSE)
    if (isTRUE(nd$type == "nest")) for (cn in nd$child_names) walk(nd$children[[cn]], paste0(path, "/", cn))
  }
  walk(fit$root, "root"); if (length(rows)) do.call(rbind, rows) else NULL
}

# --- posterior-draw prediction (coherent cut draws) --------------------------
# Returns [n x J x D] array of fine-class probabilities. Draw d threads imputation
# m(d) through the whole tree + picks one parent draw per internal node.
.ncut_predict_draw <- function(nd, X_new, share, m, P) {
  if (nd$type == "singleton") { P[, nd$fine] <- P[, nd$fine] + share; return(P) }
  if (nd$type == "even")      { for (fc in nd$fine) P[, fc] <- P[, fc] + share / length(nd$fine); return(P) }
  if (nd$type == "leaf") {
    sub <- .ncut_softmax(X_new %*% nd$beta_draws[[m]]); colnames(sub) <- nd$fine
    for (fc in nd$fine) P[, fc] <- P[, fc] + share * sub[, fc]; return(P)
  }
  Xn <- X_new
  if (nd$use_iv && length(nd$iv_children)) {
    ivm <- vapply(nd$iv_children, function(cn) .ncut_node_iv(nd$children[[cn]], X_new, m), numeric(nrow(X_new)))
    colnames(ivm) <- paste0("IV_", nd$iv_children); Xn <- cbind(X_new, ivm)
  }
  B <- nd$parent_draws[[m]]; b <- B[, , sample.int(dim(B)[3], 1)]
  cp <- .ncut_softmax(Xn %*% b); colnames(cp) <- nd$child_names
  for (cn in nd$child_names) P <- .ncut_predict_draw(nd$children[[cn]], X_new, share * cp[, cn], m, P)
  P
}
predict_nested_cut <- function(fit, X_new, D = 200, summary = c("draws", "mean")) {
  summary <- match.arg(summary); X_new <- as.matrix(X_new); J <- length(fit$fine_classes)
  ms <- if (D <= fit$M) round(seq(1, fit$M, length.out = D)) else rep_len(seq_len(fit$M), D)
  zero <- function() matrix(0, nrow(X_new), J, dimnames = list(NULL, fit$fine_classes))
  # "mean" ACCUMULATES a running sum -> O(n*J) memory (never holds all D draws); "draws" must
  # materialise the full [n x J x D] array (that IS the requested output — large at big n/D).
  if (summary == "mean") {
    S <- zero(); for (i in seq_along(ms)) S <- S + .ncut_predict_draw(fit$root, X_new, rep(1, nrow(X_new)), ms[i], zero())
    return(S / length(ms))
  }
  arr <- array(0, c(nrow(X_new), J, length(ms)), dimnames = list(NULL, fit$fine_classes, NULL))
  for (i in seq_along(ms)) arr[, , i] <- .ncut_predict_draw(fit$root, X_new, rep(1, nrow(X_new)), ms[i], zero())
  arr
}

# =============================================================================
# effective_symmetric_table(fit, X) — FINE-CLASS effective symmetric coefficients
# =============================================================================
# The nested model's native parameters are one zero-sum matrix PER NEST. For a
# flat-style heatplot / cross-model comparison, collapse them to a single [J x P]
# table of AVERAGE MARGINAL EFFECTS  mean_i d log P(j | x_i) / d x_p , centred over
# fine classes j (=> zero-sum, comparable to the flat symmetric coefficients).
# Chain rule through the per-level softmaxes and the inclusive-value derivatives;
# uncertainty comes from the M inclusive-value imputations (leaf/IV uncertainty).
.ncut_meanB <- function(nd, m) if (nd$type == "leaf") nd$beta_draws[[m]] else apply(nd$parent_draws[[m]], c(1, 2), mean)

.ncut_dIV_dx <- function(nd, X, m) {                    # d IV_nd / d x : [n x P]
  if (nd$type %in% c("singleton", "even")) return(matrix(0, nrow(X), ncol(X)))
  if (nd$type == "leaf") { B <- nd$beta_draws[[m]]; P <- .ncut_softmax(X %*% B); return(P %*% t(B)) }
  nd_d <- .ncut_node_deriv(nd, X, m)
  D <- matrix(0, nrow(X), ncol(X)); for (k in seq_along(nd_d$children)) D <- D + nd_d$P[, k] * nd_d$dv[[k]]
  D
}
.ncut_node_deriv <- function(nd, X, m) {                # node child-probs + d v_child/d x per child
  Xd <- X; dIVc <- list(); Pn <- ncol(X)
  if (nd$use_iv && length(nd$iv_children)) {
    ivm <- vapply(nd$iv_children, function(cn) .ncut_node_iv(nd$children[[cn]], X, m), numeric(nrow(X)))
    colnames(ivm) <- paste0("IV_", nd$iv_children); Xd <- cbind(X, ivm)
    for (cn in nd$iv_children) dIVc[[cn]] <- .ncut_dIV_dx(nd$children[[cn]], X, m)
  }
  B <- .ncut_meanB(nd, m); P <- .ncut_softmax(Xd %*% B)             # [n x K]
  dv <- vector("list", ncol(B))
  for (k in seq_len(ncol(B))) {
    g <- matrix(B[seq_len(Pn), k], nrow(X), Pn, byrow = TRUE)       # gamma part (const over pixels)
    if (length(dIVc)) for (cn in nd$iv_children) g <- g + B[Pn + match(cn, nd$iv_children), k] * dIVc[[cn]]
    dv[[k]] <- g
  }
  list(P = P, dv = dv, children = nd$child_names)
}
.ncut_ame_rec <- function(nd, X, m, base, env) {                    # accumulate d logP(j)/dx per fine class
  if (nd$type %in% c("singleton", "even")) {
    for (fc in nd$fine) env$AME[[fc]] <- (env$AME[[fc]] %||% 0) + colSums(base); return(invisible())
  }
  if (nd$type == "leaf") {
    B <- nd$beta_draws[[m]]; P <- .ncut_softmax(X %*% B); Ebar <- P %*% t(B)
    for (fi in seq_along(nd$fine)) {
      dev <- matrix(B[, fi], nrow(X), ncol(X), byrow = TRUE) - Ebar
      env$AME[[nd$fine[fi]]] <- (env$AME[[nd$fine[fi]]] %||% 0) + colSums(base + dev)
    }; return(invisible())
  }
  nd_d <- .ncut_node_deriv(nd, X, m)
  Ebar <- matrix(0, nrow(X), ncol(X)); for (k in seq_along(nd_d$children)) Ebar <- Ebar + nd_d$P[, k] * nd_d$dv[[k]]
  for (k in seq_along(nd_d$children))
    .ncut_ame_rec(nd$children[[nd_d$children[k]]], X, m, base + (nd_d$dv[[k]] - Ebar), env)
}
effective_symmetric_table <- function(fit, X, center = TRUE) {
  X <- as.matrix(X); Pn <- ncol(X); J <- length(fit$fine_classes)
  arr <- array(NA_real_, c(J, Pn, fit$M), dimnames = list(fit$fine_classes, colnames(X), NULL))
  for (m in seq_len(fit$M)) {
    env <- new.env(); env$AME <- list()
    .ncut_ame_rec(fit$root, X, m, matrix(0, nrow(X), Pn), env)
    A <- t(vapply(fit$fine_classes, function(cl) (env$AME[[cl]] %||% rep(0, Pn)) / nrow(X), numeric(Pn)))
    if (center) A <- sweep(A, 2, colMeans(A), "-")                 # zero-sum over fine classes
    arr[, , m] <- A
  }
  list(median = apply(arr, c(1, 2), median),
       q025   = apply(arr, c(1, 2), quantile, .025),
       q975   = apply(arr, c(1, 2), quantile, .975))              # each [J x P], zero-sum over J
}

# --- summary: per-node coef tables + lambda with 95% CI ----------------------
summary_nested_cut <- function(fit) {
  out <- list(); walk <- function(nd, path) {
    if (nd$type == "leaf") {
      bm <- Reduce(`+`, nd$beta_draws) / length(nd$beta_draws)
      out[[path]] <<- list(level = path, classes = nd$fine, coef = bm, lambda = NULL)
    } else if (nd$type == "nest") {
      bm <- apply(simplify2array(lapply(nd$parent_draws, function(a) apply(a, c(1, 2), mean))), c(1, 2), mean)
      lam <- if (length(nd$lambda_draws)) t(vapply(nd$lambda_draws, function(v)
        c(median = median(v), q025 = unname(quantile(v, .025)), q975 = unname(quantile(v, .975))), numeric(3))) else NULL
      out[[path]] <<- list(level = path, classes = nd$child_names, coef = bm, lambda = lam)
      for (cn in nd$child_names) walk(nd$children[[cn]], paste0(path, "/", cn))
    }
  }
  walk(fit$root, "root"); out
}
