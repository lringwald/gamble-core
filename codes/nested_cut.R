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

# --- per-node focal routing: IV IDENTIFICATION -------------------------------
# `IV_c = logsumexp(X'beta_leaf)` is BUILT FROM the leaf utilities, so ANY covariate that
# feeds a leaf AND also sits at that leaf's parent node is collinear with its own IV BY
# CONSTRUCTION -- an algebraic identity, not a data artifact you can clean away. Focal LU
# composition is the worst offender (neighbourhood persistence is the strongest leaf
# predictor, so IV_c is nearly LINEAR in focal_c: measured R^2(IV_Cropland ~ focal block)
# = 0.70). The parent then cannot separate lambda (the effect THROUGH the IV) from the
# direct focal effect -> lambda flips NEGATIVE and mechanically dumps probability onto the
# residual singletons in class-dominated regions (the Sweden/Finland "urban wash").
#
# Fix = route focal to the level of the choice it governs. LEAVES ALWAYS KEEP THE FULL FINE
# FOCAL BLOCK, so cross-class -> subtype signal is preserved ("this forest subtype appears
# next to cropland, that one next to urban"). At an INTERNAL node the fine focal columns of
# each IV-supplying child are replaced per `focal_rule`:
#   "macro_totals" (default): one aggregate focal_<child>_tot per IV child -> keeps direct
#                             macro-persistence at the node. lambda and that macro-focal
#                             coefficient are then JOINTLY identified: report the COMBINED
#                             macro effect, do not read lambda as a pure dissimilarity.
#   "leaf_only"             : drop them entirely -> the IV is the sole channel from the
#                             neighbourhood to the macro choice, lambda cleanly identified.
#   "flat"                  : legacy no-op (lambda unidentified -- comparison/diagnostics).
# Focal columns of singleton/degenerate children and of CONTEXT classes that head no nest
# (water, wetlands, ...) always pass through unchanged: they have no IV to collide with.
.ncut_focal_xmap <- function(x_cols, node, iv_children, prefix = "focal_", rule = "macro_totals") {
  if (rule == "flat" || !length(iv_children) || is.character(node)) return(NULL)
  agg <- list(); drop <- integer(0)
  for (cn in iv_children) {
    ix <- which(x_cols %in% paste0(prefix, .ncut_fine(node[[cn]])))
    if (!length(ix)) next
    drop <- c(drop, ix)
    if (rule == "macro_totals") agg[[paste0(prefix, cn, "_tot")]] <- ix
  }
  if (!length(drop)) return(NULL)
  keep <- setdiff(seq_along(x_cols), drop)
  list(keep = keep, agg = agg, cols = c(x_cols[keep], names(agg)), rule = rule, n_dropped = length(drop))
}
# apply a node's design map: kept columns (original order) + macro totals appended at the end
# (so the IV columns still come last, as the lambda extraction assumes).
.ncut_apply_xmap <- function(xm, X) {
  if (is.null(xm)) return(X)
  out <- X[, xm$keep, drop = FALSE]
  if (length(xm$agg)) {
    A <- vapply(xm$agg, function(ix) rowSums(X[, ix, drop = FALSE]), numeric(nrow(X)))
    if (is.null(dim(A))) A <- matrix(A, nrow = nrow(X))
    colnames(A) <- names(xm$agg); out <- cbind(out, A)
  }
  out
}
# [P_orig x P_node] map with X_node = X %*% T. Pushes a node's gradient back into the
# ORIGINAL covariate space for the marginal-effects table: a macro total's slope is shared
# by every fine focal column feeding it.
.ncut_xmap_T <- function(xm, P_orig) {
  if (is.null(xm)) return(diag(P_orig))
  Tm <- matrix(0, P_orig, length(xm$keep) + length(xm$agg))
  for (j in seq_along(xm$keep)) Tm[xm$keep[j], j] <- 1
  if (length(xm$agg)) for (a in seq_along(xm$agg)) Tm[xm$agg[[a]], length(xm$keep) + a] <- 1
  Tm
}
# re_idx is POSITIONAL in X -> remap by NAME onto the node design. A macro total inherits RE
# membership if ANY of the fine focal columns feeding it had it.
.ncut_remap_re <- function(re_idx, xm, x_cols) {
  if (is.null(re_idx) || is.null(xm)) return(re_idx)
  idx <- match(intersect(x_cols[re_idx], xm$cols), xm$cols)
  if (length(xm$agg)) { off <- length(xm$keep)
    for (a in seq_along(xm$agg)) if (any(xm$agg[[a]] %in% re_idx)) idx <- c(idx, off + a) }
  idx <- sort(unique(idx)); if (!length(idx)) NULL else idx
}

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

# in-RAM counterpart of the disk path's per-group draws: postb_total is [k, p_all, G, ndraw], already
# zero-sum expanded per group by the sampler. Returns M slices [k x p_all x G] on the SAME index grid
# as .ncut_beta_draws, so pooled draw m and per-group draw m are the same posterior sample.
.ncut_re_draws <- function(fit, M) {
  B <- fit$postb_total
  if (is.null(B) || length(dim(B)) != 4L) return(NULL)
  d <- dim(B); nd <- d[4]
  idx <- if (nd <= M) rep_len(seq_len(nd), M) else round(seq(1, nd, length.out = M))
  lapply(idx, function(s) array(B[, , , s], c(d[1], d[2], d[3])))
}
.ncut_re_levels <- function(fit) {
  dn <- dimnames(fit$postb_total)
  if (!is.null(dn) && length(dn) >= 3L && !is.null(dn[[3]])) as.character(dn[[3]]) else NULL
}

.ncut_fit_block <- function(X, Y, group_idx, use_re, re_idx, niter, nburn, thin, init_state = NULL,
                            disk_path = NULL, chain_id = NULL, prog_label = NULL) {
  mnlogit_rcpp_sym(
    X = X, Y = Y, intercept = FALSE, symmetric = TRUE, baseline = which.max(colSums(Y)),
    niter = niter, nburn = nburn, thin = thin,
    use_re = use_re, group_idx = if (use_re) group_idx else NULL,
    re_idx = if (use_re && !is.null(re_idx)) re_idx else 1:ncol(X),
    save_posterior_to_disk = !is.null(disk_path), disk_path = disk_path %||% tempdir(), calc_loo = FALSE,
    # SUPPORT PRIOR: fe_support = (n/PR)^strength multiplies the HORSESHOE precision and nothing
    # else (mnlogit_rcpp_sym.R, c_v and hs_prec_mat only). On the diagonal path that precision never
    # reaches the C++ draw, so strength=2 has always been inert there. Under symmetric_hs it IS live,
    # and on the real root design it inflates c_v by a median 106x / max 2.8e5x: measured -125 nats
    # held-out with every driver crushed to 0.0000. So it MUST be 0 whenever the symmetric kernel is
    # on. See memory: fe-horseshoe-never-reaches-draw / symmetric-hs-gate (2026-08-21).
    support_prior_strength = if (isTRUE(as.logical(Sys.getenv("NCUT_SYM_HS", "FALSE")))) 0 else 2,
    # RE VARIANCE: update_re_precision_hc_sym builds its sum of squares from CATEGORY-CENTRED
    # deviations while gibbs_step_re_ncp draws the REs UNCENTRED. Under symmetric_hs that zeroes the
    # random effects outright (RE sd 0.0000, -51.6 nats). center_ss=FALSE estimates the variance in
    # the draw's own coordinates and keeps the across-category pooling. Inert when symmetric_hs=FALSE
    # (the plain updater is used then), so the diagonal path stays bit-identical.
    re_prec_center = FALSE,
    # PER-FAMILY / PER-BLOCK GLOBAL SCALES. NCUT_HS_GROUPS="blocks" puts every non-block covariate in
    # ONE family and gives each const-sum block its OWN tau -- the variant that bought 6.7 nats on the
    # GLOBIOM root, where the fitted per-block scales spanned 450x. Splitting the NON-block covariates
    # into small families (3-6 cols) measured WORSE there: the per-family shape is too small to
    # estimate its own tau. Anything else is parsed as "name=regex;name=regex".
    hs_groups = local({
      g <- Sys.getenv("NCUT_HS_GROUPS", "")
      if (!nzchar(g)) NULL
      else if (identical(tolower(g), "blocks")) list(all = ".")
      else {
        kv <- strsplit(strsplit(g, ";", fixed = TRUE)[[1]], "=", fixed = TRUE)
        setNames(lapply(kv, function(z) paste(z[-1], collapse = "=")), vapply(kv, `[`, "", 1))
      }
    }),
    init_jitter = 0.1,
    const_sum_blocks = "auto",                           # native const-sum detection+drop+zero-sum reconstruction (full-rank fit, full-column posterior)
    # FULL BAYES ON THE RE SLAB (2026-08-13). The cap IS part of the estimation, so estimate it rather
    # than fixing it: no data split is consumed, uncertainty propagates, and it adapts per node.
    # Measured on Forests: estimating BEATS the best fixed value (held-out -725.5 vs -728.3 at c2=4),
    # and c2 is identified -- starts of 4 and 100 both converge to 4.12. collapse_slab_c2 is now just
    # the STARTING value; 4 is near the posterior so burn-in is short.
    re_regularize = TRUE, estimate_slab_c2 = TRUE, collapse_slab_c2 = 4, slab_df_re = 10,
    # FE horseshoe ON, slab FIXED (2026-08-13). Replicated on two nodes with the RE slab estimated:
    # held-out +14.2 (Forests) and +15.5 (Pasture) vs use_horseshoe = FALSE. NOTE an earlier reading
    # put this at only +2 -- that was measured in the OLD RE configuration; the two interact.
    # estimate_c2 stays FALSE: the FE slab is only weakly identified (priors 25x apart -> posteriors
    # 2.3x apart, vs 1.0x for the RE slab) and estimating it does not improve held-out.
    # FE slab ESTIMATED with a TIGHT prior, at the user's explicit request. MEASURED COST (Pasture,
    # n=4000, 600 iter): held-out -796.4 vs -795.1 FIXED (nu=50 is worse still at -798.9). Tightening
    # also made identification WORSE in relative terms -- prior centres 25x apart give posteriors 6.1x
    # apart at nu=20 vs 2.3x at nu=4 -- because a tighter prior holds the posterior nearer its own
    # centre. And the likelihood pulls c2 to 1.93, BELOW the prior mode 3.64, yet that scores worse
    # out-of-sample: the slab fitting in-sample structure that does not generalise. So this is a
    # principled-consistency choice (the slab is part of the estimation), not a fit improvement.
    use_horseshoe = TRUE, estimate_c2 = TRUE, slab_df = 20, slab_s2 = 4,
    # LAMBDA MUST NOT BE SHRUNK (2026-08-14). `horseshoe_idx` defaults to NULL, and the sampler then
    # reads that as 1:k -- EVERY column, including the `IV_*` inclusive values carried at internal
    # nodes. But lambda is derived straight from that coefficient (lambda = zero-sum IV coef *
    # K/(K-1)), so a sparsity prior pulls it toward 0, i.e. toward the RUM boundary, and a small
    # lambda then reads as strong within-nest correlation when it is really prior shrinkage. lambda is
    # a STRUCTURAL parameter, not a candidate for sparsity. This only became live when the null-index
    # fix made the horseshoe actually reach the draw -- before that hs_idx was empty and nothing was
    # shrunk, which is why it was never visible in earlier measurements.
    # SCOPE: `IV_*` ONLY. Leaf nodes carry no IV columns, so their fits stay BIT-IDENTICAL and the
    # measured FE-horseshoe gains (+14.2 Forests / +15.5 Pasture) still stand. The intercept is also
    # in hs_idx by the same NULL default and is arguably as wrong, but excluding it would move every
    # leaf fit, so it is opt-in: NCUT_HS_NO_INTERCEPT=TRUE. NCUT_HS_ALL=TRUE restores the historical
    # shrink-everything behaviour for a measurement arm.
    horseshoe_idx = local({
      cn <- colnames(X)
      if (is.null(cn) || isTRUE(as.logical(Sys.getenv("NCUT_HS_ALL", "FALSE")))) return(NULL)
      excl <- grep("^IV_", cn)
      if (isTRUE(as.logical(Sys.getenv("NCUT_HS_NO_INTERCEPT", "FALSE"))))
        excl <- union(excl, which(cn == "intercept"))
      keep <- setdiff(seq_along(cn), excl)
      # No exclusions -> return NULL, which is what the sampler already assumed: keeps leaf fits
      # bit-identical rather than merely equivalent.
      if (!length(excl) || !length(keep)) NULL else keep
    }),
    # SYMMETRIC (zero-sum/CLR) HORSESHOE, opt-in via NCUT_SYM_HS=TRUE. Principled argument FOR it: with
    # it off, the FE horseshoe is a DIAGONAL penalty in baseline-removed coords, so shrinkage depends on
    # which class is baseline -- and `baseline` below is which.max(colSums(Y)), i.e. DATA-DEPENDENT and
    # DIFFERENT PER NODE. symmetric_hs applies c_v * (I - 11'/p_all), shrinking the rotation-invariant
    # magnitude instead: baseline-invariant, and consistent with the zero-sum coding used everywhere
    # else. MEASURED COST (2026-08-17, Forests, n=6000, 9 RE covs, 4 chains x 2000 iter): held-out
    # -1190.4 vs -1061.2 diagonal, i.e. 129 nats WORSE, with ESS up ~800 -- the over-shrinkage
    # signature (one scale per covariate over its whole zero-sum vector, so a driver cannot be kept for
    # one class and shrunk for another). CAVEAT: that arm ran at Rhat 1.91, so it was not converged.
    # FIXED 2026-08-19 (two independent defects in mnlogit_rcpp_sym.R): (1) block_sym applied Mb
    # WITHIN each equation (implicitly Mb %x% I_p), a BASELINE-DEPENDENT penalty on 50 of 68
    # columns -- now kron(Msym, Mb), invariant to 7 s.f.; (2) the complement redraw conditioned on
    # the post-ASIS mu_R while Pb_lik encoded the pre-draw one. Symmetric+RE on the real root node:
    # McFadden -0.018 -> +0.216, worst nest ratio 0.23 -> 0.78, suite 35/35 with C7 at 0.94x.
    # STILL DEFAULTS FALSE: the only head-to-head vs diagonal (+0.216 vs +0.260) is IN-SAMPLE at
    # one node. The old 129/133-nat 'symmetric loses' gate is VOID -- those arms had REs, so they
    # measured the bug. A proper held-out comparison has not been run.
    symmetric_hs = isTRUE(as.logical(Sys.getenv("NCUT_SYM_HS", "FALSE"))),
    init_state = init_state,                              # HOT-START across imputations (see below)
    chain_id = chain_id,                                 # non-NULL -> no stray txtProgressBar; labels output
    progress_cb = .ncut_fit_cb(prog_label)               # throttled per-fit heartbeat (or silent)
  )
}

# read `ndraws` draws from streamed batches: the zero-sum POOLED coefficients [k x p_all] AND, when
# the fit had random effects, the matching PER-GROUP coefficients [k x p_all x G].
#
# WHY THE PER-GROUP DRAWS ARE NOT OPTIONAL. The likelihood only ever sees beta_g = mu + b_g; mu on its
# own is NOT the population-averaged effect, and because softmax is nonlinear no single coefficient
# vector reproduces a group-heterogeneous MNL (mean_g(beta_g) fails too). Keeping mu alone made the
# saved fit score McFadden -0.37 in-sample -- worse than the grand mean -- while the per-group betas
# scored +0.17 on the same data. The sampler already writes `beta` into every batch; this reads it.
#
# TWO PASSES, so memory stays bounded (the point of stream_disk): pass 1 reads only `mu` (small) to
# count the draws and choose the retained indices; pass 2 re-reads and keeps `beta` ONLY at those
# indices, dropping each batch as soon as it is scanned. Peak = one batch, independent of niter.
.ncut_draws_from_disk <- function(disk_path, ndraws) {
  meta <- qs2::qs_read(file.path(disk_path, "model_metadata.qs"))
  p_all <- length(meta$cat_names); bl <- which(meta$cat_names == (meta$baseline_name %||% ""))[1]
  if (is.na(bl)) bl <- p_all; pp <- (seq_len(p_all))[-bl]
  files <- list.files(disk_path, "posterior_batch_.*\\.qs$", full.names = TRUE)  # one fit per tempdir -> one chain
  files <- files[order(as.integer(gsub(".*batch_([0-9]+)_.*", "\\1", basename(files))))]
  # same expansion the in-RAM path applies (mnlogit_rcpp_sym: expand_to_zs / beta_c_zs)
  to_zs <- function(m) { full <- matrix(0, nrow(m), p_all); full[, pp] <- m; sweep(full, 1, rowMeans(full), "-") }
  mus <- list(); for (f in files) { b <- qs2::qs_read(f); mus <- c(mus, lapply(b, `[[`, "mu")); rm(b) }
  idx <- if (length(mus) <= ndraws) rep_len(seq_along(mus), ndraws) else round(seq(1, length(mus), length.out = ndraws))
  pooled <- lapply(idx, function(i) to_zs(mus[[i]]))
  if (!isTRUE(meta$use_re)) return(list(pooled = pooled, re = NULL, group_levels = NULL))
  keep <- sort(unique(idx)); re_by_i <- vector("list", max(keep)); off <- 0L
  for (f in files) {
    b <- qs2::qs_read(f)
    for (i in keep[keep > off & keep <= off + length(b)]) {
      B <- b[[i - off]]$beta                                   # [k x p x G], baseline-removed
      # `B[, , g]` DROPS the p dimension when p == 1, i.e. for any TWO-class node (baseline removal
      # leaves one column) -- to_zs then gets a vector, nrow() is NULL and matrix() dies with
      # "non-numeric matrix extent". Curated trees make this reachable: Waterbodies has 2 leaves.
      # Reconstruct the [k x p] matrix explicitly, as the in-RAM path already does.
      re_by_i[[i]] <- array(vapply(seq_len(dim(B)[3]),
                                   function(g) to_zs(matrix(B[, , g], dim(B)[1], dim(B)[2])),
                                   matrix(0, dim(B)[1], p_all)), c(dim(B)[1], p_all, dim(B)[3]))
    }
    off <- off + length(b); rm(b)
  }
  # group_levels are the sampler's `unique(group_idx)` APPEARANCE ORDER -- the 3rd dim is a POSITION,
  # so downstream must key via match(g, group_levels), never by raw id or sorted level.
  list(pooled = pooled, re = lapply(idx, function(i) re_by_i[[i]]), group_levels = meta$group_levels)
}

# fit a node's sub-MNL and return `ndraws` zero-sum pooled draws [k x p_all] + final_state (hot-start).
# stream_disk=TRUE keeps the sampler's big per-group array off RAM (safer at scale / high niter).
.ncut_fit_draws <- function(X, Y, group_idx, use_re, re_idx, niter, nburn, thin, ndraws,
                            init_state = NULL, stream_disk = FALSE, chain_id = NULL, prog_label = NULL,
                            persist_dir = NULL) {
  if (isTRUE(stream_disk)) {
    # persist_dir keeps the sampler's FULL streamed posterior (every retained draw of beta/mu/sigma_re/
    # horseshoe/log_lik) instead of deleting it. Default stays a tempdir because the batches are large,
    # but keeping them is what makes a fit RECOVERABLE: nested_cut itself only ever retains `ndraws`
    # thinned draws, so any later change to what is extracted (as in the 2026-08-12 per-group fix)
    # otherwise costs a full re-fit. One dir PER CHAIN -- the readers glob posterior_batch_*.qs and
    # assume a single chain per directory.
    td <- persist_dir %||% tempfile("ncut_")
    dir.create(td, recursive = TRUE, showWarnings = FALSE)
    if (is.null(persist_dir)) on.exit(unlink(td, recursive = TRUE), add = TRUE)  # clean even on crash/interrupt
    fit <- .ncut_fit_block(X, Y, group_idx, use_re, re_idx, niter, nburn, thin, init_state,
                           disk_path = td, chain_id = chain_id, prog_label = prog_label)
    d <- .ncut_draws_from_disk(td, ndraws)
    list(draws = d$pooled, re_draws = d$re, group_levels = d$group_levels, final_state = fit$final_state)
  } else {
    fit <- .ncut_fit_block(X, Y, group_idx, use_re, re_idx, niter, nburn, thin, init_state,
                           chain_id = chain_id, prog_label = prog_label)
    list(draws = .ncut_beta_draws(fit, ndraws),
         re_draws = if (isTRUE(use_re)) .ncut_re_draws(fit, ndraws) else NULL,
         group_levels = if (isTRUE(use_re)) .ncut_re_levels(fit) else NULL,
         final_state = fit$final_state)
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

# Same, for the PER-GROUP coefficients [k x p_all x G]. This is the one that matters: since the
# per-group betas are what prediction uses, pooled-mu convergence can look fine while the reportable
# quantity is stuck (exactly what happened on the count model -- mu at ESS 3000, reportable at 15).
# Cells are subsampled to `max_cells` because k*p*G runs to ~10^4 and every cell costs an rhat call.
.ncut_rhat_re <- function(chain_re, max_cells = 3000L, seed = 1L) {
  if (is.null(chain_re) || any(vapply(chain_re, is.null, TRUE))) return(NULL)
  nc <- length(chain_re); nd <- min(vapply(chain_re, length, 0L)); if (nc < 2L || nd < 2L) return(NULL)
  dm <- dim(chain_re[[1]][[1]]); if (length(dm) != 3L) return(NULL)
  k <- dm[1]; p <- dm[2]; G <- dm[3]
  arr <- array(0, c(k, p, G, nd, nc))
  for (c in seq_len(nc)) for (s in seq_len(nd)) arr[, , , s, c] <- chain_re[[c]][[s]]
  cells <- expand.grid(a = seq_len(k), b = seq_len(p), g = seq_len(G))
  if (nrow(cells) > max_cells) { set.seed(seed); cells <- cells[sort(sample(nrow(cells), max_cells)), ] }
  rh <- ess <- numeric(nrow(cells))
  for (i in seq_len(nrow(cells))) {
    m <- matrix(arr[cells$a[i], cells$b[i], cells$g[i], , ], nd, nc)
    rh[i]  <- tryCatch(posterior::rhat(m),     error = function(e) NA_real_)
    ess[i] <- tryCatch(posterior::ess_bulk(m), error = function(e) NA_real_)
  }
  list(max_rhat_re = max(rh, na.rm = TRUE), frac_rhat_re_lt_1.01 = mean(rh < 1.01, na.rm = TRUE),
       frac_rhat_re_lt_1.1 = mean(rh < 1.1, na.rm = TRUE), median_ess_re = median(ess, na.rm = TRUE),
       n_cells_re = nrow(cells), n_cells_total = k * p * G)
}

# run n_chains sub-fits (cold, different seeds) for a CONVERGENCE CHECK: pool their draws for the model
# and compute Rhat/ESS across chains. Used on leaves and on each nest evaluated at the MEAN inclusive value.
.ncut_fit_chains <- function(X, Y, group_idx, use_re, re_idx, niter, nburn, thin, ndraws, n_chains,
                             stream_disk = FALSE, conv_draws = 150L, n_cores = 1L, prog_label = NULL,
                             persist_dir = NULL) {
  one <- function(c) .ncut_fit_draws(X, Y, group_idx, use_re, re_idx, niter, nburn, thin, conv_draws,
                                     stream_disk = stream_disk, chain_id = c,   # cold, independent chain
                                     prog_label = sprintf("%s ch%d/%d", prog_label %||% "leaf", c, n_chains),
                                     persist_dir = if (is.null(persist_dir)) NULL else file.path(persist_dir, sprintf("chain%d", c)))
  nco <- max(1L, min(as.integer(n_cores), n_chains))
  # fork the independent chains when n_cores>1 (mclapply: COW-shared X/Y, no Windows fork -> serial fallback).
  # mclapply is a synchronous fork-join: it reaps ALL its forks before returning, so no worker survives into
  # the next stage. mc.cleanup=TRUE also SIGTERMs any stray child if the master exits (belt-and-suspenders for
  # an interrupted run); gc() then releases the reaped forks' copy-on-write memory promptly between stages.
  reslist <- if (nco > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(seq_len(n_chains), one, mc.cores = nco, mc.preschedule = FALSE, mc.cleanup = TRUE)
  else lapply(seq_len(n_chains), one)
  if (nco > 1L) gc(FALSE)
  bad <- vapply(reslist, function(r) inherits(r, "try-error") || is.null(r$draws), logical(1))
  if (any(bad)) stop(".ncut_fit_chains: ", sum(bad), "/", n_chains, " chain(s) failed (see mclapply warnings)")
  chains <- lapply(reslist, `[[`, "draws"); fs <- reslist[[1]]$final_state
  allp <- unlist(chains, recursive = FALSE)                       # pool all chains' draws
  sel  <- round(seq(1, length(allp), length.out = ndraws))        # subsample ndraws for the model
  # the per-group draws must ride the SAME index grid, so pooled draw m and per-group draw m stay the
  # same posterior sample (they are concatenated in identical per-chain order).
  allre <- unlist(lapply(reslist, `[[`, "re_draws"), recursive = FALSE)
  re_chains <- lapply(reslist, `[[`, "re_draws")
  cv <- .ncut_rhat(chains)
  cv_re <- .ncut_rhat_re(re_chains)          # convergence of the REPORTABLE quantity, not just mu
  list(draws = allp[sel],
       re_draws = if (length(allre) == length(allp)) allre[sel] else NULL,
       group_levels = reslist[[1]]$group_levels,
       final_state = fs, conv = c(cv, cv_re))
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
.ncut_prog_init <- function(total, on, fit_on = FALSE, fit_sec = 30) {
  .ncut_prog$on <- isTRUE(on); .ncut_prog$total <- total
  .ncut_prog$fit_on <- isTRUE(on) && isTRUE(fit_on); .ncut_prog$fit_sec <- fit_sec
  .ncut_prog$done <- 0L; .ncut_prog$t0 <- Sys.time(); .ncut_prog$tlast <- .ncut_prog$t0 }
# per-fit heartbeat: throttled relay of the sampler's "Iteration i/niter [phase]" for one chain/leaf.
.ncut_make_cb <- function(label, every_sec = 30) {
  e <- new.env(parent = emptyenv()); e$last <- as.numeric(Sys.time()) - 1e6
  function(message = NULL, ...) {                          # `message` shadows base::message -> use base::
    if (is.null(message)) return(invisible())
    now <- as.numeric(Sys.time())
    if (now - e$last >= every_sec) { e$last <- now
      base::message(sprintf("[nested_cut]     %-26s %s", label, sub("^Chain [^:]*: ", "", message))) }
    invisible()
  }
}
# progress_cb for a sub-fit: heartbeat if fit-progress is on and we have a label, else silent.
.ncut_fit_cb <- function(prog_label) {
  if (isTRUE(.ncut_prog$fit_on) && !is.null(prog_label)) .ncut_make_cb(prog_label, .ncut_prog$fit_sec %||% 30)
  else function(...) invisible(NULL)
}
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
# where a node's FULL streamed posterior lives when keep_posterior is on (NULL = use a tempdir and
# delete it, the default). One subdir per node; .ncut_fit_chains adds chainN under it.
.ncut_persist_dir <- function(store, path)
  if (is.null(store) || is.null(store$posterior_dir)) NULL else file.path(store$posterior_dir, .ncut_path_key(path))
.ncut_config_hash <- function(X, Y, tree, niter, nburn, thin, M, use_re, iv_mode, n_chains,
                              draws_per_impute, moment_rank, focal_rule = "macro_totals") {
  Xf <- X[is.finite(X)]
  # store_format 2 = nodes carry per-group draws (beta_re_draws / parent_re_draws). Bumping this
  # invalidates format-1 stores ON PURPOSE: those nodes hold pooled mu only, so resuming into a
  # format-2 run would silently mix predictable and unpredictable nodes in one tree.
  key <- list(store_format = 2L,
              cols = colnames(X), n = nrow(X), p = ncol(X), ycols = colnames(Y),
              sx = sum(Xf), sx2 = sum(Xf^2), sy = sum(Y), niter = niter, nburn = nburn, thin = thin,
              M = M, dpi = draws_per_impute, use_re = isTRUE(use_re), iv_mode = iv_mode,
              n_chains = n_chains, mrank = moment_rank, focal_rule = focal_rule,
              tree = paste(sort(.ncut_fine(tree)), collapse = "|"),
              # PRIOR SIGNATURE (2026-08-21). Without this the hash covers only the DESIGN and the
              # sampling geometry, so two runs that differ ONLY in the prior collide: the second
              # silently reloads the first's cached draws and reports "fit done in 0.0 min". That
              # makes every prior A/B run through run_nested_cut.R invalid unless the store was
              # cleared by hand -- exactly how a symmetric-vs-diagonal arm can come back identical.
              # These are the env knobs .ncut_fit_block reads at fit time.
              prior = paste(Sys.getenv("NCUT_SYM_HS", ""), Sys.getenv("NCUT_HS_GROUPS", ""),
                            Sys.getenv("NCUT_HS_ALL", ""), Sys.getenv("NCUT_HS_NO_INTERCEPT", ""),
                            sep = "|"))
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
                       stream_disk = FALSE, n_chains = 1L, n_cores = 1L,
                       focal_rule = "macro_totals", focal_prefix = "focal_", path = "root", store = NULL) {
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
                                 stream_disk, n_chains, n_cores, focal_rule = focal_rule,
                                 focal_prefix = focal_prefix, path = paste0(path, "/", cn), store = store)
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
    pdir <- .ncut_persist_dir(store, path)
    res <- if (n_chains > 1L)
      .ncut_fit_chains(X[keep, , drop = FALSE], Ysub, group_idx[keep], use_re, re_idx, niter, nburn, thin, M, n_chains, stream_disk, n_cores = n_cores, prog_label = path, persist_dir = pdir)
    else .ncut_fit_draws(X[keep, , drop = FALSE], Ysub, group_idx[keep], use_re, re_idx, niter, nburn, thin, M, stream_disk = stream_disk, chain_id = 1L, prog_label = path,
                         persist_dir = if (is.null(pdir)) NULL else file.path(pdir, "chain1"))
    .ncut_prog_tick(sprintf("leaf {%s}%s", paste(fine, collapse = ","), if (n_chains > 1L) sprintf(" x%d", n_chains) else ""))
    res_leaf <- list(type = "leaf", fine = fine, classes = fine, beta_draws = res$draws,
                     beta_re_draws = res$re_draws, group_levels = res$group_levels, conv = res$conv)
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
                                 focal_rule = focal_rule, focal_prefix = focal_prefix,
                                 path = paste0(path, "/", cn), store = store)
  }
  # aggregated child membership Y at this node
  Ynode <- vapply(child_names, function(cn) rowSums(Y[, .ncut_fine(node[[cn]]), drop = FALSE]), numeric(nrow(Y)))
  colnames(Ynode) <- child_names
  keep <- in_pixels & (rowSums(Ynode) > 0)
  # which children can supply an IV (multi-class leaf or internal, i.e. non-degenerate)
  iv_children <- child_names[vapply(children, function(c) c$type %in% c("leaf", "nest"), logical(1))]
  K <- length(child_names)

  # ---- this node's design: route focal away from its own IV columns (see .ncut_focal_xmap).
  # Only when IVs are actually used -- with use_iv=FALSE (factorized) there is nothing to collide with.
  xmap  <- if (use_iv) .ncut_focal_xmap(colnames(X), node, iv_children, focal_prefix, focal_rule) else NULL
  Xnd   <- .ncut_apply_xmap(xmap, X)
  re_nd <- .ncut_remap_re(re_idx, xmap, colnames(X))
  if (!is.null(xmap))
    message(sprintf("[nested_cut] node {%s}: focal rule '%s' -> %d fine focal col(s) %s (design %d -> %d cols)",
      paste(child_names, collapse = ","), focal_rule, xmap$n_dropped,
      if (length(xmap$agg)) paste0("aggregated to ", paste(names(xmap$agg), collapse = ", ")) else "dropped",
      ncol(X), ncol(Xnd)))

  # per-imputation IV of each iv-child on the train grid (M coherent fields)
  iv_by_m <- if (use_iv && length(iv_children))
    lapply(seq_len(M), function(m) {
      ivm <- vapply(iv_children, function(cn) .ncut_node_iv(children[[cn]], X, m), numeric(nrow(X)))
      if (is.null(dim(ivm))) ivm <- matrix(ivm, ncol = length(iv_children))
      colnames(ivm) <- paste0("IV_", iv_children); ivm }) else NULL

  # ---- IV IDENTIFICATION DIAGNOSTIC: how much of each IV is ALREADY spanned by this node's
  # own design? R^2 -> 1 means lambda_c is not separately identified from the direct effects
  # (that is the pathology the focal routing above fixes) -- reported so it is measured, not
  # assumed. Cheap: one QR on the node design, reused across children.
  iv_r2 <- NULL
  if (!is.null(iv_by_m)) {
    ivbar <- Reduce(`+`, iv_by_m) / length(iv_by_m)
    qrX <- qr(Xnd[keep, , drop = FALSE])
    iv_r2 <- vapply(seq_along(iv_children), function(cc) {
      y <- ivbar[keep, cc]; v <- var(y)
      if (!is.finite(v) || v <= 0) NA_real_ else max(0, 1 - var(qr.resid(qrX, y)) / v)
    }, numeric(1))
    names(iv_r2) <- iv_children
    message(sprintf("[nested_cut] node {%s}: IV~design R2 = %s%s", paste(child_names, collapse = ","),
      paste(sprintf("%s %.2f", iv_children, iv_r2), collapse = ", "),
      if (any(iv_r2 > 0.9, na.rm = TRUE)) "   <-- WEAK lambda identification" else ""))
  }

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
  fit_draws <- vector("list", length(fit_fields)); fit_re <- vector("list", length(fit_fields))
  node_levels <- NULL; warm <- NULL
  for (j in seq_along(fit_fields)) {
    Xnode <- if (is.null(fit_fields[[j]])) Xnd else cbind(Xnd, fit_fields[[j]])
    Yk <- Ynode[keep, , drop = FALSE]; Yk <- Yk / rowSums(Yk)
    nb_j <- if (j == 1L || is.null(warm)) nburn else nburn_warm
    res_j <- .ncut_fit_draws(Xnode[keep, , drop = FALSE], Yk, group_idx[keep], use_re, re_nd,
                             nb_j + (niter - nburn), nb_j, thin, draws_per_impute, init_state = warm, stream_disk = stream_disk,
                             chain_id = j, prog_label = sprintf("%s imp%d/%d", path, j, length(fit_fields)),
                             persist_dir = { pd <- .ncut_persist_dir(store, path)
                                             if (is.null(pd)) NULL else file.path(pd, sprintf("imp%d", j)) })
    warm <- res_j$final_state
    Bd <- res_j$draws; kk <- nrow(Bd[[1]]); Kfit <- ncol(Bd[[1]])
    pa <- array(unlist(Bd), c(kk, Kfit, length(Bd)))
    if (Kfit == length(child_names)) dimnames(pa) <- list(colnames(Xnode), child_names, NULL)
    fit_draws[[j]] <- pa
    # per-group counterpart: [k x Kfit x G x ndraw], same draw order as `pa`
    if (!is.null(res_j$re_draws)) {
      Rd <- res_j$re_draws; G <- dim(Rd[[1]])[3]
      fit_re[[j]] <- array(unlist(Rd), c(kk, Kfit, G, length(Rd)))
      node_levels <- res_j$group_levels
    }
    .ncut_prog_tick(sprintf("nest {%s} %d/%d", paste(child_names, collapse = ","), j, length(fit_fields)))
  }

  # ---- resolve to an M-length parent-draw list (resample sigma fits by weight) ----
  # the SAME resampled index must be applied to the per-group draws, or draw m's pooled and per-group
  # coefficients would come from different fields.
  sel_m <- if (length(fit_draws) == M && node_mode != "moments") seq_along(fit_draws)
    else sample(seq_along(fit_draws), M, replace = TRUE, prob = fit_w)
  parent_draws <- fit_draws[sel_m]
  parent_re_draws <- if (any(!vapply(fit_re, is.null, TRUE))) fit_re[sel_m] else NULL

  # ---- lambda per iv-child (IV rows follow the X columns, in iv_children order) ----
  # A FACTORIZED node has no lambda at all, so leave the list EMPTY rather than filling it with NULL
  # placeholders: `length(lambda_draws)` is what summary_nested_cut and nested_cut_identification test
  # to decide whether to report lambda, and a list of NULLs passes that test and then fails inside
  # vapply ("values must be length 3").
  lambda_draws <- if (use_iv) setNames(vector("list", length(iv_children)), iv_children) else list()
  if (use_iv && length(iv_children)) { Pn <- ncol(Xnd)
    for (cn in iv_children) { r <- Pn + match(cn, iv_children); cc <- which(child_names == cn)
      lambda_draws[[cn]] <- unlist(lapply(parent_draws, function(A) A[r, cc, ])) * K / (K - 1) }
  }
  # ---- convergence check (n_chains>1): fit this nest at the MEAN inclusive value across chains -> Rhat/ESS.
  # (The M imputation fits above stay single-chain; this separate multi-chain fit checks the sampler mixes.)
  node_conv <- NULL
  if (n_chains > 1L) {
    Xc <- if (!is.null(iv_by_m)) { ivm <- Reduce(`+`, iv_by_m) / length(iv_by_m)
                                   colnames(ivm) <- paste0("IV_", iv_children); cbind(Xnd, ivm) } else Xnd
    Yc <- Ynode[keep, , drop = FALSE]; Yc <- Yc / rowSums(Yc)
    rc <- .ncut_fit_chains(Xc[keep, , drop = FALSE], Yc, group_idx[keep], use_re, re_nd,
                           niter, nburn, thin, draws_per_impute, n_chains, stream_disk, n_cores = n_cores,
                           prog_label = sprintf("%s conv", path),
                           persist_dir = { pd <- .ncut_persist_dir(store, path)
                                           if (is.null(pd)) NULL else file.path(pd, "conv") })
    node_conv <- rc$conv
    .ncut_prog_tick(sprintf("nest {%s} conv x%d", paste(child_names, collapse = ","), n_chains))
  }
  res_nest <- list(type = "nest", node = node, child_names = child_names, children = children,
       iv_children = iv_children, use_iv = use_iv, K = K, iv_mode_used = node_mode, conv = node_conv,
       parent_draws = parent_draws, parent_re_draws = parent_re_draws, group_levels = node_levels,
       x_cols = colnames(Xnd), lambda_draws = lambda_draws,
       xmap = xmap, focal_rule = focal_rule,   # xmap replays the design at predict/AME time
       iv_r2 = iv_r2)
  if (!is.null(cache_file)) { to_cache <- res_nest; to_cache$children <- NULL; saveRDS(to_cache, cache_file) }
  res_nest                                                                   # children rebuilt from their caches on resume
}

# --- per-group linear predictor ---------------------------------------------
# eta with the row's OWN group coefficients where they exist. Rows whose group is unknown to the fit
# (NA after match) and fits without REs fall back to the pooled coefficients -- the honest default for
# an unseen country, and the only option when b_g was never estimated.
.ncut_at <- function(lst, m) if (is.null(lst)) NULL else lst[[m]]
.ncut_eta <- function(Xn, Bpool, Bre, gpos) {
  if (is.null(Bre) || is.null(gpos)) return(Xn %*% Bpool)
  U <- Xn %*% Bpool                                              # default: pooled (covers NA/unseen)
  G <- dim(Bre)[3]
  for (g in unique(gpos[!is.na(gpos)])) {
    if (g < 1L || g > G) next
    r <- which(gpos == g)
    U[r, ] <- Xn[r, , drop = FALSE] %*% Bre[, , g]
  }
  U
}
# one nest draw: pooled [k x K] and its per-group [k x K x G] must come from the SAME sub-draw s
.ncut_nest_draw <- function(nd, m) {
  B <- nd$parent_draws[[m]]; s <- sample.int(dim(B)[3], 1)
  Re <- nd$parent_re_draws
  list(b = B[, , s],
       bre = if (is.null(Re)) NULL else array(Re[[m]][, , , s], dim(Re[[m]])[1:3]))
}
# first group-level vector found in the tree (all nodes share one group_idx)
.ncut_group_levels <- function(nd) {
  if (!is.null(nd$group_levels)) return(nd$group_levels)
  if (identical(nd$type, "nest")) for (ch in nd$children) {
    lv <- .ncut_group_levels(ch); if (!is.null(lv)) return(lv)
  }
  NULL
}

# inclusive value of a node evaluated on X_new for imputation m (recursive) --
.ncut_node_iv <- function(nd, X_new, m, gpos = NULL) {
  if (nd$type %in% c("singleton", "even")) return(rep(0, nrow(X_new)))     # degenerate -> no IV signal
  if (nd$type == "leaf")
    return(.ncut_lse(.ncut_eta(X_new, nd$beta_draws[[m]], .ncut_at(nd$beta_re_draws, m), gpos)))
  # internal: rebuild its design (focal routing FIRST, then its children's IV for imputation
  # m -- exactly the column order it was fit with), pick a parent draw
  Xn <- .ncut_apply_xmap(nd$xmap, X_new)
  if (nd$use_iv && length(nd$iv_children)) {
    ivm <- vapply(nd$iv_children, function(cn) .ncut_node_iv(nd$children[[cn]], X_new, m, gpos), numeric(nrow(X_new)))
    colnames(ivm) <- paste0("IV_", nd$iv_children); Xn <- cbind(Xn, ivm)
  }
  d <- .ncut_nest_draw(nd, m)
  .ncut_lse(.ncut_eta(Xn, d$b, d$bre, gpos))
}

# =============================================================================
# nested_cut_fit(X, Y, tree, ...)  — X MUST include an "intercept" column
# =============================================================================
nested_cut_fit <- function(X, Y, tree, use_iv = TRUE, group_idx = NULL, use_re = FALSE, re_idx = NULL,
                           M = 40, draws_per_impute = 20, niter = 1200, nburn = 400, thin = 1L,
                           min_pixels = 50, nburn_warm = NULL,
                           iv_mode = c("draws", "moments", "auto"), moment_rank = 1L,
                           moment_skew_tol = 1.0, moment_kurt_tol = 2.0, stream_disk = FALSE, progress = TRUE,
                           n_chains = 1L, n_cores = 1L, store_dir = NULL, keep_posterior = FALSE,
                           focal_rule = c("macro_totals", "leaf_only", "flat"), focal_prefix = "focal_",
                           fit_progress = TRUE, fit_progress_sec = 30) {
  iv_mode <- match.arg(iv_mode); focal_rule <- match.arg(focal_rule)
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
    h <- .ncut_config_hash(X, Y, tree, niter, nburn, thin, M, use_re, iv_mode, n_chains, draws_per_impute,
                           moment_rank, focal_rule)
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
    store <- list(dir = nodes_dir,
                  # keep_posterior: retain the sampler's FULL streamed posterior per node/chain under
                  # <store_dir>/posterior/<node>/chainN instead of a deleted tempdir. Large (order GB
                  # at production niter) but it is the difference between "re-read the batches" and
                  # "re-fit for two days" the next time what we extract changes.
                  posterior_dir = if (isTRUE(keep_posterior)) file.path(store_dir, "posterior") else NULL)
    if (isTRUE(keep_posterior))
      message(sprintf("[nested_cut] keep_posterior: FULL per-draw posterior retained under '%s' (expect GBs).",
                      file.path(store_dir, "posterior")))
  }
  .ncut_prog_init(.ncut_count_fits(tree, M, iv_mode, moment_rank, n_chains), progress,
                  fit_on = fit_progress, fit_sec = fit_progress_sec)         # + per-chain/leaf heartbeat
  root <- .ncut_node(tree, X, Y, rep(TRUE, nrow(X)), use_iv, use_re, group_idx, re_idx,
                     M, niter, nburn, thin, min_pixels, draws_per_impute, nburn_warm,
                     iv_mode, moment_rank, moment_skew_tol, moment_kurt_tol, stream_disk, n_chains, n_cores,
                     focal_rule = focal_rule, focal_prefix = focal_prefix, path = "root", store = store)
  fit <- structure(list(root = root, tree = tree, fine_classes = colnames(Y), use_iv = use_iv,
                 X_cols = colnames(X), M = M, draws_per_impute = draws_per_impute,
                 iv_mode = iv_mode, moment_rank = moment_rank, stream_disk = stream_disk, n_chains = n_chains,
                 focal_rule = focal_rule, focal_prefix = focal_prefix, store_dir = store_dir),
            class = "nested_cut")
  if (n_chains > 1L) { cv <- nested_cut_convergence(fit)
    if (!is.null(cv)) {
      message(sprintf("[nested_cut] convergence (%d chains) POOLED mu: worst Rhat %.3f | min %%(Rhat<1.01) %.0f%% | median ESS %.0f",
                    n_chains, max(cv$max_rhat, na.rm = TRUE), min(cv$pct_rhat_lt_1.01, na.rm = TRUE), median(cv$median_ess, na.rm = TRUE)))
      # The per-group betas are what prediction uses, so THIS is the line to judge the run on --
      # pooled mu can look converged while the reportable quantity is stuck.
      if ("max_rhat_re" %in% names(cv) && any(is.finite(cv$max_rhat_re)))
        message(sprintf("[nested_cut] convergence (%d chains) PER-GROUP beta (REPORTABLE): worst Rhat %.3f | min %%(Rhat<1.01) %.0f%% | median ESS %.0f  (per-node -> fit$convergence)",
                    n_chains, max(cv$max_rhat_re, na.rm = TRUE), min(cv$pct_rhat_re_lt_1.01, na.rm = TRUE), median(cv$median_ess_re, na.rm = TRUE)))
    }
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

# --- lambda identification report -------------------------------------------
# Per internal node x IV-child: R^2 of the inclusive value on that node's OWN design (how
# much of the IV is already spanned by the direct covariates) next to the fitted lambda and
# its CI. High R^2 + a lambda outside [0,1] = the collinearity pathology; the focal routing
# (fit$focal_rule) is what keeps R^2 down. `n_focal_routed` = fine focal columns replaced.
nested_cut_identification <- function(fit) {
  rows <- list(); walk <- function(nd, path) {
    if (isTRUE(nd$type == "nest")) {
      # `iv_children` is every child regardless of use_iv (:520), so a FACTORIZED node reaches this
      # loop too -- and `list()[[cn]]` is an out-of-bounds error, not NULL. There is no lambda to
      # report without IVs, so gate on use_iv; the walk still recurses and returns NULL overall.
      if (isTRUE(nd$use_iv) && length(nd$iv_children)) for (cn in nd$iv_children) {
        lam <- nd$lambda_draws[[cn]]
        rows[[paste(path, cn)]] <<- data.frame(node = path, iv_child = cn,
          iv_r2 = unname((nd$iv_r2 %||% NA_real_)[cn]),
          lambda = if (length(lam)) median(lam) else NA_real_,
          lambda_q025 = if (length(lam)) unname(quantile(lam, .025)) else NA_real_,
          lambda_q975 = if (length(lam)) unname(quantile(lam, .975)) else NA_real_,
          in_unit_interval = if (length(lam)) median(lam) > 0 && median(lam) <= 1 else NA,
          focal_rule = nd$focal_rule %||% "flat",
          n_focal_routed = if (is.null(nd$xmap)) 0L else as.integer(nd$xmap$n_dropped),
          stringsAsFactors = FALSE)
      }
      for (cn in nd$child_names) walk(nd$children[[cn]], paste0(path, "/", cn))
    }
  }
  walk(fit$root, "root")
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

# per-node convergence table (only populated when fit with n_chains > 1): leaves + each nest
# at the mean IV. Split-Rhat / bulk-ESS across chains on the pooled coefficients.
nested_cut_convergence <- function(fit) {
  rows <- list(); walk <- function(nd, path) {
    if (!is.null(nd$conv)) rows[[path]] <<- data.frame(node = path,
        max_rhat = nd$conv$max_rhat, pct_rhat_lt_1.01 = 100 * nd$conv$frac_rhat_lt_1.01,
        pct_rhat_lt_1.1 = 100 * nd$conv$frac_rhat_lt_1.1, median_ess = nd$conv$median_ess,
        n_params = nd$conv$n_params,
        # per-group beta = the quantity prediction actually uses; NA on a no-RE fit
        max_rhat_re = nd$conv$max_rhat_re %||% NA_real_,
        pct_rhat_re_lt_1.01 = 100 * (nd$conv$frac_rhat_re_lt_1.01 %||% NA_real_),
        median_ess_re = nd$conv$median_ess_re %||% NA_real_,
        n_cells_re = nd$conv$n_cells_re %||% NA_integer_, stringsAsFactors = FALSE)
    if (isTRUE(nd$type == "nest")) for (cn in nd$child_names) walk(nd$children[[cn]], paste0(path, "/", cn))
  }
  walk(fit$root, "root"); if (length(rows)) do.call(rbind, rows) else NULL
}

# --- posterior-draw prediction (coherent cut draws) --------------------------
# Returns [n x J x D] array of fine-class probabilities. Draw d threads imputation
# m(d) through the whole tree + picks one parent draw per internal node.
.ncut_predict_draw <- function(nd, X_new, share, m, P, gpos = NULL) {
  if (nd$type == "singleton") { P[, nd$fine] <- P[, nd$fine] + share; return(P) }
  if (nd$type == "even")      { for (fc in nd$fine) P[, fc] <- P[, fc] + share / length(nd$fine); return(P) }
  if (nd$type == "leaf") {
    sub <- .ncut_softmax(.ncut_eta(X_new, nd$beta_draws[[m]], .ncut_at(nd$beta_re_draws, m), gpos))
    colnames(sub) <- nd$fine
    for (fc in nd$fine) P[, fc] <- P[, fc] + share * sub[, fc]; return(P)
  }
  Xn <- .ncut_apply_xmap(nd$xmap, X_new)                      # same focal routing as at fit time
  if (nd$use_iv && length(nd$iv_children)) {
    ivm <- vapply(nd$iv_children, function(cn) .ncut_node_iv(nd$children[[cn]], X_new, m, gpos), numeric(nrow(X_new)))
    colnames(ivm) <- paste0("IV_", nd$iv_children); Xn <- cbind(Xn, ivm)
  }
  d <- .ncut_nest_draw(nd, m)
  cp <- .ncut_softmax(.ncut_eta(Xn, d$b, d$bre, gpos)); colnames(cp) <- nd$child_names
  for (cn in nd$child_names) P <- .ncut_predict_draw(nd$children[[cn]], X_new, share * cp[, cn], m, P, gpos)
  P
}
# group_idx: the group (e.g. country) of each row of X_new, on the SAME id scale used at fit time.
# REQUIRED for a fit with random effects -- see the note on .ncut_draws_from_disk for why the pooled
# coefficients alone do not reproduce the data. Rows whose group was never seen fall back to pooled.
predict_nested_cut <- function(fit, X_new, D = 200, summary = c("draws", "mean"), group_idx = NULL) {
  summary <- match.arg(summary); X_new <- as.matrix(X_new); J <- length(fit$fine_classes)
  lev <- .ncut_group_levels(fit$root)
  gpos <- NULL
  if (!is.null(group_idx)) {
    if (is.null(lev)) warning("predict_nested_cut: fit has no per-group draws; group_idx ignored (pooled prediction).")
    else {
      if (length(group_idx) != nrow(X_new)) stop("predict_nested_cut: length(group_idx) != nrow(X_new)")
      gpos <- match(as.character(group_idx), lev)
      if (anyNA(gpos)) warning(sprintf("predict_nested_cut: %d row(s) in %d unseen group(s) -> pooled coefficients.",
                                       sum(is.na(gpos)), length(unique(group_idx[is.na(gpos)]))))
    }
  } else if (!is.null(lev)) {
    warning("predict_nested_cut: this fit HAS random effects but no group_idx was supplied. ",
            "Pooled-only prediction does not reproduce a group-heterogeneous MNL (softmax is nonlinear, ",
            "so neither mu nor mean_g(beta_g) is a valid stand-in). Pass group_idx.")
  }
  zero <- function() matrix(0, nrow(X_new), J, dimnames = list(NULL, fit$fine_classes))
  # "mean" ACCUMULATES a running sum -> O(n*J) memory (never holds all D draws); "draws" must
  # materialise the full [n x J x D] array (that IS the requested output — large at big n/D).
  ms <- if (D <= fit$M) round(seq(1, fit$M, length.out = D)) else rep_len(seq_len(fit$M), D)
  if (summary == "mean") {
    S <- zero(); for (i in seq_along(ms)) S <- S + .ncut_predict_draw(fit$root, X_new, rep(1, nrow(X_new)), ms[i], zero(), gpos)
    return(S / length(ms))
  }
  arr <- array(0, c(nrow(X_new), J, length(ms)), dimnames = list(NULL, fit$fine_classes, NULL))
  for (i in seq_along(ms)) arr[, , i] <- .ncut_predict_draw(fit$root, X_new, rep(1, nrow(X_new)), ms[i], zero(), gpos)
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
  Xnd <- .ncut_apply_xmap(nd$xmap, X)                   # node design (focal routed)
  Tm  <- .ncut_xmap_T(nd$xmap, ncol(X))                 # [P_orig x P_node]: X_node = X %*% Tm
  Xd <- Xnd; dIVc <- list(); Pn <- ncol(Xnd)
  if (nd$use_iv && length(nd$iv_children)) {
    ivm <- vapply(nd$iv_children, function(cn) .ncut_node_iv(nd$children[[cn]], X, m), numeric(nrow(X)))
    colnames(ivm) <- paste0("IV_", nd$iv_children); Xd <- cbind(Xnd, ivm)
    for (cn in nd$iv_children) dIVc[[cn]] <- .ncut_dIV_dx(nd$children[[cn]], X, m)
  }
  B <- .ncut_meanB(nd, m); P <- .ncut_softmax(Xd %*% B)             # [n x K]
  dv <- vector("list", ncol(B))
  for (k in seq_len(ncol(B))) {
    # gamma part in the ORIGINAL x-space: a macro total's slope is shared by every fine
    # focal column feeding it (chain rule through X_node = X %*% Tm). Const over pixels.
    g <- matrix(as.numeric(Tm %*% B[seq_len(Pn), k]), nrow(X), ncol(X), byrow = TRUE)
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
      # Keep only children that actually have draws: an empty/NULL entry would make vapply fail and
      # take the whole summary (i.e. the reporting for a finished fit) down with it.
      .ld <- nd$lambda_draws[vapply(nd$lambda_draws, function(v) length(v) > 0L, TRUE)]
      lam <- if (length(.ld)) t(vapply(.ld, function(v)
        c(median = median(v), q025 = unname(quantile(v, .025)), q975 = unname(quantile(v, .975))), numeric(3))) else NULL
      out[[path]] <<- list(level = path, classes = nd$child_names, coef = bm, lambda = lam)
      for (cn in nd$child_names) walk(nd$children[[cn]], paste0(path, "/", cn))
    }
  }
  walk(fit$root, "root"); out
}
