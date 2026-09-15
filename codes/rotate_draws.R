# =============================================================================
# rotate_draws.R — posterior draws -> the from/to contrast table the downscaler consumes
# =============================================================================
# The downscaling engine works in TRANSITIONS, not levels: what it needs for each covariate is the
# utility difference between a source and a destination class, per country. In a multinomial logit
# that is beta[,to,g] - beta[,from,g], which is baseline-free -- the baseline cancels in the
# difference, so the artifact does not depend on which class the sampler happened to use as
# reference. The baseline column itself is carried as an explicit zero so it participates as both a
# source and a destination like any other class.
#
# Summaries are taken ACROSS DRAWS per (covariate, from, to, country): median plus a 95% interval.
# Median rather than mean deliberately -- these are log-odds contrasts with occasional heavy tails on
# rare classes, where a mean is dragged by a handful of draws.
#
#   rot <- rotate_draws("<run>/posterior", groups_all = <all grid countries>)
#
# Countries present in the grid but not in the fit get the POOLED (fixed-effect) contrast, which is
# the correct fallback: no country deviation was estimated for them, so they take the EU-wide effect.
# They are flagged `pooled = TRUE` so a consumer can tell an estimate from a fallback.
# =============================================================================
suppressMessages({library(data.table)})

rotate_draws <- function(dirs, groups = NULL, groups_all = NULL, max_draws = 1000L, zero_cols = NULL,
                         probs = c(0.025, 0.975), chains = NULL, verbose = TRUE) {
  stopifnot(length(dirs) >= 1)
  acc_re <- NULL; acc_mu <- NULL; cn <- NULL; vn <- NULL; grp <- NULL
  for (d in dirs) {
    chs <- if (!is.null(chains)) chains else sort(unique(as.integer(
      sub(".*_chain_(\\d+)\\.qs$", "\\1", list.files(d, "^posterior_batch_.*_chain_\\d+\\.qs$")))))
    for (ci in chs) {
      f <- tryCatch(recover_mnlogit_posterior(d, chain_id = ci), error = function(e) NULL)
      if (is.null(f) || is.null(f$postb_total)) next
      if (is.null(cn)) { cn <- f$cat_names; vn <- f$var_names; grp <- dimnames(f$postb_total)[[3]] }
      # THIN PER CHAIN, on load. Accumulating every chain at full length first needs
      # k x p x G x draws x 8 bytes EACH -- at 84 x 43 x 25 x 2500 that is ~1.8 GB a chain, ~7 GB for
      # four, and the process simply dies. Thinning here bounds peak memory by max_draws regardless of
      # how many segments are pooled, which is the whole point of being able to extend a run.
      nd_i <- dim(f$postb_total)[4]
      per <- max(1L, ceiling(nd_i / max(1L, ceiling(max_draws / max(1L, length(chs) * length(dirs))))))
      idx <- seq(1L, nd_i, by = per)
      acc_re <- c(acc_re, list(f$postb_total[, , , idx, drop = FALSE]))
      acc_mu <- c(acc_mu, list(f$postb_pooled[, , idx, drop = FALSE]))
      if (verbose) cat(sprintf("  loaded %s chain %d: %d draws -> kept %d\n",
                               basename(dirname(d)), ci, nd_i, length(idx)))
      rm(f); gc(FALSE)
    }
  }
  if (!length(acc_re)) stop("no recoverable draws under: ", paste(dirs, collapse = ", "))
  B  <- do.call(function(...) abind4(...), acc_re)      # k x p x G x draws
  MU <- do.call(cbind3, acc_mu)                          # k x p x draws
  K <- dim(B)[1]; P <- dim(B)[2]; G <- dim(B)[3]; D <- dim(B)[4]
  # Thin to keep the rotation in memory: 1332 class pairs x K x G is already large, and a median does
  # not need every draw. Regular thinning preserves the chain mixture.
  if (D > max_draws) { idx <- round(seq(1, D, length.out = max_draws))
    B <- B[, , , idx, drop = FALSE]; MU <- MU[, , idx, drop = FALSE]; D <- length(idx) }
  # Group LABELS must be supplied. The streamed arrays carry only positional dimnames ("1".."26"),
  # which match no country in the grid -- so every country would fall to the pooled fallback and all
  # 26 country random effects would be silently discarded while the artifact still looked complete.
  if (!is.null(groups)) {
    if (length(groups) != G) stop(sprintf("groups has %d labels but the posterior has %d groups", length(groups), G))
    grp <- as.character(groups)
  } else if (is.null(grp) || !length(grp) || all(grepl("^[0-9]+$", grp))) {
    stop("group labels unavailable from the posterior (positional dimnames only). Pass groups = <fit$re_group_names>; ",
         "without them every group would be treated as unfitted.")
  }
  if (is.null(vn)) vn <- paste0("V", seq_len(K))
  if (is.null(cn)) cn <- paste0("C", seq_len(P))

  extra <- if (!is.null(groups_all)) setdiff(groups_all, grp) else character(0)
  Gall <- c(grp, extra)
  if (verbose) cat(sprintf("  %d covariates x %d classes x %d groups (%d fitted + %d pooled fallback) x %d draws\n",
                           K, P, length(Gall), G, length(extra), D))
  out <- vector("list", P * (P - 1L)); z <- 1L
  qn <- paste0("value_q", sub("^0\\.", "", format(probs)))
  for (fi in seq_len(P)) for (ti in seq_len(P)) {
    if (fi == ti) next
    dif <- B[, ti, , , drop = FALSE] - B[, fi, , , drop = FALSE]     # k x 1 x G x D
    dm  <- matrix(dif, nrow = K * G, ncol = D)
    md  <- matrixStats::rowMedians(dm); qq <- matrixStats::rowQuantiles(dm, probs = probs)
    dt <- data.table(ks = rep(vn, G), from_class = cn[fi], to_class = cn[ti],
                     group = rep(grp, each = K), value_median = md,
                     q1 = qq[, 1], q2 = qq[, 2], pooled = FALSE)
    if (length(extra)) {                                   # pooled fallback for unfitted countries
      dmu <- matrix(MU[, ti, , drop = FALSE] - MU[, fi, , drop = FALSE], nrow = K, ncol = D)
      mm  <- matrixStats::rowMedians(dmu); qm <- matrixStats::rowQuantiles(dmu, probs = probs)
      dt <- rbind(dt, data.table(ks = rep(vn, length(extra)), from_class = cn[fi], to_class = cn[ti],
                                 group = rep(extra, each = K), value_median = rep(mm, length(extra)),
                                 q1 = rep(qm[, 1], length(extra)), q2 = rep(qm[, 2], length(extra)),
                                 pooled = TRUE))
    }
    out[[z]] <- dt; z <- z + 1L
  }
  res <- rbindlist(out[seq_len(z - 1L)])
  setnames(res, c("q1", "q2"), qn)

  # INERT COVARIATES -> exactly zero. A covariate that is identically zero in the design cannot enter
  # any utility, so its contrast is zero by construction whatever the posterior happens to hold for
  # it. Forced rather than trusted: such a column is dropped by the sampler's rank guard and was, up
  # to 2026-09-14, handed the un-centring shift because a zero-variance column was mistaken for an
  # intercept -- which put a country-varying coefficient on a covariate that carries no information,
  # straight into this artifact and on to the downscaler. Posteriors written before that fix are
  # repaired by this step, exactly: the corruption was additive on those rows alone and never
  # re-entered the chain (verified -- 83 of 84 rows and the log-likelihood bit-identical).
  if (length(zero_cols)) {
    .hit <- unique(res$ks[res$ks %in% zero_cols])
    if (length(.hit)) {
      .vcols <- intersect(names(res), c("value_median", qn, "value"))
      for (.v in .vcols) res[ks %in% .hit, (.v) := 0]
      message(sprintf(">>> rotate_draws: zeroed %d inert covariate(s) (identically zero in the design): %s",
                      length(.hit), paste(.hit, collapse = ", ")))
    }
  }
  res[]
}

# bind a list of k x p x G x d arrays along the DRAW dimension
abind4 <- function(...) { L <- list(...); d <- dim(L[[1]])
  array(unlist(L), dim = c(d[1], d[2], d[3], sum(vapply(L, function(x) dim(x)[4], 1L)))) }
cbind3 <- function(...) { L <- list(...); d <- dim(L[[1]])
  array(unlist(L), dim = c(d[1], d[2], sum(vapply(L, function(x) dim(x)[3], 1L)))) }
