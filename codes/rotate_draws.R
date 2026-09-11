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

rotate_draws <- function(dirs, groups = NULL, groups_all = NULL, max_draws = 1000L,
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
      acc_re <- c(acc_re, list(f$postb_total)); acc_mu <- c(acc_mu, list(f$postb_pooled))
      if (verbose) cat(sprintf("  loaded %s chain %d: %s draws\n", basename(dirname(d)), ci, dim(f$postb_total)[4]))
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
  res[]
}

# bind a list of k x p x G x d arrays along the DRAW dimension
abind4 <- function(...) { L <- list(...); d <- dim(L[[1]])
  array(unlist(L), dim = c(d[1], d[2], d[3], sum(vapply(L, function(x) dim(x)[4], 1L)))) }
cbind3 <- function(...) { L <- list(...); d <- dim(L[[1]])
  array(unlist(L), dim = c(d[1], d[2], sum(vapply(L, function(x) dim(x)[3], 1L)))) }
