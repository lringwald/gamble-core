#!/usr/bin/env Rscript
# =============================================================================
# BMLEH_Los1_CAPRI — performance report + the rotated artifact for the downscaler
# =============================================================================
#   Rscript projects/BMLEH_Los1_CAPRI/gamble_model/make_report.R <posterior-or-run dir>
#
# Writes, next to the run:
#   report.html               convergence, held-out performance, per-class fit, drivers
#   beta_rotated.rds / .csv.gz  (covariate x from_class x to_class x country) median + 95% interval
#
# The report leads with convergence, not fit, because a good held-out score on an unconverged chain
# is not evidence of anything -- and this model's per-parameter medians look healthy while the joint
# chain does not mix.
# =============================================================================
suppressMessages({library(Rcpp); library(RcppArmadillo); library(data.table); library(matrixStats)})
a <- commandArgs(trailingOnly = TRUE)
RUN <- if (length(a)) a[1] else stop("usage: make_report.R <run dir>")
if (basename(RUN) == "posterior") RUN <- dirname(RUN)
D <- file.path(RUN, "posterior"); if (!dir.exists(D)) stop("no posterior/ under ", RUN)
# Generalised from projects/BMLEH_Los1_CAPRI/gamble_model/make_report.R (2026-09-17) so one report
# serves every project. The design path and the label come from the CALLER; the project name is
# derived from the run directory (results/gamble_model/<PROJECT>/prior_<ts>) when not given.
INPUT   <- Sys.getenv("REPORT_INPUT", Sys.getenv("BM_INPUT", ""))
PROJECT <- Sys.getenv("REPORT_PROJECT", basename(dirname(RUN)))
if (!nzchar(INPUT)) stop("set REPORT_INPUT=<design dump>  (the .rds the fit was built on)")
if (!file.exists(INPUT)) stop("design dump not found: ", INPUT)
source("codes/mnl_aux_func.R"); source("codes/mnlogit_rcpp_sym.R"); source("codes/rotate_draws.R")

inp <- readRDS(INPUT)
# CAPRI codes are 2-char uppercase; GLOB_country gives names. Decide from the labels themselves rather
# than trusting a flag, since the dump does not record which column produced them.
# Prefer what the dump RECORDS over what the labels look like. The shape heuristic is only a fallback
# for dumps built before re_group_col was written.
GRPCOL <- if (!is.null(inp$re_group_col)) inp$re_group_col else
          if (all(grepl("^[A-Z]{2}$", inp$re_group_names))) "CAPRI_NUTS" else "GLOB_country"
CAPRI_KEYED <- identical(GRPCOL, "CAPRI_NUTS")
GRP_SLICED <- if (!is.null(inp$re_group_sliced)) isTRUE(inp$re_group_sliced) else CAPRI_KEYED
X <- as.matrix(inp$X_mat); X[!is.finite(X)] <- 0
Y <- as.matrix(inp$Y_pixel); g <- inp$group_idx_vec
keep <- colSums(Y) > 0; Y <- Y[, keep, drop = FALSE]; Y <- Y / rowSums(Y); cats <- colnames(Y)
# HELD-OUT SET: REUSE THE FIT'S OWN, never recompute it. prior_fit.rds records `rows` -- the exact
# absolute design indices the fit scored -- so the report measures the same pixels the fit did. This
# used to hardcode set.seed(20260909) + a random 20%, which silently scored a DIFFERENT hold-out than
# the fit for any run using the hybrid spatial split (block + country-stratified random). A report
# that disagrees with the fit it describes is worse than no report.
.pf <- file.path(RUN, "prior_fit.rds")
if (file.exists(.pf) && !is.null((.z <- readRDS(.pf))$rows)) {
  te <- sort(as.integer(.z$rows))
  cat(sprintf("held-out: reusing the fit's own %d scored rows (%s split)\n",
              length(te), .z$split$mode %||% "unrecorded"))
} else {
  warning("no prior_fit.rds$rows under ", RUN, " -- falling back to a RANDOM 20% split, which will ",
          "NOT match a spatially-split fit. Scores below describe different pixels than the fit did.",
          call. = FALSE)
  set.seed(20260909); te <- sort(sample(nrow(X), round(nrow(X) * 0.2)))
}
# Buffer rows (dropped from training by the spatial split) are not recorded, so the complement
# over-counts training slightly. It is used only to decide which GROUPS were seen, which it gets right.
tr <- setdiff(seq_len(nrow(X)), te)
gl <- unique(g[tr]); gi <- match(g[te], gl); rows <- which(!is.na(gi))
Xt <- X[te, ][rows, , drop = FALSE]; Yt <- Y[te, ][rows, , drop = FALSE]; git <- gi[rows]; J <- ncol(Y)

chs <- sort(unique(as.integer(sub(".*_chain_(\\d+)\\.qs$", "\\1",
        list.files(D, "^posterior_batch_.*_chain_\\d+\\.qs$")))))
F <- lapply(chs, function(ci) recover_mnlogit_posterior(D, chain_id = ci))
nd <- vapply(F, function(f) dim(f$postb_total)[4], 1L)
pred <- function(B) { U <- t(vapply(seq_len(nrow(Xt)), function(i) as.numeric(Xt[i, ] %*% B[, , git[i]]), numeric(J)))
  U <- U - apply(U, 1, max); P <- exp(U); P / rowSums(P) }
Bpool <- Reduce(`+`, lapply(F, function(f) apply(f$postb_total, c(1,2,3), sum))) / sum(nd)
P <- pred(Bpool)
ll <- sum(Yt * log(pmax(P, 1e-12)))
nl <- sum(Yt * log(pmax(matrix(colMeans(Yt), nrow(Yt), J, byrow = TRUE), 1e-12)))

# ---- convergence, per structural block -------------------------------------------------------
rh <- function(M) { M <- M[is.finite(rowSums(M)), , drop = FALSE]; n <- nrow(M)
  if (n < 4 || ncol(M) < 2) return(NA_real_)
  W <- mean(apply(M, 2, var)); Bv <- n * var(colMeans(M)); sqrt(((n - 1)/n * W + Bv/n)/W) }
ess <- function(M) { M <- M[is.finite(rowSums(M)), , drop = FALSE]
  if (nrow(M) < 8) return(NA_real_)
  sum(apply(M, 2, function(x) { ac <- acf(x, plot = FALSE, lag.max = min(50, length(x) %/% 4))$acf[-1]
    k <- which(ac < 0.05)[1]; if (is.na(k)) k <- length(ac)
    length(x) / (1 + 2 * sum(ac[seq_len(max(1, k - 1))])) })) }
blk <- function(mats, nm) { r <- vapply(mats, rh, 1); e <- vapply(mats, ess, 1)
  data.table(block = nm, n = length(mats), rhat_med = median(r, na.rm = TRUE), rhat_max = max(r, na.rm = TRUE),
             frac_bad = mean(r > 1.05, na.rm = TRUE), ess_med = median(e, na.rm = TRUE), ess_min = min(e, na.rm = TRUE)) }
LL <- do.call(cbind, lapply(F, function(f) f$post_log_lik))
mu_l <- lapply(seq_len(dim(F[[1]]$postb_pooled)[1]), function(i) lapply(seq_len(dim(F[[1]]$postb_pooled)[2]),
  function(j) do.call(cbind, lapply(F, function(f) f$postb_pooled[i, j, ]))))
mu_l <- unlist(mu_l, recursive = FALSE)
sre <- list(do.call(cbind, lapply(F, function(f) colMeans(f$post_sigma_re))))
c2l <- if (!is.null(F[[1]]$post_slab_c2)) list(do.call(cbind, lapply(F, function(f) f$post_slab_c2))) else NULL
conv <- rbindlist(c(list(blk(list(LL), "log_lik (joint)"), blk(mu_l, "mu (fixed effects)"),
  blk(sre, "sigma_re")), if (!is.null(c2l)) list(blk(c2l, "slab_c2")) else NULL))

per <- data.table(class = cats, obs = colMeans(Yt), pred = colMeans(P))
per[, `:=`(ratio = pred / pmax(obs, 1e-12), dLL = colSums(Yt * log(pmax(P, 1e-12))))]
setorder(per, -obs)
m <- apply(Reduce(`+`, lapply(F, function(f) apply(f$postb_pooled, c(1,2), sum))) / sum(nd), c(1,2), identity)
# ---- inert covariates: identically zero in the design, therefore zero effect -------------------
# Such a column cannot enter any utility, so any value the posterior carries for it is an artifact.
# It is repaired here rather than only in the artifact, so the drivers table, the heat plots and the
# country tables all agree. See rotate_draws() for the defect this repairs and why it is exact.
.inert <- which(apply(X, 2, function(z) { v <- var(z); (is.na(v) || v < 1e-12) && !all(z == 1) })) 
INERT <- if (length(.inert)) colnames(inp$X_mat)[.inert] else character(0)
if (length(INERT)) {
  .bad <- INERT[vapply(.inert, function(i) max(abs(Bpool[i, , ])) > 1e-10 || max(abs(m[i, ])) > 1e-10, TRUE)]
  cat(sprintf(">>> inert covariates (identically zero in the design): %s%s\n",
              paste(INERT, collapse = ", "),
              if (length(.bad)) sprintf("  [REPAIRED non-zero posterior on: %s]", paste(.bad, collapse = ", ")) else ""))
  for (i in .inert) { Bpool[i, , ] <- 0; m[i, ] <- 0 }
}

# COMPARABLE SCALES. A raw coefficient is per UNIT of its covariate, so log1p_Pop and flat_share are
# not on the same footing and ranking them against each other is meaningless. Multiply by the
# covariate's SD to get a PER-SD effect: "what a one-standard-deviation move in this driver does to
# this class's utility". This is the same convention the count-model report uses (sweep by sd_x).
# Rule B (codes/mnl_aux_func.R display_sds): per SD for UNBOUNDED covariates, per UNIT for shares and
# bounded indices. This used to be a bare apply(X, 2, sd), i.e. rule A -- which scaled focal_* shares
# and soil levels by their SD too, so the heat panels and beta_rotated (which uses display_sds) were
# on DIFFERENT scales inside one report while both captions said "per SD".
SDX <- display_sds(X); SDX[!is.finite(SDX) | SDX == 0] <- 1
m_sd <- m * SDX                                            # [cov x class], per-SD
drv <- data.table(cov = if (!is.null(F[[1]]$var_names)) F[[1]]$var_names else paste0("V", seq_len(nrow(m))),
                  rms_per_sd = sqrt(rowMeans(m_sd^2)),
                  max_per_sd = apply(abs(m_sd), 1, max))[order(-rms_per_sd)]

# ---- country effects: how far each country departs from the pooled model -----------------------
# Bpool is [cov x class x group] (posterior mean of b_g); m is [cov x class] (the pooled mu). The
# RANDOM EFFECT is the difference, so a country's departure from the pooled response is Bpool[,,k]-m.
#
# Reported as an RMS over classes rather than a signed value: in a sum-to-zero MNL a country's
# deviation is spread across classes with both signs, so a mean would cancel to nearly nothing and
# read as "no country effect" when the effect is real. Magnitude is the honest summary; the signed,
# per-class numbers are in beta_rotated, which is the artifact meant for that.
vn <- if (!is.null(F[[1]]$var_names)) F[[1]]$var_names else paste0("V", seq_len(dim(Bpool)[1]))
stopifnot(length(inp$re_group_names) == dim(Bpool)[3])
re_dev <- vapply(seq_len(dim(Bpool)[3]), function(k) Bpool[, , k] - m,
                 matrix(0, dim(Bpool)[1], dim(Bpool)[2]))          # [cov x class x group]
.int_i <- which(vn == "intercept")
# intercept_share is a SHARE OF SUMMED SQUARES, not a ratio of the two RMS columns: re_rms averages
# over all 84 covariate rows and intercept_rms over one, so their ratio is not bounded by 1 and is
# not a share of anything. (It first shipped that way here and read 6.8-8.9.)
cty <- data.table(
  country = inp$re_group_names,
  pixels  = as.integer(table(factor(g, levels = seq_along(inp$re_group_names)))),
  re_rms  = apply(re_dev, 3, function(A) sqrt(mean(A^2))),
  intercept_rms = if (length(.int_i)) apply(re_dev, 3, function(A) sqrt(mean(A[.int_i, ]^2))) else NA_real_,
  intercept_share = if (length(.int_i))
    apply(re_dev, 3, function(A) sum(A[.int_i, ]^2) / pmax(sum(A^2), 1e-12)) else NA_real_)
setorder(cty, -re_rms)

# country x covariate RE magnitude, for the top drivers only: 84 covariates x 25 countries is a wall
# of numbers, and the covariates that carry no country variation carry no information either.
dev_cc <- t(apply(re_dev, 3, function(A) sqrt(rowMeans(A^2))))      # [group x cov]
dimnames(dev_cc) <- list(inp$re_group_names, vn)
top_cov <- names(sort(colMeans(dev_cc), decreasing = TRUE))[seq_len(min(12L, ncol(dev_cc)))]
cc <- as.data.table(dev_cc[, top_cov, drop = FALSE], keep.rownames = "country")
setorderv(cc, top_cov[1], -1)

# ---- the rotated artifact ---------------------------------------------------------------------
GA <- tryCatch({ suppressMessages(library(arrow))
  gw <- Sys.getenv("GAMBLE_GRIDWORK_DIR", "../LAMASUS_gridwork/output")
  gf <- sort(list.files(gw, "^one_kmID_master_mapping_.*\\.parquet$", full.names = TRUE))
  gg <- as.data.table(arrow::read_parquet(gf[length(gf)]))
  # Take the full country list from the SAME keying the fit used. Reading GLOB_country while the fit
  # is keyed on CAPRI matches nothing, so every country would land in the pooled fallback and the
  # artifact would silently carry no country effects at all.
  if (!GRPCOL %in% names(gg)) stop("the fit's RE column '", GRPCOL, "' is not in the grid mapping")
  x <- as.character(gg[[GRPCOL]])
  if (GRP_SLICED) x <- substr(x, 1, 2)          # the driver slices CAPRI_NUTS to 2 chars for the RE key
  x <- unique(x); x[!is.na(x) & nzchar(x)]
}, error = function(e) NULL)
# sds = the DISPLAY scale (rule B, codes/mnl_aux_func.R). Without it the artifact carries raw
# coefficients in raw covariate units, which span four orders of magnitude here -- ranking them raw
# puts climate last when per SD it is one of the largest families.
rot <- rotate_draws(D, groups = inp$re_group_names, groups_all = GA, max_draws = 800L,
                    zero_cols = INERT, verbose = FALSE, sds = display_sds(X))
saveRDS(rot, file.path(RUN, "beta_rotated.rds")); fwrite(rot, file.path(RUN, "beta_rotated.csv.gz"))

# POSTERIOR SIGN-PROBABILITY per (covariate, class), accumulated chain by chain rather than by
# binding every draw into one array: 84 x 43 x all-draws is large enough to matter and nothing here
# needs the draws themselves, only how often they agree on a sign.
.pos <- matrix(0, nrow(m), ncol(m)); .tot <- 0L
for (f in F) {
  pb <- f$postb_pooled                                   # [cov x class x draws]
  .pos <- .pos + apply(pb > 0, c(1, 2), sum); .tot <- .tot + dim(pb)[3]
}
p_sign <- pmax(.pos, .tot - .pos) / max(.tot, 1L)        # in [0.5, 1]
p_sign[!is.finite(p_sign)] <- 0.5
# an inert covariate has no posterior to speak of -- its "certainty" is an artifact of a degenerate
# draw, so it must not render as a confident zero
if (length(.inert)) p_sign[.inert, ] <- 0.5

# ---- (1) CALIBRATION + per-class skill ---------------------------------------------------------
# McFadden alone says nothing about whether a predicted 30% IS 30%. Pool every (pixel, class) pair,
# bin by predicted share, and compare the bin's mean prediction against its mean observation. On the
# diagonal = calibrated; below = over-predicting that range; above = under-predicting.
# RESTRICT TO CELLS THAT SAY SOMETHING. Most (pixel, class) pairs are structurally zero -- 43
# classes, few present in any one pixel -- so quantile bins over everything put 11 of 12 bins at a
# predicted share of ~0 and the curve described the sparsity, not the calibration. Keep cells where
# the model predicts, or the data shows, a non-trivial share.
.keep_c <- as.numeric(P) > 1e-3 | as.numeric(Yt) > 1e-3
.pv <- as.numeric(P)[.keep_c]; .ov <- as.numeric(Yt)[.keep_c]
calib_n_kept <- sum(.keep_c); calib_n_all <- length(.keep_c)
.bk <- unique(stats::quantile(.pv, seq(0, 1, length.out = 13), na.rm = TRUE))
.bin <- cut(.pv, breaks = .bk, include.lowest = TRUE, labels = FALSE)
calib <- data.table(bin = .bin, pred = .pv, obs = .ov)[!is.na(bin),
           .(n = .N, pred = mean(pred), obs = mean(obs)), by = bin][order(pred)]
calib[, gap := obs - pred]
# one number for the whole curve: mean |obs-pred| across bins, weighted by bin size
calib_mae <- sum(abs(calib$gap) * calib$n) / sum(calib$n)

# per-class skill: each class's own McFadden against its own null (its held-out mean share)
# PER-CLASS CONTRIBUTION, in nats. The previous version reported 1 - ll_c/ll0_c and called it a
# per-class McFadden: it is not one. ll_c = sum_i Yt[i,c] log P[i,c] is this class's share-weighted
# contribution to the joint log-score, not a likelihood for class c on its own, so the ratio has no
# McFadden interpretation and is not bounded the way a pseudo-R2 is.
# What IS well defined is how many nats the model gains on this class against the null, and what
# fraction of the total gain that represents -- both additive and directly comparable across classes.
per_skill <- data.table(class = cats,
  obs = colMeans(Yt), pred = colMeans(P),
  ll  = colSums(Yt * log(pmax(P, 1e-12))),
  ll0 = colSums(Yt * log(pmax(matrix(colMeans(Yt), nrow(Yt), J, byrow = TRUE), 1e-12))))
per_skill[, `:=`(gain_nats = ll - ll0, ratio = pred / pmax(obs, 1e-12))]
per_skill[, share_of_gain := gain_nats / sum(gain_nats)]
setorder(per_skill, -gain_nats)

# (2) SUBSTITUTION: when a class is under-predicted, where does the mass go?
# PARTIAL correlation, not raw. Shares sum to 1, so the errors sum to ZERO by construction and every
# class is mechanically negatively correlated with every other. The largest class absorbs everyone
# else's residual, so a raw correlation names it as the partner for almost everything and says
# nothing about the model: under a NULL model (column means, no fitting) the dominant class came out
# as the most-negative partner for 36 of 43 classes here. Conditioning on all other classes removes
# that shared path and leaves genuine pairwise substitution.
# One class must be dropped first -- the compositional correlation matrix is singular by
# construction (rank J-1) -- and a small ridge keeps the inverse stable for rare classes.
.E <- P - Yt
.drop_c <- which.max(colSums(Yt))                       # the largest class carries the constraint
.Es <- .E[, -.drop_c, drop = FALSE]; .cs <- cats[-.drop_c]
.C <- suppressWarnings(cor(.Es)); .C[!is.finite(.C)] <- 0
.PR <- tryCatch(solve(.C + diag(1e-6, ncol(.C))), error = function(e) NULL)
if (!is.null(.PR)) {
  .d <- sqrt(diag(.PR))
  .PC <- -.PR / outer(.d, .d); diag(.PC) <- 0            # partial correlation
  subs <- data.table(class = .cs,
    partner = .cs[apply(.PC, 1, which.min)],
    partial_r = apply(.PC, 1, min))[order(partial_r)][seq_len(min(12L, length(.cs)))]
} else {
  subs <- data.table(class = character(0), partner = character(0), partial_r = numeric(0))
}
subs_dropped <- cats[.drop_c]

# ---- (2) WHAT THE COUNTRY RE BUYS --------------------------------------------------------------
# Rescore the SAME held-out pixels with the pooled coefficients only (mu, no country deviation).
# Cheap: both arrays are already in hand. Tests, per run, the standing claim that the country RE
# carries nearly all of the log-score skill.
.pred_pooled <- function() { U <- Xt %*% m; U <- U - apply(U, 1, max); Pp <- exp(U); Pp / rowSums(Pp) }
P_pooled  <- .pred_pooled()
ll_pooled <- sum(Yt * log(pmax(P_pooled, 1e-12)))
re_gain   <- data.table(
  model    = c("pooled (mu only)", "with country RE"),
  log_lik  = c(ll_pooled, ll),
  mcfadden = c(1 - ll_pooled/nl, 1 - ll/nl))
re_gain[, gain_nats := log_lik - ll_pooled]

# ---- (3) SHRINKAGE PROFILE ---------------------------------------------------------------------
# The horseshoe does the heavy lifting and is otherwise invisible. Sorted |coefficient| over every
# (covariate, class) cell, plus how many survive at each magnitude -- the "effective" model size as
# against the nominal one.
# On the PER-SD scale, so a threshold means the same thing for every covariate.
.co <- abs(as.numeric(m_sd)); .co <- sort(.co[is.finite(.co)], decreasing = TRUE)
shrink <- data.table(
  threshold = c(1, 0.5, 0.1, 0.05, 0.01, 0.001),
  n_above   = vapply(c(1, 0.5, 0.1, 0.05, 0.01, 0.001), function(t) sum(.co > t), 0L))
shrink[, pct := round(100 * n_above / length(.co), 1)]
n_credible <- sum(p_sign >= 0.9, na.rm = TRUE)

# ---- (4) DATA / DESIGN SUMMARY -----------------------------------------------------------------
prev <- data.table(class = cats, area_share = colSums(Y) / sum(Y))[order(-area_share)]
prev[, cum := cumsum(area_share)]
n_rare <- sum(prev$area_share < 0.001)
# BASELINE: report the one the FIT used, not only the one the design recorded. estimate_prior.R
# takes `bl <- which.max(colSums(Y))` (the most prevalent class) and ignores the driver's configured
# BASELINE_CLASS, so the two can disagree -- here the design says Natural_unmanaged while the fit
# used Forests_managed. Printing only the recorded one describes a model that was not fitted.
.bl_fit    <- cats[which.max(colSums(Y))]
.bl_design <- if (!is.null(inp$baseline_class)) as.character(inp$baseline_class) else "n/a"
dsum <- data.table(
  item  = c("pixels (fit)", "pixels (held out)", "classes", "covariates", "RE groups",
            "RE-bearing covariates", "inert covariates",
            "baseline class (used by the fit)", "baseline class (recorded in the design)",
            "classes below 0.1% of area", "outcome year"),
  value = c(format(nrow(X), big.mark = ","), format(nrow(Yt), big.mark = ","), J, ncol(X),
            length(inp$re_group_names),
            sum(apply(re_dev, 1, function(A) sqrt(mean(A^2))) > 1e-10), length(INERT),
            .bl_fit,
            paste0(.bl_design, if (!identical(.bl_fit, .bl_design)) "   <- DIFFERS from the fit" else ""),
            n_rare, if (!is.null(inp$out_year)) paste(unique(inp$out_year), collapse = ",") else "n/a"))

# ---- (5) SPATIAL RESIDUALS ---------------------------------------------------------------------
# Total-variation distance per held-out pixel, binned onto a coarse grid. Structure in this map is a
# missing covariate: a well-specified model should leave residuals that look like noise in space.
.hold <- te[rows]                                        # held-out row indices in the ORIGINAL design
spat <- NULL
if (!is.null(inp$coord_X) && !is.null(inp$coord_Y) && length(inp$coord_X) >= max(.hold)) {
  tvd <- 0.5 * rowSums(abs(P - Yt))
  sp  <- data.table(x = inp$coord_X[.hold], y = inp$coord_Y[.hold], tvd = tvd)
  sp  <- sp[is.finite(x) & is.finite(y) & is.finite(tvd)]
  if (nrow(sp) > 50) {
    NB <- 44L
    sp[, `:=`(bx = cut(x, NB, labels = FALSE), by = cut(y, NB, labels = FALSE))]
    spat <- sp[, .(tvd = mean(tvd), n = .N), by = .(bx, by)][n >= 3]
    spat_rng <- range(spat$tvd)
  }
}

# ---- effect matrices for the heat plots --------------------------------------------------------
# POOLED (mu): the model's common response, covariate x class. Signed, and the sign is the point.
# TOTAL (b_g): pooled + that country's random effect. Only the RE-bearing covariates differ between
# countries -- for every other covariate b_g IS mu by construction, so a country x covariate view of
# the total would be 78 identical columns. The country panel is therefore restricted to the
# covariates that actually carry an RE, and says so.
cls <- if (!is.null(colnames(Y))) colnames(Y) else paste0("C", seq_len(dim(m)[2]))

# ALL covariates, strongest first. Previously the top 22 only, which hid exactly the rows a reader
# goes looking for: the ones that turned out to be flat.
ord_mu  <- order(sqrt(rowMeans(m_sd^2)), decreasing = TRUE)      # rank on the comparable scale
mu_heat <- m_sd[ord_mu, , drop = FALSE]; rownames(mu_heat) <- vn[ord_mu]
# certainty for the tile SIZE: rescale P(sign) from [0.5,1] onto [0,1], so a coefficient whose sign
# is a coin flip gets the minimum tile and one whose sign never changes across draws fills the cell.
# MATRIX FIRST, floor second. pmax() copies attributes from its FIRST argument
# (mostattributes(mmm) <- attributes(elts[[1L]])), so pmax(0, M) returns a plain VECTOR -- the scalar
# 0 has no dim to copy. heat() then dies on CR[i, j] with "incorrect number of dimensions", which
# reads like a shape bug in the matrices and is not: every matrix upstream is correct.
mu_cred <- pmax((p_sign[ord_mu, , drop = FALSE] - 0.5) * 2, 0)

# TOTAL effect b_g = mu + RE, as country x covariate, for EVERY covariate. A covariate with no random
# effect has b_g == mu in every country, so its column is perfectly uniform -- which is the point:
# the pooled-only covariates identify themselves visually, with no separate legend to consult. Cells
# are the RMS over classes, because a signed mean cancels across the sum-to-zero classes.
re_mag   <- apply(re_dev, 1, function(A) sqrt(mean(A^2)))
re_rows  <- which(re_mag > 1e-10)                        # covariates that actually vary by country
# COVARIATE x COUNTRY, deliberately the same orientation and row order as the pooled panel above, so
# the two can be read across: row i is the same covariate in both.
# TOTAL effect b_g = mu + RE, per SD, as covariate x country. Signed: the mean over classes of the
# per-SD coefficient would cancel across the sum-to-zero classes, so take the class where this
# covariate acts most strongly under the POOLED model and report b_g for that class. That keeps a
# sign, keeps the scale comparable, and makes the pooled-only covariates identifiable -- their row is
# flat across countries because b_g IS mu for them.
.ref_cls <- apply(abs(m_sd), 1, which.max)
tot_heat <- vapply(seq_len(dim(Bpool)[3]),
                   function(k) (Bpool[, , k] * SDX)[cbind(seq_len(nrow(m)), .ref_cls)],
                   numeric(nrow(m)))                     # [cov x group], signed, per SD
dimnames(tot_heat) <- list(vn, inp$re_group_names)
tot_heat <- tot_heat[ord_mu, , drop = FALSE]
tot_lab  <- paste0(vn, " \u2192 ", cats[.ref_cls])[ord_mu]   # which class each row refers to
# mark which covariates carry an RE, so "flat row" and "has no RE" can be checked against each other
# rather than one being inferred from the other
re_mark <- ifelse(seq_len(nrow(m)) %in% re_rows, "RE", "")[ord_mu]

# ---- run settings & sampler specification ------------------------------------------------------
# run_config.rds is written by estimate_prior.R from 2026-09-14. Runs fitted before that have only
# lineage.txt, so the switches are genuinely UNRECORDED for them -- say so rather than printing this
# script's current defaults, which would assert something about a run nobody checked.
RC <- { f <- file.path(RUN, "run_config.rds"); if (file.exists(f)) readRDS(f) else NULL }
LIN <- { f <- file.path(RUN, "lineage.txt"); if (file.exists(f)) readLines(f) else character(0) }
.kvdt <- function(l) data.table(setting = names(l),
                                value = vapply(l, function(v)
                                  paste(format(v, trim = TRUE), collapse = ", "), ""))
set_mcmc <- if (!is.null(RC)) .kvdt(RC$mcmc) else NULL
set_data <- if (!is.null(RC)) .kvdt(RC$data) else NULL
set_sw   <- if (!is.null(RC)) .kvdt(RC$switches) else NULL

# ---- HTML ---------------------------------------------------------------------------------------
esc <- function(x) gsub("<", "&lt;", as.character(x), fixed = TRUE)
tbl <- function(dt, digits = 3, hi = NULL) {
  dt <- copy(dt); for (j in names(dt)) if (is.numeric(dt[[j]])) set(dt, j = j, value = round(dt[[j]], digits))
  paste0("<table><thead><tr>", paste0("<th>", esc(names(dt)), "</th>", collapse = ""), "</tr></thead><tbody>",
    paste0(vapply(seq_len(nrow(dt)), function(i) {
      cls <- if (!is.null(hi) && hi(dt[i])) " class=\"warn\"" else ""
      paste0("<tr", cls, ">", paste0("<td>", vapply(dt[i], function(v) esc(format(v)), ""), "</td>", collapse = ""), "</tr>") }, ""),
      collapse = ""), "</tbody></table>") }
# Heat table. Diverging (signed effects, blue-white-red) or sequential (magnitudes, white-red).
# vlim is the 95th percentile of |value| rather than the max: one extreme cell would otherwise
# compress every other cell to white and the plot would show nothing but that outlier.
# ---- inline SVG helpers (no external libraries; the report must open anywhere) ------------------
.svg_calib <- function(d, w = 420, h = 300, pad = 42) {
  rng <- range(c(d$pred, d$obs), finite = TRUE); if (diff(rng) <= 0) rng <- rng + c(-1, 1) * 1e-3
  sx <- function(v) pad + (v - rng[1]) / diff(rng) * (w - pad - 8)
  sy <- function(v) h - pad - (v - rng[1]) / diff(rng) * (h - pad - 8)
  pts <- paste(sprintf("%.1f,%.1f", sx(d$pred), sy(d$obs)), collapse = " ")
  dots <- paste(sprintf('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="#c1543b" fill-opacity=".75"/>',
                  sx(d$pred), sy(d$obs), 2 + 4 * sqrt(d$n / max(d$n))), collapse = "")
  paste0('<svg viewBox="0 0 ', w, ' ', h, '" width="', w, '" height="', h, '" role="img">',
    sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#9aa4ae" stroke-dasharray="4 3"/>',
            sx(rng[1]), sy(rng[1]), sx(rng[2]), sy(rng[2])),
    sprintf('<polyline points="%s" fill="none" stroke="#c1543b" stroke-width="1.6"/>', pts), dots,
    sprintf('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#5b6672"/>', pad, h-pad, w-8, h-pad),
    sprintf('<line x1="%d" y1="%d" x2="%d" y2="%.1f" stroke="#5b6672"/>', pad, 8, pad, h-pad),
    sprintf('<text x="%.1f" y="%d" font-size="11" fill="#5b6672" text-anchor="middle">predicted share</text>', (w+pad)/2, h-10),
    sprintf('<text x="12" y="%.1f" font-size="11" fill="#5b6672" text-anchor="middle" transform="rotate(-90 12 %.1f)">observed share</text>', h/2, h/2),
    '</svg>') }

.svg_shrink <- function(v, w = 420, h = 300, pad = 46) {
  v <- v[v > 0]; n <- length(v); lo <- log10(max(min(v), 1e-6)); hi <- log10(max(v))
  sx <- function(i) pad + (i - 1) / max(n - 1, 1) * (w - pad - 8)
  sy <- function(z) h - pad - (log10(z) - lo) / max(hi - lo, 1e-9) * (h - pad - 8)
  step <- max(1L, as.integer(n / 600)); idx <- seq(1, n, by = step)
  pts <- paste(sprintf("%.1f,%.1f", sx(idx), sy(pmax(v[idx], 10^lo))), collapse = " ")
  gl <- paste(vapply(c(1, 0.1, 0.01, 0.001), function(t) if (t >= 10^lo && t <= 10^hi)
        sprintf('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e4e8ec"/><text x="%d" y="%.1f" font-size="9" fill="#9aa4ae">%g</text>',
                pad, sy(t), w-8, sy(t), 6, sy(t)+3, t) else "", ""), collapse = "")
  paste0('<svg viewBox="0 0 ', w, ' ', h, '" width="', w, '" height="', h, '" role="img">', gl,
    sprintf('<polyline points="%s" fill="none" stroke="#3b6ec1" stroke-width="1.5"/>', pts),
    sprintf('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#5b6672"/>', pad, h-pad, w-8, h-pad),
    sprintf('<text x="%.1f" y="%d" font-size="11" fill="#5b6672" text-anchor="middle">coefficients, largest first (n=%s)</text>',
            (w+pad)/2, h-10, format(n, big.mark=",")),
    '</svg>') }

.svg_spatial <- function(d, w = 460, h = 430, pad = 8) {
  if (is.null(d) || !nrow(d)) return("")
  bx <- range(d$bx); by <- range(d$by)
  cw <- (w - 2*pad) / (diff(bx) + 1); ch <- (h - 2*pad) / (diff(by) + 1)
  q <- stats::quantile(d$tvd, c(.05, .95), na.rm = TRUE); lo <- q[1]; hi <- max(q[2], lo + 1e-9)
  col <- function(z) { t <- max(0, min(1, (z - lo) / (hi - lo)))
    sprintf("rgb(%d,%d,%d)", round(255 - 40*t), round(255 - 170*t), round(255 - 190*t)) }
  cells <- paste(sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s"><title>TVD %.3f (n=%d)</title></rect>',
      pad + (d$bx - bx[1]) * cw, pad + (by[2] - d$by) * ch, cw + .5, ch + .5,
      vapply(d$tvd, col, ""), d$tvd, d$n), collapse = "")
  paste0('<svg viewBox="0 0 ', w, ' ', h, '" width="', w, '" height="', h, '" role="img">', cells, '</svg>') }

# Heat tile, matching the count-model report's convention so the two read the same way.
#   SIZE   = posterior certainty of the sign. 0 sitting mid-posterior (p=0.5) -> 12% tile; 0 outside
#            the draws (p=1) -> full cell. Certainty is a CONTINUOUS quantity and a flat dot hides
#            the difference between 0.951 and 0.999.
#   COLOUR = signed effect, red positive / blue negative, opacity by |value| relative to vlim.
# So a large pale tile is "small but certain" and a small saturated tile is "large but uncertain" --
# the two cases a reader most needs to separate, and neither is visible from colour alone.
tile <- function(val, cred, vlim, title, txt) {
  if (!is.finite(val)) return('<td class="hm"><div class="tbg"></div><span class="hval"></span></td>')
  cc <- min(max(cred, 0), 1); size <- 12 + 88 * cc
  inten <- min(abs(val) / vlim, 1) * 0.85
  col <- if (val >= 0) sprintf("rgba(178,24,43,%.3f)", 0.05 + inten)
         else          sprintf("rgba(33,102,172,%.3f)", 0.05 + inten)
  sprintf('<td class="hm" title="%s"><div class="tbg"></div><div class="tl" style="width:%.1f%%;height:%.1f%%;background:%s"></div><span class="hval">%s</span></td>',
          esc(title), size, size, col, txt)
}

# M = effect matrix, CR = matching certainty matrix in [0,1]. Values are ALREADY per-SD when the
# caller scaled them; vlim is the 95th percentile of |value| so one outlier cannot wash the rest out.
heat <- function(M, CR, rowlab, collab, digits = 2, rowmark = NULL, unit = "per SD") {
  v <- abs(as.numeric(M)); v <- v[is.finite(v)]
  vlim <- stats::quantile(v, 0.95, na.rm = TRUE); if (!is.finite(vlim) || vlim == 0) vlim <- 1
  hd <- paste0("<th></th>", paste0(sprintf('<th class="rot"><span>%s</span></th>', esc(collab)), collapse = ""))
  bd <- paste0(vapply(seq_len(nrow(M)), function(i) paste0(
      "<tr><th class=\"rl\">", if (!is.null(rowmark) && nzchar(rowmark[i]))
        sprintf('<span class="rm">%s</span>', esc(rowmark[i])) else "", esc(rowlab[i]), "</th>",
      paste0(vapply(seq_len(ncol(M)), function(j) tile(M[i,j], if (is.null(CR)) 1 else CR[i,j], vlim,
        sprintf("%s / %s: %.4g  (P(sign)=%.3f)", rowlab[i], collab[j], M[i,j],
                if (is.null(CR)) NA_real_ else CR[i,j]),
        if (is.finite(M[i,j]) && abs(M[i,j]) >= 10^(-digits)) formatC(M[i,j], format="f", digits=digits) else ""), ""),
        collapse = ""), "</tr>"), ""), collapse = "")
  paste0('<table class="heat"><thead><tr>', hd, '</tr></thead><tbody>', bd, '</tbody></table>',
         sprintf('<p class="sub">Tile SIZE = posterior certainty of sign (12%% at p=0.5, full cell at p=1). Colour = signed effect (%s), saturating at |%.3g|.</p>',
                 unit, vlim)) }

conv_bad <- conv[block == "log_lik (joint)", rhat_max] > 1.05
html <- paste0('<meta charset="utf-8"><title>', esc(PROJECT), ' prior — ', esc(basename(RUN)), '</title>
<style>
 /* House design system, shared with postprocess/build_html.R so the flat report, the nested report
    and this one read as ONE product: same green-teal accent, same surfaces, same dark mode. The
    three-way theme declaration is deliberate -- bare :root is the light palette, the media query
    follows the OS, and [data-theme] lets the toggle WIN in both directions. */
 :root{--bg:#FAFBFA;--surf:#fff;--ink:#14211E;--mut:#5E6E69;--line:#E5EAE7;--acc:#0F7A67;
       --warn:#B26B00;--warnbg:#FBF3E4;--good:#0F7A67}
 @media(prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#0E1613;--surf:#14201C;
       --ink:#E8F0EC;--mut:#93A39D;--line:#25332E;--acc:#54BCA6;--warn:#E0A050;--warnbg:#241d10;
       --good:#54BCA6}}
 :root[data-theme=dark]{--bg:#0E1613;--surf:#14201C;--ink:#E8F0EC;--mut:#93A39D;--line:#25332E;
       --acc:#54BCA6;--warn:#E0A050;--warnbg:#241d10;--good:#54BCA6}
 *{box-sizing:border-box}
 body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.6 ui-sans-serif,system-ui,
      -apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
 /* .page is the PAGE container. `.wrap` keeps its existing meaning in this report -- the
    horizontal scroll box around a wide table -- so it is left alone. */
 .page{max-width:1160px;margin:0 auto;padding:0 22px 72px}
 header.top{position:sticky;top:0;z-index:20;background:color-mix(in srgb,var(--bg) 88%,transparent);
      backdrop-filter:blur(8px);border-bottom:1px solid var(--line)}
 .hd{max-width:1160px;margin:0 auto;padding:14px 22px 10px}
 h1{font-size:22px;font-weight:640;letter-spacing:-.01em;margin:0 0 2px;text-wrap:balance}
 h2{font-size:15.5px;font-weight:620;margin:34px 0 10px;padding-bottom:.35rem;letter-spacing:-.005em;
    border-bottom:2px solid var(--line);scroll-margin-top:96px}
 .sub{color:var(--mut);font-size:13px;margin:0 0 1.3rem}
 .hd .sub{margin:0 0 10px}
 nav{display:flex;gap:3px;flex-wrap:wrap}
 nav a{color:var(--mut);text-decoration:none;font-size:11.5px;font-weight:600;padding:4px 9px;
       border:1px solid var(--line);border-radius:999px;white-space:nowrap}
 nav a:hover{color:var(--ink);border-color:var(--acc);box-shadow:0 0 0 1px var(--acc) inset}
 table{border-collapse:collapse;width:100%;margin:.6rem 0;font-size:12.8px}
 th{text-align:left;border-bottom:2px solid var(--line);padding:.42rem .55rem;color:var(--mut);
    font-weight:600;font-size:11px;text-transform:uppercase;letter-spacing:.03em}
 td{border-bottom:1px solid var(--line);padding:.36rem .55rem;font-variant-numeric:tabular-nums}
 tbody tr:hover td{background:color-mix(in srgb,var(--acc) 6%,transparent)}
 tr.warn td{background:var(--warnbg)}
 .box{border:1px solid color-mix(in srgb,var(--warn) 40%,var(--line));border-left:3px solid var(--warn);
      background:var(--warnbg);padding:13px 16px;margin:1rem 0;border-radius:10px;font-size:13.5px}
 .box.ok,.ok{border-color:color-mix(in srgb,var(--good) 35%,var(--line));border-left-color:var(--good);
      background:color-mix(in srgb,var(--good) 8%,var(--surf))}
 /* Stat cards. Same markup as before (.kv > div > b + span) -- only the presentation changed. */
 .kv{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:.7rem 0 1.4rem}
 .kv div{background:var(--surf);border:1px solid var(--line);border-radius:12px;padding:13px 14px;min-width:0}
 .kv b{display:block;font-size:23px;font-weight:650;letter-spacing:-.02em;font-variant-numeric:tabular-nums}
 .kv span{display:block;color:var(--mut);font-size:12px;margin-top:2px}
 code{background:color-mix(in srgb,var(--acc) 9%,var(--surf));padding:.1rem .32rem;border-radius:4px;
      font-size:12.5px}
 pre{background:var(--surf);border:1px solid var(--line);border-radius:10px;padding:12px 14px;
     overflow-x:auto;font-size:12.5px}
 pre code{background:none;padding:0}
 .wrap{overflow-x:auto;border:1px solid var(--line);border-radius:10px;background:var(--surf);
       padding:0 10px}
 /* FIGURES AND HEAT TILES STAY LIGHT IN BOTH THEMES, exactly as build_html.R pins its heat table.
    The diverging red/blue ramp and the sequential spatial ramp are calibrated against white; re-
    mapping them per theme would change what a colour MEANS between two readers of the same report. */
 svg{background:#fff;border:1px solid var(--line);border-radius:10px;display:block}
 table.heat{border-collapse:separate;border-spacing:1px;font-size:11px;width:auto;background:#fff;
   color:#14211E;border-radius:8px}
 table.heat th{border:0;text-transform:none;letter-spacing:0}
 table.heat th.rl{text-align:right;padding:.1rem .45rem;white-space:nowrap;border:0;font-weight:500;
   color:#14211E;font-size:11px}
 table.heat td.hm{width:28px;min-width:28px;height:22px;padding:0;border:0;position:relative;
   text-align:center;vertical-align:middle}
 table.heat tbody tr:hover td{background:none}
 table.heat td.hm .tbg{position:absolute;inset:0;background:#f7f8f9;border:1px solid #eef1f3}
 table.heat td.hm .tl{position:absolute;left:50%;top:50%;transform:translate(-50%,-50%)}
 table.heat td.hm .hval{position:relative;font-size:8.5px;font-variant-numeric:tabular-nums;
   color:#22262a;mix-blend-mode:multiply}
 table.heat th.rl .rm{display:inline-block;font-size:8px;letter-spacing:.04em;color:#0F7A67;
   font-weight:700;vertical-align:middle}
 th.rot{height:112px;vertical-align:bottom;padding:0;border:0;width:26px;min-width:26px}
 th.rot>span{writing-mode:vertical-rl;transform:rotate(180deg);white-space:nowrap;
   font-weight:500;font-size:10.5px;color:#5E6E69}
 footer{color:var(--mut);font-size:12px;margin-top:44px;border-top:1px solid var(--line);padding-top:16px}
 .themebtn{position:fixed;right:14px;bottom:14px;z-index:30;border:1px solid var(--line);
   background:var(--surf);color:var(--ink);border-radius:999px;padding:8px 12px;cursor:pointer;
   font:inherit;font-size:12px}
 @media(max-width:820px){.page,.hd{padding-left:16px;padding-right:16px}}
</style>
<header class="top"><div class="hd">
<h1>', esc(PROJECT), ' — prior module</h1>
<p class="sub">', esc(basename(RUN)), ' · ', length(chs), ' chains × ', nd[1], ' draws · ',
 format(nrow(X), big.mark = ","), ' pixels · ', J, ' classes · ', ncol(X), ' covariates · ',
 length(inp$re_group_names), ' REs on ', esc(GRPCOL), ' · generated ', format(Sys.time(), "%Y-%m-%d %H:%M"), '</p>
<nav><!--NAV--></nav></div></header>
<div class="page">

<h2>Run settings &amp; sampler specification</h2>',
 if (is.null(RC)) paste0('<div class="box"><b>Sampler switches were not recorded for this run.</b><br>',
   'It predates <code>run_config.rds</code> (2026-09-14). Only the lineage below is known; the ',
   'switches are NOT shown rather than guessed from current defaults, which would assert something ',
   'about a run nobody checked. Re-running writes the full record.</div>') else '',
 if (length(LIN)) paste0('<pre><code>', esc(paste(LIN, collapse = "\n")), '</code></pre>') else '',
 if (!is.null(RC)) paste0(
   '<p class="sub">MCMC</p><div class="wrap">', tbl(set_mcmc, 4), '</div>',
   '<p class="sub">Data &amp; design</p><div class="wrap">', tbl(set_data, 4), '</div>',
   '<p class="sub">Sampler switches &mdash; every non-default argument the sampler was called with, ',
   'including the ones hardcoded in <code>estimate_prior.R</code>: a switch being constant is not a ',
   'reason to hide it.</p><div class="wrap">', tbl(set_sw, 4), '</div>') else '', '

<h2>Convergence</h2>
<div class="box', if (!conv_bad) ' ok' else '', '"><b>',
 if (conv_bad) 'The joint chain has not mixed.' else 'Joint chain mixed.',
 '</b><br>', if (conv_bad) paste0('log_lik Rhat is ', round(conv[block=="log_lik (joint)", rhat_max],2),
 '. Note that <em>mu</em> looks healthy by its median — per-parameter medians flatter a chain whose ',
 'chains have each settled somewhere different. Point estimates of abundant classes are usable; ',
 'joint statements and every uncertainty interval are not.') else
 'Per-block Rhat is within tolerance; intervals can be read as posterior intervals.', '</div>',
 '<div class="wrap">', tbl(conv, 3, function(r) isTRUE(r$rhat_max > 1.05)), '</div>

<h2>Held-out performance</h2>
<div class="kv">
 <div><b>', round(1 - ll/nl, 4), '</b><span>McFadden R²</span></div>
 <div><b>', format(round(ll), big.mark=","), '</b><span>log-likelihood</span></div>
 <div><b>', format(round(nl), big.mark=","), '</b><span>null log-likelihood</span></div>
 <div><b>', format(nrow(Yt), big.mark=","), '</b><span>held-out pixels</span></div>
</div>

<h2>Calibration</h2>
<p class="sub">Binned by predicted share over the ', format(calib_n_kept, big.mark=","), ' of
 ', format(calib_n_all, big.mark=","), ' (pixel, class) cells where the model or the data shows a
 non-trivial share (&gt;0.001); the structurally-zero remainder would otherwise put 11 of 12 bins at
 ~0 and describe the sparsity rather than the calibration. On the dashed diagonal = calibrated.
 Mean absolute gap <b>', round(calib_mae, 4), '</b>. McFadden says how much better than the null the
 model is; this says whether a predicted 30% <em>is</em> 30%, which is the question a downscaler
 actually depends on.</p>
<div style="display:flex;gap:1.4rem;flex-wrap:wrap;align-items:flex-start">
 <div>', .svg_calib(calib), '</div>
 <div class="wrap" style="flex:1;min-width:280px">', tbl(calib[, .(bin, n, pred, obs, gap)], 4), '</div>
</div>

<h2>What the country random effect buys</h2>
<p class="sub">The same held-out pixels scored twice: once with the pooled coefficients alone, once
 with each country&rsquo;s own deviation. The difference is what the hierarchy is worth on this run,
 rather than as a general claim.</p>
<div class="wrap">', tbl(re_gain, 4), '</div>

<h2>Shrinkage</h2>
<p class="sub">All ', format(nrow(m) * ncol(m), big.mark = ","), ' (covariate &times; class)
 coefficients, largest first, on a log scale. The horseshoe does the regularising and is otherwise
 invisible in this report; this is the <em>effective</em> model size against the nominal one.
 <b>', format(n_credible, big.mark = ","), '</b> cells have posterior sign-probability &ge; 0.9.</p>
<div style="display:flex;gap:1.4rem;flex-wrap:wrap;align-items:flex-start">
 <div>', .svg_shrink(.co), '</div>
 <div class="wrap" style="flex:1;min-width:260px">', tbl(shrink, 3), '</div>
</div>

<h2>Where it misses, in space</h2>
<p class="sub">Total-variation distance between predicted and observed shares per held-out pixel,
 averaged onto a coarse grid (bins with &lt;3 pixels dropped). Structure here is a <em>missing
 covariate</em>: a well-specified model leaves residuals that look like noise in space. Darker = worse.
 Hover a cell for its value.</p>
', if (!is.null(spat)) paste0('<div>', .svg_spatial(spat), '</div>') else
   '<div class="box">No coordinates in the design dump &mdash; spatial residuals unavailable.</div>', '

<h2>Which classes get confused</h2>
<p class="sub"><b>Partial</b> correlation of held-out errors &mdash; the strongest negative partner
 after conditioning on every other class. Raw correlation cannot be used here: shares sum to 1, so
 errors sum to zero by construction and the largest class absorbs everyone&rsquo;s residual. Under a
 null model with no fitting at all, a raw correlation named the dominant class as the partner for 36
 of 43 classes &mdash; measuring the arithmetic, not the model. <code>', esc(subs_dropped), '</code>
 is dropped to make the matrix invertible (the composition is rank J&minus;1).</p>
<div class="wrap">', tbl(subs, 3), '</div>

<h2>Per-class skill</h2>
<p class="sub">Nats gained against the null on each class, and that class&rsquo;s share of the total
 gain &mdash; both additive, so they say directly where the model&rsquo;s skill comes from. (An earlier
 version reported a per-class &ldquo;McFadden&rdquo;; that ratio has no such interpretation, because a
 single class&rsquo;s share-weighted contribution is not a likelihood on its own.)
 ', n_rare, ' of ', J, ' classes hold less than 0.1% of area.</p>
<div class="wrap">', tbl(per_skill[, .(class, obs, pred, ratio, gain_nats, share_of_gain)], 4,
   function(r) isTRUE(r$ratio > 2 | r$ratio < 0.5)), '</div>

<h2>Design summary</h2>
<div class="wrap">', tbl(dsum), '</div>
<p class="sub">Class prevalence, largest first.</p>
<div class="wrap">', tbl(head(prev, 20), 5), '</div>

<h2>Per-class fit</h2>
<p class="sub">ratio = predicted ÷ observed share. Rows shaded where the model is out by more than 2×.</p>
<div class="wrap">', tbl(per[, .(class, obs, pred, ratio, dLL)], 5,
   function(r) isTRUE(r$ratio > 2 | r$ratio < 0.5)), '</div>

<h2>What drives the model</h2>
<div class="wrap">', tbl(head(drv, 20), 3), '</div>

<h2>Pooled effects &mu;</h2>
<p class="sub">All ', nrow(mu_heat), ' covariates &times; ', ncol(mu_heat), ' classes, <b>per standard
 deviation of the covariate</b> so the rows are comparable with each other, ordered by per-SD RMS.
 <span style="color:#c1543b">Red</span> raises that class&rsquo;s share, <span style="color:#3b6ec1">blue</span>
 lowers it; tile SIZE is the posterior certainty of the sign. Hover for the value and P(sign).</p>
<div class="wrap">', heat(mu_heat, mu_cred, rownames(mu_heat), cls), '</div>

<h2>Total effects by country (&mu; + RE)</h2>
<p class="sub"><code>b<sub>g</sub></code> = &mu; + that country&rsquo;s random effect, per SD, for
 <b>every</b> covariate. Each row is shown for the single class that covariate acts on most strongly
 under the pooled model (named in the row label), so the value keeps its sign and its scale.
 ', length(re_rows), ' of ', length(vn), ' covariates carry a random effect (marked
 <span class="rm">RE</span>); for the rest <code>b<sub>g</sub> = &mu;</code> exactly, so their
 <b>rows are flat across countries</b> &mdash; readable at a glance rather than from a list.</p>
<div class="wrap">', heat(tot_heat, NULL, tot_lab, colnames(tot_heat), rowmark = re_mark), '</div>

<h2>Country effects (fixed + random)</h2>
<p class="sub">How far each country departs from the pooled model, as RMS over classes of
 <code>b<sub>g</sub> &minus; &mu;</code>. Signed per-class values are in <code>beta_rotated</code>;
 a signed summary here would cancel across the sum-to-zero classes and read as no effect at all.</p>
<div class="wrap">', tbl(cty, 4), '</div>
<p class="sub">Per-driver RE magnitude, top ', length(top_cov), ' drivers by country variation.</p>
<div class="wrap">', tbl(cc, 4), '</div>

<h2>Artifact for the downscaler</h2>', if (!CAPRI_KEYED) paste0('<div class="box"><b>Group keys are NOT CAPRI codes.</b><br>This fit is grouped on <code>', esc(GRPCOL), '</code>, giving labels like <code>', esc(inp$re_group_names[1]), '</code>. The downscaler keys on CAPRI country codes, which are not ISO (BL, IR, EL, CS, KO, MO) and are <em>not</em> translatable from these after the fact — the geometries disagree at borders, so the pixel-level mapping is many-to-many. Re-fit with <code>DRIVER_RE_GROUP_COL=CAPRI_NUTS</code> before using this downstream: the table below is the right <em>shape</em> on the wrong keys.</div>') else '<div class="box ok">Group keys are CAPRI country codes, as the downscaler expects.</div>', '
<p><code>beta_rotated.rds</code> · <code>beta_rotated.csv.gz</code> — ',
 format(nrow(rot), big.mark=","), ' rows: covariate × from_class × to_class × country, as
 <code>value_median</code> with a 95% interval. The contrast is <code>beta[,to,g] − beta[,from,g]</code>,
 so it is baseline-free and antisymmetric by construction. ',
 uniqueN(rot[pooled == TRUE, group]), ' of ', uniqueN(rot$group), ' countries are not in the fit and
 carry the pooled fixed-effect contrast instead, flagged <code>pooled = TRUE</code>.</p>
<footer>', esc(PROJECT), ' · ', esc(basename(RUN)), ' · generated ',
 format(Sys.time(), "%Y-%m-%d %H:%M"), ' by <code>postprocess/model_report.R</code>. Held-out rows
 are the fit\'s own; see the convergence note before quoting any interval.</footer>
</div>
<button class="themebtn" id="themebtn">☾ dark</button>
<script>
const tb=document.getElementById("themebtn");
function cur(){return document.documentElement.getAttribute("data-theme")||
  (matchMedia("(prefers-color-scheme:dark)").matches?"dark":"light")}
tb.textContent=cur()==="dark"?"☀ light":"☾ dark";
tb.addEventListener("click",()=>{const n=cur()==="dark"?"light":"dark";
  document.documentElement.setAttribute("data-theme",n);
  tb.textContent=n==="dark"?"☀ light":"☾ dark"});
</script>')

# ---- section anchors + nav -----------------------------------------------------------------
# Ids are DERIVED FROM THE HEADINGS, so a section added later is picked up automatically and the
# nav can never drift out of step with the body -- the usual failure of a hand-maintained contents
# list. Presentation only: no heading text is altered.
.titles <- regmatches(html, gregexpr("(?<=<h2>)[^<]*(?=</h2>)", html, perl = TRUE))[[1]]
.slug <- function(x) { x <- gsub("&[#a-zA-Z0-9]+;", "", x)
                       gsub("(^-|-$)", "", tolower(gsub("[^A-Za-z0-9]+", "-", x))) }
.ids <- .slug(.titles)
for (.i in seq_along(.titles))
  html <- sub(paste0("<h2>", .titles[.i], "</h2>"),
              paste0('<h2 id="', .ids[.i], '">', .titles[.i], '</h2>'), html, fixed = TRUE)
html <- sub("<!--NAV-->",
            paste0(sprintf('<a href="#%s">%s</a>', .ids, .titles), collapse = ""), html, fixed = TRUE)
writeLines(html, file.path(RUN, "report.html"))

# ---- markdown ------------------------------------------------------------------------------------
# The same report in a form that DIFFS. The HTML is the readable artifact; this one is the one you can
# put in git and compare between runs, which is how a regression in held-out score or a country
# slipping out of the fit becomes visible without opening anything.
mdt <- function(dt, digits = 4) {
  dt <- copy(dt); for (j in names(dt)) if (is.numeric(dt[[j]])) set(dt, j = j, value = round(dt[[j]], digits))
  hdr <- paste0("| ", paste(names(dt), collapse = " | "), " |")
  sep <- paste0("|", paste(rep(" --- ", ncol(dt)), collapse = "|"), "|")
  rows <- vapply(seq_len(nrow(dt)), function(i)
    paste0("| ", paste(vapply(dt[i], function(v) format(v), ""), collapse = " | "), " |"), "")
  paste(c(hdr, sep, rows), collapse = "\n")
}
md <- paste0(
paste0("# ", PROJECT, " - prior module\n\n"),
"`", basename(RUN), "` | ", length(chs), " chains x ", nd[1], " draws | ",
  format(nrow(X), big.mark = ","), " pixels | ", J, " classes | ", ncol(X), " covariates | ",
  length(inp$re_group_names), " REs on `", GRPCOL, "` | generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n",
"## Run settings & sampler specification\n\n",
if (is.null(RC)) paste0("**Sampler switches were not recorded for this run.** It predates ",
  "`run_config.rds` (2026-09-14). Only the lineage below is known; the switches are omitted rather ",
  "than guessed from current defaults.\n\n") else "",
if (length(LIN)) paste0("```\n", paste(LIN, collapse = "\n"), "\n```\n\n") else "",
if (!is.null(RC)) paste0("### MCMC\n\n", mdt(set_mcmc, 4), "\n\n",
  "### Data & design\n\n", mdt(set_data, 4), "\n\n",
  "### Sampler switches\n\nEvery non-default argument the sampler was called with, including those ",
  "hardcoded in `estimate_prior.R`.\n\n", mdt(set_sw, 4), "\n\n") else "",
"## Convergence\n\n",
if (conv_bad) paste0("**The joint chain has not mixed.** log_lik Rhat is ",
  round(conv[block=="log_lik (joint)", rhat_max], 2),
  ". Per-parameter medians flatter a chain whose chains have each settled somewhere different. ",
  "Point estimates of abundant classes are usable; joint statements and every uncertainty interval are not.\n\n")
else "Per-block Rhat is within tolerance; intervals can be read as posterior intervals.\n\n",
mdt(conv, 3), "\n\n",
"## Held-out performance\n\n",
"| metric | value |\n| --- | --- |\n",
"| McFadden R2 | ", round(1 - ll/nl, 4), " |\n",
"| log-likelihood | ", format(round(ll), big.mark = ","), " |\n",
"| null log-likelihood | ", format(round(nl), big.mark = ","), " |\n",
"| held-out pixels | ", format(nrow(Yt), big.mark = ","), " |\n\n",
"## Calibration\n\nEvery (pixel, class) pair binned by predicted share. Mean absolute gap **",
round(calib_mae, 4), "**.\n\n", mdt(calib[, .(bin, n, pred, obs, gap)], 4), "\n\n",
"## What the country random effect buys\n\nSame held-out pixels, scored with pooled coefficients vs with country deviations.\n\n",
mdt(re_gain, 4), "\n\n",
"## Shrinkage\n\nPer-SD scale. ", format(nrow(m)*ncol(m), big.mark=","), " coefficients; ",
format(n_credible, big.mark=","), " with posterior sign-probability >= 0.9.\n\n", mdt(shrink, 3), "\n\n",
"## Which classes get confused\n\nStrongest negative PARTIAL error-correlation partner per class (raw correlation just recovers the sum-to-one constraint; `",
subs_dropped, "` dropped for invertibility).\n\n",
mdt(subs, 3), "\n\n",
"## Per-class skill\n\nNats gained against the null per class, and share of the total gain.\n\n",
mdt(per_skill[, .(class, obs, pred, ratio, gain_nats, share_of_gain)], 4), "\n\n",
"## Design summary\n\n", mdt(dsum), "\n\n",
"Class prevalence, largest first.\n\n", mdt(head(prev, 20), 5), "\n\n",
if (!is.null(spat)) paste0("## Spatial residuals\n\nMean total-variation distance per held-out pixel on a coarse grid: ",
  sprintf("min %.3f, median %.3f, max %.3f over %d bins.\n\n", min(spat$tvd), median(spat$tvd), max(spat$tvd), nrow(spat)))
else "",
"## Per-class fit\n\nratio = predicted / observed share.\n\n",
mdt(per[, .(class, obs, pred, ratio, dLL)], 5), "\n\n",
"## What drives the model\n\nPer-SD effect sizes, comparable across covariates.\n\n", mdt(head(drv, 20), 3), "\n\n",
"## Country effects (fixed + random)\n\n",
"RMS over classes of `b_g - mu`. Signed per-class values are in `beta_rotated`.\n\n",
mdt(cty, 4), "\n\n",
"Per-driver RE magnitude, top ", length(top_cov), " drivers by country variation.\n\n",
mdt(cc, 4), "\n\n",
"## Artifact for the downscaler\n\n",
if (!CAPRI_KEYED) paste0("**Group keys are NOT CAPRI codes.** This fit is grouped on `", GRPCOL,
  "`. The downscaler keys on CAPRI country codes, which are not ISO (BL, IR, EL, CS, KO, MO) and are ",
  "not translatable after the fact. Re-fit with `DRIVER_RE_GROUP_COL=CAPRI_NUTS` before using downstream.\n\n")
else "Group keys are CAPRI country codes, as the downscaler expects.\n\n",
"`beta_rotated.rds` / `beta_rotated.csv.gz` - ", format(nrow(rot), big.mark = ","),
" rows: covariate x from_class x to_class x country, as `value_median` with a 95% interval. ",
uniqueN(rot[pooled == TRUE, group]), " of ", uniqueN(rot$group),
" countries are not in the fit and carry the pooled fixed-effect contrast, flagged `pooled = TRUE`.\n")
writeLines(md, file.path(RUN, "report.md"))
cat(sprintf("\nwrote %s\n     %s\n     %s\n     %s\n", file.path(RUN, "report.html"),
    file.path(RUN, "report.md"),
    file.path(RUN, "beta_rotated.rds"), file.path(RUN, "beta_rotated.csv.gz")))
cat(sprintf("held-out McFadden %.4f | log_lik Rhat %.3f | rotated rows %s\n",
    1 - ll/nl, conv[block == "log_lik (joint)", rhat_max], format(nrow(rot), big.mark = ",")))
