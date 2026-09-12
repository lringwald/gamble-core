# =============================================================================
# Organic-vs-conventional livestock split (δ_org) — spatial PATTERN fit, NUTS3.
# Mirrors fit_composition.R but a 2-category {organic, conventional} symmetric MNL on the organic-
# grassland proxy (off the Livestock_organic 1km map: n_org = Pasture_HIO+LIO, n_conv = Pasture_HI+LI).
# Produces dOrg = sum-to-zero organic contrast per driver. We learn the spatial PATTERN here; the national
# LEVEL is anchored to Eurostat BOV_org/SGT_org in build_subclass_parameters.R.
#   Env: D_NITER (20000), D_NBURN (10000), D_CORTHRESH (0.7), D_EXCL_ACCESS (TRUE), D_COVARIATES/D_PROTECT
# =============================================================================
suppressMessages({ library(data.table); library(future); library(future.apply); source("codes/mnlogit_rcpp_sym.R") })
m <- fread("output/composition/organic_training_nuts3.csv")
NITER <- as.integer(Sys.getenv("D_NITER","20000")); NBURN <- as.integer(Sys.getenv("D_NBURN","10000"))
NCHAINS <- as.integer(Sys.getenv("D_NCHAINS","4"))
CORTHRESH <- as.numeric(Sys.getenv("D_CORTHRESH","0.7"))

drivers0 <- setdiff(names(m), c("nuts3","country","n_org","n_conv"))
Xraw <- as.matrix(m[, drivers0, with=FALSE]); Xraw[!is.finite(Xraw)] <- 0
# areas are ALREADY log1p from prepare_composition_training (same construction as totals/grid X) -> gamma valid.
Xraw <- Xraw[, apply(Xraw, 2, sd) > 1e-8, drop=FALSE]; drivers0 <- colnames(Xraw)
# exclude the same human-footprint / raw-terrain cluster the totals + D/O/F fits drop (keep biophysical).
EXCL <- c(if (as.logical(Sys.getenv("D_EXCL_ACCESS","TRUE")))
            c("GHM_HI","CISI","RAI","GHM_HI_sd","CISI_sd","Slope_rad","Elevation","Aspect_cos_mean","Aspect_sin_mean")
          else character(0),
          trimws(strsplit(Sys.getenv("D_EXCLUDE",""), ",")[[1]]))
EXCL <- EXCL[nzchar(EXCL)]
Xraw <- Xraw[, setdiff(colnames(Xraw), EXCL), drop=FALSE]; drivers0 <- colnames(Xraw)

# --- DECORRELATE (identical helper/preferences as fit_composition.R) ---
# grassland columns named by pattern: GLOBIOM calls them Pasture_HI/LI, BMLEH Grassland_intensive/
# _extensive. A hardcoded name simply drops out of keep_pref under the other classification, which
# silently changes WHICH driver survives decorrelation rather than erroring.
keep_pref <- c("flat_share","upland_share","Slope_rad_sd",
               grep("^lu_area_(Pasture|Grassland)", names(m), value = TRUE),
               "lu_area_Cropland_HI","Growing_Degree_Days_gdd5","Annual_Precipitation_bio12","GDP","allPA_share")
decorr <- function(X, thresh, prefer, protect = character(0)) {
  C <- abs(cor(X)); C[is.na(C)] <- 0; diag(C) <- 0
  repeat {
    if (max(C) < thresh) break
    ij <- which(C == max(C), arr.ind=TRUE)[1, ]; i <- ij[1]; j <- ij[2]
    ni <- rownames(C)[i]; nj <- rownames(C)[j]; ip <- ni %in% protect; jp <- nj %in% protect
    if (ip && jp) { C[i,j] <- 0; C[j,i] <- 0; next }
    d <- if (mean(C[i,]) >= mean(C[j,])) i else j
    if (ip) d <- j else if (jp) d <- i
    else { if (ni %in% prefer && !(nj %in% prefer)) d <- j
           if (nj %in% prefer && !(ni %in% prefer)) d <- i }
    C <- C[-d, -d, drop=FALSE]
  }
  colnames(C)
}
.wl <- trimws(strsplit(Sys.getenv("D_COVARIATES",""), ",")[[1]]); .wl <- .wl[nzchar(.wl)]
if (length(.wl)) {
  sel <- intersect(.wl, colnames(Xraw))
  cat(sprintf("== ORGANIC | D_COVARIATES whitelist -> %d drivers ==\n  kept: %s\n", length(sel), paste(sel, collapse=", ")))
} else {
  PROT <- intersect(c(grep("^lu_area_", colnames(Xraw), value=TRUE),
                      trimws(strsplit(Sys.getenv("D_PROTECT", "GDP,Pop,allPA_area"), ",")[[1]])), colnames(Xraw))
  sel <- decorr(Xraw, CORTHRESH, keep_pref, protect = PROT)
  cat(sprintf("== ORGANIC | decorrelate %d -> %d drivers (max|cor|<%.2f; protected %d: all LU + GDP/Pop) ==\n  kept: %s\n",
              length(drivers0), length(sel), CORTHRESH, length(PROT), paste(sel, collapse=", ")))
}

X <- cbind(intercept = 1, Xraw[, sel, drop=FALSE]); k <- ncol(X)
Yr <- as.matrix(m[, .(org = n_org, conv = n_conv)]); Yr[Yr < 0] <- 0
Yr <- Yr * (1000 / median(rowSums(Yr))); keep <- rowSums(round(Yr)) > 0     # areas -> pseudo-counts (as fit_composition)
X <- X[keep, ]; Yr <- Yr[keep, ]; grp <- as.integer(factor(m$country[keep])); Y <- round(Yr)

plan(multisession, workers = NCHAINS)
cat(sprintf(">>> Running %d chains of mnlogit_rcpp_sym (%d iterations)...\n", NCHAINS, NITER))

res_list <- future_lapply(1:NCHAINS, function(cid) {
  source("codes/mnlogit_rcpp_sym.R")
  set.seed(cid)
  mnlogit_rcpp_sym(
    X = X, Y = Y, intercept = FALSE, baseline = 2, niter = NITER, nburn = NBURN,
    use_horseshoe = TRUE, horseshoe_idx = 2:k, symmetric_hs = TRUE, equation_specific_hs = FALSE,
    use_re = TRUE, group_idx = grp, re_idx = 1L,   # INTERCEPT ONLY -- see note below
    # re_idx = 1:k was tried and reverted 2026-08-21. This table is 189 rows / 26 countries
    # (median 6 rows per country, min 1) and the decorrelation screen PROTECTS all 18 lu_area_*
    # columns, so k stays ~22-30: re_idx = 1:k is ~1144-1560 RE parameters on 189 rows, i.e.
    # 6-8 PER OBSERVATION. The count model removed a 0.30-per-obs configuration for exactly this
    # reason -- there the RE variances ran to the slab cap and absorbed the fixed effects,
    # leaving 2 of 37 pooled coefficients credible. Shrinkage cannot rescue random slopes that
    # have no information behind them. Intercept-only is already 0.28/obs, at that same edge. use_car = FALSE, use_spike_slab = FALSE,
    re_prec_center = FALSE,
    # RE VARIANCE SHRINKAGE. The Finnish cap is ESTIMATED, not fixed at 100: full Bayes on this slab
    # measured better held-out and is self-calibrating (starts of 4 and 100 both converge to ~4.1, so
    # a fixed 100 only wastes burn-in sitting above the posterior).
    re_regularize = TRUE, estimate_slab_c2 = TRUE, collapse_slab_c2 = 4, slab_df_re = 10,
    disable_separation_detection = TRUE, standardize = TRUE, method = c("center","scale"), calc_loo = FALSE)
}, future.seed = TRUE)

pb <- do.call(abind::abind, c(lapply(res_list, function(f) f$postb_pooled), list(along = 3)))
nd <- dim(pb)[3]                          # [k, 2, draws]
mu   <- (pb[,1,] + pb[,2,]) / 2
dOrg <- pb[,1,] - mu                                              # sum-to-zero organic contrast (org=+dOrg, conv=-dOrg)
ci <- function(M) t(apply(M, 1, function(z) c(mean=mean(z), lo=quantile(z,.05), hi=quantile(z,.95))))
sl <- setdiff(colnames(X), "intercept")
cc <- ci(dOrg); rownames(cc) <- colnames(X); cc <- cc[sl, , drop=FALSE]; cc <- cc[order(-abs(cc[,"mean"])), ]
vis <- sign(cc[,2]) == sign(cc[,3])
cat("\n--- δ_ORGANIC (per-SD; ✓=CI excludes 0 [visible], ·=solid 0) ---\n")
for (r in 1:min(12, nrow(cc))) cat(sprintf("  %s %-28s %+.2f  [%+.2f, %+.2f]\n",
  ifelse(vis[r],"✓","·"), rownames(cc)[r], cc[r,1], cc[r,2], cc[r,3]))
cat(sprintf("  visible drivers: %d / %d\n", sum(vis), nrow(cc)))
fit <- res_list[[1]]
S <- Y[,1] / rowSums(Y); eta <- X %*% cbind(rowMeans(dOrg), -rowMeans(dOrg)); P <- exp(eta) / rowSums(exp(eta))
cat(sprintf("\nfit: cor(pred_org, obs_org)=%.2f | country RE-var=%.3f\n",
            cor(P[,1], S), if (!is.null(fit$post_sigma_re)) mean(fit$post_sigma_re) else NA))
saveRDS(list(sel = sel, dOrg = dOrg, drivers = colnames(X)), "output/composition/organic_fit.rds")
cat("saved -> output/composition/organic_fit.rds\n")
