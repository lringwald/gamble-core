# =============================================================================
# Driver-conditioned subtype composition (δ) — production fit.
# 1) DECORRELATE the driver matrix (greedy |cor| pruning, keeping interpretable drivers) so the
#    horseshoe can give VISIBLE, identifiable effects instead of compensating collinear ±coefs.
# 2) Fit multinomial D/O/F ~ curated drivers via mnlogit_rcpp_sym (horseshoe + country RE).
# 3) Report δ per subtype with posterior intervals: a driver is either visibly nonzero
#    (CI excludes 0) or a SOLID ZERO (tight CI around 0).
#   Env: D_SPECIES (bov|sgt), D_NITER (4000), D_NBURN (2000), D_CORTHRESH (0.7)
# =============================================================================
suppressMessages({ library(data.table); source("codes/mnlogit_rcpp_sym.R") })
SP <- Sys.getenv("D_SPECIES", "bov")
m  <- fread(file.path("output/composition", paste0(SP, "_training_nuts2.csv")))
NITER <- as.integer(Sys.getenv("D_NITER","4000")); NBURN <- as.integer(Sys.getenv("D_NBURN","2000"))
CORTHRESH <- as.numeric(Sys.getenv("D_CORTHRESH","0.7"))

drivers0 <- setdiff(names(m), c("nuts2","nD","nO","nF","country"))
Xraw <- as.matrix(m[, drivers0, with=FALSE]); Xraw[!is.finite(Xraw)] <- 0
# NB: areas/GDP/Pop are ALREADY log1p(sum) from prepare_composition_training (same construction as
# the totals X + grid X), so do NOT transform again here -- that keeps gamma = beta_total + delta valid.
# drop constant / near-constant drivers (no info, break correlation)
Xraw <- Xraw[, apply(Xraw, 2, sd) > 1e-8, drop=FALSE]; drivers0 <- colnames(Xraw)
# Optionally exclude the human-development/accessibility cluster: it is partly OUTCOME-confounded
# (developed lowlands ARE dairy) and separation-inflates (GHM δ ~ +15). Keep biophysical drivers.
# Exclude the same set as the totals/grid X (so gamma = beta_total + delta applies to the grid):
# the human-footprint cluster (GHM/CISI/RAI + _sd) AND the raw terrain means (Slope/Elevation/Aspect,
# now superseded by the terrain shares). GDP/Pop kept.
EXCL <- c(if (as.logical(Sys.getenv("D_EXCL_ACCESS","TRUE")))
            c("GHM_HI","CISI","RAI","GHM_HI_sd","CISI_sd","Slope_rad","Elevation","Aspect_cos_mean","Aspect_sin_mean")
          else character(0),
          trimws(strsplit(Sys.getenv("D_EXCLUDE",""), ",")[[1]]))                  # extra exclude-list
EXCL <- EXCL[nzchar(EXCL)]
Xraw <- Xraw[, setdiff(colnames(Xraw), EXCL), drop=FALSE]; drivers0 <- colnames(Xraw)

# --- DECORRELATE: greedy drop of the most-collinear var until max|cor| < thresh ---
# prefer to KEEP these interpretable drivers (drop their collinear partners instead)
keep_pref <- c("flat_share","upland_share","Slope_rad_sd","lu_area_Pasture_HI","lu_area_Pasture_LI",
               "lu_area_Cropland_HI","Growing_Degree_Days_gdd5","Annual_Precipitation_bio12",
               "GHM_HI","GDP","allPA_share")
decorr <- function(X, thresh, prefer, protect = character(0)) {
  C <- abs(cor(X)); C[is.na(C)] <- 0; diag(C) <- 0
  repeat {
    if (max(C) < thresh) break
    ij <- which(C == max(C), arr.ind=TRUE)[1, ]; i <- ij[1]; j <- ij[2]
    ni <- rownames(C)[i]; nj <- rownames(C)[j]
    ip <- ni %in% protect; jp <- nj %in% protect
    if (ip && jp) { C[i,j] <- 0; C[j,i] <- 0; next }              # both protected -> keep both
    d <- if (mean(C[i,]) >= mean(C[j,])) i else j                 # default: drop higher mean-cor
    if (ip) d <- j else if (jp) d <- i                            # never drop a protected var
    else { if (ni %in% prefer && !(nj %in% prefer)) d <- j        # else prefer-keep
           if (nj %in% prefer && !(ni %in% prefer)) d <- i }
    C <- C[-d, -d, drop=FALSE]
  }
  colnames(C)
}
# D_COVARIATES = comma-separated whitelist (use EXACTLY these, skip decorrelation); else decorrelate.
.wl <- trimws(strsplit(Sys.getenv("D_COVARIATES",""), ",")[[1]]); .wl <- .wl[nzchar(.wl)]
if (length(.wl)) {
  sel <- intersect(.wl, colnames(Xraw))
  cat(sprintf("== %s | D_COVARIATES whitelist -> %d drivers (decorrelation skipped) ==\n  kept: %s\n",
              toupper(SP), length(sel), paste(sel, collapse=", ")))
} else {
  # PROTECT from decorrelation pruning: ALL land-use classes + GDP + Pop (always kept), + D_PROTECT.
  PROT <- intersect(c(grep("^lu_area_", colnames(Xraw), value=TRUE),
                      trimws(strsplit(Sys.getenv("D_PROTECT", "GDP,Pop,allPA_area"), ",")[[1]])), colnames(Xraw))
  sel <- decorr(Xraw, CORTHRESH, keep_pref, protect = PROT)
  cat(sprintf("== %s | decorrelate %d -> %d drivers (max|cor|<%.2f; protected %d: all LU + GDP/Pop) ==\n  kept: %s\n",
              toupper(SP), length(drivers0), length(sel), CORTHRESH, length(PROT), paste(sel, collapse=", ")))
}

X <- cbind(intercept = 1, Xraw[, sel, drop=FALSE]); k <- ncol(X)
Yr <- as.matrix(m[, .(D=nD, O=nO, F=nF)]); Yr[Yr<0] <- 0
Yr <- Yr * (1000 / median(rowSums(Yr))); keep <- rowSums(round(Yr)) > 0
X <- X[keep, ]; Yr <- Yr[keep, ]; grp <- as.integer(factor(m$country[keep])); Y <- round(Yr)

set.seed(1)
fit <- mnlogit_rcpp_sym(
  X = X, Y = Y, intercept = FALSE, baseline = 3, niter = NITER, nburn = NBURN,
  use_horseshoe = TRUE, horseshoe_idx = 2:k, symmetric_hs = TRUE, equation_specific_hs = FALSE,
  use_re = TRUE, group_idx = grp, re_idx = 1L, use_car = FALSE, use_spike_slab = FALSE,
  disable_separation_detection = TRUE, standardize = TRUE, method = c("center","scale"), calc_loo = FALSE)

pb <- fit$postb_pooled; nd <- dim(pb)[3]                          # [k,3,draws]
mu <- (pb[,1,] + pb[,2,]) / 3
dD <- pb[,1,] - mu; dO <- pb[,2,] - mu; dF <- -mu                 # sum-to-zero δ per draw
ci <- function(M) t(apply(M, 1, function(z) c(mean=mean(z), lo=quantile(z,.05), hi=quantile(z,.95))))
sl <- setdiff(colnames(X), "intercept")
rep_sub <- function(M, lab) {
  cc <- ci(M); rownames(cc) <- colnames(X); cc <- cc[sl, , drop=FALSE]
  cc <- cc[order(-abs(cc[,"mean"])), ]
  vis <- sign(cc[,2]) == sign(cc[,3])                            # CI excludes 0 -> visible
  cat(sprintf("\n--- δ_%s (per-SD; ✓=CI excludes 0 [visible], ·=solid 0) ---\n", lab))
  for (r in 1:min(10,nrow(cc))) cat(sprintf("  %s %-28s %+.2f  [%+.2f, %+.2f]\n",
    ifelse(vis[r],"✓","·"), rownames(cc)[r], cc[r,1], cc[r,2], cc[r,3]))
  cat(sprintf("  visible drivers: %d / %d\n", sum(vis), nrow(cc)))
}
rep_sub(dD, "DAIRY"); rep_sub(dO, "MEAT")

# RE-aware fit
tt <- fit$postb_total; S <- Y/rowSums(Y)
if (length(dim(tt))==4) { tm <- apply(tt,c(1,2,3),mean); eta <- t(sapply(1:nrow(X), function(i) X[i,] %*% tm[,,grp[i]])) } else eta <- X %*% rbind(rowMeans(dD),rowMeans(dO),rowMeans(dF))
P <- exp(eta)/rowSums(exp(eta))
cat(sprintf("\nfit (RE-aware): cor(pred,obs) D=%.2f O=%.2f F=%.2f | country RE-var=%.3f\n",
            cor(P[,1],S[,1]),cor(P[,2],S[,2]),cor(P[,3],S[,3]), if(!is.null(fit$post_sigma_re)) mean(fit$post_sigma_re) else NA))
saveRDS(list(sel=sel, dD=dD, dO=dO, dF=dF, drivers=colnames(X)),
        file.path("output/composition", paste0(SP, "_composition_fit.rds")))
cat(sprintf("saved -> output/composition/%s_composition_fit.rds\n", SP))
