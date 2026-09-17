#!/usr/bin/env Rscript
# Reconstruction engine (slopes-exact + conditional intercepts) for one model branch.
# Produces a cached artifact with: coords, obs/pred shares, per-draw class + per-country
# totals (with CI), exact slope summaries, convergence, and fit metrics.
suppressMessages({library(qs2); library(data.table)})
`%||%` <- function(a,b) if(is.null(a)) b else a
source("codes/mnlogit_rcpp_sym.R")
args <- commandArgs(trailingOnly=TRUE)
BR <- args[1]                                   # AGMIP | GLOBIOM | BIOCLIMA
NDRAW_PER_CHAIN <- as.integer(args[2] %||% "60")
OUT <- "output/report_work/cache"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

RECON_MIN_GROUP_PIXELS <- 100L   # fallback group filter when metadata lacks group_levels
# discover the newest saved batch dir + matching design for the branch (no hardcoded paths)
.dir <- Sys.glob(sprintf("output/saved_model_outputs/*%s*", BR))
.dat <- Sys.glob(sprintf("output/dat_pixel_FULL_*%s*.rds", BR))
if (!length(.dir) || !length(.dat)) stop(sprintf("engine: no saved batches / dat_pixel_FULL for %s under output/.", BR))
cfg <- list(dir = .dir[which.max(file.mtime(.dir))], dat = .dat[which.max(file.mtime(.dat))],
            scheme = c(AGMIP="AgMIP crop-type", GLOBIOM="GLOBIOM land-use", BIOCLIMA="BIOCLIMA intensity")[BR] %||% BR,
            re = if (grepl("CAPRI", .dir[1])) "CAPRI_NUTS(country)" else "GLOB_country")
cat(sprintf(">>> [%s] batches: %s | design: %s\n", BR, basename(cfg$dir), basename(cfg$dat)))
meta <- qs_read(file.path(cfg$dir,"model_metadata.qs")); covs<-meta$cov_names; cats<-meta$cat_names
J <- length(cats); covx <- setdiff(covs,"intercept")
bpos <- which(cats==(meta$baseline_name %||% "")); if(length(bpos)!=1) bpos<-J; nb <- (1:J)[-bpos]
icr <- which(covs=="intercept"); sr <- setdiff(seq_len(length(covs)), icr)
# CLEAN-BATCH DETECTION: if the batches were written by the fixed sampler (correct intercept
# back-transform), use the RECOVERED intercepts directly -> genuine posterior totals + CIs.
# Otherwise (pre-fix corrupt intercept) fall back to conditional per-draw intercept refit.
use_true_int <- isTRUE(meta$batches_back_transformed)
intercept_mode <- if (use_true_int) "recovered (genuine posterior)" else "conditional refit (slopes exact, intercept re-solved)"
cat(sprintf("[%s] intercept mode: %s\n", BR, intercept_mode))

# ---- build fitted design (dat + filters + kept groups), exact fit rows/order ----
d <- as.data.table(readRDS(cfg$dat))
for(v in c("RAI","Pop","GDP")){lv<-paste0("log1p_",v); if(lv%in%covs && v%in%names(d)) d[[lv]]<-log1p(d[[v]])}
Ym<-as.matrix(d[,..cats]); d<-d[rowSums(Ym)>0]; Ym<-as.matrix(d[,..cats]); d<-d[rowSums(Ym)>=1]
# group set: metadata's group_levels if recorded (clean re-fits), else data-driven pixel-count filter.
# The recovered betas have a FIXED group dimension -> assert the reconstructed set reproduces it.
n_grp_fit <- dim(qs_read(list.files(cfg$dir, "posterior_batch_1_chain_1.qs", full.names=TRUE)[1])[[1]]$beta)[3]
keep_g <- if (!is.null(meta$group_levels)) meta$group_levels else { gt<-table(d$Grouping_Key); names(gt)[gt>=RECON_MIN_GROUP_PIXELS] }
d <- d[Grouping_Key %in% keep_g]
if (length(unique(d$Grouping_Key)) != n_grp_fit)
  stop(sprintf("engine: reconstructed group set (%d) != fitted batches' groups (%d) — the auto filter can't reproduce the fit's coverage. Supply the exact group list, or use a clean re-fit whose metadata carries group_levels.",
               length(unique(d$Grouping_Key)), n_grp_fit))
X <- cbind(intercept=1, as.matrix(d[,..covx])); X[!is.finite(X)]<-0
Yraw <- as.matrix(d[,..cats]); parea <- rowSums(Yraw); Y <- Yraw/parea
grp <- as.character(d$Grouping_Key); gf<-as.factor(grp); gidx<-as.integer(gf); glev<-levels(gf); uq<-unique(gidx)
coordX <- d$X; coordY <- d$Y
n <- nrow(X)
cat(sprintf("[%s] n=%d J=%d groups=%d baseline=%s\n", BR, n, J, length(uq), cats[bpos]))

lse <- function(U){m<-apply(U,1,max); m+log(rowSums(exp(U-m)))}
fit_int <- function(O, Yg, a0){
  a<-a0; nll<-function(a){U<-O; U[,nb]<-sweep(U[,nb,drop=FALSE],2,a,"+"); -sum(rowSums(Yg*U)-lse(U))}
  cur<-nll(a)
  for(it in 1:60){ U<-O; U[,nb]<-sweep(U[,nb,drop=FALSE],2,a,"+"); E<-exp(U-apply(U,1,max)); P<-E/rowSums(E); Pnb<-P[,nb,drop=FALSE]
    g<-colSums(Pnb-Yg[,nb,drop=FALSE]); H<-crossprod(Pnb); diag(H)<-diag(H)-colSums(Pnb); H<- -H
    step<-tryCatch(solve(H+diag(1e-6,length(nb)),g),error=function(e) g*1e-3)
    t<-1; repeat{an<-a-t*step; nn<-nll(an); if(nn<=cur-1e-9||t<1e-6)break; t<-t/2}; if(nn>=cur-1e-9)break; a<-an; cur<-nn }
  a
}

# precompute per-group row indices + Xg (slope part)
gi <- lapply(uq, function(g) which(gidx==g))
Xslp <- lapply(gi, function(rr) X[rr, sr, drop=FALSE])
Xall <- lapply(gi, function(rr) X[rr, , drop=FALSE])

sumP <- matrix(0, n, J)
tot_draws <- list()                                  # per-draw class totals [J]
totc_draws <- list()                                 # per-draw per-country totals [J x G]
pooled_slopes <- vector("list", 4)                   # exact slope draws for convergence/betas
sumBt <- NULL                                        # accumulation of total effects (FE + RE) [P x J x G]
nchain <- 4L; ndraw_used <- 0L
a_warm <- lapply(uq, function(x) rep(0, length(nb)))  # warm-start intercepts per group

t0 <- Sys.time()
for(ch in 1:nchain){
  r <- recover_mnlogit_posterior(cfg$dir, chain_id=ch)
  nd <- dim(r$postb_total)[4]
  keep <- unique(round(seq(1, nd, length.out=min(NDRAW_PER_CHAIN, nd))))
  pooled_slopes[[ch]] <- r$postb_pooled[, , keep, drop=FALSE]     # [P,J,K] exact slopes(+corrupt intercept row, dropped later)
  for(s in keep){
    Bt <- r$postb_total[,,,s]                                     # [P,J,G] (intercept row corrupt, slopes exact)
    if (is.null(sumBt)) sumBt <- Bt else sumBt <- sumBt + Bt
    Ps <- matrix(0, n, J); totc <- matrix(0, J, length(uq))
    for(gi2 in seq_along(uq)){
      B25 <- Bt[,,gi2]                                            # [P,J] already J cols (zero-sum expanded)
      if (use_true_int) {
        U <- Xall[[gi2]] %*% B25                                 # full coeffs incl correct intercept row
      } else {
        O <- Xslp[[gi2]] %*% B25[sr,,drop=FALSE]                 # slopes only
        a <- fit_int(O, Y[gi[[gi2]],,drop=FALSE], a_warm[[gi2]]); a_warm[[gi2]]<-a
        U <- O; U[,nb]<-sweep(U[,nb,drop=FALSE],2,a,"+")         # + conditionally re-solved intercepts
      }
      E<-exp(U-apply(U,1,max)); Pg<-E/rowSums(E)
      Ps[gi[[gi2]],] <- Pg
      totc[,gi2] <- colSums(Pg * parea[gi[[gi2]]])
    }
    sumP <- sumP + Ps; ndraw_used <- ndraw_used + 1L
    tot_draws[[length(tot_draws)+1]] <- rowSums(totc)
    totc_draws[[length(totc_draws)+1]] <- totc
  }
  rm(r); gc()
  cat(sprintf("  chain %d done (%.1f min elapsed)\n", ch, as.numeric(difftime(Sys.time(),t0,units="mins"))))
}
pointP <- sumP / ndraw_used
meanBt <- sumBt / ndraw_used  # [P, J, G] posterior mean of total effects (FE + RE)

# ---- exact slope summary (pooled symmetric), convergence on slopes ----
PS <- abind::abind(pooled_slopes, along=4)             # [P,J,K,chain]
dimnames(PS) <- list(covs, cats, NULL, NULL)
sl_rows <- sr                                          # slope covariate rows (exclude intercept)
# rhat/ess per slope param via posterior
flat <- matrix(aperm(PS[sl_rows,,,,drop=FALSE], c(1,2,3,4)), nrow=length(sl_rows)*J, ncol=dim(PS)[3]*dim(PS)[4])
rh <- apply(array(flat, c(nrow(flat), dim(PS)[3], dim(PS)[4])), 1, function(m) tryCatch(posterior::rhat(m), error=function(e) NA))
ess <- apply(array(flat, c(nrow(flat), dim(PS)[3], dim(PS)[4])), 1, function(m) tryCatch(posterior::ess_bulk(m), error=function(e) NA))
conv <- list(rhat=rh, ess=ess, frac_rhat_lt_1.01=mean(rh<1.01,na.rm=TRUE), frac_rhat_lt_1.1=mean(rh<1.1,na.rm=TRUE),
             median_ess=median(ess,na.rm=TRUE))
# rotated betas: pairwise category contrasts of pooled slopes, median + 95% CI
draws3 <- matrix(PS[,,,], nrow=dim(PS)[1]*dim(PS)[2])  # not used; compute contrasts directly
PSflat <- array(PS, c(dim(PS)[1], dim(PS)[2], dim(PS)[3]*dim(PS)[4]))  # [P,J,draws]
rot <- list()
for(v in sl_rows){
  med <- apply(PSflat[v,,],1,median); lo<-apply(PSflat[v,,],1,quantile,.025); hi<-apply(PSflat[v,,],1,quantile,.975)
  rot[[covs[v]]] <- data.frame(cov=covs[v], class=cats, median=med, q025=lo, q975=hi)
}
rot <- rbindlist(rot)
# DISPLAY SCALE, computed here because this is the only point in the report pipeline where the
# design X is in scope. render_heatplot.R and build_html.R read this cache and never see a design,
# so without it they can only plot RAW coefficients -- which span four orders of magnitude in this
# design and rank climate last when per SD it is among the largest effects. Rule B lives in
# codes/mnl_aux_func.R (display_sds): per SD for unbounded covariates, per unit for shares.
source("codes/mnl_aux_func.R")
sds_display <- display_sds(X)
rot[, sd_x := sds_display[cov]]
rot[is.na(sd_x), sd_x := 1]
rot[, `:=`(eff_median = median * sd_x, eff_q025 = q025 * sd_x, eff_q975 = q975 * sd_x)]

# ---- totals with CI (class + per-country) ----
TD <- do.call(rbind, tot_draws)                       # [S x J]
class_tot <- data.table(class=cats, actual=colSums(Yraw),
                        est_median=apply(TD,2,median), est_q025=apply(TD,2,quantile,.025), est_q975=apply(TD,2,quantile,.975))
# per country
TCc <- simplify2array(totc_draws)                     # [J x G x S]
ctry_tot <- rbindlist(lapply(seq_along(uq), function(gi2){
  actual_g <- colSums(Yraw[gi[[gi2]],,drop=FALSE])
  md <- apply(TCc[,gi2,],1,median); lo<-apply(TCc[,gi2,],1,quantile,.025); hi<-apply(TCc[,gi2,],1,quantile,.975)
  data.table(group=glev[uq[gi2]], class=cats, actual=actual_g, est_median=md, est_q025=lo, est_q975=hi)
}))

# ---- fit metrics ----
obs_share <- colMeans(Y); pred_share <- colMeans(pointP)
pixcorr <- sapply(seq_len(J), function(j) if(sd(pointP[,j])>0 && sd(Y[,j])>0) cor(pointP[,j],Y[,j]) else NA)
# stored trustworthy McFadden
b1 <- qs_read(list.files(cfg$dir,"posterior_batch_1_chain_1.qs",full.names=TRUE))
pll <- sapply(b1, function(x) x$log_lik)              # subset of draws' stored ll
pbar <- colMeans(Y); LLn <- n*sum(pbar[pbar>0]*log(pbar[pbar>0])); LLp <- sum(Yraw[Yraw>0]/parea[row(Yraw)][Yraw>0])   # approx
LLp <- sum(Y[Y>0]*log(Y[Y>0]))
McF_recon <- 1 - sum(Y*log(pmax(pointP,1e-12)))/LLn
fitm <- data.table(class=cats, obs_share=round(obs_share,5), pred_share=round(pred_share,5), pixel_corr=round(pixcorr,3))

saveRDS(list(branch=BR, cfg=cfg, meta=meta, n=n, J=J, cats=cats, baseline=cats[bpos],
             coordX=coordX, coordY=coordY, grp=grp, glev=glev[uq],
             obs=Y, pred=pointP, parea=parea,
             class_tot=class_tot, ctry_tot=ctry_tot, fitm=fitm, rot=rot, conv=conv,
             sds_display=sds_display, scale_note=display_scale_note(),
             McF_recon=McF_recon, LLn=LLn, LLp=LLp, ndraw_used=ndraw_used,
             use_true_int=use_true_int, intercept_mode=intercept_mode,
             meanBt=meanBt),
        file.path(OUT, paste0(BR, ".rds")))
cat(sprintf("[%s] DONE. draws used=%d | McFadden(recon)=%.3f | rhat<1.01: %.1f%% | median ESS(slopes)=%.0f\n",
    BR, ndraw_used, McF_recon, 100*conv$frac_rhat_lt_1.01, conv$median_ess))
print(fitm[order(-obs_share)][1:min(8,J)])
