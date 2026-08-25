# =============================================================================
# Step 2: Bayesian driver-conditioned composition (δ) fit, via mnlogit_rcpp_sym.
# Multinomial D/O/F ~ drivers at NUTS2, horseshoe (symmetric) on the slopes, baseline = F.
# Output is baseline-relative coefficients -> centered to the sum-to-zero δ_j contrasts.
#   Env: D_SPECIES (bov|sgt), D_NITER (3000), D_NBURN (1500)
# =============================================================================
suppressMessages({ library(data.table); source("codes/mnlogit_rcpp_sym.R") })
SP <- Sys.getenv("D_SPECIES", "bov")
m  <- fread(file.path("output/composition", paste0(SP, "_training_nuts2.csv")))
NITER <- as.integer(Sys.getenv("D_NITER","3000")); NBURN <- as.integer(Sys.getenv("D_NBURN","1500"))

drivers <- setdiff(names(m), c("nuts2","nD","nO","nF","country"))
Yr <- as.matrix(m[, .(D=nD, O=nO, F=nF)]); Yr[Yr < 0] <- 0
# scale counts so the median region total ≈ 1000 (moderate multinomial size: preserves relative
# region weight, avoids tiny-region round-to-zero AND huge-count PG strain), then round + drop zeros.
Yr <- Yr * (1000 / median(rowSums(Yr)))
keep <- rowSums(round(Yr)) > 0
m <- m[keep]; Yr <- Yr[keep, , drop=FALSE]
X <- cbind(intercept = 1, as.matrix(m[, drivers, with=FALSE])); X[!is.finite(X)] <- 0
Y <- round(Yr)
k <- ncol(X)
cat(sprintf("== δ-fit %s | %d NUTS2 x %d cov x 3 cats (baseline=F) | %d/%d ==\n",
            toupper(SP), nrow(X), k, NITER, NBURN))

grp <- as.integer(factor(m$country))                       # country index for the random intercept
set.seed(1)
fit <- mnlogit_rcpp_sym(
  X = X, Y = Y, intercept = FALSE, baseline = 3,            # F (followers) = reference
  niter = NITER, nburn = NBURN,
  use_horseshoe = TRUE, horseshoe_idx = 2:k, symmetric_hs = TRUE, equation_specific_hs = FALSE,
  use_re = TRUE, group_idx = grp, re_idx = 1L,             # country random INTERCEPT (national baseline)
  use_car = FALSE, use_spike_slab = FALSE, disable_separation_detection = TRUE,
  standardize = TRUE, method = c("center","scale"), calc_loo = FALSE)

pb <- fit$postb_pooled                                      # [k, 3, draws] (col3 = baseline F ~ 0)
nd <- dim(pb)[3]
# δ posterior per draw: center each draw across {D,O,F} -> sum-to-zero contrasts
deltaD <- pb[, 1, ]; deltaO <- pb[, 2, ]; deltaF <- matrix(0, k, nd)
mu <- (deltaD + deltaO + deltaF) / 3
dD <- deltaD - mu; dO <- deltaO - mu; dF <- deltaF - mu     # each [k x draws]
mn <- function(M) rowMeans(M); sdf <- function(M) apply(M, 1, sd)
B  <- rbind(D = mn(dD), O = mn(dO), F = mn(dF)); colnames(B) <- colnames(X)
sl <- setdiff(colnames(X), "intercept")
# rank by the dairy-suckler contrast strength, report mean ± sd + P(δ_D>0)
pgt <- rowMeans(dD > 0)[match(sl, colnames(X))]
tab <- data.table(driver = sl, dD = mn(dD)[match(sl,colnames(X))], sd_dD = sdf(dD)[match(sl,colnames(X))],
                  P_dairy_up = pgt)
tab <- tab[order(-abs(dD))][1:12]
cat("\n--- δ_DAIRY (vs followers, sum-to-zero) with posterior uncertainty, top 12 ---\n")
print(tab[, .(driver, dD=round(dD,2), sd=round(sd_dD,2), P_dairy_up=round(P_dairy_up,2))])

# RE-aware fit: per region use its COUNTRY's total coefficients (FE + country RE)
S <- Y/rowSums(Y)
tt <- fit$postb_total                                       # [k, 3, n_groups, draws]
if (length(dim(tt)) == 4) {
  totm <- apply(tt, c(1,2,3), mean)                          # [k,3,groups]
  eta <- t(sapply(1:nrow(X), function(i) X[i,] %*% totm[, , grp[i]]))   # [n x 3]
} else eta <- X %*% t(B)
P <- exp(eta)/rowSums(exp(eta))
cat(sprintf("\nfit (RE-aware): cor(pred,obs) D=%.2f O=%.2f F=%.2f | RE-var(country)=%.3f\n",
            cor(P[,1],S[,1]), cor(P[,2],S[,2]), cor(P[,3],S[,3]),
            if(!is.null(fit$post_sigma_re)) mean(fit$post_sigma_re) else NA))
saveRDS(list(delta_mean = B, delta_sd = rbind(D=sdf(dD),O=sdf(dO),F=sdf(dF)),
             postb = pb, drivers = colnames(X)), file.path("output/composition", paste0(SP, "_delta_fit.rds")))
cat(sprintf("saved -> output/composition/%s_delta_fit.rds\n", SP))
