# joint_fe_re_shrink: one kappa per covariate gating BOTH mu and b_g.
# x1 = signal in both blocks; x2 = pure noise in both (the covariate the gate should retire);
# x3 = signal ONLY in the RE block (must SURVIVE -- the gate must not kill group-varying effects).
suppressMessages({library(Rcpp); library(RcppArmadillo)})
sys.source("codes/mnlogit_rcpp_sym.R", environment())
set.seed(7)
n <- 2400; G <- 20; K <- 3; kx <- 4
X <- cbind(1, matrix(rnorm(n*(kx-1)), n, kx-1)); colnames(X) <- c("(Intercept)","x1","x2","x3")
g <- sample(G, n, TRUE)
mu_true <- rbind(c(0,0), c(1.2,-0.8), c(0,0), c(0,0))          # x2, x3 have NO pooled effect
sd_re   <- c(0.5, 0.6, 0.0, 0.9)                                # x2 has NO group variation; x3 has a lot
b <- array(0, c(kx, K-1, G)); for (v in 1:kx) b[v,,] <- rnorm((K-1)*G, 0, sd_re[v])
eta <- matrix(0, n, K-1)
for (j in 1:(K-1)) eta[,j] <- X %*% mu_true[,j] + rowSums(X * t(b[,j,g]))
eta <- cbind(0, eta); P <- exp(eta)/rowSums(exp(eta))
Y <- t(apply(P,1,function(p) rmultinom(1,1,p))); colnames(Y) <- paste0("c",1:K)

fit2 <- function(joint, sd) { set.seed(sd); capture.output(f <- mnlogit_rcpp_sym(
  Y=Y, X=X, group_idx=g, niter=700, nburn=300, use_re=TRUE, use_horseshoe=TRUE,
  re_idx=1:kx, joint_fe_re_shrink=joint)); f }
fit <- function(joint) fit2(joint, 99)
off <- fit(FALSE); on <- fit(TRUE)
# OFF at extra seeds: without a noise floor, ordinary MCMC drift reads as "the gate worked".
off2 <- fit2(FALSE, 100); off3 <- fit2(FALSE, 101)

sm <- function(f, v) { m <- apply(f$postb_pooled[v,,,drop=FALSE],1:2,mean); max(abs(m)) }
sr <- function(f, v) max(abs(f$sigma_beta_pooled[v,]))
cat("\n=== joint FE/RE gate: |mu| and RE sd by covariate ===\n")
cat(sprintf("%-6s %18s %18s | %18s %18s\n","cov","|mu| OFF","re_sd OFF","|mu| ON","re_sd ON"))
for (v in 2:kx) cat(sprintf("%-6s %18.4f %18.4f | %18.4f %18.4f\n",
  colnames(X)[v], sm(off,v), sr(off,v), sm(on,v), sr(on,v)))
if (!is.null(on$post_kappa)) { cat("\nkappa (posterior mean):\n")
  km <- rowMeans(on$post_kappa); for (v in 2:kx) cat(sprintf("  %-6s %.4f\n", colnames(X)[v], km[v])) }

ok1 <- sm(on,2) > 0.4                       # real pooled effect survives
# x2 is noise in BOTH blocks: |mu| must stay ~0 (it already is, so asserting a strict
# decrease there would test Monte Carlo noise) and its RE sd must be shrunk toward the true 0.
# x2 must be gated BELOW the seed-noise floor, not merely below one OFF run.
floor2 <- min(sr(off,3), sr(off2,3), sr(off3,3))
ok2 <- sm(on,3) < 0.05 && sr(on,3) < floor2
ok3 <- sr(on,4) > 0.3                       # RE-only signal survives the gate
cat(sprintf("\n[%s] x1 pooled effect survives the gate (|mu|=%.3f)\n", if(ok1)"PASS" else "FAIL", sm(on,2)))
cat(sprintf("[%s] x2 gated below the seed-noise floor: |mu|=%.4f ~ 0, re_sd ON=%.4f < min OFF=%.4f\n",
  if(ok2)"PASS" else "FAIL", sm(on,3), sr(on,3), floor2))
cat(sprintf("\nSEED-NOISE FLOOR  x1 OFF over 3 seeds: %.3f %.3f %.3f  (ON=%.3f)\n",
  sr(off,2), sr(off2,2), sr(off3,2), sr(on,2)))
if (!is.null(on$post_kappa)) cat(sprintf(
  "KAPPA IDENTIFICATION  min over draws = %.3f (a gate that truly retires a covariate reaches ~1e-2;\n",
  min(on$post_kappa[2:kx, ])),
  "                      kappa is aliased with sigma in the RE block, so it sits near 1 -- see docs.\n")
cat(sprintf("[%s] x3 RE-only signal survives (re_sd=%.3f, true 0.9)\n", if(ok3)"PASS" else "FAIL", sr(on,4)))
# aliasing probe: if kappa shrinks and sigma inflates to compensate, the RE block is unidentified
cat(sprintf("\nALIASING PROBE  sum(re_sd) OFF=%.3f  ON=%.3f  (blow-up => kappa/sigma aliased)\n",
  sum(sapply(2:kx, function(v) sr(off,v))), sum(sapply(2:kx, function(v) sr(on,v)))))
if (!all(ok1,ok2,ok3)) quit(status=1)
