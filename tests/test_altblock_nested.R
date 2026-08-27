# ---------------------------------------------------------------------------
# NATIVE ALT-SPEC BLOCKS through nested_cut: does NCUT_ALT_BLOCKS build an n x K_node matrix
# aligned to the node's ALTERNATIVES, remove the source columns from the design, and reach the
# sampler? Covers BOTH node types via one resolver (.ncut_alt_fine_of): at a leaf an alternative is
# a fine class, at a nest it is the set of fine classes under that child.
# ---------------------------------------------------------------------------
suppressMessages(suppressWarnings({source("codes/mnlogit_rcpp_sym.R"); source("codes/nested_cut.R")}))
ok <- TRUE
chk <- function(nm, cond, det = "") { ok <<- ok && isTRUE(cond)
  cat(sprintf("  [%s] %-52s %s\n", if (isTRUE(cond)) "PASS" else "FAIL", nm, det)) }

cats <- c("Wheat", "Maize", "Forest", "Urban")
set.seed(5); n <- 1200
Fsh <- matrix(rgamma(n * 4, 1.2), n, 4); Fsh <- Fsh / rowSums(Fsh)
colnames(Fsh) <- paste0("focal_", cats)
X <- cbind(intercept = 1, x1 = rnorm(n), x2 = rnorm(n), Fsh)
cat("\n--- unit: block construction ---\n")
Sys.setenv(NCUT_ALT_BLOCKS = "spatial:focal_:shared")
sp <- .ncut_alt_specs()
chk("spec parsed", length(sp) == 1 && sp[[1]]$name == "spatial" && sp[[1]]$coef == "shared")

# LEAF: alternatives are the fine classes
ab_leaf <- .ncut_alt_blocks(X, cats, .ncut_alt_fine_of(cats), sp)
chk("leaf block is n x K", nrow(ab_leaf$blocks$spatial$Z) == n && ncol(ab_leaf$blocks$spatial$Z) == 4L)
chk("leaf block columns are the alternatives", identical(colnames(ab_leaf$blocks$spatial$Z), cats))
chk("leaf consumed the 4 focal cols", length(ab_leaf$drop) == 4L)
chk("leaf block rows sum to 1", max(abs(rowSums(ab_leaf$blocks$spatial$Z) - 1)) < 1e-10)

# NEST: alternatives are children; a child's column must be the SUM over its fine classes
node <- list(Crop = c("Wheat", "Maize"), Nature = c("Forest"), Built = c("Urban"))
ab_nest <- .ncut_alt_blocks(X, names(node), .ncut_alt_fine_of(node), sp)
Zn <- ab_nest$blocks$spatial$Z
chk("nest block is n x n_children", nrow(Zn) == n && ncol(Zn) == 3L)
chk("nest Crop col == focal_Wheat + focal_Maize (renormalised)",
    max(abs(Zn[, "Crop"] - (Fsh[, "focal_Wheat"] + Fsh[, "focal_Maize"]))) < 1e-10)

cat("\n--- integration: through nested_cut_fit (flat tree) ---\n")
b <- matrix(0, ncol(X), length(cats))
b[2, ] <- c(.8, -.5, .3, 0); b[3, ] <- c(-.4, .6, .1, 0)
U <- X %*% b + 1.1 * Fsh                       # true spatial-lag delta = 1.1
P <- exp(U - apply(U, 1, max)); P <- P / rowSums(P)
Y <- t(apply(P, 1, function(pr) rmultinom(1, 30, pr))); colnames(Y) <- cats
tree <- setNames(as.list(cats), cats)
fit <- suppressWarnings(suppressMessages(nested_cut_fit(
  X = X, Y = Y, tree = tree, use_iv = FALSE, use_re = FALSE,
  niter = 400, nburn = 200, thin = 1L, M = 1L)))
chk("nested_cut_fit runs with the block on", is.list(fit) && !is.null(fit$root),
    if (is.list(fit)) paste("top:", paste(names(fit), collapse=",")) else "not a list")
Sys.unsetenv("NCUT_ALT_BLOCKS")
cat(sprintf("\n  [%s] native alt-spec blocks wire through nested_cut\n", if (ok) "PASS" else "FAIL"))
if (!ok) quit(status = 1L)
