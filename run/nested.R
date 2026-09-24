# =============================================================================
# run/nested.R — fit the NESTED land-use model (RStudio-friendly)
# =============================================================================
# HOW TO USE
#   1. Open gamble-core.Rproj (sets the working directory to the repo root).
#   2. Edit the CONTROL PANEL below.
#   3. Source the file (Ctrl/Cmd-Shift-S). The fit lands in `fit`; the saved path
#      is printed at the end. INSPECT snippets are at the bottom.
#
# This is a thin wrapper over run_nested_cut.R -- the SAME code path as the shell
# runners, so an interactive fit and a batch fit are identical. Do not reimplement
# the fit here: run_nested_cut.R owns the curated-tree lookup, the resumable store
# and the provenance that gets written into the saved object.
# =============================================================================

## ============================ CONTROL PANEL ============================== ##
DESIGN     <- Sys.getenv("DESIGN_PATH", "output/designs/pixel_model_inputs.rds")
BRANCH     <- "GLOBIOM"        # "GLOBIOM" | "AGMIP" | "BIOCLIMA" -- only used when the dump carries
                               # no curated `class_nest` column; otherwise the curated tree wins.
VARIANT    <- "factorized"     # "factorized" = levels fit independently (VALIDATED CHOICE)
                               # "iv"         = inclusive-value coupling + lambda. Much slower (the
                               #                root is refit per imputation) and on the 27-class
                               #                GLOBIOM tree lambda fell outside (0,1] on 4 of 6 nests.
RE_BLOCK   <- "intercept"      # "intercept"       = country random intercept (best held-out + converges)
                               # "intercept+socio" = + 7 random slopes (ties on skill, RE Rhat 1.5 vs 1.1)
SYMMETRIC_HS <- TRUE           # STANDARD as of 2026-09-16 (BMLEH_Los1). Zero-sum/CLR horseshoe:
                               # shrinkage invariant to which class is the baseline -- which matters
                               # because `baseline` is which.max(colSums(Y)), i.e. data-dependent and
                               # DIFFERENT PER NODE.
                               # The old "~120 nats worse held-out" gate is VOID: those arms ran the
                               # pre-2026-08-19 bug, at Rhat 1.91. Set FALSE for a diagonal
                               # measurement arm. See docs/sampler_model_specification.md sec. 8.0.
SUBSAMPLE  <- 6000L            # 0 = all 64,173 pixels (production). 6000 ~ a 10-15 min smoke.
NITER      <- 200L             # production = 1000
M          <- 25L              # posterior draws carried per node (also IV imputations when VARIANT="iv")
N_CHAINS   <- 1L               # 4 -> per-node Rhat/ESS in fit$convergence, at ~4x the cost
STORE_DIR  <- NULL             # NULL = auto (tagged by variant/RE/HS). Kill-safe + resumable.
                               # !! A store fingerprints data/niter/M but NOT the sampler version --
                               # !! never point a run at a store built by an older sampler.
## ========================================================================= ##

if (!file.exists("codes/nested_cut.R"))
  stop("Working directory is not the repo root. Open gamble-core.Rproj, or setwd() to gamble-core/.")

if (!file.exists(DESIGN)) {
  cands <- c(list.files("output/designs", pattern = "^pixel_model_inputs.*\\.rds$", full.names = TRUE),
             list.files("output", pattern = "^pixel_model_inputs.*\\.rds$", full.names = TRUE))
  cands <- unique(cands[file.exists(cands)])
  if (length(cands) > 0) {
    cands <- cands[order(file.mtime(cands), decreasing = TRUE)]
    DESIGN <- cands[1]
    message(">>> Auto-detected newest design dump: ", DESIGN)
  } else {
    stop("Design dump not found: ", DESIGN, "\n  Build one with run/flat.R (BUILD_DESIGN_ONLY <- TRUE) or TASK=flat_design.")
  }
}

USE_IV <- switch(VARIANT, factorized = "FALSE", iv = "TRUE",
                 stop("VARIANT must be 'factorized' or 'iv'"))
if (!RE_BLOCK %in% c("intercept", "intercept+socio")) stop("RE_BLOCK must be 'intercept' or 'intercept+socio'")
if (isTRUE(SYMMETRIC_HS) && !any(grepl("kronecker(Msym", readLines("codes/mnlogit_rcpp_sym.R"), fixed = TRUE)))
  stop("SYMMETRIC_HS=TRUE but codes/mnlogit_rcpp_sym.R lacks the 2026-08-19 fix -- results would be invalid.")

if (is.null(STORE_DIR))
  STORE_DIR <- sprintf("output/ncut_store_%s_%s_%s%s", BRANCH, VARIANT,
                       gsub("[+]", "_", RE_BLOCK), if (isTRUE(SYMMETRIC_HS)) "_symHS" else "")

Sys.setenv(NCUT_USE_IV = USE_IV, NCUT_RE_COLS = RE_BLOCK, NCUT_SUB = as.character(SUBSAMPLE),
           NCUT_SYM_HS = if (isTRUE(SYMMETRIC_HS)) "TRUE" else "FALSE", NCUT_STORE_DIR = STORE_DIR)

# run_nested_cut.R reads NCUT_ARGS when it exists (see its `args <- ...` line):
#   BRANCH  design.rds  M  NITER  USE_RE  IVMODE  STREAM  N_CHAINS  N_CORES
NCUT_ARGS <- c(BRANCH, DESIGN, M, NITER, "TRUE", "auto", "TRUE", N_CHAINS,
               max(1L, min(N_CHAINS, parallel::detectCores() - 1L)))
message(sprintf(">>> %s | %s | RE=%s | symHS=%s | n=%s | niter=%d",
                BRANCH, VARIANT, RE_BLOCK, SYMMETRIC_HS,
                if (SUBSAMPLE == 0) "all" else SUBSAMPLE, NITER))
source("drivers/run_nested_cut.R", echo = FALSE)

## ------------------------------- INSPECT ---------------------------------- ##
# s <- summary_nested_cut(fit); names(s)                    # per-node coefficient tables
# s[["root/Forests"]]$coef                                  # one node's coefficients
# for (nm in names(s)) print(s[[nm]]$lambda)                # lambdas (VARIANT = "iv" only)
# nested_cut_identification(fit)                            # IV~design R^2 vs lambda (iv only)
# fit$convergence                                           # Rhat / ESS  (N_CHAINS > 1 only)
#
# Score it (in-sample) against any other fit on the same design:
#   system2("Rscript", c("postprocess/score_nested_cut.R", "output/nested_cut_<...>.rds"))
