# =============================================================================
# run/score.R — score fitted nested models on the design each one recorded
# =============================================================================
# Reports held-in log-likelihood, McFadden, MAD and per-class predicted/observed
# composition for one or more fits, side by side.
#
# THIS IS IN-SAMPLE. Every production fit trained on all rows, so these numbers
# measure reproduction, not generalisation, and they favour the richer random-
# effects block. For an honest RE comparison use experiments/mixing/re_idx_tradeoff.R,
# which holds out a country-stratified 25%.
# =============================================================================

## ============================ CONTROL PANEL ============================== ##
FITS <- c(
  # "output/nested_cut_GLOBIOM_factorized_2026-08-18_1039.rds"
)
if (!length(FITS)) FITS <- head(sort(Sys.glob("output/nested_cut_*.rds"), decreasing = TRUE), 3)
SCORE_N <- 0L     # 0 = all rows; a positive number subsamples for speed
SCORE_D <- 0L     # 0 = use every posterior draw the fit carries
## ========================================================================= ##

if (!file.exists("postprocess/score_nested_cut.R"))
  stop("Working directory is not the repo root. Open gamble-core.Rproj, or setwd() to gamble-core/.")
Sys.setenv(SCORE_N = as.character(SCORE_N), SCORE_D = as.character(SCORE_D))
message(">>> scoring:\n", paste0("    ", basename(FITS), collapse = "\n"))
system2("Rscript", c("postprocess/score_nested_cut.R", shQuote(FITS)))
