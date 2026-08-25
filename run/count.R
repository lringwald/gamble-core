# =============================================================================
# run/count.R — fit the LIVESTOCK COUNT model (NB / Poisson), RStudio-friendly
# =============================================================================
# Open gamble-core.Rproj, edit the CONTROL PANEL, source the file.
#
# NOTE: one run fits BOTH livestock aggregates -- BOV (cattle) and SGT (sheep+goats).
# There is no species switch; they are estimated in the same pass.
# =============================================================================

## ============================ CONTROL PANEL ============================== ##
# Fix macOS Accelerate/OpenMP fork bug
Sys.setenv(OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)
RUN_MODE  <- "production"     # "test" = quick settings | "production" = full run
NITER     <- 40000L     # the count Gibbs mixes SLOWLY: dispersion r and the RE variance need long
                        # chains. For anything publication-grade use
                        # experiments/count/long_run_segmented.R, which is RESUMABLE -- long single
                        # jobs get killed on this machine.
YEARS     <- ""         # "" = driver default; else e.g. "2000,2010,2020"
MASTER_PARQUET <- "/Users/leopoldringwald/gamble_local_data/prior_model_1km_master_inputs.parquet"
                        # local copy; reading the parquet off Google Drive intermittently fails with
                        # `IOError ... [errno 60]`. "" = let the driver resolve it.
CLASS_COLS <- "GLOBIOM_UNFCCC,GLOBIOM_mngmt" # GLOBIOM settings
DOWNSCALE_GRID <- "" # "" runs the MCMC; "10km" builds the grid and stops
## ========================================================================= ##

if (!file.exists("run_ls_count_model.R"))
  stop("Working directory is not the repo root. Open gamble-core.Rproj, or setwd() to gamble-core/.")
if (!RUN_MODE %in% c("production", "test")) stop("RUN_MODE must be 'production' or 'test'")

vars <- c(DRIVER_RUN_MODE = RUN_MODE, DRIVER_NITER = as.character(NITER), 
          DRIVER_CLASS_COLS = CLASS_COLS, DRIVER_DOWNSCALE_GRID = DOWNSCALE_GRID,
          DRIVER_N_CHAINS = "4")
if (nzchar(YEARS))          vars <- c(vars, DRIVER_YEARS = YEARS)
if (nzchar(MASTER_PARQUET)) vars <- c(vars, GAMBLE_MASTER_PARQUET = MASTER_PARQUET)
do.call(Sys.setenv, as.list(vars))

message(sprintf(">>> count model (BOV + SGT) | mode=%s | niter=%d", RUN_MODE, NITER))
source("run_ls_count_model.R", echo = FALSE)
