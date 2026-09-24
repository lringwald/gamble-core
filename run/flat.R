# =============================================================================
# run/flat.R — build the DESIGN and/or fit the FLAT pixel land-use model
# =============================================================================
# HOW TO USE
#   1. Open gamble-core.Rproj (working directory = repo root).
#   2. Edit the CONTROL PANEL. 3. Source the file.
#
# TWO JOBS, one switch:
#   BUILD_DESIGN_ONLY = TRUE   assemble X/Y and STOP, writing output/designs/pixel_model_inputs.rds.
#                              This is the input run/nested.R needs. ~30 min.
#   BUILD_DESIGN_ONLY = FALSE  assemble the design AND fit the flat (single-level) MNL over all
#                              classes at once -- no nesting.
#
# CLASSIFICATION is the whole point of this script: the class column(s) in the LUM mapping decide
# what the model predicts. `GLOBIOM_subclass` also carries a curated `_nest` column, which is what
# gives run/nested.R its 6-nest tree.
# =============================================================================

## ============================ CONTROL PANEL ============================== ##
CLASSIFICATION <- "GLOBIOM_subclass"
  # "GLOBIOM_subclass"            27 classes + a curated 6-nest tree (Forests/Cropland/Pasture/
  #                               Urban/Natural/Waterbodies). The production taxonomy.
  # "GLOBIOM_UNFCCC,GLOBIOM_mngmt" legacy 19-class GLOBIOM (no curated tree -> prefix nesting)
  # "AgMIP_label"                 crop-type target (auto-enables the HRL crop split)
  # "BIOCLIMA_DS_reporting"       BIOCLIMA reporting classes

BUILD_DESIGN_ONLY <- TRUE   # TRUE = write design dump (output/designs/pixel_model_inputs.rds) and stop.
                            # FALSE = fit the flat MNL model (using existing design dump or building it).
REBUILD_DESIGN    <- FALSE  # When FALSE and fitting: reuse existing design dump if available.
DESIGN_PATH       <- "output/designs/pixel_model_inputs.rds" # path to design dump

PROMOTE_NATURAL_OTHER <- TRUE   # move Natural_other out of the non-choosable residual
MODEL_YEARS  <- "2018"      # target year(s) of the land-use map
FOCAL_YEARS  <- "2010"      # year the focal (neighbourhood composition) covariate is taken from
COV_YEARS    <- "2020"      # exogenous covariate year (spei48_2018 is renamed _2020 upstream)
MASTER_PARQUET <- "/Users/leopoldringwald/gamble_local_data/prior_model_1km_master_inputs.parquet"
  # Local copy of the 1 km covariate parquet. Set to "" to use the newest copy the driver resolves.

# --- MCMC Sampling Settings (used when BUILD_DESIGN_ONLY = FALSE) ---
NITER   <- 1000L        # total sweeps
NBURN   <- 250L         # discarded sweeps; MUST be < NITER
NCHAINS <- 1L           # parallel chains
THIN    <- 1L           # retain every k-th draw
SUBSAMPLE <- 0L         # 0 = all pixels; >0 subsamples for rapid testing

# --- BART on the terrain surface (used when BUILD_DESIGN_ONLY = FALSE) ---
USE_BART   <- FALSE     # TRUE = terrain enters as a tree ensemble instead of linear terms
BART_COLS  <- "topo"    # "topo" = Slope_rad, Elevation, Aspect_cos/sin_mean (+ lon/lat when ADD_COORDS)
ADD_COORDS <- FALSE     # TRUE = add lon/lat as design columns
COORD_TYPE <- "lonlat"  # "lonlat" = EPSG:4326 degrees | "laea" = native metres
## ========================================================================= ##

if (!file.exists("codes/mnl_aux_func.R"))
  stop("Working directory is not the repo root. Open gamble-core.Rproj, or setwd() to gamble-core/.")

if (!isTRUE(BUILD_DESIGN_ONLY) && NBURN >= NITER)
  stop(sprintf("NBURN (%d) must be < NITER (%d)", NBURN, NITER))

# Determine whether design dump needs to be assembled
need_build <- isTRUE(BUILD_DESIGN_ONLY) || isTRUE(REBUILD_DESIGN) || !file.exists(DESIGN_PATH)

if (need_build) {
  message(sprintf(">>> Assembling DESIGN dump: classification=%s | years %s (focal %s, cov %s)",
                  CLASSIFICATION, MODEL_YEARS, FOCAL_YEARS, COV_YEARS))
  vars <- c(DRIVER_CLASS_COLS = CLASSIFICATION,
            DRIVER_PROMOTE_NATURAL_OTHER = if (isTRUE(PROMOTE_NATURAL_OTHER)) "TRUE" else "FALSE",
            DRIVER_MODEL_YEARS = MODEL_YEARS, DRIVER_FOCAL_YEARS = FOCAL_YEARS,
            DRIVER_COV_YEARS = COV_YEARS,
            DRIVER_DUMP_INPUTS = "TRUE",
            DRIVER_DUMP_EXIT = "TRUE",
            DRIVER_DUMP_PATH = DESIGN_PATH)
  if (isTRUE(ADD_COORDS)) vars <- c(vars, DRIVER_ADD_COORDS = "TRUE", DRIVER_COORD_TYPE = COORD_TYPE)
  if (nzchar(MASTER_PARQUET) && file.exists(MASTER_PARQUET)) vars <- c(vars, GAMBLE_MASTER_PARQUET = MASTER_PARQUET)
  do.call(Sys.setenv, as.list(vars))
  source("drivers/run_lu_pixel_model.R", echo = FALSE)
  message(">>> Design dump built -> ", DESIGN_PATH)
}

if (!isTRUE(BUILD_DESIGN_ONLY)) {
  message(sprintf(">>> Fitting FLAT MNL model: %d sweeps (burn %d, thin %d) x %d chain(s) | BART %s",
                  NITER, NBURN, THIN, NCHAINS, if (isTRUE(USE_BART)) "ON" else "OFF"))
  fit_vars <- c(DESIGN_PATH = DESIGN_PATH,
                NITER = as.character(NITER),
                NBURN = as.character(NBURN),
                THIN = as.character(THIN),
                N_CHAINS = as.character(NCHAINS),
                SUBSAMPLE = as.character(SUBSAMPLE),
                USE_BART = if (isTRUE(USE_BART)) "TRUE" else "FALSE")
  do.call(Sys.setenv, as.list(fit_vars))
  source("drivers/run_flat_fit.R", echo = FALSE)
}

## ------------------------------- NEXT ------------------------------------- ##
# Design written to output/designs/pixel_model_inputs.rds -> now run run/nested.R
# Check what you built:
#   inp <- readRDS("output/designs/pixel_model_inputs.rds")
#   dim(inp$X_mat); colnames(inp$Y_pixel); inp$class_nest
