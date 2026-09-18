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

BUILD_DESIGN_ONLY <- TRUE   # TRUE = write the design dump and stop (what run/nested.R consumes)
PROMOTE_NATURAL_OTHER <- TRUE   # move Natural_other out of the non-choosable residual
MODEL_YEARS  <- "2018"      # target year(s) of the land-use map
FOCAL_YEARS  <- "2010"      # year the focal (neighbourhood composition) covariate is taken from
COV_YEARS    <- "2020"      # exogenous covariate year (spei48_2018 is renamed _2020 upstream)
MASTER_PARQUET <- "/Users/leopoldringwald/gamble_local_data/prior_model_1km_master_inputs.parquet"
  # Local copy of the 1 km covariate parquet. Reading it straight off Google Drive intermittently
  # fails with `IOError ... [errno 60] Operation timed out` under arrow's parallel reads.
  # Set to "" to use the newest copy the driver can resolve itself.

# --- only used when BUILD_DESIGN_ONLY = FALSE ---
NITER   <- 1000L
NCHAINS <- 1L
## ========================================================================= ##

if (!file.exists("drivers/run_lu_pixel_model.R"))
  stop("Working directory is not the repo root. Open gamble-core.Rproj, or setwd() to gamble-core/.")

vars <- c(DRIVER_CLASS_COLS = CLASSIFICATION,
          DRIVER_PROMOTE_NATURAL_OTHER = if (isTRUE(PROMOTE_NATURAL_OTHER)) "TRUE" else "FALSE",
          DRIVER_MODEL_YEARS = MODEL_YEARS, DRIVER_FOCAL_YEARS = FOCAL_YEARS,
          DRIVER_COV_YEARS = COV_YEARS, DRIVER_NITER = as.character(NITER),
          DRIVER_NCHAINS = as.character(NCHAINS))
# BOTH dump flags are required: DRIVER_DUMP_EXIT alone is nested INSIDE the DUMP_INPUTS block, so on
# its own it neither dumps nor exits -- it silently runs a full flat fit instead.
if (isTRUE(BUILD_DESIGN_ONLY)) vars <- c(vars, DRIVER_DUMP_INPUTS = "TRUE", DRIVER_DUMP_EXIT = "TRUE")
if (nzchar(MASTER_PARQUET))    vars <- c(vars, GAMBLE_MASTER_PARQUET = MASTER_PARQUET)
do.call(Sys.setenv, as.list(vars))

message(sprintf(">>> classification=%s | %s | years %s (focal %s, cov %s)",
                CLASSIFICATION, if (isTRUE(BUILD_DESIGN_ONLY)) "DESIGN ONLY" else "FLAT FIT",
                MODEL_YEARS, FOCAL_YEARS, COV_YEARS))
source("drivers/run_lu_pixel_model.R", echo = FALSE)

## ------------------------------- NEXT ------------------------------------- ##
# Design written to output/designs/pixel_model_inputs.rds -> now run run/nested.R.
# Check what you built:
#   inp <- readRDS("output/designs/pixel_model_inputs.rds")
#   dim(inp$X_mat); colnames(inp$Y_pixel); inp$class_nest   # class_nest = the curated tree
