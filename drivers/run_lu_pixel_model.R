# =============================================================================
# run_prior_module_pixel_level_model.R
# =============================================================================
# Self-Contained Pixel Level Model (Multi-Year Panel)
# =============================================================================

rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(data.table)
library(arrow)
library(future)
library(future.apply)
require(progressr)
library(qs2)

handlers(global = TRUE)
handlers("cli")

# Helper function to find the most recently modified file matching a pattern
get_latest_file <- function(dir_path, pattern) {
  files <- list.files(path = dir_path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    return(NULL)
  }
  info <- file.info(files)
  latest_file <- rownames(info[order(info$mtime, decreasing = TRUE), ])[1]
  return(basename(latest_file))
}

# =========================================================================
# 1. PATHS & CONFIGURATION
# =========================================================================
# Upstream data locations (env-overridable so gamble-core can point anywhere). The raster data-prep
# reads the big 1km land-use maps / grid mapping / crop types from GRIDWORK_DIR (~5 GB, produced by the
# upstream LAMASUS_gridwork pipeline) — NOT from this repo's input/. Defaults assume the sibling projects
# sit next to this repo in the same parent folder. Set GAMBLE_* to relocate.
LS_DIR       <- Sys.getenv("GAMBLE_LS_DIR",       "/Users/leopoldringwald/Library/CloudStorage/OneDrive-IIASA/R/LAND_SUPPLY_ELASTICITY")
DS_DIR       <- Sys.getenv("GAMBLE_DS_DIR",       "../LAMASUS_downscaling/")
# CANONICAL DATA ROOT. Everything under ../LAMASUS_* is superseded: those are frozen older exports.
# The maintained copies live in cascadinggamble-core/data, and the differences are not cosmetic --
# LAMASUS_downscaling/input has no 2010 LUM map at all (so the count model silently used 2018 for the
# 2010 tier), and the gridwork grid mapping is the 2026-05-28 build against 2026-08-12 here, which
# differ on 3.6% of cells in NUTS2 -- the column the Eurostat crop shares join on. CAPRI_NUTS is
# identical between them, so RE keys are unaffected.
CASCADE_DATA <- Sys.getenv("GAMBLE_CASCADE_DATA", "../cascadinggamble-core/data")
GRIDWORK_DIR <- Sys.getenv("GAMBLE_GRIDWORK_DIR", file.path(CASCADE_DATA, "02_intermediate"))
AUXDATA_DIR  <- Sys.getenv("GAMBLE_AUXDATA_DIR",  file.path(CASCADE_DATA, "aux_files"))

# Dynamically pick the latest mapping file
mapping_file <- get_latest_file(AUXDATA_DIR, "^one_kmID_master_mapping_.*\\.parquet$")
if (is.null(mapping_file)) stop("No mapping file found in ", AUXDATA_DIR)
cat(sprintf("Using latest mapping file: %s\n", mapping_file))

PIXEL_RES <- as.integer(Sys.getenv("DRIVER_PIXEL_RES", "10"))  # grid resolution in km (5 or 10); EU 10km ~ 48-64k pixels, 5km ~ 180k. CANONICAL default 10.

# --- Spatial Intersection Configuration ---
# Set to a column name (e.g., "NUTS3") to define the observation unit as the intersection
# of the grid cell and this region (e.g., ID becomes "12345_AT111"). Set to NULL for native grid cells.
PIXEL_INTERSECT_COL <- { .piv <- Sys.getenv("DRIVER_PIXEL_INTERSECT", "NUTS3")  # flexible: "" / "NULL" => native grid cells (no intersect); CANONICAL default native. Set e.g. "CAPRI_NUTS"/"NUTS3" to split by region.
                         if (nzchar(.piv) && !toupper(.piv) %in% c("NULL","NONE")) .piv else NULL }

# --- Sampler selection ------------------------------------------------------
# Which posterior sampler to run. SAMPLER is woven into MODEL_LABEL so every output
# (the saved-model directory, posterior/dat_pixel files, figures) is tagged with it
# and different samplers never collide.
#   "mnlogit_rcpp"     : exact MNL/PG, baseline + post-hoc zero-sum projection
#   "mnlogit_rcpp_sym" : symmetric MNL/PG via the Msym coupling prior (no Helmert), baseline-free
#   "lnm_gibbs"        : logistic-normal multinomial (overdispersed)  [own post-processing]
#   "mvclr_gibbs"      : CLR Gaussian, symmetric by construction      [own post-processing]
# The MNL family shares this driver's full downstream pipeline (disk recovery ->
# zero-sum projection -> convergence -> heatplot). lnm_gibbs / mvclr_gibbs return
# different shapes and are not wired into that pipeline yet (see fit section).
SAMPLER <- Sys.getenv("DRIVER_SAMPLER", "mnlogit_rcpp_sym")
.sampler_allowed <- c("mnlogit_rcpp", "mnlogit_rcpp_sym", "lnm_gibbs", "mvclr_gibbs")
if (!SAMPLER %in% .sampler_allowed) {
  stop(sprintf("Unknown SAMPLER '%s'. Options: %s", SAMPLER, paste(.sampler_allowed, collapse = ", ")))
}
use_adaptive_sampler <- (SAMPLER == "mnlogit_rcpp_sym") # the symmetric/adaptive core

MODEL_LABEL <- paste0(PIXEL_RES, "km_", SAMPLER, if (!is.null(PIXEL_INTERSECT_COL)) paste0("_by_", PIXEL_INTERSECT_COL))

# Model Settings
N_CHAINS <- as.integer(Sys.getenv("DRIVER_NCHAINS", "4"))
use_bart <- isTRUE(as.logical(Sys.getenv("DRIVER_USE_BART", "FALSE")))  # validated non-focal BART design (topo/climate/soil/socioecon -> BART, focal lags linear); OFF by default
n_trees_bart <- as.integer(Sys.getenv("DRIVER_NTREES_BART", "25"))
niter <- as.integer(Sys.getenv("DRIVER_NITER", "16000"))   # 16000-4000 = 12000 retained: the validated regime where the RE-variance cells settle (median RE-var ESS ~173, vs ~9 at 1000 kept). Override via env for a quick test.
nburn <- as.integer(Sys.getenv("DRIVER_NBURN", "4000"))
# thin: store every k-th of the 12000 retained draws -> bounds the POST-FIT recover/combine RAM (the big
# [k x J x G x nretain] posterior recovered from disk) + shrinks the disk footprint. k<<autocorr-time keeps
# ESS (RE-var ESS ~173 at 12000 -> autocorr-time ~70; thin=4 -> 3000 stored loses ~no ESS). Keeps iterations.
thin_keep <- max(1L, as.integer(Sys.getenv("DRIVER_THIN", "4")))
DO_POSTERIOR_ROTATION <- TRUE

# Adaptive Sampler Settings (apply when SAMPLER == "mnlogit_rcpp_sym")
adaptive_use_wls_init <- TRUE
adaptive_use_precision_hs <- TRUE
adaptive_use_spike_slab <- FALSE # RE-inclusion spike-slab RETIRED: experiments/mixing/exp3 shows the deterministic screen (A1/A2/A3) + continuous support prior match its recovery and mix ~1.3-1.5x better (its delta fires only ~3% once the separation penalty is on -> de-facto deterministic). Reversible (set TRUE to restore). VALIDATE LOCALLY: confirm recovery on the real multinomial, esp. PARTIAL-separation cells (exp3 was binary/clean separation; borderline cells were the delta's only real niche, now covered continuously by the support prior).
adaptive_store_delta <- TRUE
adaptive_use_car <- FALSE # Set to TRUE if country_adjacency is provided
adaptive_car_rho <- 0.95

# Separation handling (see memory: separation-handling-design)
adaptive_separation_as_prior <- TRUE # severity-scaled inclusion prior on flagged cells
adaptive_separation_soft <- FALSE # FALSE = observe-only (logs n_soft); TRUE = activate soft tier
adaptive_sep_overlap_hard <- 0.0 # support overlap <= this => HARD mask (severe separation)
adaptive_sep_overlap_neutral <- 0.5 # support overlap >= this => healthy cell (pinned delta=1)

# Parallel Settings (Windows-compatible multisession)
N_CORES <- min(parallel::detectCores() - 1, N_CHAINS)
plan(multisession, workers = N_CORES)
# Increase global size limit (e.g., 10GB) for massive pixel-level datasets
options(future.globals.maxSize = 10 * 1024^3)
cat(sprintf("Parallel Plan: %s with %d workers (Limit: 10GB)\n", "multisession", N_CORES))

# Stabilization Settings
STABILIZE_SPARSE_Y <- FALSE
ADAPTIVE_EPSILON_MAX <- 1e-4 # Max 1e-4% mass added to empty classes

# Hierarchical / RE Settings
use_re <- TRUE
RE_GROUP_COL <- Sys.getenv("DRIVER_RE_GROUP_COL", "GLOB_country") # "GLOB_country", "NUTS2", "CAPRI_NUTS", ...
# Run versioning: bucket ALL outputs (saved-model dir, dat_pixel, fits, plots) into a
# stable "production" or a throwaway "test" namespace via MODEL_LABEL — replaces the old
# manual REvNN bumping. A "test" run starts fresh (its saved-model dir is cleared at the
# start) so it never hot-starts or recovers from stale state; "production" persists, so a
# re-run resumes/extends the canonical result.
RUN_MODE <- Sys.getenv("DRIVER_RUN_MODE", "production") # "production" | "test"
if (!RUN_MODE %in% c("production", "test")) stop("RUN_MODE must be 'production' or 'test'")
# Promote a finished test run to the production bucket: at the end of the run, copy
# every artifact tagged with this test MODEL_LABEL (saved-model dir + dat_pixel + fits
# + plots) to its "_production" equivalent, overwriting production. Only acts in test mode.
PROMOTE_TO_PRODUCTION <- isTRUE(as.logical(Sys.getenv("DRIVER_PROMOTE", "FALSE")))

# Estimation Scale
ESTIMATE_ON_LEVELS <- FALSE # If TRUE, uses absolute area (km2/ha) scale instead of normalized shares
LEVEL_SCALE_FACTOR <- 1 # e.g., 1 for km2, 100 for ha (if input is km2)
# Note: Higher counts lead to tighter posterior distributions (higher precision).

# --- Year Configuration ---
# Panel tiers. MODEL_YEARS = the OUTCOME (Y) year per tier; COV_YEARS = the matching
# exogenous-covariate year; FOCAL_YEARS = the composition-source year for the focal (X_input)
# lag (see below). All three run PARALLEL (one entry per tier) and are env-overridable, so the
# Y years can be chosen without editing the file, e.g. for the GLOBIOM 2010+2018 panel:
#   DRIVER_MODEL_YEARS=2010,2018  DRIVER_COV_YEARS=2010,2020  DRIVER_FOCAL_YEARS=2000,2010
# Default here = the single BIOCLIMA 2018 tier with a 2010 (t-1) focal lag.
.years_env <- function(var, default_vec) {
  v <- Sys.getenv(var, "")
  if (!nzchar(v)) return(default_vec)
  as.integer(trimws(strsplit(v, ",")[[1]]))
}
MODEL_YEARS <- .years_env("DRIVER_MODEL_YEARS", c(2018))
COV_YEARS   <- .years_env("DRIVER_COV_YEARS",   c(2020))  # spei48_2018 is renamed to _2020 upstream

# --- Focal (X_input) lag configuration -------------------------------------
# The focal neighborhood LU-composition (the `focal_*` autoregressive "LU-lag" predictors,
# i.e. the model's X_input) can be drawn from an EARLIER year than the outcome, giving a
# true t-1 autoregressive lag instead of the contemporaneous (endogenous) same-year map.
# FOCAL_YEARS[i] is the composition-source year for tier i's focal input. Set it == MODEL_YEARS[i]
# (or NA) for legacy contemporaneous behavior. Any lagged year must be in LUM_SOURCE_REGISTRY.
# NOTE: the AGMIP scheme overrides this default to CONTEMPORANEOUS (focal_year = out_year) below,
# unless DRIVER_FOCAL_YEARS was set explicitly — HRL crop types exist only for 2018, so a lag focal
# would carry no crop-type classes (it would fall back to generic Cropland_*_other).
# --- OWN t-1 STATE (`prev_*`) ------------------------------------------------
# DRIVER_PREV_STATE=TRUE additionally emits the pixel's OWN composition at focal_year as `prev_<class>`
# share columns. The focal is a NEIGHBOURHOOD statistic (compute_focal_coord excludes the centre), so
# the pixel's own predecessor state is NOT otherwise in the design. It is the transition-agnostic
# ("net") half of the transition ladder: with a shared coefficient it is scalar inertia, per-class it
# is per-class inertia. Route it into a conditional-logit block via
# NCUT_ALT_BLOCKS="temporal:prev_:shared" -- the share of class j at t-1 is an attribute OF
# alternative j, not of the pixel.
# HARD ERROR if the focal year equals the outcome year: `prev_*` would then BE the outcome and the fit
# is circular. Not a warning, because AGMIP deliberately defaults to a contemporaneous focal (HRL crop
# types exist only for 2018), so the unsafe case is the DEFAULT on that branch.
PREV_STATE <- isTRUE(as.logical(Sys.getenv("DRIVER_PREV_STATE", "FALSE")))
.focal_years_explicit <- nzchar(Sys.getenv("DRIVER_FOCAL_YEARS", ""))
FOCAL_YEARS <- .years_env("DRIVER_FOCAL_YEARS", c(2010))
if (length(FOCAL_YEARS) != length(MODEL_YEARS) || length(COV_YEARS) != length(MODEL_YEARS)) {
  stop("MODEL_YEARS, COV_YEARS and FOCAL_YEARS must all have the same length (one entry per tier).")
}
# NOTE: LUM_SOURCE_REGISTRY (the per-year granular-LUM source files) is lineage-aware and is
# defined below, once CLASS_SCHEME is known (BIOCLIMA vs GLOBIOM data lineages).

# --- Source Functions ---
source("codes/mnl_aux_func.R")
source("postprocess/MNL_parameter_heatplot.R")
source("codes/mnlogit_rcpp.R")
source("codes/spatial_utils.R")

# =========================================================================
# 2. SHARED DATA: CLASS MAPPING
# =========================================================================
cat("\nLoading Class Mapping...\n")

# ---- Classification configuration (single source of truth) ----------------
# Edit these to reshape the outcome classes; the map->aggregate step adapts.
OUTCOME_SOURCE <- "LUM" # "LUM" | "CLC"  -> selects raw outcome + mapping + join key
# Mapping column(s) whose combined (left->right) value defines a model class.
#   LUM mapping offers: GLOBIOM_UNFCCC, GLOBIOM_mngmt, BIOCLIMA3_LU, BIOCLIMA3_mngmt
#   CLC mapping offers: GLOBIOM_UNFCCC, CLC_level_1/2/3
# e.g. c("GLOBIOM_UNFCCC"); c("GLOBIOM_UNFCCC","GLOBIOM_mngmt"); c("BIOCLIMA3_LU","BIOCLIMA3_mngmt")
CLASS_COLS <-  strsplit(Sys.getenv("DRIVER_CLASS_COLS", "AgMIP_label"), ",")[[1]]  # DEFAULT = latest BIOCLIMA target (single combined LU_mngmt column); flexible: e.g. DRIVER_CLASS_COLS="GLOBIOM_UNFCCC,GLOBIOM_mngmt"
# Tag MODEL_LABEL with the classification family (BIOCLIMA / GLOBIOM / CLC), derived
# from the first CLASS_COL's prefix (trailing version digits stripped), so runs on
# different target classifications never collide. MODEL_LABEL isn't used before here.
CLASS_SCHEME <- toupper(sub("[0-9]+$", "", sub("_.*$", "", CLASS_COLS[1])))
MODEL_LABEL <- paste0(MODEL_LABEL, "_", CLASS_SCHEME, "_", RUN_MODE)
# AGMIP focal default = CONTEMPORANEOUS (focal_year = out_year). The HRL crop types only exist for
# 2018, so a t-1 lag focal would carry no crop-type classes; a same-year focal keeps the focal_*
# neighborhood composition crop-typed and consistent with the outcome. Honoured unless the user set
# DRIVER_FOCAL_YEARS explicitly (then their choice stands, e.g. a deliberate generic-cropland lag).
# BMLEH is in the same position as AGMIP and for the same reason -- it is crop-typed from the same
# 2018-only HRL product -- so it must take the same default. Left out, it fell to the 2010 lag and
# silently lost every crop focal column: focal_softwheat, focal_barley, focal_maize and the rest
# vanished while the residual-derived ones (focal_tobacco, focal_other_crop) survived, because only
# the generic cropland reached the cascade. In the AgMIP fit focal_softwheat was among the strongest
# predictors in the model (RMS 1.84).
if (CLASS_SCHEME %in% c("AGMIP", "BMLEH") && !.focal_years_explicit) {
  FOCAL_YEARS <- MODEL_YEARS
  cat(sprintf(">>> %s: focal defaulted to CONTEMPORANEOUS (focal_year = out_year); set DRIVER_FOCAL_YEARS to override.\n", CLASS_SCHEME))
}
# Tag MODEL_LABEL with the focal (X_input) lag when any tier uses an EARLIER composition
# year than its outcome, so lagged artifacts never collide with contemporaneous ones. Gap
# is the max out_year - focal_year across lagged tiers (e.g. 2018 outcome, 2010 focal -> _focalLag8y).
.focal_lag_gaps <- mapply(function(oy, fy) if (!is.na(oy) && !is.na(fy) && fy != oy) oy - fy else 0L,
                          MODEL_YEARS, FOCAL_YEARS)
if (any(.focal_lag_gaps != 0L)) MODEL_LABEL <- paste0(MODEL_LABEL, "_focalLag", max(.focal_lag_gaps), "y")

# --- Granular LUM source registry (lineage-aware; single source of truth per year) ---------
# Maps a year -> the granular 1km LUM file + id/code column names, used by BOTH the outcome
# loader and the (possibly lagged) focal loader. Two data lineages exist for the SAME years:
#   bioclima : the BIOCLIMA/EWM product (ibiom_* / *_wetlands_enhanced)  -> for BIOCLIMA targets
#   globiom  : the GLOBIOM downscaling fit (LUM_fit_with_energy_levels_*) -> for GLOBIOM targets
# All files share the layout INSPIRE_Europe_buffer_1kmID + LUM_code_clc (== mapping LUM_Code, verified;
# 100% grid linkage). The lineage is auto-selected from CLASS_SCHEME but can be forced via
# DRIVER_LUM_LINEAGE. Add a year here to make it available as an outcome or a focal-lag source.
.LUM_REGISTRIES <- list(
  bioclima = list(
    `2018` = list(file = file.path(GRIDWORK_DIR, "ibiom_LUM_EWM_2018.rds"),       id_col = "INSPIRE_Europe_buffer_1kmID", code_col = "LUM_code_clc"),
    `2010` = list(file = file.path(GRIDWORK_DIR, "LUM_2010_wetlands_enhanced.rds"), id_col = "INSPIRE_Europe_buffer_1kmID", code_col = "LUM_code_clc")
  ),
  globiom = list(
    # keyed EEA_1kmID with a per-year code column, which is how the canonical export is shaped. The
    # two 1km keys are exactly 1:1 (6,387,084 each way, max cardinality 1), so the swap moves no area.
    `2000` = list(file = file.path(GRIDWORK_DIR, "LUM_fit_with_energy_levels_and_new_FM_2000_EEA_1kmID.rds"), id_col = "EEA_1kmID", code_col = "LUM_fit_with_energy_levels_and_new_FM_2000"),
    `2010` = list(file = file.path(GRIDWORK_DIR, "LUM_fit_with_energy_levels_and_new_FM_2010_EEA_1kmID.rds"), id_col = "EEA_1kmID", code_col = "LUM_fit_with_energy_levels_and_new_FM_2010"),
    `2018` = list(file = file.path(GRIDWORK_DIR, "LUM_fit_with_energy_levels_and_new_FM_2018_EEA_1kmID.rds"), id_col = "EEA_1kmID", code_col = "LUM_fit_with_energy_levels_and_new_FM_2018")
  )
)
LUM_DATA_LINEAGE <- { .lin <- tolower(Sys.getenv("DRIVER_LUM_LINEAGE", ""))
  # GLOBIOM and AGMIP both take their land-use AREA basis from the GLOBIOM downscaling fit
  # (LUM_fit_*, keyed LUM_code_clc); AgMIP only re-labels those codes + splits cropland by crop type.
  # BMLEH belongs here too: BMLEH_Los1_label is a relabel of the SAME LUM codes AgMIP uses (plus the
  # short-rotation and grassland-intensity splits), so it takes the same GLOBIOM area basis. Without
  # it the scheme fell through to bioclima and looked for a file that does not exist.
  if (!nzchar(.lin)) .lin <- if (CLASS_SCHEME %in% c("GLOBIOM", "AGMIP", "BMLEH")) "globiom" else "bioclima"
  if (!.lin %in% names(.LUM_REGISTRIES)) stop(sprintf("Unknown DRIVER_LUM_LINEAGE '%s' (options: %s).", .lin, paste(names(.LUM_REGISTRIES), collapse = ", ")))
  .lin }
LUM_SOURCE_REGISTRY <- .LUM_REGISTRIES[[LUM_DATA_LINEAGE]]
cat(sprintf("LUM data lineage: %s  (registered years: %s)\n", LUM_DATA_LINEAGE, paste(names(LUM_SOURCE_REGISTRY), collapse = ", ")))

# --- Crop-type split (AgMIP target): replace LUM cropland with HRL Crop Types -----------------
# When on, load_granular_lum splits the LUM arable/permanent cropland AREA into specific crop-type
# classes (wheat, maize, grapes, ...) using the Copernicus HRL Crop Types product (per-1km crop
# composition), area-conserving to LUM. The LUM cropland area with no HRL crop overlap stays in the
# residual classes Cropland_arable_other / Cropland_permanent_other (both defined in AgMIP_label).
# Auto-on for the AGMIP scheme; env-overridable via DRIVER_CROP_SPLIT.
# BMLEH is crop-typed from the same HRL product, so it defaults on too. It was relying on the caller
# passing DRIVER_CROP_SPLIT=TRUE; forgetting that would have produced a design with no crop classes
# at all and no error -- the same silent-degradation shape as the focal-year default.
DO_CROP_SPLIT <- as.logical(Sys.getenv("DRIVER_CROP_SPLIT", if (CLASS_SCHEME %in% c("AGMIP", "BMLEH")) "TRUE" else "FALSE"))
# HRL Crop Types source per year (same 1km key as LUM). HRL is a 2017-2019 average -> maps to 2018.
# YEAR-MATCHED by default, not the 2017-2019 average. The average is exactly sum/3 with an unmapped
# year counted as ZERO, and HRL's coverage varies enormously by year: BE/IE/LU/NL/CH appear only in
# 2018 (so the average keeps a third of their crop area), DE and FR have 2017/2019 totals around a
# sixth of their 2018 one (~55% lost), while stably-mapped countries (ES/BG/CZ/SE/DK) lose 0-5%.
# EU crop area is 552k/770k/412k km2 by year against 578k for the average -- a 47% swing between
# years is coverage, not rotation, since rotation moves area BETWEEN crops and leaves the national
# total flat. Averaging therefore shrinks H, and every km2 of LUM cropland that H fails to account
# for falls through to Cropland_arable_other and gets split across the residual crops (~75% fodder).
# Set DRIVER_CROP_AVG=TRUE to restore the old averaged source.
CROP_SOURCE_REGISTRY <- if (isTRUE(as.logical(Sys.getenv("DRIVER_CROP_AVG", "FALSE")))) list(
  `2018` = list(file = file.path(GRIDWORK_DIR, "Crop_Types/Crop_Types_Avg_2017_2019_1km.rds"),
                id_col = "INSPIRE_Europe_buffer_1kmID", code_col = "Crop_Type_Code", area_col = "avg_area_km2")
) else list(
  `2017` = list(file = file.path(GRIDWORK_DIR, "Crop_Types/Crop_Types_2017_1km.rds"),
                id_col = "INSPIRE_Europe_buffer_1kmID", code_col = "Crop_Type_Code", area_col = "area_km2"),
  `2018` = list(file = file.path(GRIDWORK_DIR, "Crop_Types/Crop_Types_2018_1km.rds"),
                id_col = "INSPIRE_Europe_buffer_1kmID", code_col = "Crop_Type_Code", area_col = "area_km2"),
  `2019` = list(file = file.path(GRIDWORK_DIR, "Crop_Types/Crop_Types_2019_1km.rds"),
                id_col = "INSPIRE_Europe_buffer_1kmID", code_col = "Crop_Type_Code", area_col = "area_km2")
)
# HRL Crop_Type_Code -> AgMIP crop class name + cropland group (arable|permanent). 3100/3200
# ("Unclassified arable/permanent") and 65535 ("Outside area") are NOT listed -> they fold into the
# matching Cropland_*_other residual (no HRL overlap identified as a specific crop).
CROP_LEGEND <- data.table::data.table(
  Crop_Type_Code = c(1110,1120,1130,1140,1150,1210,1220,1310,1320,1410,1420,1430,1440,  2100,2200,2310,2320),
  crop_class     = c("Wheat","Barley","Maize","Rice","Other_cereals","Fresh_vegetables","Dry_pulses",
                     "Potatoes","Sugar_beet","Sunflower","Soybeans","Rapeseed","Flax_cotton_hemp",
                     "Grapes","Olives","Fruits","Nuts"),
  crop_group     = c(rep("arable", 13), rep("permanent", 4))
)

INCLUDE_IRRIGATION <- as.logical(Sys.getenv("DRIVER_INCLUDE_IRRIGATION", "FALSE")) # switch: carve irrigated cropland into its own ..._IR class
# Organic default is SCHEME-AWARE: TRUE for GLOBIOM (canonical; GLOBIOM_mngmt encodes organic as
# HIO/LIO/IRO/O, so the carve-out yields e.g. Cropland_HIO), FALSE for BIOCLIMA (BIOCLIMA_DS_intermediate
# is blank on organic rows -> carve-out would halt). Always env-overridable.
INCLUDE_ORGANIC <- as.logical(Sys.getenv("DRIVER_INCLUDE_ORGANIC", if (CLASS_SCHEME == "GLOBIOM") "TRUE" else "FALSE"))
# Every setting whose DEFAULT depends on CLASS_SCHEME, printed together. Adding a scheme means taking
# a position on each of these, and three of them were missed for BMLEH in a single day -- the LUM
# lineage (errored loudly), the crop split (caller happened to set it), and the focal year (produced a
# plausible design with 26 focal columns instead of 44 and no error at all). Latent defaults are the
# problem; showing them is the cheapest fix.
.scheme_report <- function() cat(sprintf(
  "\nSCHEME-SENSITIVE DEFAULTS for CLASS_SCHEME=%s\n  LUM lineage   %s\n  crop split    %s\n  focal years   %s%s\n  organic       %s\n  baseline      %s\n  no_choice     %s\n\n",
  CLASS_SCHEME, LUM_DATA_LINEAGE, DO_CROP_SPLIT,
  paste(FOCAL_YEARS, collapse = ","),
  if (all(FOCAL_YEARS == MODEL_YEARS)) " (contemporaneous)" else " (LAGGED -- crop types exist only for 2018)",
  INCLUDE_ORGANIC, BASELINE_CLASS, paste(NO_CHOICE_LU, collapse = ", ")))
# Historical organic-area assumption (ported from run_prior_module_count_model.R): the organic master
# map is a single ~present-day layer applied to every year, so it overstates historical organic area.
# Downweight the organic FRACTION by the EU organic-share ratio vs the master's reference year; the
# un-carved remainder stays conventional (area-conserving). Keyed on the year of the LUM map being loaded.
ORGANIC_AREA_REF_YEAR <- 2020
.eu_org_share <- c("2000" = 3.0, "2010" = 5.2, "2018" = 7.7, "2020" = 9.1) # % of UAA (Eurostat)
organic_area_weight <- function(yr) {
  s <- .eu_org_share[as.character(yr)]
  r <- .eu_org_share[as.character(ORGANIC_AREA_REF_YEAR)]
  if (is.na(s) || is.na(r)) 1 else min(s / r, 1)
}
# Base land-use values folded into a single non-modelled "no_choice" class:
# ---- NATURE / NO-CHOICE STRUCTURE --------------------------------------------------------------
# The natural split now lives in the MAPPING, not here: `GLOBIOM_subclass` (built by
# prep/build_globiom_subclass.R) resolves the old catch-all Natural_unmanaged into
# Natural_grassland / Natural_shrubland using the BIOCLIMA_DS_reporting taxonomy, while
# reproducing all 54 managed class names EXACTLY. Select it with
#   DRIVER_CLASS_COLS="GLOBIOM_subclass"
# Two consequences follow automatically and need no code here:
#   * `nest_tree_by_prefix` groups on the first token, so >=2 Natural_* classes form a
#     "Natural" nest by themselves -- Natural_unmanaged used to be a ROOT SINGLETON.
#   * Natural_other stops being swallowed by no_choice (below), giving that nest 3 leaves.
# no_choice keeps only what is genuinely non-choosable: water, wetlands, NODATA.
# Natural_other (beaches, rock, sparse, burnt, ice) is a real unmanaged cover, so it is
# promoted to its own class and nested with the naturals.
PROMOTE_NATURAL_OTHER <- isTRUE(as.logical(Sys.getenv("DRIVER_PROMOTE_NATURAL_OTHER", "FALSE")))
NO_CHOICE_LU <- c("Waterbodies_marine", "Waterbodies_inland", "Wetlands_natural", "Natural_other", "NODATA")
if (PROMOTE_NATURAL_OTHER) NO_CHOICE_LU <- setdiff(NO_CHOICE_LU, "Natural_other")
# DRIVER_NO_CHOICE lets a project set the sink explicitly, because "non-choosable" is a property of
# the DOWNSTREAM model, not of the land. BMLEH_Los1_CAPRI carries Waterbodies_inland (INLW),
# Wetlands_natural (TWET) and Natural_other (OLND) as real target classes with CAPRI codes: folded
# into no_choice the prior emits no transitions for them and the downscaler has nothing to allocate
# them with. Give the list (comma-separated) to override, or "" to keep only NODATA non-choosable.
if (nzchar(Sys.getenv("DRIVER_NO_CHOICE", ""))) {
  .nc <- trimws(strsplit(Sys.getenv("DRIVER_NO_CHOICE"), ",")[[1]])
  NO_CHOICE_LU <- unique(c(.nc[nzchar(.nc)], "NODATA"))
  cat(sprintf("DRIVER_NO_CHOICE: non-choosable set to %s\n", paste(NO_CHOICE_LU, collapse = ", ")))
}

# The baseline must EXIST in the fitted classes. Natural_unmanaged is gone under
# GLOBIOM_subclass, so default to the largest class overall (Forests_MI, 16.4% of area);
# the driver additionally self-heals to the most prevalent class if this is absent.
BASELINE_CLASS <- Sys.getenv("DRIVER_BASELINE",
                             if (any(grepl("GLOBIOM_subclass", CLASS_COLS))) "Forests_MI" else "Natural_unmanaged")

.scheme_report()

# Where each input's mapping table + join key live.
SOURCE_REGISTRY <- list(
  LUM = list(
    # Prefer the CURATED copy inside this repo (aux_files/) -- gamble-core is meant to be
    # self-contained, and the hand-maintained taxonomy (GLOBIOM_subclass + _nest) lives there.
    # Falls back to the shared LAMASUS_downscaling copy so older checkouts keep working.
    mapping_file = local({
      .c <- c("aux_files/LUM_Code_to_macro_model_mapping.csv",
              "../LAMASUS_downscaling/aux_files/LAMASUS_LUM_thematic_mapping.csv")
      .c <- c(Sys.getenv("LUM_MAPPING_FILE", ""), .c); .c <- .c[nzchar(.c)]
      .h <- .c[file.exists(.c)]
      if (!length(.h)) stop("no LUM mapping found; looked in: ", paste(.c, collapse = ", "))
      .h[1]
    }),
    join_key = "LUM_Code", base_lu = "GLOBIOM_UNFCCC"
  ),
  CLC = list(
    mapping_file = file.path(Sys.getenv("GAMBLE_INPUT_DIR", "input"), "CLC_Code1_LEVEL123_ETL2_mapping.csv"),
    join_key = "Code1", base_lu = "GLOBIOM_UNFCCC"
  )
)
.src <- SOURCE_REGISTRY[[OUTCOME_SOURCE]]
if (is.null(.src)) stop("Unknown OUTCOME_SOURCE: ", OUTCOME_SOURCE)
cat(sprintf("Class mapping: %s\n", .src$mapping_file))

# ---- Single class-derivation function (target set, outcome, focal all use it) ----
# Adds `model_class` to a mapping table. The base class is the make.names() of the
# non-blank CLASS_COLS joined by "_". Organic/irrigation are detected from the
# descriptive LUM_label (robust to GLOBIOM/BIOCLIMA tagging quirks) and only split
# off when the matching switch is on. NO_CHOICE_LU values collapse to "no_choice".
derive_model_class <- function(map_dt, class_cols = CLASS_COLS, base_lu = .src$base_lu) {
  m <- as.data.table(copy(map_dt))
  miss <- setdiff(class_cols, names(m))
  if (length(miss)) {
    stop(sprintf(
      "CLASS_COLS missing from %s mapping: %s",
      OUTCOME_SOURCE, paste(miss, collapse = ", ")
    ))
  }
  mat <- sapply(class_cols, function(cc) trimws(as.character(m[[cc]])))
  mat <- matrix(mat, nrow = nrow(m))
  base <- apply(mat, 1, function(r) paste(r[nzchar(r) & !is.na(r)], collapse = "_"))
  m[, model_class := fifelse(base == "" | is.na(base), "NODATA", make.names(base))]
  lbl <- if ("LUM_label" %in% names(m)) tolower(as.character(m$LUM_label)) else rep("", nrow(m))
  # Suffixes attach only to a resolved base class (never to NODATA). If a switch is on
  # but the chosen CLASS_COLS are blank on those rows (e.g. BIOCLIMA3_* is empty for
  # organic rows), fail loud rather than mint a meaningless "NODATA_O" class.
  if (INCLUDE_ORGANIC) {
    if (nrow(m[grepl("organic", lbl) & model_class == "NODATA"]) > 0) {
      stop(sprintf("INCLUDE_ORGANIC=TRUE but CLASS_COLS [%s] are blank on organic rows (cannot encode organic). Use GLOBIOM_* columns or set INCLUDE_ORGANIC=FALSE.", paste(class_cols, collapse = "+")))
    }
    # Append _O ONLY when the class doesn't already carry an organic marker. GLOBIOM_mngmt already
    # encodes organic (HIO/LIO/IRO/O), so blindly appending would double-label (-> Cropland_HIO_O);
    # schemes whose mngmt does NOT encode organic still get the _O split. (Mirrors the count model.)
    m[grepl("organic", lbl) & model_class != "NODATA" & !grepl("O$", gsub("_", "", model_class)),
      model_class := paste0(model_class, "_O")]
  }
  if (INCLUDE_IRRIGATION) {
    if (nrow(m[grepl("irrigat", lbl) & model_class == "NODATA"]) > 0) {
      stop(sprintf("INCLUDE_IRRIGATION=TRUE but CLASS_COLS [%s] are blank on irrigated rows. Use GLOBIOM_* columns or set INCLUDE_IRRIGATION=FALSE.", paste(class_cols, collapse = "+")))
    }
    # Same guard as organic: skip when GLOBIOM_mngmt already encodes irrigation (IR/IRO).
    m[grepl("irrigat", lbl) & model_class != "NODATA" & !grepl("IR", gsub("_", "", model_class)),
      model_class := paste0(model_class, "_IR")]
  }
  # ---- Natural split: key on LUM_Code (robust) rather than the free-text label ----
  m[, focal_class := model_class]
  # Match on the DERIVED class as well as the raw base LU. Historically only base_lu was checked, so a
  # hand-curated leaf name in CLASS_COLS could never be routed to no_choice however it was spelled.
  # No-op on the current mapping (GLOBIOM_subclass copies GLOBIOM_UNFCCC on exactly the no_choice rows).
  m[model_class %in% NO_CHOICE_LU, model_class := "no_choice"]
  if (base_lu %in% names(m)) m[get(base_lu) %in% NO_CHOICE_LU, model_class := "no_choice"]
  m[]
}

# ---- Load mapping + derive the target class set ---------------------------
mapping_thematic <- as.data.table(fread(.src$mapping_file))

# A CURATED nest column supersedes the no_choice COLLAPSE. NO_CHOICE_LU rewrites those classes to a
# single "no_choice" column in derive_model_class -- i.e. BEFORE the tree is built -- so a curated
# taxonomy that gives Waterbodies/NODATA their own leaves (or moves Wetlands_natural under Natural)
# would be silently erased: the leaves simply would not exist in Y. When the curated column is present
# the TREE does the grouping instead, so the collapse must stand down. DRIVER_KEEP_NO_CHOICE=TRUE
# restores the historical behaviour.
if ((Sys.getenv("DRIVER_NEST_COL", paste0(CLASS_COLS[1], "_nest")) %in% names(mapping_thematic)) &&
    !isTRUE(as.logical(Sys.getenv("DRIVER_KEEP_NO_CHOICE", "FALSE")))) {
  if (length(NO_CHOICE_LU))
    cat(sprintf(">>> curated nests present -> no_choice COLLAPSE disabled (was: %s); these become real leaves.\n",
                paste(NO_CHOICE_LU, collapse = ", ")))
  NO_CHOICE_LU <- character(0)
}
# Drop irrigated/organic source rows when the switch is off, so their codes can
# never appear as classes (their area is folded back by the carve-out logic).
if ("LUM_label" %in% names(mapping_thematic)) {
  if (!INCLUDE_IRRIGATION) mapping_thematic <- mapping_thematic[!grepl("irrigat", tolower(LUM_label))]
  if (!INCLUDE_ORGANIC) mapping_thematic <- mapping_thematic[!grepl("organic", tolower(LUM_label))]
}
mapping_thematic <- derive_model_class(mapping_thematic)

target_classes <- setdiff(unique(mapping_thematic$model_class), "NODATA")
# The crop-type split mints classes (Wheat, Maize, ...) that don't exist in the mapping; register
# them so they survive the final_cats_pixel = intersect(target_classes, ...) outcome filter. Empty
# crop classes (none in the data) are still dropped later by the zero-area filter.
if (DO_CROP_SPLIT) target_classes <- union(target_classes, CROP_LEGEND$crop_class)

# --- Project target classification (optional) --------------------------------------------------
# A project may need its OWN class list rather than the scheme's: BMLEH_Los1_CAPRI wants the 43
# CAPRI-style classes in data/BMLEH_Los1_thematic_mapping.csv, which is a finer cut than any
# CLASS_SCHEME provides. codes/target_class_split.R applies an ordered rename/split cascade to the
# 1km rows immediately after the crop split -- i.e. still at the resolution where regional crop
# statistics apply, and before the pixel split. Off unless a rules file is named.
TARGET_RULES_FILE  <- Sys.getenv("DRIVER_TARGET_RULES", "")
TARGET_SHARES_FILE <- Sys.getenv("DRIVER_TARGET_SHARES", "")
DO_TARGET_RULES <- nzchar(TARGET_RULES_FILE)
TARGET_RULES <- NULL; TARGET_SHARES <- NULL
if (DO_TARGET_RULES) {
  if (!file.exists(TARGET_RULES_FILE)) stop(sprintf("DRIVER_TARGET_RULES file not found: %s", TARGET_RULES_FILE))
  source("codes/target_class_split.R")
  .re <- new.env(); sys.source(TARGET_RULES_FILE, .re)
  .is_rules <- function(n) { o <- get(n, .re); is.data.frame(o) && all(c("order","action","from_class","to") %in% names(o)) }
  .cand <- Filter(.is_rules, ls(.re))
  if (length(.cand) != 1L) stop(sprintf(
    "%s must define exactly ONE rules table (order/action/from_class/to); found %d: %s",
    basename(TARGET_RULES_FILE), length(.cand), paste(.cand, collapse = ", ")))
  TARGET_RULES <- as.data.table(get(.cand[[1]], .re))
  if (any(TARGET_RULES$action == "split")) {
    if (!nzchar(TARGET_SHARES_FILE) || !file.exists(TARGET_SHARES_FILE))
      stop(sprintf("%s has split rules, so DRIVER_TARGET_SHARES must point at a share table (got '%s').",
                   .cand[[1]], TARGET_SHARES_FILE))
    TARGET_SHARES <- fread(TARGET_SHARES_FILE)
    .missg <- setdiff(TARGET_RULES[action == "split", to], unique(TARGET_SHARES$group))
    if (length(.missg)) stop(sprintf("share table has no rows for split group(s): %s",
                                     paste(.missg, collapse = ", ")))
  }
  # Register the minted classes so the final_cats_pixel filter keeps them, and retire the classes
  # the cascade consumes. Order matters: drop the sources FIRST, because a class can be both
  # consumed and re-minted (Cropland_permanent_other_fruit is split at order 16, then renamed
  # back into existence at order 20-21 as the permanent catch-all).
  target_classes <- setdiff(target_classes, TARGET_RULES$from_class)
  target_classes <- union(target_classes, TARGET_RULES[action == "rename", to])
  if (!is.null(TARGET_SHARES)) target_classes <- union(target_classes, unique(TARGET_SHARES$target))
  cat(sprintf("Target classification: %s (%d rules: %d rename, %d split) -> %d target classes\n",
    basename(TARGET_RULES_FILE), nrow(TARGET_RULES), sum(TARGET_RULES$action == "rename"),
    sum(TARGET_RULES$action == "split"), length(target_classes)))
}

# Raw-code -> class lookup for the outcome join (keyed by the source's own join key).
class_lookup <- unique(mapping_thematic[, c(.src$join_key, "model_class", "focal_class"), with = FALSE])

# CLC fallback map (Code1 -> aggregate GLOBIOM class), used only for non-LUM panel years.
thematic_dt <- tryCatch(
  {
    clc <- as.data.table(fread(SOURCE_REGISTRY$CLC$mapping_file))
    clc[, model_class := fifelse(GLOBIOM_UNFCCC %in% NO_CHOICE_LU, "no_choice", make.names(GLOBIOM_UNFCCC))]
    unique(clc[, .(Code1, model_class)])
  },
  error = function(e) NULL
)

cat(sprintf(
  "Choice set [%s | %s%s%s]: %d categories.\n", OUTCOME_SOURCE,
  paste(CLASS_COLS, collapse = "+"),
  if (INCLUDE_IRRIGATION) " +IR" else "", if (INCLUDE_ORGANIC) " +O" else "",
  length(target_classes)
))

# =========================================================================
# 3. 1KM SOURCE DATA LOADING
# =========================================================================
cat("\nLoading 1km Source Data...\n")

COL_ID <- paste0("EEA_", PIXEL_RES, "kmID")
COL_X <- paste0("x_", PIXEL_RES, "kmID")
COL_Y <- paste0("y_", PIXEL_RES, "kmID")

cat("  Loading prior covariates...\n")
# COLUMN-PROJECTED read (performance only): read a SUPERSET of the columns the model can use -- terrain/soil/
# yield + ALL year-suffixed (_YYYY, covers the grep('_2000$') climate discovery, GDP/Pop/GHM_HI/spei48/bio per
# year) + base bio/pr/tas/spei48/socioecon. Drops only columns outside that union (never used) -> result-
# identical, but avoids loading the full 1.3GB parquet (memory + read time). Verified bit-identical X_mat.
# The master 1km parquet and the one_kmID mapping DO NOT live in the same folder: gridwork moved to
# cascadinggamble-core, and ../LAMASUS_gridwork/output (the GRIDWORK_DIR default, which still holds the
# mapping) keeps a STALE copy. Resolving the parquet via GRIDWORK_DIR therefore silently picked the June
# copy, which has GHM_HI/GHM_Ovr but NO GHM_TI -- so the GHM_VARS gate below dropped GHM_TI from the
# design without failing. Pick the NEWEST parquet across both roots (same resolver the count driver
# uses, GAMBLE_MASTER_PARQUET to override) and leave GRIDWORK_DIR for the mapping file.
PRIOR_1KM_PARQUET <- Sys.getenv("GAMBLE_MASTER_PARQUET", "")
if (!nzchar(PRIOR_1KM_PARQUET)) {
  .cands <- c("../cascadinggamble-core/data/02_intermediate/prior_model_1km_master_inputs.parquet",
              file.path(GRIDWORK_DIR, "prior_model_1km_master_inputs.parquet"))
  .cands <- .cands[file.exists(.cands)]
  if (!length(.cands)) stop("master 1km parquet not found in cascadinggamble-core or ", GRIDWORK_DIR)
  PRIOR_1KM_PARQUET <- .cands[order(file.mtime(.cands), decreasing = TRUE)][1]
}
cat(sprintf(">>> master 1km parquet: %s  (modified %s)\n", PRIOR_1KM_PARQUET,
            format(file.mtime(PRIOR_1KM_PARQUET), "%Y-%m-%d %H:%M")))
.pq_nm <- names(arrow::open_dataset(PRIOR_1KM_PARQUET))
.keep_1km <- unique(c("INSPIRE_Europe_buffer_1kmID",
  # RAI removed 2026-08-04 (see the count driver): accessibility/human-footprint is carried by the
  # GHM v3 threat groups GHM_HI + GHM_TI instead. GHM_TI is optional until the gridwork pull lands.
  "Slope_rad", "Elevation", "Aspect_cos_mean", "Aspect_sin_mean", "allPA_share", "CISI",
  "OC_TOP", "ROO", "AWC_TOP", "VS", "yield_index_Low", "yield_index_Med", "yield_index_High",
  "Growing_Degree_Days_gdd5", "Precipitation_Seasonality_bio15", "Annual_Precipitation_bio12",  # static CHELSA (descriptive names, no year suffix, don't match ^bio/pr/tas)
  grep("_[0-9]{4}$", .pq_nm, value = TRUE),
  grep("^(bio|pr|tas|spei48|GHM_HI|GHM_TI|GDP|Pop)", .pq_nm, value = TRUE)))
prior_1km <- as.data.table(dplyr::collect(dplyr::select(
  arrow::open_dataset(PRIOR_1KM_PARQUET),
  dplyr::any_of(intersect(.keep_1km, .pq_nm)))))
cat(sprintf("  prior_1km: column-projected read %d of %d cols x %d rows\n", ncol(prior_1km), length(.pq_nm), nrow(prior_1km)))
setnames(prior_1km, "spei48_2018", "spei48_2020")
prior_1km[, INSPIRE_Europe_buffer_1kmID := as.integer(INSPIRE_Europe_buffer_1kmID)]
for (.col in names(prior_1km)[sapply(prior_1km, is.double)]) {
  set(prior_1km, j = .col, value = fifelse(is.nan(prior_1km[[.col]]), NA_real_, prior_1km[[.col]]))
}
yield_cols <- c("yield_index_Low", "yield_index_Med", "yield_index_High")
# GHM threat groups present for every modelled year (GHM_HI always; GHM_TI once gridwork pulls it).
# They are SOCIO-ECONOMIC year-varying covariates, NOT climate -> excluded from the climate discovery.
GHM_VARS <- c("GHM_HI", "GHM_TI")[vapply(c("GHM_HI", "GHM_TI"),
  function(v) all(paste0(v, "_", COV_YEARS) %in% names(prior_1km)), logical(1))]
.ghm_missing <- setdiff(c("GHM_HI", "GHM_TI"), GHM_VARS)
cat(sprintf("  GHM threat groups available for all COV_YEARS: %s%s\n",
  if (length(GHM_VARS)) paste(GHM_VARS, collapse = ", ") else "(none)",
  if (length(.ghm_missing)) sprintf("   [MISSING: %s -> not modelled; run gridwork download_ghm.R]",
                                    paste(.ghm_missing, collapse = ", ")) else ""))
# Also raise it as a WARNING, not just a cat: silently dropping a requested driver variable is how
# GHM_TI vanished from the 2026-08-12 design (stale parquet had no GHM_TI_2020) and was only caught by
# diffing column names against the previous design afterwards.
if (length(.ghm_missing))
  warning(sprintf("GHM variable(s) %s NOT in %s for COV_YEARS=%s -> SILENTLY EXCLUDED from the design. ",
                  paste(.ghm_missing, collapse = ", "), basename(PRIOR_1KM_PARQUET),
                  paste(COV_YEARS, collapse = ",")),
          "Check you are on the LIVE parquet (GAMBLE_MASTER_PARQUET) before using this design.",
          call. = FALSE, immediate. = TRUE)
climate_cols_raw_2000 <- grep("_2000$", names(prior_1km), value = TRUE)
climate_cols <- sub("_2000$", "", setdiff(climate_cols_raw_2000,
  c(paste0(c("GHM_HI", "GHM_TI", "GHM_Ovr"), "_2000"), "GDP_2000", "Pop_2000")))

cat(paste0("  Loading 1km -> ", PIXEL_RES, "km grid mapping...\n"))
grid_map_pixel <- arrow::read_parquet(file.path(AUXDATA_DIR, mapping_file)) %>% as.data.table()
grid_map_pixel[, `:=`(
  INSPIRE_Europe_buffer_1kmID = as.integer(INSPIRE_Europe_buffer_1kmID),
  EEA_1kmID = as.integer(EEA_1kmID),
  ID = if (!is.null(PIXEL_INTERSECT_COL)) paste0(get(COL_ID), "_", get(PIXEL_INTERSECT_COL)) else as.integer(get(COL_ID)),
  X = as.numeric(get(COL_X)),
  Y = as.numeric(get(COL_Y)),
  Grouping_Key = as.character(get(RE_GROUP_COL))
)]

# geo_region: the region code a PROJECT's regional statistics are keyed on, kept independently of
# the RE grouping key. These are not the same thing and must not be conflated -- the RE block wants
# few, well-populated groups (country), while a target classification wants the finest region its
# statistics support (NUTS2 for the Eurostat crop shares). Note CAPRI_NUTS is the 8-char padded form
# ("NO070000"), so it matches NEITHER a NUTS2 code nor a country code without slicing; DRIVER_GEO_COL
# defaults to the plain NUTS2 column, whose first 2 chars are also the country fallback the share
# cascade uses.
GEO_REGION_COL <- Sys.getenv("DRIVER_GEO_COL", if ("NUTS2" %in% names(grid_map_pixel)) "NUTS2" else RE_GROUP_COL)
if (!GEO_REGION_COL %in% names(grid_map_pixel))
  stop(sprintf("DRIVER_GEO_COL '%s' is not a column of the grid mapping (have: %s)",
               GEO_REGION_COL, paste(grep("NUTS|country", names(grid_map_pixel), value = TRUE), collapse = ", ")))
grid_map_pixel[, geo_region := as.character(get(GEO_REGION_COL))]
if (RE_GROUP_COL == "CAPRI_NUTS") grid_map_pixel[, Grouping_Key := substr(Grouping_Key, 1, 2)]
grid_map_pixel[, pixel_weight := if ("GLOB_5arcminID_area_km2" %in% names(grid_map_pixel)) as.numeric(GLOB_5arcminID_area_km2) else 1.0]
grid_map_pixel[is.na(pixel_weight) | pixel_weight == 0, pixel_weight := 1.0]

# New Observation Unit: ID (Grid Cell) + Grouping Key
# Hold the region key in its OWN 1km lookup rather than as a column of grid_map_pixel. That table is
# area-weighted and aggregated in several places (every value column gets multiplied by
# pixel_weight), so a character key riding along inside it is a bug waiting to happen.
GEO_LOOKUP_1KM <- unique(grid_map_pixel[, .(INSPIRE_Europe_buffer_1kmID, EEA_1kmID, geo = geo_region)])
# id crosswalk for sources keyed differently from the LUM export (see .load_crop_source)
grid_map_pixel_raw_ids <- unique(grid_map_pixel[, .(INSPIRE_Europe_buffer_1kmID, EEA_1kmID)])
grid_map_pixel <- unique(grid_map_pixel[, .(INSPIRE_Europe_buffer_1kmID, EEA_1kmID, ID, X, Y, Grouping_Key, pixel_weight)])

# =========================================================================
# 4. STATIC COVARIATES
# =========================================================================
cat("\nAggregating static covariates...\n")
# Soil
spatial_cat_vars <- intersect(c("OC_TOP", "ROO", "AWC_TOP", "VS"), colnames(prior_1km))
x_pixel_soil_list <- lapply(spatial_cat_vars, function(v) {
  res <- dcast(prior_1km[!is.na(get(v)), .N, by = .(INSPIRE_Europe_buffer_1kmID, Val = get(v))],
    INSPIRE_Europe_buffer_1kmID ~ Val,
    value.var = "N", fill = 0
  )
  setDT(res)
  res <- res[grid_map_pixel, on = "INSPIRE_Europe_buffer_1kmID", nomatch = 0L]

  # Weight by pixel_weight for area-accurate aggregation
  val_cols <- setdiff(colnames(res), c("INSPIRE_Europe_buffer_1kmID", "EEA_1kmID", "ID", "X", "Y", "Grouping_Key", "pixel_weight"))
  res[, (val_cols) := lapply(.SD, function(x) x * pixel_weight), .SDcols = val_cols]

  res <- res[, lapply(.SD, sum), by = .(ID, Grouping_Key), .SDcols = val_cols]
  cls_cols <- setdiff(colnames(res), c("ID", "Grouping_Key"))
  res[, (cls_cols) := .SD / pmax(1e-9, rowSums(.SD)), .SDcols = cls_cols]
  setnames(res, cls_cols, paste0(v, "_s", cls_cols))

  # DROP share classes that no modelled pixel occupies. dcast creates a column for every value
  # present in the 1km source, but the inner join to grid_map_pixel keeps only pixels in the model
  # grid -- so a class living entirely outside it survives as an ALL-ZERO column. That column
  # carries no information and cannot be identified (X %*% beta is unchanged by its coefficient),
  # and it does real damage downstream: the sampler's rank guard drops it, and until 2026-09-14 a
  # zero-variance column was then mistaken for an intercept and handed the un-centring shift,
  # putting a phantom country-varying coefficient into beta_rotated and on to the downscaler.
  # ROO_s6 was exactly this, in the BMLEH and both GLOBIOM designs.
  .sh <- paste0(v, "_s", cls_cols)
  .dead <- .sh[vapply(.sh, function(cc) { z <- res[[cc]]; all(!is.finite(z) | z == 0) }, TRUE)]
  if (length(.dead)) {
    res[, (.dead) := NULL]
    cat(sprintf("    %s: dropped %d empty share class(es) -- no modelled pixel occupies them: %s\n",
                v, length(.dead), paste(.dead, collapse = ", ")))
  }
  return(res)
})
x_pixel_soil_list <- lapply(x_pixel_soil_list, function(dt) setkey(as.data.table(dt), ID, Grouping_Key))
x_pixel_soil <- Reduce(function(x, y) merge(x, y, all = TRUE), x_pixel_soil_list)

# Soil contrasts are now computed dynamically alongside focal contrasts in section 5c

# =========================================================================
# 5. MULTI-YEAR PANEL ASSEMBLY
# =========================================================================

# ---- Organic / irrigation carve-out (applied by load_granular_lum when the switches are on) ----
# The base LUM maps carry no organic (<code>999) or irrigated (2006) twins; the carve-out
# synthesizes them from present-day auxiliary layers, exactly as the legacy 2018 branch /
# count model do, plus the count model's YEAR-SPECIFIC organic-area downweight. Auxiliary
# layers are cached (loaded once) since load_granular_lum runs up to 2x per tier (outcome + focal).
.organic_master_cache <- NULL
.irrigation_cache <- NULL
.org_cols_carve <- c("Cropland_organic", "Livestock_organic", "Mixed_organic", "All_organic")
.load_organic_master <- function() {
  if (is.null(.organic_master_cache)) {
    om <- as.data.table(readRDS(file.path(GRIDWORK_DIR, "organic_certificaties_final_calibrated_master.rds")))
    om <- om[, .SD, .SDcols = c("EEA_1kmID", .org_cols_carve)]
    om[, EEA_1kmID := as.integer(EEA_1kmID)]
    for (j in .org_cols_carve) set(om, which(is.na(om[[j]])), j, 0)
    .organic_master_cache <<- om
  }
  .organic_master_cache
}
.load_irrigation <- function() {
  if (is.null(.irrigation_cache)) {
    ir <- as.data.table(readRDS(paste0(DS_DIR, "input/irrigation_binary_EEA_1kmID.rds")))
    ir[, EEA_1kmID := as.integer(EEA_1kmID)]
    .irrigation_cache <<- ir[irrigation_binary == 1, .(EEA_1kmID, irrigation_binary)]
  }
  .irrigation_cache
}

# Carve a 1km WIDE table (cols: EEA_1kmID + LUM<code> areas). `year` sets the organic downweight.
# Area-conserving: organic/irrigated area moves between codes, the un-carved remainder stays
# conventional. Returns the carved wide table (EEA_1kmID dropped).
apply_lum_carveouts <- function(wide, year) {
  if (INCLUDE_IRRIGATION) {
    wide <- merge(wide, .load_irrigation(), by = "EEA_1kmID", all.x = TRUE)
    wide[is.na(irrigation_binary), irrigation_binary := 0]
    non_irr_codes <- mapping_thematic[GLOBIOM_UNFCCC == "Cropland" & GLOBIOM_mngmt != "IR" & !grepl("O", GLOBIOM_mngmt), LUM_Code]
    non_irr_cols <- intersect(paste0("LUM", non_irr_codes), names(wide))
    wide[, LUM2006 := irrigation_binary * rowSums(.SD, na.rm = TRUE), .SDcols = non_irr_cols]
    wide[, (non_irr_cols) := {
      dt_sub <- .SD; rs <- rowSums(dt_sub, na.rm = TRUE); w <- dt_sub / rs; w[is.na(w)] <- 0
      pmax(dt_sub - (LUM2006 * w), 0)
    }, .SDcols = non_irr_cols]
    wide[, irrigation_binary := NULL]
  }
  if (INCLUDE_ORGANIC) {
    wide <- merge(wide, .load_organic_master(), by = "EEA_1kmID", all.x = TRUE)
    for (j in .org_cols_carve) set(wide, which(is.na(wide[[j]])), j, 0)
    arable_cols    <- intersect(c("LUM2001", "LUM2002", "LUM2003", "LUM2004", "LUM2005", if (INCLUDE_IRRIGATION) "LUM2006"), names(wide))
    permanent_cols <- intersect(c("LUM3001", "LUM3002", "LUM3003", "LUM3004", "LUM3005"), names(wide))
    pasture_cols   <- intersect(c(paste0("LUM", 14:26), "LUM4012", "LUM4013"), names(wide))
    wide[, total_arable_lum    := rowSums(.SD, na.rm = TRUE), .SDcols = arable_cols]
    wide[, total_permanent_lum := rowSums(.SD, na.rm = TRUE), .SDcols = permanent_cols]
    wide[, total_pasture_lum   := rowSums(.SD, na.rm = TRUE), .SDcols = pasture_cols]
    wide[, total_agri_lum      := total_arable_lum + total_permanent_lum + total_pasture_lum]
    # Present-day organic AREA -> per-cell FRACTION (Cropland organic shared across arable+permanent).
    wide[, Cropland_organic  := ifelse(total_arable_lum + total_permanent_lum > 0, pmin(Cropland_organic / (total_arable_lum + total_permanent_lum), 1.0), 0)]
    wide[, Livestock_organic := ifelse(total_pasture_lum > 0, pmin(Livestock_organic / total_pasture_lum, 1.0), 0)]
    wide[, Mixed_organic     := ifelse(total_agri_lum > 0, pmin(Mixed_organic / total_agri_lum, 1.0), 0)]
    wide[, All_organic       := ifelse(total_agri_lum > 0, pmin(All_organic / total_agri_lum, 1.0), 0)]
    # Historical organic-area downweight (count-model assumption), keyed on the map's year.
    ow <- organic_area_weight(year)
    if (ow < 1) {
      for (.j in .org_cols_carve) wide[, (.j) := get(.j) * ow]
      cat(sprintf("    Organic-area downweight for %s: x%.2f (EU organic share %.1f%% vs %d ref %.1f%%)\n",
        year, ow, .eu_org_share[as.character(year)], ORGANIC_AREA_REF_YEAR, .eu_org_share[as.character(ORGANIC_AREA_REF_YEAR)]))
    }
    for (col in c(arable_cols, permanent_cols, pasture_cols)) {
      wide[, primary_org_area := 0]
      if (col %in% c(arable_cols, permanent_cols)) wide[, primary_org_area := get(col) * Cropland_organic]
      else if (col %in% pasture_cols)              wide[, primary_org_area := get(col) * Livestock_organic]
      wide[, mixed_org_area := pmax(get(col) - primary_org_area, 0) * Mixed_organic]
      wide[, all_org_area   := pmax(get(col) - primary_org_area - mixed_org_area, 0) * All_organic]
      wide[, (col) := pmax(get(col) - primary_org_area - mixed_org_area - all_org_area, 0)]
      wide[, paste0(col, "999") := primary_org_area + mixed_org_area + all_org_area]
      wide[, c("primary_org_area", "mixed_org_area", "all_org_area") := NULL]
    }
    wide[, c("total_arable_lum", "total_permanent_lum", "total_pasture_lum", "total_agri_lum", .org_cols_carve) := NULL]
  }
  if ("EEA_1kmID" %in% names(wide)) wide[, EEA_1kmID := NULL]
  wide
}

# ---- General granular-LUM loader (registered years) -----------------------
# Reads the granular 1km LUM file registered for `year` in LUM_SOURCE_REGISTRY, optionally applies
# the organic/irrigation carve-out (INCLUDE_ORGANIC/INCLUDE_IRRIGATION), maps raw codes ->
# model_class + focal_class via class_lookup, aggregates area to the model grid (weighted by
# pixel_weight), and returns a compiled_long data.table with the same shape the outcome/focal
# blocks expect: ID, Grouping_Key, X, Y, model_class, focal_class, area_km2. Used for BOTH the
# outcome (Y) of a tier AND its (possibly lagged) focal input.
# ---- HRL crop-type split (AgMIP): replace LUM cropland area with specific crop classes ----------
# Operates on a per-1km-cell (join_id) class table dt (cols: join_id, model_class, focal_class, area).
# Splits the Cropland_arable_other / Cropland_permanent_other area into HRL crop classes using each
# cell's crop composition. Area-conserving per cell: crops get min(H, A) distributed by HRL shares,
# the non-overlap remainder max(A - H, 0) stays in the residual class (H = HRL crop area of that
# group in the cell, A = LUM cropland area of that group). HRL file cached (read once per year).
.crop_cache <- new.env(parent = emptyenv())
.load_crop_source <- function(year, lum_id_col = NULL) {
  key <- paste0(as.character(year), "|", if (is.null(lum_id_col)) "" else lum_id_col)
  if (is.null(.crop_cache[[key]])) {
    reg <- CROP_SOURCE_REGISTRY[[key0 <- as.character(year)]]
    cr <- as.data.table(readRDS(reg$file))[, .SD, .SDcols = c(reg$id_col, reg$code_col, reg$area_col)]
    setnames(cr, c(reg$id_col, reg$code_col, reg$area_col), c("join_id", "Crop_Type_Code", "carea"))
    cr[, join_id := as.integer(join_id)]
    # TRANSLATE the crop ids into the LUM's id space when the two sources are keyed differently.
    # HRL Crop Types are keyed INSPIRE_Europe_buffer_1kmID; the canonical LUM export is keyed EEA_1kmID.
    # Joining one against the other is not an error -- both are integers and 14.2% of the values
    # coincide by chance -- so it silently dropped 86% of the crop data and the area assert still
    # passed, because unattributed cropland simply stays in the residual. Durum wheat disappeared
    # from the design entirely before this was caught.
    if (!is.null(lum_id_col) && !identical(lum_id_col, reg$id_col)) {
      xw <- unique(grid_map_pixel_raw_ids[, .(from = as.integer(get(reg$id_col)),
                                              to   = as.integer(get(lum_id_col)))])
      xw <- xw[!is.na(from) & !is.na(to) & !duplicated(from)]
      n0 <- uniqueN(cr$join_id)
      cr <- merge(cr, xw, by.x = "join_id", by.y = "from")
      cr[, join_id := to][, to := NULL]
      .keep_frac <- uniqueN(cr$join_id) / max(n0, 1L)
      cat(sprintf("    crop ids translated %s -> %s (%s of %s kept, %.1f%%)\n", reg$id_col, lum_id_col,
                  format(uniqueN(cr$join_id), big.mark = ","), format(n0, big.mark = ","), 100 * .keep_frac))
      if (.keep_frac < as.numeric(Sys.getenv("DRIVER_CROP_ID_MIN_KEEP", "0.95")))
        stop(sprintf("crop-id translation kept only %.1f%% of ids (%s -> %s). The crosswalk is not covering this source; a partial join here silently drops crop area and the run still looks healthy.",
                     100 * .keep_frac, reg$id_col, lum_id_col))
    }
    cr <- merge(cr, CROP_LEGEND, by = "Crop_Type_Code")  # keep only legend crops (drops 65535/3100/3200)
    .crop_cache[[key]] <- cr[carea > 0]
  }
  .crop_cache[[key]]
}
apply_crop_split <- function(dt, year, lum_id_col = NULL) {
  cr <- .load_crop_source(year, lum_id_col)
  new_rows <- list()
  for (grp in c("arable", "permanent")) {
    resid_class <- if (grp == "arable") "Cropland_arable_other" else "Cropland_permanent_other"
    A <- dt[model_class == resid_class, .(A = sum(area)), by = join_id]
    if (nrow(A) == 0) next
    Hc <- cr[crop_group == grp, .(H_c = sum(carea)), by = .(join_id, crop_class)]  # HRL area per crop/cell
    Htot <- Hc[, .(H = sum(H_c)), by = join_id]
    j <- merge(A, Htot, by = "join_id", all.x = TRUE)
    j[is.na(H), H := 0]
    j[, `:=`(scale = fifelse(H > 0, pmin(1, A / H), 0), resid = pmax(A - H, 0))]
    crop_rows <- merge(Hc, j[, .(join_id, scale)], by = "join_id")[scale > 0,
                   .(join_id, model_class = crop_class, focal_class = crop_class, area = H_c * scale)]
    resid_rows <- j[resid > 1e-9, .(join_id, model_class = resid_class, focal_class = resid_class, area = resid)]
    new_rows[[grp]] <- rbind(crop_rows, resid_rows)
  }
  # Replace the residual-class rows with the crop-split rows; keep everything else untouched.
  dt_keep <- dt[!model_class %in% c("Cropland_arable_other", "Cropland_permanent_other")]
  out <- rbind(dt_keep, rbindlist(new_rows), fill = TRUE)

  # HOW MUCH CROPLAND ACTUALLY GOT TYPED. This is the check that matters, and it is deliberately
  # independent of the id mechanism: area conservation cannot see this failure at all, because
  # cropland HRL fails to match simply stays in the residual and the totals still balance. When the
  # crop ids were joined against the wrong key, 86% of the crop data vanished, the assert passed, and
  # the only symptom was one class quietly missing from the design.
  .crop_cls <- CROP_LEGEND$crop_class
  .typed <- out[model_class %in% .crop_cls, sum(area)]
  .resid <- out[model_class %in% c("Cropland_arable_other", "Cropland_permanent_other"), sum(area)]
  .frac  <- .typed / max(.typed + .resid, 1e-9)
  cat(sprintf("    crop split: %.1f%% of cropland typed (%.0f of %.0f km2)\n",
              100 * .frac, .typed, .typed + .resid))
  .minf <- as.numeric(Sys.getenv("DRIVER_CROP_MIN_TYPED", "0.25"))
  if (.frac < .minf)
    stop(sprintf("the crop split typed only %.1f%% of cropland (expected well above %.0f%%). Area is still conserved -- untyped cropland stays in the residual -- so this will NOT show up as lost area. Check that the crop source and the LUM export share an id space.",
                 100 * .frac, 100 * .minf))
  out
}

load_granular_lum <- function(year) {
  reg <- LUM_SOURCE_REGISTRY[[as.character(year)]]
  if (is.null(reg)) stop(sprintf("No LUM_SOURCE_REGISTRY entry for year %s.", year))
  if (!file.exists(reg$file)) stop(sprintf("Registered LUM file missing for %s: %s", year, reg$file))

  raw <- as.data.table(readRDS(reg$file))
  miss <- setdiff(c(reg$id_col, reg$code_col, "area_km2"), names(raw))
  if (length(miss)) stop(sprintf("LUM file %s missing column(s): %s", basename(reg$file), paste(miss, collapse = ", ")))
  raw <- raw[, .SD, .SDcols = c(reg$id_col, reg$code_col, "area_km2")]
  setnames(raw, c(reg$id_col, reg$code_col), c("join_id", "LUM_Code"))
  raw[, join_id := as.integer(join_id)]
  raw[, LUM_Code := as.integer(LUM_Code)]
  raw <- raw[area_km2 > 0]

  # Optional 1km carve-out: synthesize organic (<code>999) / irrigated (2006) twins on the raw
  # 1km cells (join_id) BEFORE the grid/pixel split, joining the EEA_1kmID-keyed organic &
  # irrigation layers via the grid's id mapping (one EEA per 1km cell). Area-conserving.
  if (INCLUDE_ORGANIC || INCLUDE_IRRIGATION) {
    j2e <- unique(grid_map_pixel[, .(join_id = as.integer(get(reg$id_col)), EEA_1kmID = as.integer(EEA_1kmID))])
    j2e <- j2e[!duplicated(join_id)]
    wide <- dcast(raw, join_id ~ LUM_Code, value.var = "area_km2", fill = 0, fun.aggregate = sum)
    setnames(wide, setdiff(names(wide), "join_id"), paste0("LUM", setdiff(names(wide), "join_id")))
    wide <- merge(wide, j2e, by = "join_id", all.x = TRUE)
    a_pre <- sum(as.matrix(wide[, .SD, .SDcols = patterns("^LUM")]), na.rm = TRUE)
    wide <- apply_lum_carveouts(wide, year)
    a_post <- sum(as.matrix(wide[, .SD, .SDcols = patterns("^LUM")]), na.rm = TRUE)
    if (abs(a_pre - a_post) > 1) stop(sprintf("load_granular_lum(%s): carve-out changed 1km area by %.1f km2 (should conserve).", year, a_pre - a_post))
    raw <- melt(wide, id.vars = "join_id", measure.vars = grep("^LUM", names(wide), value = TRUE),
                variable.name = "LUM_Code", value.name = "area_km2")
    raw[, LUM_Code := as.integer(sub("LUM", "", as.character(LUM_Code)))]
    raw <- raw[area_km2 > 0]
  }

  # Map raw LUM codes -> model_class / focal_class at the 1km (join_id) level, BEFORE the pixel
  # split, so the crop-type split can act on each cell's cropland composition. Fail loud on
  # unmapped area (never drop silently).
  area_before <- sum(raw$area_km2, na.rm = TRUE)
  raw <- merge(raw, class_lookup, by = "LUM_Code", all.x = TRUE)
  unmapped <- raw[is.na(model_class) & area_km2 > 0]
  if (nrow(unmapped) > 0) {
    bad <- unmapped[, .(area = sum(area_km2)), by = LUM_Code][order(-area)]
    stop(sprintf(
      "load_granular_lum(%s): %d LUM codes carrying %.1f km2 have no class under scheme [%s]. Codes: %s",
      year, nrow(bad), sum(bad$area), paste(CLASS_COLS, collapse = "+"),
      paste(head(bad$LUM_Code, 12), collapse = ", ")
    ))
  }
  raw[is.na(model_class) | model_class == "NODATA", model_class := "no_choice"]
  raw[is.na(focal_class), focal_class := "NODATA"]
  rawc <- raw[, .(area = sum(area_km2, na.rm = TRUE)), by = .(join_id, model_class, focal_class)]
  if (abs(area_before - sum(rawc$area)) > 1) {
    stop(sprintf("load_granular_lum(%s) lost %.1f km2 on the LUM_Code join.", year, area_before - sum(rawc$area)))
  }

  # Optional crop-type split (AgMIP): replace the Cropland_arable_other / Cropland_permanent_other
  # area with specific HRL crop-type classes per 1km cell (area-conserving). Years with no HRL
  # source keep the cropland in the residual classes (logged, not an error).
  if (DO_CROP_SPLIT) {
    if (!is.null(CROP_SOURCE_REGISTRY[[as.character(year)]])) {
      a_pre <- sum(rawc$area)
      rawc <- apply_crop_split(rawc, year, reg$id_col)
      if (abs(a_pre - sum(rawc$area)) > 1) {
        stop(sprintf("load_granular_lum(%s): crop split changed 1km area by %.1f km2 (should conserve).", year, a_pre - sum(rawc$area)))
      }
    } else {
      cat(sprintf("    (no HRL crop source for %s -> cropland kept as Cropland_*_other)\n", year))
    }
  }

  # Project target classification: rename/split the model classes onto the project's own target
  # list (see DO_TARGET_RULES above). Runs AFTER the crop split so the HRL crop classes it mints
  # (Wheat, Fruits, ...) are available as split sources, and still at 1km so geo_region carries
  # the untruncated NUTS2 code the shares are keyed on.
  if (DO_TARGET_RULES) {
    a_pre <- sum(rawc$area)
    .geo <- unique(GEO_LOOKUP_1KM[, .(join_id = as.integer(get(reg$id_col)), geo)])
    .geo <- .geo[!is.na(join_id) & !duplicated(join_id)]
    rawc <- apply_target_classification(rawc, TARGET_RULES, shares = TARGET_SHARES,
                                        geo_lookup = .geo, verbose = TRUE)
    rawc <- rawc[, .(area = sum(area, na.rm = TRUE)), by = .(join_id, model_class, focal_class)]
    if (abs(a_pre - sum(rawc$area)) > 1) {
      stop(sprintf("load_granular_lum(%s): target classification changed 1km area by %.1f km2 (should conserve).",
                   year, a_pre - sum(rawc$area)))
    }
  }

  # Grid/pixel split: distribute each 1km cell across its grid pixels (pixel_weight), then aggregate.
  gm <- grid_map_pixel[, .(join_id = as.integer(get(reg$id_col)), ID, Grouping_Key, X, Y, pixel_weight)]
  # The fan-out here is intended and inherently many-to-many: a 1km cell carries several classes and
  # is split across several grid pixels. With a fine target classification (39 classes) the product
  # exceeds data.table's default cartesian guard, which is a size heuristic rather than a
  # correctness check. Unlike the steps above, this one is NOT area-conserving by construction --
  # all.x = FALSE drops 1km cells outside the grid, and pixel_weight redistributes the rest -- so
  # there is deliberately no conservation assert here; the panel-level area totals are the check.
  merged <- merge(rawc, gm, by = "join_id", all.x = FALSE, allow.cartesian = TRUE)
  merged[is.na(pixel_weight), pixel_weight := 1.0]
  merged[, area := area * pixel_weight]

  merged[, .(area_km2 = sum(area, na.rm = TRUE)),
         by = .(ID, Grouping_Key, X, Y, model_class, focal_class)]
}

cat(paste0("\nBuilding ", PIXEL_RES, "km Level Model Panel...\n"))
T_PAIRS <- lapply(seq_along(MODEL_YEARS), function(i) {
  list(out_year = MODEL_YEARS[i], cov_year = COV_YEARS[i], focal_year = FOCAL_YEARS[i])
})

dat_pixel_list <- lapply(T_PAIRS, function(tp) {
  .fy <- if (is.na(tp$focal_year)) tp$out_year else tp$focal_year
  cat(sprintf("  Processing Tier: Outcome = %s | Covariates = %s | Focal%s = %s\n",
    tp$out_year, tp$cov_year,
    if (!is.na(tp$out_year) && !is.na(.fy) && .fy != tp$out_year) "(t-1 lag)" else "", .fy))

  # 1. OUTCOME (Y)
  temp_y_pixel_yr_wide <- NULL
  # Registry path: a granular LUM file registered for this outcome year -> use the general loader,
  # which applies the organic/irrigation carve-out internally when those switches are on. The legacy
  # 2018 EEA-file organic branch below is kept only as a fallback for years NOT in the registry.
  .use_registry_outcome <- !is.na(tp$out_year) &&
    !is.null(LUM_SOURCE_REGISTRY[[as.character(tp$out_year)]])
  if (.use_registry_outcome) {
    cat(sprintf("  Loading granular %s outcome via LUM_SOURCE_REGISTRY: %s\n",
      tp$out_year, basename(LUM_SOURCE_REGISTRY[[as.character(tp$out_year)]]$file)))
    temp_compiled_long <- load_granular_lum(tp$out_year)
    temp_y_pixel_yr_wide <- dcast(temp_compiled_long, ID + Grouping_Key ~ model_class,
      value.var = "area_km2", fill = 0, fun.aggregate = sum)
    setDT(temp_y_pixel_yr_wide)
    cat(sprintf(
      "  >>> Final Outcome Area Check (Year %s): Total area = %.1f km2\n",
      tp$out_year, sum(temp_y_pixel_yr_wide[, .SD, .SDcols = setdiff(names(temp_y_pixel_yr_wide), c("ID", "Grouping_Key"))], na.rm = TRUE)
    ))
  } else if (!is.na(tp$out_year) && tp$out_year == 2018) {
    cat(sprintf("  Processing granular 2018 LUM with Organic and Irrigation logic...\n"))

    # A. Base LUM Map
    temp_map_1km_raw <- as.data.table(readRDS(paste0(DS_DIR, "input/LUM_fit_with_energy_levels_and_new_FM_2018_EEA_1kmID.rds")))
    temp_map_1km_raw[, EEA_1kmID := as.integer(EEA_1kmID)]

    initial_area <- sum(temp_map_1km_raw$area_km2, na.rm = TRUE)
    temp_map_1km_raw <- temp_map_1km_raw[area_km2 > 0 & EEA_1kmID %in% unique(grid_map_pixel$EEA_1kmID)]
    after_grid_filter_area <- sum(temp_map_1km_raw$area_km2, na.rm = TRUE)

    cat(sprintf(
      "    Area Check: Initial = %.1f | After Grid Filter = %.1f (Loss = %.1f)\n",
      initial_area, after_grid_filter_area, initial_area - after_grid_filter_area
    ))

    # B. Irrigation layer (loaded only when the switch is on)
    if (INCLUDE_IRRIGATION) {
      temp_binary_irrigation_map <- as.data.table(readRDS(paste0(DS_DIR, "input/irrigation_binary_EEA_1kmID.rds")))
      temp_binary_irrigation_map[, EEA_1kmID := as.integer(EEA_1kmID)]
      temp_binary_irrigation_map <- temp_binary_irrigation_map[irrigation_binary == 1 & EEA_1kmID %in% unique(grid_map_pixel$EEA_1kmID)]
    }

    # C. Organic layer (loaded only when the switch is on)
    if (INCLUDE_ORGANIC) {
      temp_organic_map <- readRDS(paste0(GRIDWORK_DIR, "/organic_certificaties_final_calibrated_master.rds"))
      org_cols <- c("Cropland_organic", "Livestock_organic", "Mixed_organic", "All_organic")
      temp_organic_map <- temp_organic_map[, .SD, .SDcols = c("EEA_1kmID", org_cols)]
      for (j in org_cols) set(temp_organic_map, which(is.na(temp_organic_map[[j]])), j, 0)
    }

    # D. Wide LUM by raw code
    temp_lum_wide_1km <- dcast(temp_map_1km_raw, EEA_1kmID ~ LUM_fit_with_energy_levels_and_new_FM_2018, value.var = "area_km2", fill = 0)

    if (ncol(temp_lum_wide_1km) < 2) {
      cat(sprintf("  WARNING: Granular 2018 LUM map is empty for current pixel subset. Falling back to CLC.\n"))
      temp_y_pixel_yr_wide <- NULL # This will trigger the CLC fallback or dummy logic below
    } else {
      setnames(temp_lum_wide_1km, old = names(temp_lum_wide_1km)[-1], new = paste0("LUM", names(temp_lum_wide_1km)[-1]))

      # --- Irrigation intersection (SWITCH): carve irrigated arable into LUM2006 ---
      # Off: irrigated area stays in the rainfed intensity classes (no …_IR class).
      if (INCLUDE_IRRIGATION) {
        temp_lum_wide_1km <- merge(temp_lum_wide_1km, temp_binary_irrigation_map, by = "EEA_1kmID", all.x = TRUE)
        temp_lum_wide_1km[is.na(irrigation_binary), irrigation_binary := 0]
        temp_lum_wide_1km[is.na(area_km2), area_km2 := 0]
        temp_non_irrigated_codes <- mapping_thematic[GLOBIOM_UNFCCC == "Cropland" & GLOBIOM_mngmt != "IR" & !grepl("O", GLOBIOM_mngmt), LUM_Code]
        temp_non_irrigated_cols <- intersect(paste0("LUM", temp_non_irrigated_codes), names(temp_lum_wide_1km))
        temp_lum_wide_1km[, LUM2006 := irrigation_binary * rowSums(.SD, na.rm = TRUE), .SDcols = temp_non_irrigated_cols]
        temp_lum_wide_1km[, (temp_non_irrigated_cols) := {
          dt_sub <- .SD
          row_sums <- rowSums(dt_sub, na.rm = TRUE)
          weights <- dt_sub / row_sums
          weights[is.na(weights)] <- 0
          pmax(dt_sub - (LUM2006 * weights), 0)
        }, .SDcols = temp_non_irrigated_cols]
        temp_lum_wide_1km[, c("irrigation_binary", "area_km2") := NULL]
      }

      # --- Organic intersection (SWITCH): carve organic twins into <code>999 ---
      # Off: organic area stays in its conventional class (no …_O class).
      if (INCLUDE_ORGANIC) {
        temp_final_dt_1km <- merge(temp_lum_wide_1km, temp_organic_map, by = "EEA_1kmID", all.x = TRUE)
        for (j in org_cols) set(temp_final_dt_1km, which(is.na(temp_final_dt_1km[[j]])), j, 0)

        arable_cols <- intersect(c(
          "LUM2001", "LUM2002", "LUM2003", "LUM2004", "LUM2005",
          if (INCLUDE_IRRIGATION) "LUM2006"
        ), names(temp_final_dt_1km))
        permanent_cols <- intersect(c("LUM3001", "LUM3002", "LUM3003", "LUM3004", "LUM3005"), names(temp_final_dt_1km))
        pasture_cols <- intersect(c(paste0("LUM", 14:26), "LUM4012", "LUM4013"), names(temp_final_dt_1km))

        temp_final_dt_1km[, total_arable_lum := rowSums(.SD, na.rm = TRUE), .SDcols = arable_cols]
        temp_final_dt_1km[, total_permanent_lum := rowSums(.SD, na.rm = TRUE), .SDcols = permanent_cols]
        temp_final_dt_1km[, total_pasture_lum := rowSums(.SD, na.rm = TRUE), .SDcols = pasture_cols]
        temp_final_dt_1km[, total_agri_lum := total_arable_lum + total_permanent_lum + total_pasture_lum]

        # Cropland organic is shared across Arable AND Permanent
        temp_final_dt_1km[, Cropland_organic := ifelse(total_arable_lum + total_permanent_lum > 0, pmin(Cropland_organic / (total_arable_lum + total_permanent_lum), 1.0), 0)]
        temp_final_dt_1km[, Livestock_organic := ifelse(total_pasture_lum > 0, pmin(Livestock_organic / total_pasture_lum, 1.0), 0)]
        temp_final_dt_1km[, Mixed_organic := ifelse(total_agri_lum > 0, pmin(Mixed_organic / total_agri_lum, 1.0), 0)]
        temp_final_dt_1km[, All_organic := ifelse(total_agri_lum > 0, pmin(All_organic / total_agri_lum, 1.0), 0)]

        all_agri_lum_cols <- c(arable_cols, permanent_cols, pasture_cols)
        for (col in all_agri_lum_cols) {
          temp_final_dt_1km[, primary_org_area := 0]
          if (col %in% c(arable_cols, permanent_cols)) {
            temp_final_dt_1km[, primary_org_area := get(col) * Cropland_organic]
          } else if (col %in% pasture_cols) {
            temp_final_dt_1km[, primary_org_area := get(col) * Livestock_organic]
          }
          temp_final_dt_1km[, mixed_org_area := pmax(get(col) - primary_org_area, 0) * Mixed_organic]
          temp_final_dt_1km[, all_org_area := pmax(get(col) - primary_org_area - mixed_org_area, 0) * All_organic]

          # Update conventional and create organic
          temp_final_dt_1km[, (col) := pmax(get(col) - primary_org_area - mixed_org_area - all_org_area, 0)]
          temp_final_dt_1km[, paste0(col, "999") := primary_org_area + mixed_org_area + all_org_area]
          temp_final_dt_1km[, c("primary_org_area", "mixed_org_area", "all_org_area") := NULL]
        }
        temp_final_dt_1km[, c("total_arable_lum", "total_permanent_lum", "total_pasture_lum", "total_agri_lum", "Cropland_organic", "Livestock_organic", "Mixed_organic", "All_organic") := NULL]
      } else {
        temp_final_dt_1km <- temp_lum_wide_1km
      }

      # Negative area validation
      lum_cols <- grep("^LUM", names(temp_final_dt_1km), value = TRUE)
      neg_check <- temp_final_dt_1km[, lapply(.SD, function(x) any(x < -1e-6)), .SDcols = lum_cols]
      if (any(unlist(neg_check))) {
        neg_names <- names(neg_check)[unlist(neg_check)]
        cat(sprintf("    WARNING: Negative areas detected after carve-out in: %s. Clipping to 0.\n", paste(neg_names, collapse = ", ")))
        for (col in neg_names) temp_final_dt_1km[get(col) < 0, (col) := 0]
      }

      # Diagnostic: check area before aggregation to find the 17.6k loss
      area_pre_agg <- sum(rowSums(temp_final_dt_1km[, .SD, .SDcols = patterns("^LUM")], na.rm = TRUE))
      cat(sprintf("    Area before Step F (Aggregation) = %.1f km2\n", area_pre_agg))

      # F. Aggregate to Model Resolution (Area-Weighted)
      temp_final_dt_1km <- merge(temp_final_dt_1km, grid_map_pixel[, .(EEA_1kmID, ID, Grouping_Key, X, Y, pixel_weight)], by = "EEA_1kmID", all.x = TRUE)

      # Handle potential NAs (should be minimal after the grid filter)
      temp_final_dt_1km[is.na(pixel_weight), pixel_weight := 1.0]

      # Weight land use areas by the intersection share (pixel_weight)
      lum_cols <- grep("^LUM", colnames(temp_final_dt_1km), value = TRUE)
      temp_final_dt_1km[, (lum_cols) := lapply(.SD, function(x) x * pixel_weight), .SDcols = lum_cols]

      temp_agg_wide <- temp_final_dt_1km[, lapply(.SD, sum, na.rm = TRUE), by = .(ID, Grouping_Key, X, Y), .SDcols = patterns("^LUM")]

      temp_compiled_long <- melt(temp_agg_wide, id.vars = c("ID", "Grouping_Key", "X", "Y"), variable.name = "LUM_Code", value.name = "area_km2")
      temp_compiled_long[, LUM_Code := as.integer(sub("LUM", "", LUM_Code))]

      # Map raw LUM codes -> model_class via the single class lookup, then aggregate.
      temp_compiled_long[, LUM_Code := as.integer(LUM_Code)]
      area_before_mapping <- sum(temp_compiled_long$area_km2, na.rm = TRUE)
      temp_compiled_long <- merge(temp_compiled_long, class_lookup, by = "LUM_Code", all.x = TRUE)

      # Rock-solid guard: every area-bearing code must map. Fail loud, never drop silently.
      unmapped <- temp_compiled_long[is.na(model_class) & area_km2 > 0]
      if (nrow(unmapped) > 0) {
        bad <- unmapped[, .(area = sum(area_km2)), by = LUM_Code][order(-area)]
        stop(sprintf(
          "Outcome mapping: %d LUM codes carrying %.1f km2 have no class under scheme [%s]. Codes: %s",
          nrow(bad), sum(bad$area), paste(CLASS_COLS, collapse = "+"),
          paste(head(bad$LUM_Code, 12), collapse = ", ")
        ))
      }
      # Zero-area leftovers and NODATA fold into the non-modelled no_choice class.
      temp_compiled_long[is.na(model_class) | model_class == "NODATA", model_class := "no_choice"]
      temp_compiled_long[is.na(focal_class), focal_class := "NODATA"]
      area_after_mapping <- sum(temp_compiled_long$area_km2, na.rm = TRUE)
      if (abs(area_before_mapping - area_after_mapping) > 1) {
        stop(sprintf("Outcome mapping lost %.1f km2 on the LUM_Code join.", area_before_mapping - area_after_mapping))
      }

      temp_y_pixel_yr_wide <- dcast(temp_compiled_long, ID + Grouping_Key ~ model_class,
        value.var = "area_km2", fill = 0, fun.aggregate = sum
      )
      setDT(temp_y_pixel_yr_wide)

      # Final Area Check
      cat(sprintf(
        "  >>> Final Outcome Area Check (Year %d): Total area = %.1f km2\n",
        tp$out_year, sum(temp_y_pixel_yr_wide[, .SD, .SDcols = setdiff(names(temp_y_pixel_yr_wide), c("ID", "Grouping_Key"))], na.rm = TRUE)
      ))

      # Diagnostic: Check pixel area range
      pixel_totals <- rowSums(temp_y_pixel_yr_wide[, .SD, .SDcols = setdiff(names(temp_y_pixel_yr_wide), c("ID", "Grouping_Key"))], na.rm = TRUE)
      cat(sprintf(
        "  >>> Pixel Area Range: min = %.4f | max = %.1f km2\n",
        min(pixel_totals[pixel_totals > 0]), max(pixel_totals)
      ))
    }
  }


  if (!is.na(tp$out_year) && is.null(temp_y_pixel_yr_wide)) {
    # Standard logic for CLC years or fallback for 2018
    path_clc_yr <- file.path(GRIDWORK_DIR, paste0("CLC_Annual_TS_", tp$out_year, ".rds"))
    ls_path_clc <- file.path(LS_DIR, "input/CLC_Annual_TS_2000_2018.rds")

    clc_data_exists <- file.exists(path_clc_yr) || file.exists(ls_path_clc)

    if (clc_data_exists) {
      if (file.exists(path_clc_yr)) {
        temp_clc_yr_raw <- as.data.table(readRDS(path_clc_yr))
      } else {
        temp_clc_yr_raw <- as.data.table(readRDS(ls_path_clc))[year == tp$out_year]
      }

      if (nrow(temp_clc_yr_raw) > 0) {
        temp_clc_yr_raw[, INSPIRE_Europe_buffer_1kmID := as.integer(INSPIRE_Europe_buffer_1kmID)]
        temp_y_1km_yr <- temp_clc_yr_raw[thematic_dt, on = "Code1", nomatch = 0L]

        # Create temp_compiled_long for CLC so focal_df can use the model_class scheme
        temp_compiled_long <- temp_y_1km_yr[grid_map_pixel, on = "INSPIRE_Europe_buffer_1kmID", nomatch = 0L][
          , .(area_km2 = sum(area_km2, na.rm = TRUE)),
          by = .(ID, Grouping_Key, X, Y, model_class)
        ]

        temp_y_pixel_yr_wide <- dcast(
          temp_compiled_long,
          ID + Grouping_Key ~ model_class,
          value.var = "area_km2", fill = 0, fun.aggregate = sum
        )
        setDT(temp_y_pixel_yr_wide)

        # Ensure all expected target_classes are present
        missing_y <- setdiff(target_classes, names(temp_y_pixel_yr_wide))
        if (length(missing_y) > 0) {
          temp_y_pixel_yr_wide[, (missing_y) := 0]
        }

        # Remove columns that are entirely zero/NA
        id_cols <- c("ID", "Grouping_Key")
        val_cols <- setdiff(colnames(temp_y_pixel_yr_wide), id_cols)
        empty_cols <- val_cols[sapply(temp_y_pixel_yr_wide[, ..val_cols], function(x) all(x == 0 | is.na(x)))]
        if (length(empty_cols) > 0) {
          temp_y_pixel_yr_wide[, (empty_cols) := NULL]
        }
      }
    }
  }

  # Placeholder if still null (LUM/CLC missing for this year)
  if (is.null(temp_y_pixel_yr_wide)) {
    cat(sprintf("  WARNING: No LUM/CLC data found for year %d. Creating dummy outcome.\n", tp$out_year))
    temp_y_pixel_yr_wide <- unique(grid_map_pixel[, .(ID, Grouping_Key)])
    for (.cat in target_classes) temp_y_pixel_yr_wide[, (.cat) := 0]
  }

  # 2. COVARIATES (X-Temporal)
  # yield_cols is now globally defined
  climate_cols <- c("spei48", grep("^(bio|pr|tas)", names(prior_1km), value = TRUE))
  climate_cols <- setdiff(climate_cols, c(GHM_VARS, "GHM_Ovr", "GDP", "Pop"))

  gdp_col <- paste0("GDP_", tp$cov_year)
  pop_col <- paste0("Pop_", tp$cov_year)

  climate_cols_raw_expected <- c(paste0("spei48_", tp$cov_year), setdiff(climate_cols, "spei48"))
  climate_cols_raw <- intersect(climate_cols_raw_expected, names(prior_1km))
  if (length(climate_cols_raw) < length(climate_cols_raw_expected)) {
    cat(sprintf("  WARNING: Missing %d climate covariates for year %d!\n", length(climate_cols_raw_expected) - length(climate_cols_raw), tp$cov_year))
  }
  climate_cols_found <- sub(paste0("_", tp$cov_year, "$"), "", climate_cols_raw)

  static_chelsa_cols <- intersect(c("Growing_Degree_Days_gdd5", "Precipitation_Seasonality_bio15", "Annual_Precipitation_bio12"), names(prior_1km))

  temp_x_1km_yr <- prior_1km[, .SD, .SDcols = c("INSPIRE_Europe_buffer_1kmID", "Slope_rad", "Elevation", "Aspect_cos_mean", "Aspect_sin_mean", "allPA_share", "CISI", yield_cols, static_chelsa_cols, climate_cols_raw)]

  # Rename climate columns to drop the year suffix so they are panel-consistent
  if (length(climate_cols_raw) > 0) {
    setnames(temp_x_1km_yr, old = climate_cols_raw, new = climate_cols_found)
  }

  temp_x_1km_yr[, Pop := prior_1km[[pop_col]]]
  temp_x_1km_yr[, GDP := prior_1km[[gdp_col]]]
  for (.g in GHM_VARS) temp_x_1km_yr[[.g]] <- prior_1km[[paste0(.g, "_", tp$cov_year)]]
  cont_vars <- setdiff(colnames(temp_x_1km_yr), "INSPIRE_Europe_buffer_1kmID")
  # Distinguish aggregation type per variable to handle intensive vs extensive quantities
  mean_vars <- setdiff(cont_vars, c("Pop", "GDP"))
  sum_vars <- intersect(cont_vars, c("Pop", "GDP"))

  # Intensive variables (Slope, Elevation, etc.) -> Weighted Mean
  temp_x_mean <- temp_x_1km_yr[grid_map_pixel, on = "INSPIRE_Europe_buffer_1kmID", nomatch = 0L][
    , lapply(.SD, weighted.mean, w = pixel_weight, na.rm = TRUE),
    by = .(ID, Grouping_Key), .SDcols = mean_vars
  ]

  # Extensive variables (Pop, GDP totals) -> Weighted Sum
  # We multiply by pixel_weight to account for 1km cells split across grid units
  temp_x_sum <- temp_x_1km_yr[grid_map_pixel, on = "INSPIRE_Europe_buffer_1kmID", nomatch = 0L][
    , lapply(.SD, function(x) sum(x * pixel_weight, na.rm = TRUE)),
    by = .(ID, Grouping_Key), .SDcols = sum_vars
  ]

  temp_x_pixel_yr_cont <- merge(temp_x_mean, temp_x_sum, by = c("ID", "Grouping_Key"))

  # Re-order and clean
  setcolorder(temp_x_pixel_yr_cont, c("ID", "Grouping_Key", cont_vars))

  # 3. FOCAL (X_input neighborhood shares) -- sourced from the focal_year composition.
  #    focal_year == out_year (or NA) => contemporaneous same-year map (legacy behavior);
  #    focal_year < out_year => TRUE t-1 autoregressive lag: the neighborhood LU-composition
  #    predictors are the pixel's PREDECESSOR state, loaded from the registered earlier-year map.
  coords_pixel <- unique(grid_map_pixel[, .(ID, X, Y, Grouping_Key)])
  focal_year <- if (is.na(tp$focal_year)) tp$out_year else tp$focal_year
  is_lagged_focal <- !is.na(focal_year) && !is.na(tp$out_year) && focal_year != tp$out_year

  if (is_lagged_focal) {
    cat(sprintf("  Calculating LAGGED focal context (X_input) from %s composition (outcome %s)...\n", focal_year, tp$out_year))
    focal_compiled_long <- load_granular_lum(focal_year)
  } else if (exists("temp_compiled_long")) {
    cat("  Calculating focal context from granular thematic map...\n")
    focal_compiled_long <- temp_compiled_long
  } else {
    focal_compiled_long <- NULL
  }

  # Neighborhood composition by the granular focal_class scheme.
  if (!is.null(focal_compiled_long)) {
    temp_fcl_input_granular <- dcast(focal_compiled_long, ID + Grouping_Key + X + Y ~ focal_class, value.var = "area_km2", fill = 0, fun.aggregate = sum)
  } else {
    # Fallback to a single dummy variable so focal_df can run cleanly
    temp_fcl_input_granular <- data.table(ID = coords_pixel$ID, Grouping_Key = coords_pixel$Grouping_Key, dummy_fcl_class = 0)
    temp_fcl_input_granular <- merge(temp_fcl_input_granular, coords_pixel, by = c("ID", "Grouping_Key"), all.x = TRUE)
  }

  # Identify focal classes and EXCLUDE any aggregate/noise classes
  fcl_classes <- setdiff(colnames(temp_fcl_input_granular), c("ID", "Grouping_Key", "X", "Y", "no_choice", "NODATA"))

  # F (focal-source carry): focal_class keeps NO_CHOICE_LU classes (Waterbodies/Wetlands/Natural_other)
  # SEPARATE (derive_model_class L209 sets focal_class BEFORE the no_choice fold), but the model OUTCOME
  # folds them into no_choice -> their AREA cols never reach dat_pixel, so the predict recipe can't recompute
  # focal. Snapshot the granular AREAS (pre-normalization) so we can carry the EXTRA (non-outcome) ones below.
  .focal_src_all <- copy(temp_fcl_input_granular[, c("ID", "Grouping_Key", fcl_classes), with = FALSE])

  # Normalize to shares (ensuring no NAs if rowSum is 0)
  temp_fcl_input_granular[, (fcl_classes) := {
    rs <- rowSums(.SD, na.rm = TRUE)
    lapply(.SD, function(x) ifelse(rs > 0, x / rs, 0))
  }, .SDcols = fcl_classes]

  # OWN t-1 STATE: temp_fcl_input_granular currently holds the PIXEL'S OWN shares at focal_year --
  # the last point before compute_focal_coord spreads them over the neighbourhood. Snapshot, do not
  # recompute.
  prev_state_tier <- NULL
  if (PREV_STATE) {
    if (is.na(tp$focal_year) || is.na(tp$out_year) || tp$focal_year == tp$out_year)
      stop(sprintf(paste0("DRIVER_PREV_STATE=TRUE requires a LAGGED focal (focal_year < out_year); ",
                          "got focal_year=%s, out_year=%s. With a contemporaneous focal the `prev_*` ",
                          "columns ARE the outcome and the model is circular. Set DRIVER_FOCAL_YEARS ",
                          "earlier than DRIVER_MODEL_YEARS."), tp$focal_year, tp$out_year))
    prev_state_tier <- copy(temp_fcl_input_granular[, c("ID", "Grouping_Key", fcl_classes), with = FALSE])
    setnames(prev_state_tier, fcl_classes, paste0("prev_", fcl_classes))
    cat(sprintf("  OWN t-1 STATE: %d prev_* share column(s) from %s composition (outcome %s)\n",
                length(fcl_classes), focal_year, tp$out_year))
  }

  # Diagnostic: Check for NA coordinates (common cause of row loss in focal_df)
  na_coords <- sum(is.na(temp_fcl_input_granular$X) | is.na(temp_fcl_input_granular$Y))
  if (na_coords > 0) {
    cat(sprintf("    WARNING: %d pixels have NA coordinates and will be ignored by focal_df.\n", na_coords))
  }

  # Focal via compute_focal_coord = the EXACT routine the predict RECIPE uses -> one focal code path,
  # 0 fit/predict drift (Option B); resolution-agnostic + fixes terra's arbitrary duplicate-collapse.
  # Replaces the former focal_df/terra 3x3. complete_only=TRUE mirrors .add_focal (incomplete nbhd -> 0).
  if (!exists("compute_focal_coord")) source("experiments/focal/compute_focal_coord.R")
  temp_fcl_input_granular[, .fslice := 1L]
  temp_lu_focal_stats_yr <- compute_focal_coord(temp_fcl_input_granular, classes = fcl_classes,
    coord = c("X", "Y"), slice = ".fslice", res = PIXEL_RES * 1000, complete_only = TRUE, agg = "mean")
  setDT(temp_lu_focal_stats_yr)
  focal_cols_new <- grep("^focal_", names(temp_lu_focal_stats_yr), value = TRUE)
  for (j in focal_cols_new) set(temp_lu_focal_stats_yr, which(is.na(temp_lu_focal_stats_yr[[j]])), j, 0)
  temp_lu_focal_stats_yr[, X_merge := round(X, 0)]
  temp_lu_focal_stats_yr[, Y_merge := round(Y, 0)]
  tier_dat <- merge(temp_y_pixel_yr_wide, temp_x_pixel_yr_cont, by = c("ID", "Grouping_Key"), all.x = TRUE)
  tier_dat <- merge(tier_dat, x_pixel_soil, by = c("ID", "Grouping_Key"), all.x = TRUE)
  tier_dat <- merge(tier_dat, coords_pixel, by = c("ID", "Grouping_Key"), all.x = TRUE)
  # F: carry the EXTRA focal-source areas (focal classes NOT in the model outcome, e.g. Waterbodies/Wetlands/
  # Natural_other) so the recipe can recompute focal. Kept as SOURCE cols -> excluded from X_mat (selection
  # lists) and Y (final_cats_pixel); the recipe's .add_focal reads them alongside the model-class Y areas.
  .extra_fcl <- setdiff(fcl_classes, setdiff(names(temp_y_pixel_yr_wide), c("ID", "Grouping_Key")))
  if (length(.extra_fcl)) tier_dat <- merge(tier_dat, .focal_src_all[, c("ID", "Grouping_Key", .extra_fcl), with = FALSE], by = c("ID", "Grouping_Key"), all.x = TRUE)
  if (!is.null(prev_state_tier)) tier_dat <- merge(tier_dat, prev_state_tier, by = c("ID", "Grouping_Key"), all.x = TRUE)
  tier_dat[, X_merge := round(X, 0)]
  tier_dat[, Y_merge := round(Y, 0)]

  # Final Join for this Tier (Left Joins to preserve all land-use area)
  tier_dat <- merge(tier_dat, temp_lu_focal_stats_yr[, .SD, .SDcols = c("X_merge", "Y_merge", focal_cols_new)], by = c("X_merge", "Y_merge"), all.x = TRUE)
  tier_dat[, c("X_merge", "Y_merge") := NULL]

  # Identify focal columns and fill NAs introduced by the merge
  focal_cols_final <- grep("^focal_", names(tier_dat), value = TRUE)
  for (j in focal_cols_final) {
    set(tier_dat, which(is.na(tier_dat[[j]])), j, 0)
  }

  # Fill gaps for pixels dropped by terra focal (e.g. boundary pixels) with 0
  for (j in focal_cols_final) set(tier_dat, which(is.na(tier_dat[[j]])), j, 0)

  # focal_NODATA = the UNCOVERED share of the neighbourhood, 1 - sum(class shares).
  # WHY A COMPLEMENT, NOT THE OLD BINARY is.na() INDICATOR (which sat here commented out and, placed
  # after the NA->0 fill above, would have been all-zero anyway): it makes the focal block sum to
  # EXACTLY 1 on every row. That is the precondition for the zero-sum reconstruction -- it shifts each
  # class utility by -mean_k * rowSum(block), which an intercept can absorb ONLY if rowSum is constant.
  # Measured on the delivered design, the class columns alone sum to 1 on just 5% of rows (median
  # 0.957, 8.25% entirely zero), so without this column the block cannot be treated as compositional.
  # It also carries real information: partial neighbourhood coverage at coast/border pixels.
  # NOT an LU class -- it is DERIVED, so it must stay out of the recipe's `lu_classes` and be
  # recomputed as the complement at predict time.
  tier_dat[, focal_NODATA := pmax(0, 1 - rowSums(.SD)), .SDcols = focal_cols_final]

  tier_dat[, out_year := tp$out_year]
  tier_dat[, cov_year := tp$cov_year]

  # Area Check after joins
  y_cols <- intersect(names(tier_dat), setdiff(names(temp_y_pixel_yr_wide), c("ID", "Grouping_Key")))
  cat(sprintf(
    "  >>> Tier Data Area Check (Year %d): Total area preserved = %.1f km2\n",
    tp$out_year, sum(rowSums(tier_dat[, .SD, .SDcols = y_cols], na.rm = TRUE))
  ))
  return(tier_dat)
})

dat_pixel <- rbindlist(dat_pixel_list, fill = TRUE)

# Remove the dummy focal column if it was generated
if ("focal_dummy_fcl_class" %in% names(dat_pixel)) {
  dat_pixel[, focal_dummy_fcl_class := NULL]
}

# Explicit Filtering
skewed_vars <- c("GDP", "Pop")   # RAI removed 2026-08-04 (replaced by the GHM threat groups)
static_chelsa_cols <- intersect(c("Growing_Degree_Days_gdd5", "Precipitation_Seasonality_bio15", "Annual_Precipitation_bio12"), names(dat_pixel))

# ---- YIELD REPARAMETERISATION: level + input response (2026-09-17) -----------------------------
# The three yield_index columns are ONE suitability surface entered three times. Measured on the
# GLOBIOM design: pairwise r 0.95-0.99, PC1 carries 98.31% of the variance with near-equal loadings
# (0.575/0.582/0.575), and the ratios are near-constant (Med/Low median 1.46, High/Low 2.12). That
# leaves VIF 312.8 / 107.6 / 82.6 -- i.e. 99.7% of yield_index_Med is explained by the other two, and
# its standard error is inflated ~17.7x.
#
# The consequence is specific: the JOINT yield effect is identified, the SPLIT between the three is
# not. Only ~1.7% of residual variation distinguishes them, so individual coefficients can be large
# and opposite-signed, they move between chains and seeds, and the horseshoe shrinks an essentially
# arbitrary one of them. Prediction is largely unharmed (collinear predictors that keep the same
# relationship out of sample predict fine); ATTRIBUTION is not, and a near-collinear block is a prime
# suspect for the slow mixing seen on the 2026-09-16 GLOBIOM fit (21% of mu above Rhat 1.05).
#
# const_sum_blocks="auto" cannot catch this: it detects blocks summing to a CONSTANT, while these sum
# to a varying total (0 - 3.90, sd 1.07), and mnlogit_rcpp_sym.R states that continuous indices are
# never flagged. So nothing upstream handles it.
#
# Replaced by two quantities, formed AFTER the 1km -> pixel weighted-mean aggregation above (a mean
# of ratios is not the ratio of means):
#   yield_level    = mean(Low, Med, High)   general suitability     (~ PC1, 98.3% of the variance)
#   yield_response = High - Low             response to input intensification
# Measured effect: VIF 312/108/83 -> 9.09 / 8.53. Still correlated (r = 0.926), so improved rather
# than clean -- report them as a pair, not as independent drivers.
#
# OPEN DATA QUESTION: 5,551 of 64,472 GLOBIOM pixels (8.6%) have all three indices exactly 0. If that
# is a missing-data sentinel rather than true zero suitability, both derived columns inherit it and
# the fix is upstream in gridwork, not here.
YIELD_PARAM <- Sys.getenv("DRIVER_YIELD_PARAM", "level_contrast")   # "level_contrast" | "raw"
yield_model_cols <- yield_cols
if (identical(YIELD_PARAM, "level_contrast") && all(yield_cols %in% names(dat_pixel))) {
  dat_pixel[, yield_level    := rowMeans(.SD, na.rm = TRUE), .SDcols = yield_cols]
  dat_pixel[, yield_response := yield_index_High - yield_index_Low]
  dat_pixel[, (yield_cols) := NULL]
  yield_model_cols <- c("yield_level", "yield_response")
  cat(">>> yield: level + response contrast (collinearity VIF ~312/108/83 -> ~9). DRIVER_YIELD_PARAM=raw restores the three raw indices.\n")
} else if (!identical(YIELD_PARAM, "level_contrast")) {
  cat(">>> yield: RAW three indices (DRIVER_YIELD_PARAM=raw) -- near-collinear, VIF ~312/108/83; individual coefficients are not interpretable.\n")
}

spatial_cont_cols <- c("Slope_rad", "Elevation", "Aspect_cos_mean", "Aspect_sin_mean", "allPA_share", GHM_VARS, "CISI", skewed_vars, yield_model_cols, static_chelsa_cols, climate_cols)
soil_cols <- grep("_s[0-9]+$", colnames(dat_pixel), value = TRUE)
focal_cols <- grep("^focal_", colnames(dat_pixel), value = TRUE)

# Fill missing focal values caused by panel joining (e.g. classes present in 2018 but missing in 2000)
for (j in focal_cols) {
  set(dat_pixel, which(is.na(dat_pixel[[j]])), j, 0)
}

# Gate ONLY on covariates the model actually uses as predictors (see X_mat ~L855):
# the continuous predictors + soil + focal. Deliberately NOT climate_cols/yield_cols,
# which are unused and may be entirely NA in a given data version (e.g. spei48) — those
# must never delete the dataset.
model_cont_cols <- c(
  "Slope_rad", "Elevation", "Aspect_cos_mean", "Aspect_sin_mean",
  "allPA_share", GHM_VARS, "CISI", "spei48", skewed_vars, static_chelsa_cols
)
essential_cols <- intersect(c(model_cont_cols, soil_cols, focal_cols), colnames(dat_pixel))
final_cats_pixel <- intersect(target_classes, colnames(dat_pixel))
# Drop EMPTY outcome classes (zero total area) -- mapping ARTIFACTS: defined as a model_class in the thematic
# mapping (target_classes) but no LU code with area maps to them (e.g. Natural_protected). An all-zero softmax
# category only adds a separation burden, so remove it from the outcome. General: catches any such artifact.
.cls_area <- vapply(final_cats_pixel, function(cc) sum(dat_pixel[[cc]], na.rm = TRUE), numeric(1))
.empty_cls <- final_cats_pixel[.cls_area <= 0]
if (length(.empty_cls)) {
  cat(sprintf(">>> Dropping %d EMPTY outcome class(es) [zero area, mapping artifact]: %s\n",
              length(.empty_cls), paste(.empty_cls, collapse = ", ")))
  final_cats_pixel <- setdiff(final_cats_pixel, .empty_cls)
}

area_before_cov_filter <- sum(rowSums(dat_pixel[, ..final_cats_pixel], na.rm = TRUE))
rows_before <- nrow(dat_pixel)

# Analyze which covariates are causing a loss in y
cat("\n--- Covariate Missingness Impact Report (Area Lost) ---\n")
for (.col in essential_cols) {
  if (.col %in% colnames(dat_pixel)) {
    missing_mask <- is.na(dat_pixel[[.col]])
    if (any(missing_mask)) {
      lost_area <- sum(rowSums(dat_pixel[missing_mask, ..final_cats_pixel], na.rm = TRUE))
      if (lost_area > 1) {
        cat(sprintf("  %s: %.1f km2 lost due to NAs\n", .col, lost_area))
      }
    }
  }
}

.naf <- sort(sapply(intersect(essential_cols, names(dat_pixel)), function(.c) mean(is.na(dat_pixel[[.c]]))), decreasing = TRUE)
cat("\n--- essential cols by NA fraction (top 15) ---\n")
print(round(head(.naf, 15), 3))
if (isTRUE(as.logical(Sys.getenv("DRIVER_DIAG_EXIT", "FALSE")))) stop("Early exit for diagnostics completed!")

dat_pixel <- dat_pixel[complete.cases(dat_pixel[, ..essential_cols])]
area_after_cov_filter <- sum(rowSums(dat_pixel[, ..final_cats_pixel], na.rm = TRUE))
rows_after <- nrow(dat_pixel)

# Rock-solid guard: an empty estimation set means an essential covariate is (near-)
# fully NA. Fail loud here with the offenders instead of crashing later in the sampler.
if (rows_after == 0L) {
  stop(sprintf(
    "All %d rows dropped by complete.cases on essential covariates. Worst offenders (NA fraction): %s",
    rows_before,
    paste(sprintf("%s=%.2f", names(head(.naf, 5)), head(.naf, 5)), collapse = ", ")
  ))
}

cat(sprintf(
  "\n>>> Covariate Filtering Check:\n    Dropped %d rows and %.1f area due to missing essential covariates.\n",
  rows_before - rows_after, area_before_cov_filter - area_after_cov_filter
))

# Final Dataset Pixel Area Check (By Year)
cat("\n>>> Final Dataset Pixel Area Check (By Year):\n")
dat_pixel[,
  {
    p_totals <- rowSums(.SD, na.rm = TRUE)
    cat(sprintf(
      "    Year %d: Min Area = %.6f | Max Area = %.1f km2\n",
      .BY$out_year, min(p_totals[p_totals > 0]), max(p_totals)
    ))
    NULL
  },
  by = out_year,
  .SDcols = final_cats_pixel
]

# Final NA Summary
# Final NA Summary
any_nas <- colSums(is.na(dat_pixel))
if (any(any_nas > 0)) {
  cat("\nWARNING: Remaining NAs in dat_pixel (to be zero-filled):\n")
  print(any_nas[any_nas > 0])
}

dat_pixel[is.na(dat_pixel)] <- 0

# Save FULL dat_pixel for debugging and post-estimation analysis (including dummy years)
timestamp_str <- Sys.Date()
# INTERMEDIATES live under output/intermediate/ (2026-09-17). They are large, regenerable, and were
# cluttering the top of output/ alongside the artifacts people actually look for. Consumers glob BOTH
# locations, so files written before the move are still found.
dir.create("output/intermediate", recursive = TRUE, showWarnings = FALSE)
dat_pixel_file <- paste0("output/intermediate/dat_pixel_FULL_", MODEL_LABEL, "_", timestamp_str, ".rds")
saveRDS(dat_pixel, dat_pixel_file)
cat(sprintf("\n>>> Saved FULL dat_pixel (including dummy years) to %s\n", dat_pixel_file))

cat("\n>>> Identified Outcome (Y) Classes:\n")
print(final_cats_pixel)

cat("\n>>> Total Area per Outcome (Y) Class (Outcome Year | Covariate Year):\n")
print(dat_pixel[, lapply(.SD, sum, na.rm = TRUE), by = .(out_year, cov_year), .SDcols = final_cats_pixel])

cat("\n>>> Grand Total Area per Year (All Classes Combined):\n")
print(dat_pixel[, .(Grand_Total_Area = sum(.SD, na.rm = TRUE)), by = .(out_year, cov_year), .SDcols = final_cats_pixel])

dat_pixel <- dat_pixel[rowSums(dat_pixel[, ..final_cats_pixel]) > 0]
cat(sprintf(">>> Filtered for estimation: %d rows with valid outcomes remaining.\n", nrow(dat_pixel)))

# Sliver Filter: Remove pixels with negligible total area to avoid numerical noise
MIN_PIXEL_AREA <- 1 # km2
p_totals_temp <- rowSums(dat_pixel[, ..final_cats_pixel])
rows_before_sliver <- nrow(dat_pixel)
dat_pixel <- dat_pixel[p_totals_temp >= MIN_PIXEL_AREA]
cat(sprintf(
  ">>> Sliver Filter: Dropped %d tiny pixels (< %.4f km2)\n",
  rows_before_sliver - nrow(dat_pixel), MIN_PIXEL_AREA
))

# =========================================================================
# GROUP COVERAGE FILTER
# =========================================================================
# Controls for border cases where large grouping aggregates (e.g. countries)
# have very low spatial coverage in the modeling grid (due to border clipping, missing data, etc.),
# which can result in highly unstable group-specific random effects.
# If coverage is below 15% or total modeled area is below 100 km2, we filter the group out.

MIN_GROUP_COVERAGE <- 0.15 # Minimum 15% coverage of original physical extent
MIN_GROUP_AREA_KM2 <- 100.0 # Minimum 100 km2 of total modeled physical area

cat("\n>>> Applying Group Coverage Filter...\n")
# 1. Calculate original physical area per group from grid_map_pixel
orig_group_area <- grid_map_pixel[, .(orig_area = sum(pixel_weight, na.rm = TRUE)), by = .(Grouping_Key)]

# 3. Calculate modeled physical area per group in dat_pixel (sum of outcome category areas)
modeled_group_area <- dat_pixel[, .(modeled_area = sum(rowSums(.SD, na.rm = TRUE))),
  by = .(Grouping_Key),
  .SDcols = final_cats_pixel
]

# 4. Merge and compute coverage percent
group_coverage <- merge(modeled_group_area, orig_group_area, by = "Grouping_Key", all.y = TRUE)
group_coverage[is.na(modeled_area), modeled_area := 0.0]
group_coverage[, coverage_pct := modeled_area / orig_area]

# Identify low-coverage groups
groups_to_drop <- group_coverage[coverage_pct < MIN_GROUP_COVERAGE | modeled_area < MIN_GROUP_AREA_KM2]

if (nrow(groups_to_drop) > 0) {
  cat("    Dropped low-coverage or small-area grouping aggregates:\n")
  print(groups_to_drop[, .(Grouping_Key,
    modeled_area_km2 = round(modeled_area, 1),
    orig_area_km2 = round(orig_area, 1),
    coverage_pct = round(coverage_pct * 100, 1)
  )])

  # Filter dat_pixel to drop these groups
  rows_before_group_filter <- nrow(dat_pixel)
  dat_pixel <- dat_pixel[!Grouping_Key %in% groups_to_drop$Grouping_Key]
  cat(sprintf(
    "    Group Coverage Filter: Dropped %d rows belonging to low-coverage groups.\n",
    rows_before_group_filter - nrow(dat_pixel)
  ))
} else {
  cat("    All grouping aggregates passed coverage and area thresholds.\n")
}

# Outcome (Y) and Weights
Y_pixel_raw <- as.matrix(dat_pixel[, ..final_cats_pixel])
y_pixel_totals <- rowSums(Y_pixel_raw)

if (ESTIMATE_ON_LEVELS) {
  cat(sprintf("\n>>> Estimating on LEVELS (Factor = %g) <<<\n", LEVEL_SCALE_FACTOR))
  # Keep Y as actual area/counts, scale by factor
  Y_pixel <- Y_pixel_raw * LEVEL_SCALE_FACTOR
  weights_pixel <- rep(1, nrow(Y_pixel)) # weights are baked into Y counts
} else {
  # Standard: Normalized shares with area weights
  Y_pixel <- Y_pixel_raw / y_pixel_totals
  weights_pixel <- rep(1, nrow(Y_pixel)) # y_pixel_totals / mean(y_pixel_totals)
}

if (STABILIZE_SPARSE_Y) {
  cat("\n>>> Applying Global Adaptive Stabilization...\n")
  sparsity_fractions <- colMeans(Y_pixel_raw == 0)
  # Define epsilon for 0-share stabilization
  # We use a fixed tiny constant instead of scaling by mean area to minimize area distortion
  eps_base <- if (ESTIMATE_ON_LEVELS) LEVEL_SCALE_FACTOR * ADAPTIVE_EPSILON_MAX else ADAPTIVE_EPSILON_MAX
  eps_vec <- eps_base * (1 - sparsity_fractions)

  for (k in seq_along(final_cats_pixel)) {
    Y_pixel[, k] <- Y_pixel[, k] + eps_vec[k]
  }
  # Only re-normalize to shares if NOT on levels
  if (!ESTIMATE_ON_LEVELS) {
    Y_pixel <- Y_pixel / rowSums(Y_pixel)
  }
}


# Covariate Matrix
for (.v in skewed_vars) {
  dat_pixel[[.v]] <- log1p(dat_pixel[[.v]])
  setnames(dat_pixel, .v, paste0("log1p_", .v))
}

spatial_cont_cols_trans <- c(
  "Slope_rad", "Elevation", "Aspect_cos_mean", "Aspect_sin_mean", "allPA_share", GHM_VARS, "CISI",
  paste0("log1p_", skewed_vars), static_chelsa_cols, climate_cols, yield_model_cols
)
focal_cov_cols <- grep("^focal_", colnames(dat_pixel), value = TRUE)
prev_cov_cols  <- grep("^prev_",  colnames(dat_pixel), value = TRUE)   # own t-1 state (DRIVER_PREV_STATE)

# # Filter out dummy years for estimation: Keep only rows where at least one outcome category > 0
X_mat_full <- cbind(
  dat_pixel[, !c(final_cats_pixel), with = FALSE],
  intercept = 1
)

# Identify numeric columns
cols <- names(X_mat_full)[sapply(X_mat_full, is.numeric)]

# Replace non-finite values with 0
for (j in cols) {
  set(X_mat_full, i = which(!is.finite(X_mat_full[[j]])), j = j, value = 0)
}
dir.create("output/intermediate", recursive = TRUE, showWarnings = FALSE)
Xmat_pixel_file <- paste0("output/intermediate/Xmat_FULL_", MODEL_LABEL, "_", timestamp_str, ".rds")
saveRDS(X_mat_full, Xmat_pixel_file)

# ── APPROACH B: Zero-Sum Carving via Centered Predictors ──────────────────
# For share variables that sum to 1, centering Z_k = X_k - 1/K makes the
# likelihood algebraically equivalent to projecting beta onto a zero-sum
# space inside the MCMC (Approach B), without modifying the C++ sampler.
#
# Proof: sum(beta_k * Z_k) = sum(beta_k * (X_k - 1/K))
#      = sum(beta_k * X_k) - mean(beta) * sum(X_k)
#      = sum(beta_k * X_k) - mean(beta)         [since sum(X_k) = 1]
#      = sum((beta_k - mean(beta)) * X_k)        [= sum(beta_tilde_k * X_k)]
#
# The prior acts on the unconstrained beta (as Approach B prescribes).
# Post-processing: beta_tilde_k = beta_k - mean(beta) gives zero-sum effects.

dat_pixel_transformed <- copy(dat_pixel)

all_soil_cols <- c()
soil_groups <- list()
for (v in c("OC_TOP", "ROO", "AWC_TOP", "VS")) {
  cols <- sort(grep(paste0("^", v, "_s"), colnames(dat_pixel_transformed), value = TRUE))
  all_soil_cols <- c(all_soil_cols, cols)
  if (length(cols) >= 2) soil_groups[[v]] <- cols
}

X_mat <- cbind(
  intercept = 1,
  as.matrix(dat_pixel_transformed[, intersect(names(dat_pixel_transformed), c(spatial_cont_cols_trans, focal_cov_cols, prev_cov_cols, all_soil_cols)), with = FALSE])
)
X_mat[!is.finite(X_mat)] <- 0

# --- Mundlak device (DRIVER_MUNDLAK) --------------------------------------------------------
# b_c ~ N(mu + gamma' xbar_c, sigma^2), added in the equivalent EXPLICIT form (group means as
# design columns) so the sampler itself is untouched. Purges between-group confounding from beta
# and makes mu the population-averaged effect GIVEN group composition.
# Measured on GLOBIOM init-2000 (8k px, 4 chains, 3 splits, chains pooled): removes the FE-RE
# correlation it targets (median |r| 0.258 -> 0.139, eliminated in 2 of 3 splits) at zero
# predictive cost (-1.1 +/- 4.2 nats). It is an INTERPRETATION fix, not an accuracy fix.
# Inserted HERE (before linear_cols / horseshoe_idx_pixel are derived) so the new columns are
# covered by both; they therefore sit in the horseshoe pool and gamma can be shrunk to zero,
# which degrades gracefully back to plain RE for covariates whose group means carry nothing.
# --- Coordinate covariates (DRIVER_ADD_COORDS) ----------------------------------------------
# A smooth spatial surface: with BART these give the trees something to split space on, absorbing
# residual spatial structure the named drivers miss. dat_pixel X/Y are EPSG:3035 (ETRS89-LAEA)
# METRES, not degrees -- confirmed from their range (2.6-6.0e6 easting, 1.4-5.3e6 northing).
#   DRIVER_COORD_TYPE=lonlat (default) reprojects to EPSG:4326 -> lon/lat degrees, interpretable.
#   DRIVER_COORD_TYPE=laea keeps native metres -- equal-area, so distances are undistorted; for
#   tree SPLITS the two are near-equivalent (both monotone in position), so this is a readability
#   choice more than a modelling one. Scaled to ~unit variance either way so BART's split grid and
#   the linear prior are not dominated by raw magnitude.
if (isTRUE(as.logical(Sys.getenv("DRIVER_ADD_COORDS", "FALSE")))) {
  .ct <- tolower(Sys.getenv("DRIVER_COORD_TYPE", "lonlat"))
  .xy <- data.frame(x = dat_pixel$X, y = dat_pixel$Y)
  if (identical(.ct, "lonlat")) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("DRIVER_COORD_TYPE=lonlat needs the sf package; use DRIVER_COORD_TYPE=laea instead.")
    .pts <- sf::st_transform(sf::st_as_sf(.xy, coords = c("x","y"), crs = 3035), 4326)
    .cc  <- sf::st_coordinates(.pts); .nm <- c("lon", "lat")
  } else { .cc <- as.matrix(.xy); .nm <- c("coord_X", "coord_Y") }
  .cc <- scale(.cc); colnames(.cc) <- .nm
  X_mat <- cbind(X_mat, .cc)
  cat(sprintf(">>> Coordinates ON (%s): +%s  [range %.2f..%.2f, %.2f..%.2f after scaling]\n",
              .ct, paste(.nm, collapse=","), min(.cc[,1]), max(.cc[,1]), min(.cc[,2]), max(.cc[,2])))
}

MUNDLAK     <- isTRUE(as.logical(Sys.getenv("DRIVER_MUNDLAK", "FALSE")))
mundlak_def <- NULL
if (MUNDLAK && use_re) {
  source("codes/mundlak.R")
  .mdl_mc <- trimws(strsplit(Sys.getenv("DRIVER_MUNDLAK_COLS", "CISI,log1p_Pop,log1p_GDP,GHM_HI"), ",")[[1]])
  # DEFAULT INCLUDES THE SLOPE. Measured (3 splits): intercept-only can DISPLACE the between-group
  # correlation into the random slopes rather than remove it -- in one split it drove log1p_Pop's
  # slope correlation 0.340 -> 0.481 (43% -> 72% above the p=.05 threshold). intercept+slope
  # reduced it on BOTH rows in ALL three splits. Names absent from X are filtered out below, so
  # this degrades to "intercept" on a design without log1p_Pop.
  .mdl_rc <- trimws(strsplit(Sys.getenv("DRIVER_MUNDLAK_RE", "intercept,log1p_Pop"), ",")[[1]])
  .mdl_mc <- intersect(.mdl_mc, colnames(X_mat))
  .mdl_rc <- .mdl_rc[.mdl_rc == "intercept" | .mdl_rc %in% colnames(X_mat)]
  if (!length(.mdl_mc) || !length(.mdl_rc)) {
    warning("DRIVER_MUNDLAK set, but no usable DRIVER_MUNDLAK_COLS / DRIVER_MUNDLAK_RE in X_mat; skipping.")
  } else {
    # the SAME group vector the sampler receives (cf. group_idx_vec below), so the stored
    # group_levels line up and apply_mundlak_design() reproduces these columns at predict time
    .mdl_g <- as.integer(as.factor(dat_pixel$Grouping_Key))
    mundlak_def <- build_mundlak_design(X_mat, .mdl_g, mean_cols = .mdl_mc, re_cols = .mdl_rc)
    X_mat <- mundlak_def$X
    cat(sprintf(">>> Mundlak ON: +%d group-mean column(s) [means: %s | RE rows: %s]\n",
                length(mundlak_def$names), paste(.mdl_mc, collapse = ","), paste(.mdl_rc, collapse = ",")))
  }
}

# Detect sum-to-constant (compositional) blocks. Robust (median/MAD) detector so the focal shares
# qualify: they ARE constant-sum (median rowsum 1.0, MAD 0 -> 99.4% of valid pixels sum to exactly
# 1), but a plain mean/sd test was fooled by ~0.6% partial-coverage border pixels (lone 0.5 rowsums)
# -> the old strict 1e-4 detector missed focal entirely, driving the high focal-share Rhats.
# DIAGNOSTIC ONLY -- the samplers run their own detection (const_sum_blocks = "auto"). Printed here
# with verbose=TRUE so the accept/reject reason per block is visible in the driver log, and so a block
# that silently stops qualifying (e.g. focal_NODATA missing) is obvious rather than a mystery Rhat.
const_sum_blocks <- detect_constant_sum_blocks(X_mat, verbose = TRUE)
cat(sprintf("Const-sum blocks the samplers will treat as compositional: %s\n",
            if (length(const_sum_blocks))
              paste(sprintf("%s(%d cols)", names(const_sum_blocks), lengths(const_sum_blocks)), collapse = ", ")
            else "NONE (all blocks fall back to plain contrast coding via the rank guard)"))

# DRIVER-LEVEL drop-one identification: REMOVED 2026-08-12. The samplers already own this -- every
# caller passes const_sum_blocks = "auto" (e.g. .ncut_fit_block), and the sampler's path is the more
# complete one: it drops one column for the fit AND reconstructs the full block as zero-sum WITH the
# intercept compensation (mnlogit_rcpp_sym L3013-3016) that keeps the exported coefficients
# softmax-equivalent to the fit. Doing it here as well left the two mechanisms half-applied: the
# driver removed focal_Forests_HI, after which the remaining 21 focal columns no longer summed to a
# constant, so the sampler correctly declined to treat them as a block and focal ended up plain
# contrast-coded -- with no coefficient at all for the dropped class.
#
# With focal_NODATA now in the design the focal block sums to exactly 1 on every row, so the sampler
# admits it and every focal class gets a symmetric, comparable coefficient. `dropped_blocks` stays an
# empty list; build_fit_object() already no-ops on that.
dropped_blocks <- list()

# Prepare for MNL
# Baseline / reference class for the MNL (configured at the top). A baseline that is
# present in ~all pixels is the most numerically stable (e.g. Natural_unmanaged ~99%).
baseline_class <- BASELINE_CLASS
baseline_idx <- which(final_cats_pixel == baseline_class)
if (length(baseline_idx) != 1L) {
  # Guard: never let a missing/ambiguous baseline silently become integer(0)
  # (that produced `empirical_log_odds[-integer(0)]` -> length-0 prior_mu crash).
  fallback <- names(which.max(colMeans(Y_pixel_raw > 0))) # most-prevalent class
  warning(sprintf(
    "baseline_class '%s' not found among %d outcome classes; falling back to most-prevalent '%s'.",
    baseline_class, length(final_cats_pixel), fallback
  ))
  baseline_class <- fallback
  baseline_idx <- which(final_cats_pixel == baseline_class)
}

linear_cols <- 1:ncol(X_mat)
# Exclude the intercept BY NAME, not by position (2026-09-16). `linear_cols[-1]` assumed the
# intercept is column 1; when it is not, the horseshoe shrinks a class's baseline level and leaves
# some other covariate unshrunk, silently. The BART path below already did this by name. Matches the
# BMLEH_Los1 standard: horseshoe_idx = setdiff(seq_len(ncol(X)), icpt).
horseshoe_idx_pixel <- setdiff(linear_cols, which(colnames(X_mat) == "intercept"))

# --- BART covariate partition + AUTOMATIC support screen (gated; validated non-focal design) ---
# When use_bart: the non-focal continuous/share drivers (topo/climate/soil/socioecon/accessibility/
# yields) go to BART; linear keeps intercept + focal LU-lags (autoregressive -> must stay linear).
# screen_bart_design() auto-prunes near-constant / low-participation-ratio / collinear BART covariates
# (DRIVER_BART_SCREEN, default on; floor DRIVER_BART_MIN_PR_FRAC, default 0.5% effective support).
bart_cols <- NULL
if (isTRUE(use_bart)) {
  .icpt <- which(colnames(X_mat) == "intercept"); .focal <- grep("^focal_", colnames(X_mat))
  # DRIVER_BART_COLS selects WHICH covariates BART gets; everything else stays linear.
  #   ""      (default) = every non-focal, non-intercept column  [previous behaviour, unchanged]
  #   "topo"            = the terrain block (+ coordinates when DRIVER_ADD_COORDS), so the smooth
  #                       spatial/terrain surface is non-parametric while socio/climate stay linear
  #                       and interpretable -- the "linear spec + BART on topo" design.
  #   "a,b,c"           = an explicit comma-separated column list (regex-free, exact names).
  # focal LU-lags and the intercept ALWAYS stay linear: focal is autoregressive and must keep a
  # coefficient, and BART carries its own centred intercept.
  .bsel <- trimws(Sys.getenv("DRIVER_BART_COLS", ""))
  bart_cols <- if (!nzchar(.bsel)) {
    setdiff(seq_len(ncol(X_mat)), c(.icpt, .focal))
  } else if (identical(tolower(.bsel), "topo")) {
    .topo <- c("Slope_rad", "Elevation", "Aspect_cos_mean", "Aspect_sin_mean", "lon", "lat", "coord_X", "coord_Y")
    .hit  <- which(colnames(X_mat) %in% .topo)
    if (!length(.hit)) stop("DRIVER_BART_COLS=topo but none of ", paste(.topo, collapse=", "), " are in the design.")
    setdiff(.hit, c(.icpt, .focal))
  } else {
    .want <- trimws(strsplit(.bsel, ",")[[1]])
    .miss <- setdiff(.want, colnames(X_mat))
    if (length(.miss)) stop("DRIVER_BART_COLS names column(s) not in the design: ", paste(.miss, collapse=", "))
    setdiff(which(colnames(X_mat) %in% .want), c(.icpt, .focal))
  }
  cat(sprintf(">>> BART covariates (%s): %d of %d -> %s\n",
              if (nzchar(.bsel)) .bsel else "all non-focal", length(bart_cols), ncol(X_mat),
              paste(head(colnames(X_mat)[bart_cols], 10), collapse=", ")))
  if (isTRUE(as.logical(Sys.getenv("DRIVER_BART_SCREEN", "TRUE")))) {
    .scr <- screen_bart_design(X_mat, bart_cols,
              min_pr_frac = as.numeric(Sys.getenv("DRIVER_BART_MIN_PR_FRAC", "0.005")))
    if (length(.scr$drop_idx)) cat(sprintf(
      ">>> BART auto-screen: pruned %d low-support covariate(s) -> %s\n",
      length(.scr$drop_idx), paste(colnames(X_mat)[.scr$drop_idx], collapse = ", ")))
    bart_cols <- .scr$keep_idx
  }
  linear_cols <- setdiff(seq_len(ncol(X_mat)), bart_cols)        # NO overlap (focal+intercept+screened -> linear)
  horseshoe_idx_pixel <- setdiff(linear_cols, .icpt)
  cat(sprintf(">>> BART partition: %d BART covariates | %d linear (intercept + %d focal + %d screened-out)\n",
      length(bart_cols), length(linear_cols), length(.focal), length(linear_cols) - 1L - length(.focal)))
}

cat(sprintf("\nEstimation Sample: N = %d | Baseline: %s\n", nrow(Y_pixel), baseline_class))

# Helper: chains are combined using the comprehensive version in codes/mnl_aux_func.R

# =========================================================================
# 6. MODEL ESTIMATION
# =========================================================================
cat(sprintf("\nStep 3: Fitting %s (%d chains in parallel)...\n", SAMPLER, N_CHAINS))

# Prepare Grouping for Random Effects
if (use_re) {
  re_groups <- as.factor(dat_pixel$Grouping_Key)
  group_idx_vec <- as.integer(re_groups)
  re_group_names <- levels(re_groups)
  cat(sprintf(">>> RE Model: Detected %d groups using key: %s\n", length(re_group_names), RE_GROUP_COL))
} else {
  group_idx_vec <- NULL
}

# Dump the assembled REAL model inputs so the BART-hybrid validation can iterate standalone
# (full LUM classes + complete driver set), without re-running the heavy data build each time.
if (isTRUE(as.logical(Sys.getenv("DRIVER_DUMP_INPUTS", "FALSE")))) {
  # Project runs need their own dump: a BMLEH design and a GLOBIOM design cannot share one path
  # without one silently overwriting the other.
  # DESIGN DUMPS live under output/designs/ (2026-09-17). Everything that consumes a dump takes an
  # explicit path, so the only thing this default changes is where an unqualified build lands.
  .dump_path <- Sys.getenv("DRIVER_DUMP_PATH", "output/designs/pixel_model_inputs.rds")
  dir.create(dirname(.dump_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(X_mat = X_mat, Y_pixel = Y_pixel, weights_pixel = weights_pixel,
               # WHAT THIS DESIGN IS. A dump previously recorded nothing about the scheme that
               # produced it, so anything choosing between several staged dumps could only read the
               # FILENAME -- and a filename is a convention, not a fact. entrypoint.sh has to select
               # by PROJECT/CLASSIFICATION substring for exactly this reason; with these fields a
               # consumer can VERIFY the dump it picked instead of trusting what it is called.
               class_cols     = CLASS_COLS,      # DRIVER_CLASS_COLS, verbatim
               class_scheme   = CLASS_SCHEME,    # GLOBIOM / BMLEH / AGMIP / BIOCLIMA
               model_label    = MODEL_LABEL,
               pixel_res_km   = PIXEL_RES,
               pixel_intersect = PIXEL_INTERSECT_COL,
               built_at       = Sys.time(),
               mundlak_def = mundlak_def,   # NULL unless DRIVER_MUNDLAK; needed to rebuild the
                                            # same group-mean columns on new data at predict time
               group_idx_vec = group_idx_vec, re_group_names = if (use_re) re_group_names else NULL,
               # WHICH COLUMN produced those labels, and which produced geo_region. Without this a
               # consumer has to guess the keying from the label shapes, and guessing wrong is silent:
               # matching a CAPRI-keyed fit against GLOB_country names matches nothing, so every
               # country falls to the pooled fallback and the artifact looks complete with no country
               # effects in it at all.
               re_group_col = RE_GROUP_COL,
               re_group_sliced = identical(RE_GROUP_COL, "CAPRI_NUTS"),   # sliced to 2 chars for the RE key
               geo_col = if (exists("GEO_REGION_COL")) GEO_REGION_COL else NA_character_,
               baseline_class = baseline_class, col_names = colnames(X_mat),
               focal_cov_cols = focal_cov_cols, spatial_cont_cols_trans = spatial_cont_cols_trans,
               # WHICH yield parameterisation built this design. A dump with yield_level/
               # yield_response and one with the three raw indices are not interchangeable, and a
               # consumer cannot tell them apart from the column names alone without knowing the
               # convention -- so record it rather than infer it.
               yield_param = YIELD_PARAM, yield_model_cols = yield_model_cols,
               # coords/ID aligned 1:1 with X_mat rows (dat_pixel row order) -> enables spatial tile-block CV
               coord_X = if ("X" %in% names(dat_pixel)) dat_pixel$X else NULL,
               coord_Y = if ("Y" %in% names(dat_pixel)) dat_pixel$Y else NULL,
               row_ID  = if ("ID" %in% names(dat_pixel)) dat_pixel$ID else NULL,
               out_year = if ("out_year" %in% names(dat_pixel)) dat_pixel$out_year else NULL,
               # Hand-curated nesting, carried alongside the classes so the nested runner never has to
               # locate the mapping CSV itself. Two columns: the LEAF class (CLASS_COLS[1]) and its
               # ancestor path e.g. "Cropland/permanent" (blank = root leaf). Absent column -> NULL,
               # and the runner falls back to the naming-prefix scheme exactly as before.
               class_nest = local({
                 .nc <- Sys.getenv("DRIVER_NEST_COL", paste0(CLASS_COLS[1], "_nest"))
                 if (.nc %in% names(mapping_thematic) && CLASS_COLS[1] %in% names(mapping_thematic)) {
                   cat(sprintf(">>> curated nests: %s -> %s\n", CLASS_COLS[1], .nc))
                   unique(mapping_thematic[, c(CLASS_COLS[1], .nc), with = FALSE])
                 } else {
                   # NEVER fall back silently: a typo in the column name would otherwise look like a
                   # successful run that quietly used the naming-prefix tree instead of the curated one.
                   .near <- grep("nest|nets", names(mapping_thematic), value = TRUE, ignore.case = TRUE)
                   warning(sprintf("no curated nest column '%s' in the mapping -> the nested runner will fall back to a NAMING-PREFIX tree.%s",
                           .nc, if (length(.near)) paste0(" Did you mean: ", paste(.near, collapse = ", "), "?") else ""),
                           call. = FALSE, immediate. = TRUE)
                   NULL
                 }
               })),
          .dump_path)
  cat(sprintf(">>> DRIVER_DUMP_INPUTS: saved pixel inputs (X %dx%d, %d LUM classes, %d groups) -> %s\n",
              nrow(X_mat), ncol(X_mat), ncol(Y_pixel), length(unique(group_idx_vec)), .dump_path))

  # ---- GeoTIFF of the LU layer that actually enters the model ---------------------------------
  # Not a re-derivation from source: this rasterises Y_pixel_raw, the exact per-pixel class AREAS
  # the sampler is handed, so what you open in QGIS is what the model sees -- crop splits, target
  # cascade, coverage filtering and all. Areas in km2, one band per class, plus a `dominant` band
  # coding the argmax class (band levels written as the raster's categories).
  if (isTRUE(as.logical(Sys.getenv("DRIVER_DUMP_TIF", "FALSE")))) {
    if (!requireNamespace("terra", quietly = TRUE)) {
      warning("DRIVER_DUMP_TIF set but the terra package is not installed -- skipping.", call. = FALSE)
    } else if (!all(c("X","Y") %in% names(dat_pixel))) {
      warning("DRIVER_DUMP_TIF set but dat_pixel carries no X/Y -- skipping.", call. = FALSE)
    } else {
      .tif_base <- Sys.getenv("DRIVER_DUMP_TIF_PATH", sub("\\.rds$", "", .dump_path))
      .yr <- if ("out_year" %in% names(dat_pixel)) as.character(dat_pixel$out_year) else rep("all", nrow(dat_pixel))
      for (.y in unique(.yr)) {
        .k <- which(.yr == .y)
        # AGGREGATE FRAGMENTS TO THE GRID CELL FIRST.
        # PIXEL_INTERSECT_COL splits one grid cell into one row PER REGION, and those rows all
        # carry the PARENT cell's X/Y. terra::rast(type = "xyz") keeps a single row per
        # coordinate, so rasterising the rows directly stamped ONE region's composition across
        # the whole cell and silently dropped the others -- on the BMLEH 10 km design, 19,791 of
        # 66,861 rows (29.6%), with up to 9 fragments in a cell. The visible symptom was crisp
        # borders along the INTERSECT geometry (NUTS3) that the model itself never produces:
        # the design keeps the fragments separate, only this picture collapsed them.
        # Y_pixel_raw is class AREA, so the cell's composition is the SUM over its fragments --
        # exact, not an approximation. `dominant` is then taken from the AGGREGATE, since an
        # argmax cannot be averaged after the fact.
        .agg <- cbind(data.table(x = dat_pixel$X[.k], y = dat_pixel$Y[.k]),
                      as.data.table(Y_pixel_raw[.k, , drop = FALSE]))[
                        , lapply(.SD, sum), by = .(x, y)]
        .A <- as.matrix(.agg[, !c("x", "y"), with = FALSE])
        if (nrow(.agg) < length(.k))
          cat(sprintf(">>> DRIVER_DUMP_TIF: merged %d region fragments into %d grid cells\n",
                      length(.k), nrow(.agg)))
        .dom <- max.col(.A, ties.method = "first")
        .dom[rowSums(.A) <= 0] <- NA_integer_
        .xyz <- data.frame(x = .agg$x, y = .agg$y, .A, dominant = .dom, check.names = FALSE)
        .r <- terra::rast(.xyz, type = "xyz", crs = "EPSG:3035")   # dat_pixel X/Y are ETRS89-LAEA
        levels(.r[["dominant"]]) <- data.frame(value = seq_along(final_cats_pixel), class = final_cats_pixel)
        .out <- sprintf("%s_%s.tif", .tif_base, .y)
        terra::writeRaster(.r, .out, overwrite = TRUE,
                           gdal = c("COMPRESS=DEFLATE", "TILED=YES", "PREDICTOR=2"))
        cat(sprintf(">>> DRIVER_DUMP_TIF: %s  (%d pixels, %dm, %d class bands + dominant)\n",
                    .out, length(.k), PIXEL_RES * 1000L, length(final_cats_pixel)))
      }
    }
  }

  # Early exit after the inputs dump: data-prep only (minutes), no flat fit. The dumped
  # output/pixel_model_inputs.rds feeds run_nested_cut.R (FROM_INPUTS mode) or any standalone model.
  if (isTRUE(as.logical(Sys.getenv("DRIVER_DUMP_EXIT", "FALSE")))) {
    cat(">>> DRIVER_DUMP_EXIT: inputs dumped; exiting before model fit.\n"); quit(save = "no", status = 0)
  }

  # ---- F: build + self-check the RECIPE so fit & predict share one assembly (predict_prior from raw). ----
  # SELF-DIAGNOSING: recipe_selfcheck reports max|dX| vs this exact X_mat. nonfocal ~0 => design reproduced;
  # focal ~4e-4 => the KNOWN terra(focal_df at fit) vs compute_focal_coord(at predict) drift (see predict-
  # wrapper memory). Errors visibly (tryCatch) if dat_pixel lacks the granular LU-area columns focal needs.
  suppressWarnings(tryCatch({
    source("codes/prior_model_predict.R"); source("experiments/focal/compute_focal_coord.R")
    .recipe <- list(
      feature_cols = setdiff(colnames(X_mat), c("intercept", focal_cov_cols)),
      # focal_NODATA is DERIVED (1 - sum of the class shares), not an LU layer -- it must stay out of
      # lu_classes or compute_focal_coord would look for a "NODATA" class; add_nodata rebuilds it.
      lu_classes = sub("^focal_", "", setdiff(focal_cov_cols, "focal_NODATA")),
      add_nodata = "focal_NODATA" %in% focal_cov_cols, outcome_classes = colnames(Y_pixel),
      coord = c("X", "Y"), slice = "out_year", group_col = "Grouping_Key",
      transforms = list(list(fn = "log1p", cols = skewed_vars, prefix = "log1p_")),
      res = PIXEL_RES * 1000, col_order = colnames(X_mat), focal_cols = focal_cov_cols,
      group_levels_factor = if (use_re) re_group_names else NULL,
      group_levels_appear = if (use_re) unique(group_idx_vec) else NULL)
    .sc <- recipe_selfcheck(.recipe, dat_pixel, X_mat)
    cat(sprintf(">>> RECIPE self-check: ncol match=%s | nonfocal max|dX|=%.2e | FOCAL: mean=%.2e p99=%.2e max=%.2e frac>0.01=%.4f (terra-vs-coord; small mean + few boundary spikes = benign)\n",
                .sc$ncol_match, .sc$max_dX_nonfocal, .sc$mean_dX_focal, .sc$p99_dX_focal, .sc$max_dX_focal, .sc$frac_focal_gt01))
    saveRDS(.recipe, "output/pixel_model_recipe.rds")
    cat(">>> saved recipe -> output/pixel_model_recipe.rds (build_prior_model(fit, recipe, ...) after fitting to ship a self-contained predictor)\n")
  }, error = function(e) cat(sprintf(">>> RECIPE build skipped (%s) -- ensure dat_pixel carries the granular LU-area cols for focal.\n", conditionMessage(e)))))
  quit(save = "no")
}

# Determine steps for progress bar and identify hot-starting chains
model_disk_path <- file.path("output/saved_model_outputs", paste0(MODEL_LABEL, if (use_re) "_RE_" else "_pooled_", RE_GROUP_COL))
# Test runs start clean so stale batches/state never poison hot-start or recovery
# (this is what the old REvNN bumping was for). Production persists to allow resume.
if (RUN_MODE == "test" && dir.exists(model_disk_path)) {
  .stale <- list.files(model_disk_path, full.names = TRUE)
  if (length(.stale) > 0) {
    cat(sprintf(">>> RUN_MODE=test: clearing %d stale file(s) in %s for a fresh run\n", length(.stale), model_disk_path))
    unlink(.stale, recursive = TRUE)
  }
}

# Copy every artifact tagged with this test MODEL_LABEL to its "_production" twin
# (saved-model dir + flat output files + plots), overwriting the production bucket.
promote_test_to_production <- function() {
  prod_label <- sub("_test$", "_production", MODEL_LABEL)
  # 1. saved-model state directory
  if (dir.exists(model_disk_path)) {
    dst_dir <- sub(MODEL_LABEL, prod_label, model_disk_path, fixed = TRUE)
    if (dir.exists(dst_dir)) unlink(dst_dir, recursive = TRUE)
    dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(list.files(model_disk_path, full.names = TRUE), dst_dir, recursive = TRUE, overwrite = TRUE)
  }
  # 2. flat output files + plots whose names carry this test MODEL_LABEL
  cand <- list.files(c("output", "output/plots"), full.names = TRUE)
  cand <- cand[!dir.exists(cand) & grepl(MODEL_LABEL, basename(cand), fixed = TRUE)]
  n <- 0L
  for (f in cand) {
    dst <- file.path(dirname(f), sub(MODEL_LABEL, prod_label, basename(f), fixed = TRUE))
    if (file.copy(f, dst, overwrite = TRUE)) n <- n + 1L
  }
  cat(sprintf(">>> Promoted to production: saved-model dir + %d file(s)  [%s -> %s]\n", n, MODEL_LABEL, prod_label))
}

# -------------------------------------------------------------------------
# Standalone path for the logistic-normal (LNM) and CLR Gaussian samplers.
# They return per-draw SYMMETRIC (baseline-free) coefficient arrays
# [cov, cat, draws] directly, so we fit, summarize, plot and save here and
# skip the MNL-specific recovery/zero-sum/convergence downstream entirely.
# LNM is multinomial over area/counts; CLR is Gaussian over shares.
# -------------------------------------------------------------------------
total_steps <- N_CHAINS * niter

# Unified dispatch: EVERY model variant is fit by the same future_lapply over chains,
# called with one shared `common_args` (matched to the MNL setup) plus the variant's own
# `sampler_extra`. Each sampler takes the exact same input data (X_mat, Y_pixel,
# group_idx_vec, ...) and handles its own internals: CLR normalises Y to shares inside the
# function; LNM uses rowSums(Y) as N; both ignore the MNL-only args via their `...`.
sampler_extra <- switch(SAMPLER,
  mnlogit_rcpp = list(use_horseshoe = TRUE, equation_specific_hs = TRUE, support_prior_strength = 2, re_asis = FALSE), # RE-scale ASIS implemented in core.cpp but DISABLED: validates on the symmetric core (cor 0.994) but the baseline-coded base core shows a residual bias (cor 0.87) not yet root-caused. Production uses mnlogit_rcpp_sym.
  mnlogit_rcpp_sym = list(
    bart_symmetric = TRUE, bart_base = 0.90, bart_power = 3.0, bart_k = 2.0,   # validated BART: symmetric (CLR) + tight depth prior (surface cor 0.74 vs 0.56); ignored when use_bart=FALSE
    # PER-CLASS leaf shrinkage k_j = bart_k * sqrt(p_max/p_j). dbarts rescales the response every
    # setResponse, so the leaf prior is relative to the WORKING-RESPONSE RANGE -- and Polya-Gamma
    # standardises that range across classes (measured 31.4-74.9 over 26 ensembles). Every class
    # therefore gets the same prior latitude while signal-to-noise varies with prevalence, so rare
    # classes fill it with noise. Measured held-out on 27 GLOBIOM classes, WITH bart_symmetric,
    # against a linear-topography arm: total +38.4 -> +75.4, rare classes -51.6 -> -9.0, common
    # +90.0 -> +84.4. The two fixes are complementary: symmetric alone helps only the common
    # classes, k~prev alone only the rare ones. DRIVER_BART_K_PREV=FALSE restores flat k.
    bart_k_prevalence = isTRUE(as.logical(Sys.getenv("DRIVER_BART_K_PREV", "TRUE"))),
    bart_k_prev_cap   = as.numeric(Sys.getenv("DRIVER_BART_K_PREV_CAP", "10")),
    store_bart_trees = TRUE, do_slim_trees = TRUE,                              # calibrated slim trees -> out-of-sample prediction via reconstruct_bart_f_mean
    use_horseshoe = TRUE, equation_specific_hs = FALSE, symmetric_hs = TRUE,   # HS ON 2026-08-13:
    # +14.2 (Forests) / +15.5 (Pasture) held-out with the RE slab estimated. estimate_c2 stays FALSE
    # (FE slab only weakly identified; estimating it gains nothing).
    use_wls_init = adaptive_use_wls_init, use_precision_hs = adaptive_use_precision_hs,
    use_spike_slab = adaptive_use_spike_slab, store_delta = adaptive_store_delta,
    use_car = adaptive_use_car, country_adjacency = NULL, car_rho = adaptive_car_rho,
    estimate_c2 = TRUE, slab_df = 20, slab_s2 = 4,   # FE slab estimated with a TIGHT prior (requested).
                                 # Costs ~1.3 held-out LL vs FIXED on Pasture; kept for consistency with
                                 # the RE slab, not because it improves fit. See nested_cut.R for detail.
    estimate_slab_c2 = TRUE, collapse_slab_c2 = 4, slab_df_re = 10,  # full Bayes on the RE slab: the
                                 # cap is estimated, not fixed (beats the best fixed value on held-out;
                                 # c2 identified -- starts 4 and 100 both converge to ~4.1).
    const_sum_blocks = "auto",   # 2026-08-13: the DRIVER-LEVEL drop-one was removed (the samplers own
                                 # this, and their path also does the zero-sum reconstruction WITH the
                                 # intercept compensation). Without passing this the flat driver got NO
                                 # const-sum treatment at all -- the sampler default is NULL, so the
                                 # redundant columns fell through to the rank guard instead. Matches
                                 # .ncut_fit_block.
    # SUPPORT STRENGTHS SET EXPLICITLY (2026-09-16). These used to be inherited from the single
    # `support_prior_strength = 2`, which nobody chose: it drove BOTH the FE horseshoe multiplier and
    # the RE sparse-group gate. On the FE side under symmetric_hs that factor IS live (it scales c_v,
    # mnlogit_rcpp_sym.R ~3638), and nested_cut measures (n/PR)^2 there as a median 106x / max 2.8e5x
    # precision inflation costing ~125 nats. BMLEH_Los1 and nested_cut both set FE to 0; only this
    # caller inherited 2.
    # The RE side is a SEPARATE decision: 1 is nested_cut's reasoned value ("re: 1 always" --
    # deterministic, identified, stops a country with no within-country variation in a covariate
    # contributing a free RE for it). BMLEH does NOT pass re_support_strength at all, so it runs at
    # the sampler default 0 with the gate OFF -- an inheritance, not a choice, and the same accident
    # being fixed here. Keeping 1; see docs/sampler_model_specification.md sec. 8.1b.
    fe_support_strength = 0, re_support_strength = 1,
    # RE-scale ASIS ON: see codes/nested_cut.R for the controlled comparison and the counter-
    # measurement. Standard as of 2026-09-16.
    re_asis = TRUE, re_regularize = TRUE, init_jitter = 0.1, # const-sum blocks now handled BY THE SAMPLER (see above) -> not passed to any sampler. participation-ratio RE prior + per-chain overdispersed init. re_regularize=TRUE: RE variance = PROPER regularised horseshoe via a 1-D SLICE on log tau_raw (NOT the old post-draw cap), effective precision tau_eff = tau_raw + 1/collapse_slab_c2 (default c2=100 -> SD<=10) used consistently in the C++ standard (update_re_precision_hc_sym) so the heavy half-Cauchy tail can't run off; stable + exact (exp9 + the oracle gate at c2=100). Also support-CONSISTENT variance (C++ weights ss by 1/rs^2, matching the support-scaled RE draw - fixes a pre-existing inconsistency). RE-scale ASIS OFF: measured to HURT MNL RE-var ESS (73->27). collapse_re_var OFF but REWORKED + validated standalone (per-covariate zero-sum-POOLED slice collapse, exp7/8/9) -> enable (collapse_re_var=TRUE + collapse_re_var_validated=TRUE) only AFTER experiments/mixing/gate_collapse_production.R PASSES on the real design. TUNE collapse_slab_c2 to the largest plausible group-slope SD (raise if sigma~10-12 cells are real; lower to regularise them harder). Retained draws now default to 12000 (niter 16000 - nburn 4000).
    separation_as_prior = adaptive_separation_as_prior, separation_soft = adaptive_separation_soft,
    sep_overlap_hard = adaptive_sep_overlap_hard, sep_overlap_neutral = adaptive_sep_overlap_neutral
  ),
  lnm_gibbs = list(symmetric = TRUE, use_cpp_core = TRUE, re_screen = TRUE, support_prior_strength = 2, re_asis = TRUE, init_jitter = 0.1), # baseline-free symcore + RE-identifiability screen (A2+A3) + participation-ratio RE prior + per-category RE-scale ASIS interweave (guarded against PG overflow at full slopes) + per-chain jittered init
  mvclr_gibbs = list(zero_prevalence_scaling = TRUE, re_asis = TRUE, re_screen = TRUE, support_prior_strength = 2, init_jitter = 0.1) # RE-ASIS + screen via IN-DRAW precision pin (post-hoc pin runs CLR away; in-draw is stable) + support-aware prior (participation ratio n/PR) + per-chain jittered init
)
sampler_extra$use_tempering <- TRUE
sampler_extra$tempering_T0 <- 0.05
sampler_file <- file.path("codes", paste0(SAMPLER, ".R"))
# Every sampler takes IDENTICAL inputs, including re_idx = random slopes on all covariates.
# (LNM/CLR ridge their per-group RE precision internally so full random slopes stay PD.)
# RE structure: attach random effects only to covariates whose EFFECT on land use plausibly varies BY COUNTRY
# (policy/development-driven), + the country INTERCEPT (baseline), instead of random slopes on ALL covariates.
# Default RE set: SOCIO-ECONOMIC (GDP/Pop/human-modification GHM_HI + GHM_TI/infrastructure CISI) +
# PROTECTED-AREA share (allPA_share: LU-constraint varies by national enforcement) + TERRAIN (Slope_rad,
# Elevation: the marginal-land/abandonment threshold varies by national agricultural intensity). Climate/soil
# (biophysical, ~universal) and the focal block (mixing) stay POOLED. Far fewer REs (better mixing + parsimony)
# while keeping the essential country-baseline intercept RE (carries most of the skill). re_regularize self-
# prunes any RE whose effect doesn't actually vary by country. Env DRIVER_RE_VARS overrides; "all" = full slopes.
# re_idx MUST be a subset of linear_cols (BART covariates are pooled/aspatial -> if a var sits in BART, auto-excluded).
.re_spec <- Sys.getenv("DRIVER_RE_VARS", "log1p_GDP,log1p_Pop,GHM_HI,GHM_TI,CISI,allPA_share,Slope_rad,Elevation")
if (identical(tolower(trimws(.re_spec)), "all")) {
  re_idx_use <- linear_cols
} else {
  .socio_re <- which(colnames(X_mat) %in% trimws(strsplit(.re_spec, ",")[[1]]))
  re_idx_use <- sort(intersect(c(which(colnames(X_mat) == "intercept"), .socio_re), linear_cols))
  cat(sprintf(">>> RE on %d covariate(s) [intercept + socio-economic]: %s\n",
              length(re_idx_use), paste(colnames(X_mat)[re_idx_use], collapse = ", ")))
}

res_list <- with_progress({
  p_bar <- progressor(steps = total_steps)
  future_lapply(seq_len(N_CHAINS), function(i) {
    source("codes/mnl_aux_func.R") # shared helpers (%||%, screen_re_design, ...) — needed in the worker;
    source(sampler_file) # MNL samplers also source it internally, CLR/LNM do not.
    fs_file <- file.path(model_disk_path, sprintf("final_state_chain_%d.qs", i))
    chain_init_state <- NULL
    if (file.exists(fs_file)) {
      cat(sprintf("Chain %d: Found final state on disk. Hot-starting...\n", i))
      chain_init_state <- qs2::qs_read(fs_file)
    }

    common_args <- list(
      X = X_mat, Y = Y_pixel, intercept = FALSE, # intercept handled in X
      empirical_intercept_prior = TRUE, baseline = baseline_idx,
      use_re = use_re, group_idx = group_idx_vec, re_idx = re_idx_use,
      use_bart = use_bart, linear_idx = linear_cols, y_weight = weights_pixel,
      bart_idx = if (use_bart) bart_cols else NULL, n_trees_bart = n_trees_bart,
      save_posterior_to_disk = TRUE, disk_path = model_disk_path,
      niter = niter, nburn = nburn, thin = thin_keep, standardize = TRUE, method = c("center", "scale"),
      calc_loo = FALSE, horseshoe_idx = horseshoe_idx_pixel,
      p0_mu = round(0.6 * length(horseshoe_idx_pixel)), gamma_matched_pg = FALSE,
      re_scale_A = 1, chain_id = i, progress_cb = function(...) p_bar(),
      init_state = chain_init_state
    )
    do.call(get(SAMPLER), c(common_args, sampler_extra))
  }, future.seed = TRUE, future.stdout = NA)
})

# -------------------------------------------------------------------------
# LNM / CLR post-processing: each chain returns a per-draw SYMMETRIC (baseline-free)
# coefficient array [cov, cat, draws], so we summarize, plot and save here and skip the
# MNL-specific recovery / zero-sum / convergence downstream.
# -------------------------------------------------------------------------
if (SAMPLER %in% c("lnm_gibbs", "mvclr_gibbs")) {
  getd <- function(f) if (!is.null(f$postB)) f$postB else f$Bsym_draws
  if (is.null(getd(res_list[[1]]))) stop(sprintf("%s returned no per-draw coefficient array.", SAMPLER))
  d1 <- getd(res_list[[1]])
  P <- dim(d1)[1]
  K <- dim(d1)[2]
  nd <- dim(d1)[3]
  cov_names_processed <- dimnames(d1)[[1]]
  cat_all <- dimnames(d1)[[2]]
  arr <- array(0, c(P, K, nd, N_CHAINS), dimnames = list(cov_names_processed, cat_all, NULL, NULL))
  for (ch in seq_len(N_CHAINS)) arr[, , , ch] <- getd(res_list[[ch]])

  B_mean <- apply(arr, c(1, 2), mean)
  # Per-coefficient Gelman Rhat across chains (NA for a single chain).
  rhat <- if (N_CHAINS > 1) {
    apply(arr, c(1, 2), function(M) {
      n <- nrow(M)
      W <- mean(apply(M, 2, var))
      B <- n * var(colMeans(M))
      if (!is.finite(W) || W <= 0) NA_real_ else sqrt(((n - 1) / n * W + B / n) / W)
    })
  } else {
    matrix(NA_real_, P, K)
  }
  cat(sprintf(
    ">>> %s fit: max|B|=%.2f | median Rhat=%s | worst Rhat=%s\n", SAMPLER, max(abs(B_mean)),
    ifelse(all(is.na(rhat)), "NA (1 chain)", sprintf("%.3f", median(rhat, na.rm = TRUE))),
    ifelse(all(is.na(rhat)), "NA", sprintf("%.3f", max(rhat, na.rm = TRUE)))
  ))

  saveRDS(
    build_fit_object(SAMPLER, arr, cov_names_processed, cat_all, baseline_class, res_list[[1]]$Sigma_mean, dropped_blocks),
    paste0("output/", SAMPLER, "_pixel_model_fit_", MODEL_LABEL, "_", timestamp_str, ".rds")
  )

  cat("\n>>> Generating Parameter Heatplot...\n")
  # Display scale (rule B) -- see codes/mnl_aux_func.R. X_mat is the fitted design, so the SDs are
  # the ones the sampler actually standardised on.
  .sd_disp <- display_sds(X_mat, cov_names_processed)
  understandable_labels <- format_mnl_labels(cov_names_processed, .sd_disp)
  cat_sym <- c(setdiff(cat_all, baseline_class), intersect(cat_all, baseline_class))
  reordered <- reorder_mnl_covariates(cov_names_processed, cat_names = cat_sym)
  ord <- match(reordered, cov_names_processed)
  flat_draws <- array(arr[ord, cat_sym, , , drop = FALSE], c(P, K, nd * N_CHAINS),
    dimnames = list(cov_names_processed[ord], cat_sym, NULL)
  )
  p_heat <- MNL_parameter_heatplot(
    draws = flat_draws, cov_names = understandable_labels[ord], cat_names = cat_sym,
    exclude_intercept = TRUE,
    # Was unscaled ("Raw symmetric coefficients"), which is not comparable across drivers: raw units
    # span four orders of magnitude here, so climate rendered as a hairline. Same rule B as the MNL path.
    scaling_vec = .sd_disp[reordered],
    title = paste0(MODEL_LABEL, " ", SAMPLER, " Coefficient Heatplot (baseline-free symmetric)"),
    subtitle = "Tile size = Credibility, Color = Median. Symmetric (zero-sum) coefficients."
  ) + theme(axis.text.y = element_text(size = 9, hjust = 1))
  dir.create("output/plots", recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0("output/plots/", SAMPLER, "_pixel_", MODEL_LABEL, "_heatplot.png"), p_heat, width = 12, height = 20)
  p_heat
  # --- Sampler-agnostic downstream artifacts (mirror the MNL pipeline) -------------------
  # (b1) Convergence check: pooled per-coefficient Gelman Rhat (group-level Rhat needs per-group
  #      draws, which LNM/CLR do not store — pooled only).
  conv_df <- data.table::as.data.table(expand.grid(ks = cov_names_processed, lu_to = cat_all, stringsAsFactors = FALSE))
  conv_df[, `:=`(rhat = as.vector(rhat), B_mean = as.vector(B_mean))]
  saveRDS(as.data.frame(conv_df), paste0("output/", SAMPLER, "_pixel_convergence_check_", MODEL_LABEL, "_", timestamp_str, ".rds"))

  # (b2) Rotated betas: ALL pairwise category contrasts (to - from) per covariate, pooled across
  #      groups — same tidy schema as the MNL mnl_pixel_beta_rotated artifact. Uses the symmetric
  #      (zero-sum) pooled draws directly; per-group contrasts would need stored RE draws.
  draws3 <- array(arr, c(P, K, nd * N_CHAINS), dimnames = list(cov_names_processed, cat_all, NULL))
  .rq <- if (requireNamespace("matrixStats", quietly = TRUE)) {
    function(M) cbind(matrixStats::rowMedians(M), matrixStats::rowQuantiles(M, probs = c(0.025, 0.975)))
  } else {
    function(M) t(apply(M, 1, stats::quantile, probs = c(0.5, 0.025, 0.975)))
  }
  rot_list <- vector("list", K * (K - 1))
  ci <- 1L
  for (fi in seq_len(K)) {
    for (ti in seq_len(K)) {
      if (fi == ti) next
      qq <- .rq(draws3[, ti, ] - draws3[, fi, ])
      # sd_x carries the DISPLAY scale with the artifact (rule B, codes/mnl_aux_func.R), so any
      # consumer can put coefficients on a comparable footing without needing the design: multiply
      # value_* by sd_x. Without it a reader ranks raw units and concludes climate does nothing.
      rot_list[[ci]] <- data.table::data.table(
        ks = cov_names_processed, from_class = cat_all[fi], to_class = cat_all[ti],
        value_median = qq[, 1], value_q025 = qq[, 2], value_q975 = qq[, 3],
        sd_x = as.numeric(display_sds(X_mat, cov_names_processed)[cov_names_processed])
      )
      ci <- ci + 1L
    }
  }
  saveRDS(
    as.data.frame(data.table::rbindlist(rot_list)),
    paste0("output/", SAMPLER, "_pixel_beta_rotated_", MODEL_LABEL, "_", timestamp_str, ".rds")
  )
  cat(sprintf(">>> Saved convergence-check + rotated-beta (pooled, %d pairwise contrasts) artifacts.\n", K * (K - 1)))

  cat(sprintf("\nEstimation Complete (%s). Saved fit + heatplot to output/.\n", SAMPLER))
  if (RUN_MODE == "test" && PROMOTE_TO_PRODUCTION) promote_test_to_production()
  quit(save = "no", status = 0)
}

# 7. RECOVER POSTERIOR FROM DISK
# Change this to a specific directory path if you want to recover an existing run from disk
# res_list <- "output/saved_model_outputs/10km_mnlogit_rcpp_sym_BIOCLIMA_test_RE_CAPRI_NUTS"

cat("\n>>> Recovering full posterior from disk batches...\n")
source("codes/mnlogit_rcpp_sym.R")
if (is.character(res_list)) {
  dir_path <- res_list
  res_list <- lapply(seq_len(N_CHAINS), function(ch) {
    recover_mnlogit_posterior(dir_path, chain_id = ch)
  })
} else {
  res_list <- lapply(res_list, function(r) {
    if (!is.null(r$posterior_store)) {
      return(recover_mnlogit_posterior(r$posterior_store))
    }
    return(r)
  })
}

# 8. ZERO-SUM POST-PROCESSING (Approach B extraction)
# For each posterior draw, project the share coefficients onto the zero-sum
# subspace: beta_tilde_k = beta_k - mean(beta) for each share group.
cov_names_processed <- colnames(X_mat)

# Identify covariate indices for each share group
focal_idx <- which(cov_names_processed %in% focal_cov_cols)
soil_group_idx <- lapply(soil_groups, function(cols) which(cov_names_processed %in% cols))

all_share_groups <- c(list(focal = focal_idx), soil_group_idx)
cat(sprintf(
  "\n>>> Zero-sum projection on %d share groups (%s)...\n",
  length(all_share_groups),
  paste(names(all_share_groups), sapply(all_share_groups, length), sep = ":", collapse = ", ")
))

for (ch in seq_along(res_list)) {
  # A. postb_pooled [K_cov, P_cat, draws]
  if (!is.null(res_list[[ch]]$postb_pooled)) {
    A <- res_list[[ch]]$postb_pooled
    for (grp_name in names(all_share_groups)) {
      idx <- all_share_groups[[grp_name]]
      if (length(idx) >= 2) {
        # For each category and draw, subtract the mean across the group
        grp_mean <- colMeans(A[idx, , , drop = FALSE]) # [P_cat, draws]
        for (i in idx) {
          A[i, , ] <- A[i, , ] - grp_mean
        }
      }
    }
    res_list[[ch]]$postb_pooled <- A
  }

  # B. postb_total [K_cov, P_cat, draws] or [K_cov, P_cat, G, draws]
  if (!is.null(res_list[[ch]]$postb_total)) {
    A <- res_list[[ch]]$postb_total
    d <- length(dim(A))
    for (grp_name in names(all_share_groups)) {
      idx <- all_share_groups[[grp_name]]
      if (length(idx) >= 2) {
        if (d == 3) {
          grp_mean <- colMeans(A[idx, , , drop = FALSE])
          for (i in idx) A[i, , ] <- A[i, , ] - grp_mean
        } else if (d == 4) {
          grp_mean <- colMeans(A[idx, , , , drop = FALSE])
          for (i in idx) A[i, , , ] <- A[i, , , ] - grp_mean
        }
      }
    }
    res_list[[ch]]$postb_total <- A
  }
}
cat("    Zero-sum projection complete.\n")

# 9. ASSESS CONVERGENCE
assess_convergence(res_list,
  name = paste0(MODEL_LABEL, "_", if (use_re) "RE" else "pooled"),
  categories = final_cats_pixel,
  baseline_idx = baseline_idx,
  type = if (use_re) "re" else "pooled"
)

res_full <- combine_chains(res_list, keep_chains = TRUE)

# # =========================================================================
# # POST-PROCESSING RANDOM EFFECTS CLEANUP
# # =========================================================================
# # Forces country-specific random deviations to zero (total = pooled) if:
# # The covariate has negligible or zero standard deviation within a country (group_sd < 1e-4),
# # meaning there is no variation within that group to estimate a country-specific deviation.
# if (use_re && length(dim(res_full$postb_total)) == 4) {
#   cat("\n>>> Applying Post-Processing RE Cleanup based on within-group variation...\n")
#
#   group_names_vec <- if (exists("re_group_names") && !is.null(re_group_names)) {
#     re_group_names
#   } else if (!is.null(res_full$nuts0_names)) {
#     res_full$nuts0_names
#   } else {
#     paste0("group_", seq_len(dim(res_full$postb_total)[3]))
#   }
#
#   cleaned_count <- 0
#   for (v in seq_len(dim(res_full$postb_total)[1])) {
#     cov_name <- colnames(X_mat)[v]
#     if (cov_name == "intercept") next
#
#     for (m in seq_len(dim(res_full$postb_total)[3])) {
#       group_rows <- which(group_idx_vec == m)
#
#       # Calculate standard deviation of covariate within this country
#       group_sd <- if (length(group_rows) > 1) sd(X_mat[group_rows, v], na.rm = TRUE) else 0
#
#       if (is.na(group_sd) || group_sd < 1e-4) {
#         # Force deviation to 0: set postb_total = postb_pooled for this covariate, country
#         res_full$postb_total[v, , m, ] <- res_full$postb_pooled[v, , ]
#         cleaned_count <- cleaned_count + 1
#       }
#     }
#   }
#   cat(sprintf("    Successfully forced %d country-covariate deviations with no within-group variation back to the pooled mean.\n", cleaned_count))
# }
#
# # (Back-transformation is now handled automatically by recover_mnlogit_posterior)

saveRDS(res_full, paste0("output/mnl_pixel_model_fit_", paste0(MODEL_LABEL, "_", if (use_re) "RE" else "pooled"), "_", timestamp_str, ".rds"))

# =========================================================================
# IN-PIPELINE FIT METRICS  (correct-by-construction: uses the SAME X_mat +
# group index the sampler was given, so predicted shares match the model's
# own utilities). Validated against the stored post_log_lik first -- if a
# reconstructed draw's log_lik does NOT match, the metrics are NOT reported
# (guards against the off-pipeline reconstruction drift we hit in analysis).
# McFadden pseudo-R2 + per-class calibration -> printed + saved to CSV.
# =========================================================================
if (!is.null(res_full$postb_total) && !is.null(res_full$post_log_lik)) {
  cat("\n>>> Computing in-pipeline fit metrics (validated against stored log_lik)...\n")
  .Bt  <- res_full$postb_total                 # RE: [cov,cat,group,draw,chain]; pooled: [cov,cat,draw,chain]
  .pll <- res_full$post_log_lik                 # [draw, chain]
  .Xl  <- as.matrix(X_mat[, linear_cols, drop = FALSE])
  .J   <- dim(.Bt)[2]
  .isRE <- (length(dim(.Bt)) == 5L) && use_re
  .glev <- if (.isRE) unique(group_idx_vec) else NULL      # appearance-order group columns (predict_shares convention)
  .lse  <- function(U) { m <- apply(U, 1, max); m + log(rowSums(exp(U - m))) }
  .Uof  <- function(Bg) {                                   # linear utilities for one draw (BART-free path)
    if (!.isRE) return(.Xl %*% Bg)
    U <- matrix(0, nrow(.Xl), .J)
    for (g in unique(group_idx_vec)) { r <- which(group_idx_vec == g); U[r, ] <- .Xl[r, , drop = FALSE] %*% Bg[, , match(g, .glev)] }
    U
  }
  # ---- TRUSTWORTHY metric (no reconstruction): McFadden from the model's OWN stored log_lik ----
  # post_log_lik is what the sampler itself computed for each draw -> comparable to the same-scale
  # null (predict the global mean shares). No beta->utility replication, so always valid.
  .pbar <- colMeans(Y_pixel); .n <- nrow(Y_pixel)
  .LLn <- .n * sum(.pbar[.pbar > 0] * log(.pbar[.pbar > 0]))
  .LLp <- sum(Y_pixel[Y_pixel > 0] * log(Y_pixel[Y_pixel > 0]))
  .llv <- as.numeric(.pll)
  cat(sprintf("    STORED-log_lik fit (trustworthy): null LL=%.0f | perfect LL=%.0f\n", .LLn, .LLp))
  cat(sprintf("    McFadden R2 [mean draw]=%.3f  [median]=%.3f  [best]=%.3f   (in-sample; conservative — per-draw not posterior-predictive)\n",
              1 - mean(.llv)/.LLn, 1 - median(.llv)/.LLn, 1 - max(.llv)/.LLn))
  cat(sprintf("    normalized fit [median draw]=%.3f  [0=null,1=perfect]  | %% of draws beating null = %.1f%%\n",
              (median(.llv)-.LLn)/(.LLp-.LLn), 100*mean(.llv > .LLn)))
  saveRDS(data.frame(metric=c("LL_null","LL_perfect","LL_mean_draw","LL_median_draw","LL_best_draw",
                              "McF_mean","McF_median","McF_best","frac_draws_beat_null"),
                     value=c(.LLn,.LLp,mean(.llv),median(.llv),max(.llv),
                             1-mean(.llv)/.LLn,1-median(.llv)/.LLn,1-max(.llv)/.LLn,mean(.llv>.LLn))),
          paste0("output/mnl_pixel_fit_stored_", MODEL_LABEL, "_", timestamp_str, ".rds"))

  if (isTRUE(use_bart)) cat("    (note: BART is ON -> reconstructed per-class utilities below OMIT the BART term)\n")
  # ---- validation guard: reconstructed draw-1 log_lik must match the stored value ----
  .b1 <- if (.isRE) .Bt[, , , 1, 1] else .Bt[, , 1, 1]
  .U1 <- .Uof(.b1); .ll1 <- sum(Y_pixel * (.U1 - .lse(.U1)))
  .stored1 <- if (!is.null(dim(.pll))) .pll[1, 1] else .pll[1]
  cat(sprintf("    [validate] reconstructed draw-1 log_lik = %.0f | stored = %.0f | diff = %.1f\n", .ll1, .stored1, .ll1 - .stored1))
  if (!isTRUE(use_bart) && abs(.ll1 - .stored1) < max(50, 1e-3 * abs(.stored1))) {
    .ndr <- dim(.Bt)[if (.isRE) 4 else 3]; .nch <- dim(.Bt)[if (.isRE) 5 else 4]
    .ss <- unique(round(seq(1, .ndr, length.out = min(40L, .ndr))))
    .P <- matrix(0, nrow(.Xl), .J); .cnt <- 0L
    for (ch in seq_len(.nch)) for (s in .ss) {
      .Bg <- if (.isRE) .Bt[, , , s, ch] else .Bt[, , s, ch]
      .U <- .Uof(.Bg); .mx <- apply(.U, 1, max); .E <- exp(.U - .mx); .P <- .P + .E / rowSums(.E); .cnt <- .cnt + 1L
    }
    .P <- .P / .cnt; colnames(.P) <- final_cats_pixel
    .pbar <- colMeans(Y_pixel); .n <- nrow(Y_pixel)
    .LLn <- .n * sum(.pbar[.pbar > 0] * log(.pbar[.pbar > 0]))
    .LLp <- sum(Y_pixel[Y_pixel > 0] * log(Y_pixel[Y_pixel > 0]))
    .LLm <- sum(Y_pixel * log(pmax(.P, 1e-12)))
    cat(sprintf("    >>> McFadden pseudo-R2 = %.3f | normalized fit = %.3f  [0=null, 1=perfect]  (%d draws)\n",
                1 - .LLm / .LLn, (.LLm - .LLn) / (.LLp - .LLn), .cnt))
    .fm <- data.frame(class = final_cats_pixel, obs_share = round(.pbar, 5), pred_share = round(colMeans(.P), 5),
                      pixel_corr = round(vapply(seq_len(.J), function(j) if (sd(.P[, j]) > 0 && sd(Y_pixel[, j]) > 0) cor(.P[, j], Y_pixel[, j]) else NA_real_, numeric(1)), 3))
    .fm <- .fm[order(-.fm$obs_share), ]; print(.fm, row.names = FALSE)
    write.csv(.fm, paste0("output/mnl_pixel_fit_metrics_", MODEL_LABEL, "_", timestamp_str, ".csv"), row.names = FALSE)
    cat(sprintf("    >>> saved output/mnl_pixel_fit_metrics_%s_%s.csv\n", MODEL_LABEL, timestamp_str))
  } else if (!isTRUE(use_bart)) {
    cat("    !! VALIDATION FAILED (reconstructed log_lik != stored) -> fit metrics NOT reported (would be untrustworthy).\n")
  }
}

# Aligned cross-sampler fit object (same schema as LNM/CLR): symmetric K-category beta.
.ppm <- res_full$postb_pooled # [cov, cat, draws, chains]
.catm <- dimnames(.ppm)[[2]]
if (dim(.ppm)[2] == length(final_cats_pixel) - 1L) { # baseline-removed -> expand to symmetric K
  .ppm <- abind::abind(.ppm, array(0, c(dim(.ppm)[1], 1, dim(.ppm)[3], dim(.ppm)[4])), along = 2)
  .ppm <- sweep(.ppm, c(1, 3, 4), apply(.ppm, c(1, 3, 4), mean), "-") # center across categories -> zero-sum
  .catm <- c(final_cats_pixel[-baseline_idx], final_cats_pixel[baseline_idx])
} else if (is.null(.catm)) .catm <- final_cats_pixel
saveRDS(
  build_fit_object(SAMPLER, .ppm, cov_names_processed, .catm, baseline_class, NULL, dropped_blocks),
  paste0("output/", SAMPLER, "_pixel_model_fit_", MODEL_LABEL, "_", timestamp_str, ".rds")
)
cat(sprintf(
  ">>> Saved aligned cross-sampler fit object: B[%d cov x %d cat] -> output/%s_pixel_model_fit_...\n",
  dim(.ppm)[1], dim(.ppm)[2], SAMPLER
))

# 9. GENERATE PARAMETER HEATPLOT
cat("\n>>> Generating Parameter Heatplot...\n")
# Create understandable labels and scaling info
# Display scale = rule B, now defined ONCE in codes/mnl_aux_func.R (display_sds) so every plot,
# table and artifact in the repo uses the same convention: per 1 SD for unbounded covariates, per
# 1 unit for shares and bounded indices (which are already zero-sum projected, so the raw
# coefficient IS the per-unit effect). This block used to define that rule inline, here only, which
# is why the other heatplot, the report heatplot and the rotated artifacts all showed raw
# coefficients -- and why climate looked inert when per SD it is one of the largest families.
sd_final <- display_sds(X_mat, cov_names_processed)

understandable_labels <- format_mnl_labels(cov_names_processed, sd_final)

# Reorder covariates for better grouping in the plot
# Pass cat_names so focal shares are sorted to match the Y-category order (diagonal)
cat_names_sym <- c(final_cats_pixel[-baseline_idx], final_cats_pixel[baseline_idx])
reordered_covs <- reorder_mnl_covariates(cov_names_processed, cat_names = cat_names_sym)
cov_order_idx <- match(reordered_covs, cov_names_processed)

# Calculate Full Symmetric Effects
cat("  -> Computing symmetric effects (centered across all categories including baseline)...\n")
plot_draws <- res_full$postb_pooled[cov_order_idx, , , , drop = FALSE]
dim_p <- dim(plot_draws)

if (dim_p[2] == length(final_cats_pixel)) {
  # It's already symmetrically expanded by recover_mnlogit_posterior
  sym_draws_centered <- plot_draws
  cat_names_sym <- if (!is.null(dimnames(plot_draws)[[2]])) dimnames(plot_draws)[[2]] else c(final_cats_pixel[-baseline_idx], final_cats_pixel[baseline_idx])
  dimnames(sym_draws_centered)[[2]] <- cat_names_sym
} else {
  # Needs manual expansion (for older pooled non-symmetric batches)
  baseline_zeros <- array(0, dim = c(dim_p[1], 1, dim_p[3], dim_p[4]))
  sym_draws <- abind::abind(plot_draws, baseline_zeros, along = 2)
  mean_across_cats <- apply(sym_draws, c(1, 3, 4), mean)
  sym_draws_centered <- sweep(sym_draws, c(1, 3, 4), mean_across_cats, "-")
  cat_names_sym <- c(final_cats_pixel[-baseline_idx], final_cats_pixel[baseline_idx])
  dimnames(sym_draws_centered)[[2]] <- cat_names_sym
}

p_heat <- MNL_parameter_heatplot(
  draws             = sym_draws_centered,
  cov_names         = understandable_labels[cov_order_idx],
  cat_names         = cat_names_sym,
  exclude_intercept = TRUE,
  scaling_vec       = sd_final[reordered_covs],
  title             = paste0(MODEL_LABEL, "km MNL Coefficient Heatplot (Symmetric)"),
  subtitle          = "Tile size = Credibility, Color = Median Estimate. Symmetric effects centered across all categories."
) +
  theme(axis.text.y = element_text(size = 9, hjust = 1))

ggsave(paste0("output/plots/mnl_pixel_", MODEL_LABEL, "_heatplot.png"), p_heat, width = 12, height = 20)

cat("\nEstimation Complete. Diagnostics and plots saved to output/.\n")

# 7. EXTRACT POSTERIOR BETA MEANS
# =========================================================================
cat("\nStep 4: Extracting Posterior Beta Means...\n")
# Calculate the posterior median across MCMC draws
if (length(dim(res_full$postb_total)) == 5) {
  # Random effects model: [covariates, categories, groups, draws, chains]
  # GROUP LABELLING (fixed 2026-08-05). postb_total[, , k, , ] is the k-th group in APPEARANCE
  # order: the sampler builds its per-group arrays from `groups <- unique(group_idx)`
  # (mnlogit_rcpp_sym.R ~L595), NOT from the sorted factor levels. group_idx_vec =
  # as.integer(factor) carries SORTED-level ids, so slice k is re_group_names[unique(group_idx_vec)[k]].
  # Labelling with the sorted levels directly mislabelled ALL 26 countries here (the pixel rows are
  # ordered spatially, not by country: slice 1 "Austria" was actually Portugal, "Belgium" was Spain).
  # Affects the exported per-country betas / RE plots only -- predict_shares is passed
  # group_levels = unique(group_idx) explicitly (recipe$group_levels_appear) and was never wrong.
  group_names_vec <- if (exists("re_group_names") && !is.null(re_group_names)) {
    .appear <- if (exists("group_idx_vec") && !is.null(group_idx_vec)) unique(as.integer(group_idx_vec)) else seq_along(re_group_names)
    if (length(.appear) != length(re_group_names))
      stop(sprintf("group labelling: %d appearance-order groups vs %d level names -- cannot key postb_total safely.",
                   length(.appear), length(re_group_names)))
    re_group_names[.appear]
  } else if (!is.null(res_full$nuts0_names)) {
    res_full$nuts0_names
  } else {
    paste0("group_", seq_len(dim(res_full$postb_total)[3]))
  }

  # Extract dimensions
  d_ks <- dim(res_full$postb_total)[1]
  d_lu <- dim(res_full$postb_total)[2]
  d_Ns <- dim(res_full$postb_total)[3]
  d_draws <- dim(res_full$postb_total)[4]
  d_chains <- dim(res_full$postb_total)[5]
  N_params <- d_ks * d_lu * d_Ns

  # 1. Create the final coordinate data frame directly using expand.grid.
  # By listing variables in dimension order (ks, lu_to, Ns), it perfectly
  # matches the column-major memory layout of the flattened array.
  check_point_rhat_df <- data.table::as.data.table(expand.grid(
    ks = cov_names_processed,
    lu_to = if (d_lu == ncol(Y_pixel)) colnames(Y_pixel) else colnames(Y_pixel)[-baseline_idx],
    Ns = group_names_vec,
    stringsAsFactors = FALSE
  ))
  setnames(check_point_rhat_df, "Ns", RE_GROUP_COL)

  # 2. Flatten the 5D array to a 2D Matrix [N_params, draws * chains]
  # This doesn't duplicate the memory if done right, and is much faster to traverse
  flat_mat <- matrix(res_full$postb_total, nrow = N_params, ncol = d_draws * d_chains)

  cat("Calculating posterior medians across MCMC draws for each group...\n")
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    check_point_rhat_df$value_median <- matrixStats::rowMedians(flat_mat)
    check_point_rhat_df$value_sd <- matrixStats::rowSds(flat_mat)
    quantiles <- matrixStats::rowQuantiles(flat_mat, probs = c(0.025, 0.975))
    check_point_rhat_df$value_q025 <- quantiles[, 1]
    check_point_rhat_df$value_q975 <- quantiles[, 2]
    prob_pos <- matrixStats::rowMeans2(flat_mat > 0)
    check_point_rhat_df$prob_dir <- pmax(prob_pos, 1 - prob_pos)
  } else {
    check_point_rhat_df$value_median <- apply(flat_mat, 1, median)
    check_point_rhat_df$value_sd <- apply(flat_mat, 1, sd)
    check_point_rhat_df$value_q025 <- apply(flat_mat, 1, quantile, probs = 0.025)
    check_point_rhat_df$value_q975 <- apply(flat_mat, 1, quantile, probs = 0.975)
    prob_pos <- apply(flat_mat > 0, 1, mean)
    check_point_rhat_df$prob_dir <- pmax(prob_pos, 1 - prob_pos)
  }
  check_point_rhat_df$is_sig_95 <- sign(check_point_rhat_df$value_q025) == sign(check_point_rhat_df$value_q975)

  cat("Calculating posterior rhat across MCMC chains for each group...\n")
  # 3. Calculate rhat using a pre-allocated loop (O(1) memory footprint)
  rhat_vals <- numeric(N_params)
  for (i in seq_len(N_params)) {
    # Reshape the flat row back into a [draws, chains] matrix for posterior::rhat
    chain_mat <- matrix(flat_mat[i, ], nrow = d_draws, ncol = d_chains)
    rhat_vals[i] <- posterior::rhat(chain_mat)
  }
  check_point_rhat_df$value_rhat <- rhat_vals

  # Clean up memory
  rm(flat_mat)
  gc()

  saveRDS(as.data.frame(check_point_rhat_df), paste0("output/mnl_pixel_convergence_check_", MODEL_LABEL, "_", timestamp_str, ".rds"))

  # Save to CSV
  if (use_re) {
    beta_point_df <- as.data.frame(check_point_rhat_df)[, c(RE_GROUP_COL, "lu_to", "ks", "value_median")] %>%
      dplyr::mutate(value = value_median) %>%
      dplyr::select(!value_median)
  } else {
    beta_point_df <- as.data.frame(check_point_rhat_df)[, c("lu_to", "ks", "value_median")] %>%
      dplyr::mutate(value = value_median) %>%
      dplyr::select(!value_median)
  }
  beta_point_df$lu_to <- as.character(beta_point_df$lu_to)
  beta_point_df$ks <- as.character(beta_point_df$ks)

  out_file <- paste0("output/mnl_pixel_beta_median_", paste0(MODEL_LABEL, if (use_re) "_RE_" else "_pooled_", if (use_re) RE_GROUP_COL else ""), "_", timestamp_str, ".csv")
  write.csv(beta_point_df, out_file, row.names = FALSE)
  cat(sprintf("Saved posterior beta means to %s\n", out_file))

  # =========================================================================
  # 7.5 EXTRACT ROTATED PARAMETERS FOR ALL BASELINE COMBINATIONS
  # =========================================================================
  if (DO_POSTERIOR_ROTATION) {
    cat("\nStep 4.5: Extracting rotated parameters for all baseline combinations and extending to missing grid countries...\n")

    categories_all <- final_cats_pixel
    d_lu_all <- length(categories_all)

    # 1. Identify missing grid countries
    if (exists("grid_map_pixel") && "Grouping_Key" %in% names(grid_map_pixel)) {
      all_grid_countries <- unique(grid_map_pixel[["Grouping_Key"]])
      all_grid_countries <- all_grid_countries[!is.na(all_grid_countries)]
    } else {
      all_grid_countries <- group_names_vec
    }
    missing_countries <- setdiff(all_grid_countries, group_names_vec)
    G_all <- length(group_names_vec) + length(missing_countries)
    all_countries <- c(group_names_vec, missing_countries)

    cat(sprintf("Found %d countries in sample, %d missing in grid.\n", length(group_names_vec), length(missing_countries)))

    # 2. Reconstruct full array with baseline as 0
    # res_full$postb_total: [K, P, G, draws, chains] -> flatten to [K, P, G, draws * chains]
    beta_flat_re <- matrix(res_full$postb_total, nrow = d_ks, ncol = d_lu * d_Ns * d_draws * d_chains)
    dim(beta_flat_re) <- c(d_ks, d_lu, d_Ns, d_draws * d_chains)

    # 3. Reconstruct pooled mean array with baseline as 0
    # res_full$postb_pooled: [K, P, draws, chains]
    mu_flat <- matrix(res_full$postb_pooled, nrow = d_ks, ncol = d_lu * d_draws * d_chains)
    dim(mu_flat) <- c(d_ks, d_lu, d_draws * d_chains)

    if (d_lu == d_lu_all) {
      beta_full <- beta_flat_re
      mu_full <- mu_flat
    } else {
      beta_full <- array(0.0, dim = c(d_ks, d_lu_all, d_Ns, d_draws * d_chains))
      non_base_idx <- (1:d_lu_all)[-baseline_idx]
      beta_full[, non_base_idx, , ] <- beta_flat_re

      mu_full <- array(0.0, dim = c(d_ks, d_lu_all, d_draws * d_chains))
      mu_full[, non_base_idx, ] <- mu_flat
    }

    # 4. Extend to missing countries using pooled mean
    beta_extended <- array(0.0, dim = c(d_ks, d_lu_all, G_all, d_draws * d_chains))
    beta_extended[, , 1:d_Ns, ] <- beta_full
    if (length(missing_countries) > 0) {
      for (i in seq_along(missing_countries)) {
        beta_extended[, , d_Ns + i, ] <- mu_full
      }
    }

    # Clean up memory
    rm(beta_flat_re, beta_full, mu_flat, mu_full)
    gc()

    # 5. Rotate and summarize
    total_pairs <- d_lu_all * (d_lu_all - 1)
    cat(sprintf("Rotating parameters across %d combinations...\n", total_pairs))

    results_list <- vector("list", total_pairs)
    counter <- 1

    pb <- utils::txtProgressBar(min = 0, max = total_pairs, style = 3)
    for (f_idx in seq_along(categories_all)) {
      from_class <- categories_all[f_idx]
      for (t_idx in seq_along(categories_all)) {
        if (f_idx == t_idx) next
        to_class <- categories_all[t_idx]

        # Difference: to - from. Dimension: [K, G_all, draws * chains]
        diff_array <- beta_extended[, t_idx, , , drop = FALSE] - beta_extended[, f_idx, , , drop = FALSE]
        # Flatten to [K * G_all, draws * chains] for fast quantiles
        diff_mat <- matrix(diff_array, nrow = d_ks * G_all, ncol = d_draws * d_chains)

        if (requireNamespace("matrixStats", quietly = TRUE)) {
          meds <- matrixStats::rowMedians(diff_mat)
          quants <- matrixStats::rowQuantiles(diff_mat, probs = c(0.025, 0.975))
        } else {
          meds <- apply(diff_mat, 1, median)
          quants <- apply(diff_mat, 1, quantile, probs = c(0.025, 0.975))
          if (is.null(dim(quants))) quants <- matrix(quants, ncol = 2, byrow = TRUE) else quants <- t(quants)
        }

        dt <- data.table::data.table(
          ks = rep(cov_names_processed, G_all),
          from_class = from_class,
          to_class = to_class,
          group = rep(all_countries, each = d_ks),
          value_median = meds,
          value_q025 = quants[, 1],
          value_q975 = quants[, 2]
        )

        results_list[[counter]] <- dt
        counter <- counter + 1
        utils::setTxtProgressBar(pb, counter - 1)
      }
    }
    close(pb)

    cat("Binding and saving rotated parameters...\n")
    final_rotated_df <- data.table::rbindlist(results_list)
    data.table::setnames(final_rotated_df, "group", RE_GROUP_COL)

    saveRDS(as.data.frame(final_rotated_df), paste0("output/mnl_pixel_beta_rotated_FULL_", MODEL_LABEL, "_", timestamp_str, ".rds"))
  }
} else {
  # Pooled model: [covariates, categories, draws, chains]
  d_ks <- dim(res_full$postb_pooled)[1]
  d_lu <- dim(res_full$postb_pooled)[2]
  d_draws <- dim(res_full$postb_pooled)[3]
  d_chains <- dim(res_full$postb_pooled)[4]
  N_params <- d_ks * d_lu

  check_point_rhat_df <- data.table::as.data.table(expand.grid(
    ks = cov_names_processed,
    lu_to = if (d_lu == ncol(Y_pixel)) colnames(Y_pixel) else colnames(Y_pixel)[-baseline_idx],
    stringsAsFactors = FALSE
  ))

  flat_mat <- matrix(res_full$postb_pooled, nrow = N_params, ncol = d_draws * d_chains)

  cat("Calculating posterior medians across MCMC draws...\n")
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    check_point_rhat_df$value_median <- matrixStats::rowMedians(flat_mat)
    check_point_rhat_df$value_sd <- matrixStats::rowSds(flat_mat)
    quantiles <- matrixStats::rowQuantiles(flat_mat, probs = c(0.025, 0.975))
    check_point_rhat_df$value_q025 <- quantiles[, 1]
    check_point_rhat_df$value_q975 <- quantiles[, 2]
    prob_pos <- matrixStats::rowMeans2(flat_mat > 0)
    check_point_rhat_df$prob_dir <- pmax(prob_pos, 1 - prob_pos)
  } else {
    check_point_rhat_df$value_median <- apply(flat_mat, 1, median)
    check_point_rhat_df$value_sd <- apply(flat_mat, 1, sd)
    check_point_rhat_df$value_q025 <- apply(flat_mat, 1, quantile, probs = 0.025)
    check_point_rhat_df$value_q975 <- apply(flat_mat, 1, quantile, probs = 0.975)
    prob_pos <- apply(flat_mat > 0, 1, mean)
    check_point_rhat_df$prob_dir <- pmax(prob_pos, 1 - prob_pos)
  }
  check_point_rhat_df$is_sig_95 <- sign(check_point_rhat_df$value_q025) == sign(check_point_rhat_df$value_q975)

  cat("Calculating posterior rhat across MCMC chains...\n")
  rhat_vals <- numeric(N_params)
  for (i in seq_len(N_params)) {
    chain_mat <- matrix(flat_mat[i, ], nrow = d_draws, ncol = d_chains)
    rhat_vals[i] <- posterior::rhat(chain_mat)
  }
  check_point_rhat_df$value_rhat <- rhat_vals

  rm(flat_mat)
  gc()

  saveRDS(as.data.frame(check_point_rhat_df), paste0("output/mnl_pixel_convergence_check_pooled_", MODEL_LABEL, "_", timestamp_str, ".rds"))

  # Save to CSV
  beta_point_df <- as.data.frame(check_point_rhat_df)[, c("lu_to", "ks", "value_median")] %>%
    dplyr::mutate(value = value_median) %>%
    dplyr::select(!value_median)
  beta_point_df$lu_to <- as.character(beta_point_df$lu_to)
  beta_point_df$ks <- as.character(beta_point_df$ks)

  out_file <- paste0("output/mnl_pixel_beta_median_pooled_", MODEL_LABEL, "_", timestamp_str, ".csv")
  write.csv(beta_point_df, out_file, row.names = FALSE)
  cat(sprintf("Saved posterior beta means to %s\n", out_file))
}

# =========================================================================
# 10. PROMOTE TEST -> PRODUCTION (optional, DRIVER_PROMOTE=TRUE)
# =========================================================================
if (RUN_MODE == "test" && PROMOTE_TO_PRODUCTION) promote_test_to_production()

cat("\nBuilding Diagnostic HTML Report...\n")
source("postprocess/nested_report.R")
