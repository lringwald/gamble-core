# =============================================================================
# run_prior_module_count_model.R
# =============================================================================
# Self-Contained Count Model for Livestock (Multi-Year Panel)
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

# (Diagnostic output directories are pre-created below, MODEL_LABEL-derived, after RUN_MODE is set.)

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
LS_DIR <- Sys.getenv("GAMBLE_LS_DIR", "../LAND_SUPPLY_ELASTICITY")   # relative (sibling of gamble-core) + env-overridable; was a hardcoded OneDrive path, broke on the Google Drive migration
# Canonical data root -- ../LAMASUS_* paths are superseded frozen exports. This also gives the count
# model the 2010 LUM map, which LAMASUS_downscaling/input does not have: the year rule below fell back
# to 2018 for the 2010 tier, so 2010 livestock counts were explained by land use from EIGHT YEARS
# LATER purely because the file was missing from the directory it was reading.
CASCADE_DATA <- Sys.getenv("GAMBLE_CASCADE_DATA", "../cascadinggamble-core/data")
DS_DIR <- Sys.getenv("GAMBLE_DS_DIR", "../LAMASUS_downscaling/")
# Same split as the pixel driver: bulk rasters/LUM in 02_intermediate, lookup tables in aux_files.
# DS_DIR survives only for irrigation_binary_EEA_1kmID.rds, which has no canonical copy and is read
# only when INCLUDE_IRRIGATION=TRUE (off by default in both models).
GRIDWORK_DIR <- Sys.getenv("GAMBLE_GRIDWORK_DIR", file.path(CASCADE_DATA, "02_intermediate"))
AUXDATA_DIR  <- Sys.getenv("GAMBLE_AUXDATA_DIR",  file.path(CASCADE_DATA, "aux_files"))
# Master 1km covariate parquet: the gridwork pipeline MIGRATED to ../cascadinggamble-core (Snakemake),
# whose data/02_intermediate copy is the LIVE one -- it carries the new GHM v3 threat groups (TI/NS/AG)
# that the legacy LAMASUS_gridwork copy does not. Pick the NEWEST available; env GAMBLE_MASTER_PARQUET
# overrides. The other inputs (1km mapping, LSU counts, organic) still come from GRIDWORK_DIR.
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

# GHM v3 threat groups arrive under one of TWO namings: the merged form GHM_<CODE>_<year> (HI, Ovr)
# or the RAW Earth-Engine export name GHM_v3_<CODE>_<year>_EU_UK_<timestamp>_aggregated (TI, NS, AG --
# merge_all_variables.R's get_clean_name only renamed HI|Overall). Resolve either, then rename to the
# clean form on read so everything downstream sees GHM_<CODE>_<year>.
.ghm_resolve <- function(schema_names, code, year) {
  cn <- paste0("GHM_", code, "_", year)
  if (cn %in% schema_names) return(cn)
  raw <- grep(sprintf("^GHM_v3_%s_%s_", code, year), schema_names, value = TRUE)
  if (length(raw)) raw[1] else NA_character_
}

# Run versioning (matches the pixel driver): bucket ALL outputs into a stable "production" or a
# throwaway "test" namespace via MODEL_LABEL. A "test" run starts fresh (its saved-model dir is
# cleared); "production" persists so a re-run resumes. DRIVER_PROMOTE copies a finished test run's
# artifacts to the "_production" twin.
RUN_MODE <- Sys.getenv("DRIVER_RUN_MODE", "test") # "production" | "test"
if (!RUN_MODE %in% c("production", "test")) stop("RUN_MODE must be 'production' or 'test'")
PROMOTE_TO_PRODUCTION <- isTRUE(as.logical(Sys.getenv("DRIVER_PROMOTE", "FALSE")))
# Label reflects the OBSERVATION RESOLUTION (the admin units are NUTS3), not a generic "Admin".
RESOLUTION_LABEL <- Sys.getenv("DRIVER_RESOLUTION_LABEL", "NUTS3")
# Sampler token woven into MODEL_LABEL so every artifact is traceable — mirrors the pixel driver, whose
# label is <res>km_<sampler>_by_<intersect>_<scheme>_<mode> (e.g. 10km_mnlogit_rcpp_sym_by_NUTS3_GLOBIOM_production).
# Count analogue: <units>_<sampler>_<scheme>_<mode> (units = NUTS3; sampler = count_rcpp; no separate intersect).
SAMPLER <- Sys.getenv("DRIVER_SAMPLER", "count_rcpp")
MODEL_LABEL <- paste0(RESOLUTION_LABEL, "_", SAMPLER)   # scheme + mode appended below (after CLASS_SCHEME is known)

# =========================================================================
# CLASSIFICATION CONFIGURATION (from pixel model)
# =========================================================================
# The classification scheme defines how the source map is parsed into `model_class` categories.
CLASS_COLS <- strsplit(Sys.getenv("DRIVER_CLASS_COLS", "GLOBIOM_UNFCCC,GLOBIOM_mngmt"), ",")[[1]]
OUTCOME_SOURCE <- Sys.getenv("DRIVER_OUTCOME_SOURCE", "LUM")
CLASS_SCHEME <- toupper(sub("[0-9]+$", "", sub("_.*$", "", CLASS_COLS[1])))
MODEL_LABEL <- paste0(MODEL_LABEL, "_", CLASS_SCHEME, "_", RUN_MODE)   # -> e.g. NUTS3_count_rcpp_GLOBIOM_production (pixel-style ordering)

INCLUDE_IRRIGATION <- FALSE # switch: carve irrigated cropland into its own ..._IR class
INCLUDE_ORGANIC <- TRUE # switch: carve organic twins into their own ..._O class
# Organic AREA was far below today's in the past (EU organic share of UAA: ~3% 2000, ~5% 2010,
# ~8% 2018, ~9% 2020), but organic_certificaties_..._master.rds is a single ~present-day map applied
# to every year -> it overstates historical organic area. Downweight the organic FRACTION by the
# EU organic-share ratio vs the master's reference year; the area not carved as organic stays
# conventional (area-conserving). Mirrors the organic-livestock 2000=marginal / 2010=backcast rule.
ORGANIC_AREA_REF_YEAR <- 2020
.eu_org_share <- c("2000" = 3.0, "2010" = 5.2, "2018" = 7.7, "2020" = 9.1) # % of UAA (Eurostat)
organic_area_weight <- function(yr) {
  s <- .eu_org_share[as.character(yr)]
  r <- .eu_org_share[as.character(ORGANIC_AREA_REF_YEAR)]
  if (is.na(s) || is.na(r)) 1 else min(s / r, 1)
}
# Within-admin-unit HETEROGENEITY drivers: alongside the weighted MEAN of each indicator, also
# aggregate its weighted SD over the 1km pixels (a unit mixing lowland+upland grazes differently
# than a uniform one). Added as <var>_sd covariates; the horseshoe selects the informative ones.
INCLUDE_SD_DRIVERS <- as.logical(Sys.getenv("DRIVER_INCLUDE_SD", "TRUE"))
SD_VARS <- strsplit(Sys.getenv(
  "DRIVER_SD_VARS",
  "Slope_rad,Elevation,GHM_HI,GHM_TI,CISI,Growing_Degree_Days_gdd5,Annual_Precipitation_bio12"
), ",")[[1]]   # intersected with the vars actually present, so GHM_TI is a no-op until it lands

# Derived terrain SHARE features: nonlinear summaries the mean Elevation/Aspect can't capture --
# the fraction of the unit that is flat / steep / lowland / upland (e.g. "flat grazeable area").
# Computed per 1km pixel by threshold, then weighted-averaged to the admin unit.
INCLUDE_TERRAIN_FEATURES <- as.logical(Sys.getenv("DRIVER_INCLUDE_TERRAIN", "TRUE"))
SLOPE_FLAT_RAD <- as.numeric(Sys.getenv("DRIVER_SLOPE_FLAT_RAD", "0.087")) # < ~5deg  -> flat
SLOPE_STEEP_RAD <- as.numeric(Sys.getenv("DRIVER_SLOPE_STEEP_RAD", "0.21")) # > ~12deg -> steep
ELEV_LOW_M <- as.numeric(Sys.getenv("DRIVER_ELEV_LOW_M", "300")) # < 300m   -> lowland
ELEV_HIGH_M <- as.numeric(Sys.getenv("DRIVER_ELEV_HIGH_M", "1000")) # > 1000m  -> upland

NO_CHOICE_LU <- c("Waterbodies_marine", "Waterbodies_inland", "Wetlands_natural", "Natural_other", "NODATA")

SOURCE_REGISTRY <- list(
  LUM = list(
    # The maintained copy lives in THIS repo. The count model was reading an external one in
    # ../LAMASUS_downscaling that is a frozen older export: same 92 LUM codes, and every shared column
    # byte-identical except GLOBIOM_subclass, where the external file collapses the five urban classes
    # into a single "Urban". It does NOT carry BMLEH_Los1_label, so pointing the count model at a
    # project class column failed with "CLASS_COLS missing" -- correct, but for a confusing reason.
    # Verified before switching: LUM_label, GLOBIOM_UNFCCC, GLOBIOM_mngmt, BIOCLIMA3_*, ETL2_*,
    # BIOCLIMA_DS_*, AgMIP_label all identical across the two.
    mapping_file = Sys.getenv("DRIVER_MAPPING_FILE", "aux_files/LUM_Code_to_macro_model_mapping.csv"),
    join_key = "LUM_Code", base_lu = "GLOBIOM_UNFCCC"
  )
)
.src <- SOURCE_REGISTRY[[OUTCOME_SOURCE]]
if (is.null(.src)) stop("Unknown OUTCOME_SOURCE: ", OUTCOME_SOURCE)

derive_model_class <- function(map_dt, class_cols = CLASS_COLS, base_lu = .src$base_lu) {
  m <- as.data.table(copy(map_dt))
  miss <- setdiff(class_cols, names(m))
  if (length(miss)) stop(sprintf("CLASS_COLS missing: %s", paste(miss, collapse = ", ")))

  mat <- sapply(class_cols, function(cc) trimws(as.character(m[[cc]])))
  mat <- matrix(mat, nrow = nrow(m))
  base <- apply(mat, 1, function(r) paste(r[nzchar(r) & !is.na(r)], collapse = "_"))
  m[, model_class := fifelse(base == "" | is.na(base), "NODATA", make.names(base))]

  lbl <- if ("LUM_label" %in% names(m)) tolower(as.character(m$LUM_label)) else rep("", nrow(m))
  if (INCLUDE_ORGANIC) {
    if (nrow(m[grepl("organic", lbl) & model_class == "NODATA"]) > 0) stop("INCLUDE_ORGANIC=TRUE but CLASS_COLS are blank.")
    # GLOBIOM_mngmt already encodes organic (HIO/LIO/other_O), so appending _O here double-labels
    # (-> Cropland_HIO_O). Only append when the class doesn't already carry the organic marker.
    m[
      grepl("organic", lbl) & model_class != "NODATA" & !grepl("O$", gsub("_", "", model_class)),
      model_class := paste0(model_class, "_O")
    ]
  }
  if (INCLUDE_IRRIGATION) {
    if (nrow(m[grepl("irrigat", lbl) & model_class == "NODATA"]) > 0) stop("INCLUDE_IRRIGATION=TRUE but CLASS_COLS are blank.")
    m[grepl("irrigat", lbl) & model_class != "NODATA", model_class := paste0(model_class, "_IR")]
  }
  m[, focal_class := model_class]
  if (base_lu %in% names(m)) m[get(base_lu) %in% NO_CHOICE_LU, model_class := "no_choice"]
  m[]
}

# ---- Load mapping + derive the target class set ----
mapping_thematic <- as.data.table(fread(.src$mapping_file))
if ("LUM_label" %in% names(mapping_thematic)) {
  if (!INCLUDE_IRRIGATION) mapping_thematic <- mapping_thematic[!grepl("irrigat", tolower(LUM_label))]
  if (!INCLUDE_ORGANIC) mapping_thematic <- mapping_thematic[!grepl("organic", tolower(LUM_label))]
}
mapping_thematic <- derive_model_class(mapping_thematic)
class_lookup <- unique(mapping_thematic[, c(.src$join_key, "model_class", "focal_class"), with = FALSE])


# Model Settings (env-driven, matching run_prior_module_pixel_level_model.R conventions)
N_CHAINS <- as.integer(Sys.getenv("DRIVER_NCHAINS", "4"))
use_re <- TRUE # Random Effects enabled
RE_GROUP_REQUEST <- Sys.getenv("DRIVER_RE_GROUP_COL", "GLOB_country")  # "GLOB_country", "CAPRI_NUTS", ...
# CAPRI_NUTS is a REGION code ("DE11", ~246 of them), but the RE group here has to be a COUNTRY: the
# downstream parameter table (build_livestock_capri.R) is keyed by country, and an RE over ~246
# regions is a different model, not a re-keying of this one. So asking for CAPRI_NUTS means "group on
# the CAPRI country", i.e. its 2-char prefix -- the same slice the pixel driver applies and records as
# re_group_sliced. Resolve it to the real column name here so every artifact, path and diagnostic
# records the key that was ACTUALLY used rather than the one that was asked for.
RE_GROUP_COL <- if (identical(RE_GROUP_REQUEST, "CAPRI_NUTS")) "CAPRI_country" else RE_GROUP_REQUEST
if (!identical(RE_GROUP_COL, RE_GROUP_REQUEST))
  cat(sprintf(">>> RE group: %s requested -> using %s (CAPRI_NUTS sliced to its 2-char country)\n",
              RE_GROUP_REQUEST, RE_GROUP_COL))
# POOLED RE VARIANCE across the outcome columns (one scale per predictor instead of one per
# (predictor, column)). Ported from the MNL's category-pooled RE block, but UNCENTRED: the count
# columns are independent rates with no baseline and no zero-sum constraint, so the category-centring
# that zeroed the MNL's random effects has no meaning here. Off by default (bit-identical: verified
# max|beta diff| = max|sigma diff| = 0.000e+00 vs the pre-port sampler).
COUNT_RE_PREC_POOLED <- as.logical(Sys.getenv("COUNT_RE_PREC_POOLED", "TRUE"))
COUNT_USE_COUNTRY_SHRINKAGE <- as.logical(Sys.getenv("COUNT_USE_COUNTRY_SHRINKAGE", "TRUE"))
COUNT_SLAB_C2_COUNTRY <- as.numeric(Sys.getenv("COUNT_SLAB_C2_COUNTRY", "4.0"))
use_bart <- FALSE
n_trees_bart <- 20
target_livestock <- c("BOV", "SGT") # BOV = bovines (cattle); SGT = small grazers (sheep + goats). Fit separately. (defined early for path setup)
# After the totals, also fit the subtype COMPOSITION (D/O/F δ per species) so one run produces the
# full pipeline. Run as separate R processes (the composition uses a different sampler / C++ core).
FIT_COMPOSITION <- as.logical(Sys.getenv("DRIVER_FIT_COMPOSITION", "TRUE"))
# DOWNSCALING prediction grid: build the covariate space X at EEA_10kmID x NUTS3 (the 10km x NUTS3
# allocation units), reusing the EXACT covariate aggregation so the fitted coefficients apply. There
# is no outcome at this resolution -> covariate-only, no fitting. "" = normal NUTS3 run.
DOWNSCALE_GRID <- Sys.getenv("DRIVER_DOWNSCALE_GRID", "") # "" = off; "10km" = build the grid
DOWNSCALE_MAP <- Sys.getenv(
  "DRIVER_DOWNSCALE_MAP",
  mapping_file <- get_latest_file(AUXDATA_DIR, "^one_kmID_master_mapping_.*\\.parquet$")
)

# --- Count unit (LSU per modeled count) ---
# Raw LSU totals reach ~1.3M. At those magnitudes the NB dispersion r is statistically
# UNIDENTIFIED (the variance can't pin it), so r is prior-driven and rides a flat ridge with
# the intercept -> poor mixing (total-loglik ESS ~11, GDP ESS ~14, intercept Rhat ~1.11).
# Modeling livestock in COARSER units (round(LSU / LSU_PER_COUNT)) moves the counts into the
# regime where r IS identified: every weak axis mixes 3-48x better (validated: GDP ESS 14->671,
# total-loglik ESS 11->428, intercept Rhat ->1.01) and the fit is preserved/improved. This is a
# units/aggregation choice -- still ABSOLUTE counts with the log-area offset (NOT shares); the
# slope relationships (what feed the prior) are unchanged, only the intercept shifts by -log(unit).
# Default 100 is safe for both BOV (mean ~620) and the smaller SGT (mean ~105); raise toward 1000
# for even better BOV mixing, lower to 1 to recover raw-LSU counts. See count-offset-bug memory.
LSU_PER_COUNT <- as.numeric(Sys.getenv("DRIVER_LSU_PER_COUNT", "100"))

# Pre-create the MODEL_LABEL-derived diagnostic dirs (avoids runtime cloud-lock issues). Names match
# assess_convergence's tolower(gsub(" ","_", name)), name = paste0(MODEL_LABEL,"_",tgt,"_",RE/pooled).
for (.tgt in target_livestock) {
  .nm <- tolower(gsub(" ", "_", paste0(MODEL_LABEL, "_", .tgt, "_", if (use_re) "RE" else "pooled")))
  dir.create(file.path("output/diagnostics", .nm), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("output/plots/diagnostics", .nm), recursive = TRUE, showWarnings = FALSE)
}

# Promote a finished test run to the production bucket (matches the pixel driver): copy every artifact
# tagged with this test MODEL_LABEL -- per-target saved-model dirs + flat output files + diagnostics
# dirs (which carry the lowercased label) -- to the "_production" twin, overwriting. Only in test mode.
promote_test_to_production <- function() {
  prod_label <- sub("_test$", "_production", MODEL_LABEL)
  # 1. per-target saved-model state dirs
  for (tgt in target_livestock) {
    src <- file.path("output/saved_model_outputs", paste0(MODEL_LABEL, "_", tgt, if (use_re) "_RE_" else "_pooled_", RE_GROUP_COL))
    if (dir.exists(src)) {
      dst <- sub(MODEL_LABEL, prod_label, src, fixed = TRUE)
      if (dir.exists(dst)) unlink(dst, recursive = TRUE)
      dir.create(dst, recursive = TRUE, showWarnings = FALSE)
      file.copy(list.files(src, full.names = TRUE), dst, recursive = TRUE, overwrite = TRUE)
    }
  }
  # 2. flat output files carrying this MODEL_LABEL (dat_admin, checkpoints-as-files, CSVs)
  cand <- list.files("output", full.names = TRUE, recursive = FALSE)
  cand <- cand[!dir.exists(cand) & grepl(MODEL_LABEL, basename(cand), fixed = TRUE)]
  nf <- 0L
  for (f in cand) {
    dst <- file.path(dirname(f), sub(MODEL_LABEL, prod_label, basename(f), fixed = TRUE))
    if (!identical(dst, f) && file.copy(f, dst, overwrite = TRUE)) nf <- nf + 1L
  }
  # 3. diagnostics + plot dirs (named with the lowercased label)
  lo <- tolower(MODEL_LABEL)
  lo_prod <- tolower(prod_label)
  for (base in c("output/diagnostics", "output/plots/diagnostics")) {
    for (d in list.dirs(base, recursive = FALSE)) {
      if (grepl(lo, basename(d), fixed = TRUE)) {
        dst <- file.path(base, sub(lo, lo_prod, basename(d), fixed = TRUE))
        if (dir.exists(dst)) unlink(dst, recursive = TRUE)
        dir.create(dst, recursive = TRUE, showWarnings = FALSE)
        file.copy(list.files(d, full.names = TRUE), dst, recursive = TRUE, overwrite = TRUE)
      }
    }
  }
  cat(sprintf(">>> Promoted to production: saved-model dirs + %d flat file(s) + diagnostics  [%s -> %s]\n", nf, MODEL_LABEL, prod_label))
}
niter <- as.integer(Sys.getenv("DRIVER_NITER", "12000")) # 12000-2000 = 10000 retained: the dispersion r and RE variances need it (audit: r ESS was 4 at 500 retained). Override via env for a quick test.
nburn <- as.integer(Sys.getenv("DRIVER_NBURN", "4000"))
# RAM controls (the fit phase held ALL ~8000 retained draws x N_CHAINS in RAM -> crashes). thin stores every
# k-th draw (posterior arrays shrink ~k x; k<<autocorr-time keeps ESS); calc_loo OFF drops the [nretain x n]
# pointwise matrix + loo PSIS (a heavy diagnostic, not needed for the prior).
thin_keep <- max(1L, as.integer(Sys.getenv("DRIVER_THIN", "5")))
calc_loo_flag <- as.logical(Sys.getenv("DRIVER_CALC_LOO", "FALSE"))

# Parallel Settings
N_CORES <- min(parallel::detectCores() - 1, N_CHAINS)
if (.Platform$OS.type == "windows") {
  plan(multisession, workers = N_CORES)
  cat(sprintf("Parallel Plan: %s with %d workers (Limit: 10GB)\n", "multisession", N_CORES))
} else {
  plan(multicore, workers = N_CORES)
  cat(sprintf("Parallel Plan: %s with %d workers (Limit: 10GB)\n", "multicore", N_CORES))
}
options(future.globals.maxSize = 10 * 1024^3)

# --- Year Configuration ---
# Panel years (env-driven). Default is a SINGLE timestep (2010) for diagnosing: 2010 has BOTH
# the LSU outcome AND CLC land-use covariates (the CLC panel only runs to 2018, so 2020 carries
# no LU). Set DRIVER_YEARS="2000,2010,2020" to restore the full multi-year panel.
MODEL_YEARS <- as.integer(strsplit(Sys.getenv("DRIVER_YEARS", "2000,2010,2020"), ",")[[1]])
COV_YEARS <- MODEL_YEARS # covariate years match outcome years

# --- Source Functions ---
source("codes/mnl_aux_func.R")
source("codes/count_rcpp.R")

# =========================================================================
# 2. SOURCE DATA LOADING
# =========================================================================
cat("\nLoading 1km Source Data & LSU Dataset...\n")

cat("  Loading prior covariates...\n")
# COLUMN-PROJECTED read: the master parquet is 1.3GB with many columns, but the count model uses only ~25
# (terrain/soil/climate + year-varying Pop/GDP/GHM_HI). Reading ONLY those cuts peak memory drastically AND
# avoids streaming the whole file from OneDrive (the OOM + mmap/read-timeout culprit). any_of() = drop missing
# gracefully. Zero effect on results: the dropped columns were never touched downstream.
# RAI (rural accessibility) REMOVED 2026-08-04 per spec; the accessibility/human-footprint signal is
# carried by the GHM v3 threat groups instead: GHM_HI (human intrusion) + GHM_TI (transport
# infrastructure). GHM_TI needs the upstream gridwork pull (download_ghm.R with selected_threat="TI"
# -> aggregate_ee_data.R -> merge_all_variables.R); any_of() means the driver runs fine before it
# lands and picks it up automatically afterwards -- the presence check below reports which arrived.
.pq_nm <- names(arrow::open_dataset(PRIOR_1KM_PARQUET)$schema)
# resolve each wanted GHM threat group x year to its physical column (clean or raw EE name)
GHM_CODES <- trimws(strsplit(Sys.getenv("DRIVER_GHM_CODES", "HI,TI"), ",")[[1]])
.ghm_map <- do.call(rbind, lapply(GHM_CODES, function(cd) data.frame(
  code = cd, year = COV_YEARS,
  phys = vapply(COV_YEARS, function(y) .ghm_resolve(.pq_nm, cd, y), character(1)),
  clean = paste0("GHM_", cd, "_", COV_YEARS), stringsAsFactors = FALSE)))
.ghm_map <- .ghm_map[!is.na(.ghm_map$phys), , drop = FALSE]

.needed_1km <- unique(c("INSPIRE_Europe_buffer_1kmID",
  "Slope_rad", "Elevation", "Aspect_cos_mean", "Aspect_sin_mean", "allPA_share", "CISI",
  "OC_TOP", "ROO", "AWC_TOP", "VS",
  "Growing_Degree_Days_gdd5", "Precipitation_Seasonality_bio15", "Annual_Precipitation_bio12",
  paste0("Pop_", COV_YEARS), paste0("GDP_", COV_YEARS), .ghm_map$phys))
prior_1km <- as.data.table(dplyr::collect(dplyr::select(
  arrow::open_dataset(PRIOR_1KM_PARQUET),
  dplyr::any_of(.needed_1km))))
# normalize raw EE export names -> GHM_<CODE>_<year>
.ren <- .ghm_map[.ghm_map$phys != .ghm_map$clean & .ghm_map$phys %in% names(prior_1km), , drop = FALSE]
if (nrow(.ren)) { setnames(prior_1km, .ren$phys, .ren$clean)
  cat(sprintf(">>> renamed %d raw GHM export column(s) -> GHM_<CODE>_<year> (e.g. %s -> %s)\n",
              nrow(.ren), .ren$phys[1], .ren$clean[1])) }
cat(sprintf(">>> prior_1km: read %d of the parquet's columns (column-projected) x %d rows\n", ncol(prior_1km), nrow(prior_1km)))
# Which GHM threat groups actually arrived? GHM_TI is optional until the gridwork pull is done.
.ghm_want <- paste0("GHM_", GHM_CODES)
GHM_VARS <- .ghm_want[vapply(.ghm_want,
  function(v) all(paste0(v, "_", COV_YEARS) %in% names(prior_1km)), logical(1))]
.ghm_missing <- setdiff(.ghm_want, GHM_VARS)
cat(sprintf(">>> GHM threat groups available for all COV_YEARS: %s%s\n",
  if (length(GHM_VARS)) paste(GHM_VARS, collapse = ", ") else "(none)",
  if (length(.ghm_missing)) sprintf("   [MISSING: %s -> not modelled; run gridwork download_ghm.R for it]",
                                    paste(.ghm_missing, collapse = ", ")) else ""))
if (!length(GHM_VARS)) stop("No GHM threat group has full year coverage -- check the master parquet.")
prior_1km[, INSPIRE_Europe_buffer_1kmID := as.integer(INSPIRE_Europe_buffer_1kmID)]
for (.col in names(prior_1km)[sapply(prior_1km, is.double)]) {
  set(prior_1km, j = .col, value = fifelse(is.nan(prior_1km[[.col]]), NA_real_, prior_1km[[.col]]))
}

cat("  Loading grid mapping...\n")
mapping_file <- get_latest_file(AUXDATA_DIR, "^one_kmID_master_mapping_.*\\.parquet$")
grid_map_raw <- arrow::read_parquet(file.path(AUXDATA_DIR, mapping_file)) %>% as.data.table()

cat("  Loading LSU processed counts...\n")
lsu_data <- readRDS(file.path(GRIDWORK_DIR, "LSU_data_processed_with_POL_20250701.rds")) %>% as.data.table()
lsu_data[, INSPIRE_Europe_buffer_1kmID := as.integer(INSPIRE_Europe_buffer_1kmID)]

# Merge LSU with mapping to get NUTS3 and NUTS0
# Enforce a strictly nested unique mapping from NUTS3 to NUTS0 and GLOB_country
# by assigning NUTS3 regions to their dominant country pixel count. This eliminates
# multi-country overlap anomalies on borders, ensuring perfect multi-level nesting.
nuts3_lookup <- grid_map_raw[!is.na(NUTS3), .N, by = .(NUTS3, NUTS0, GLOB_country)][
  order(-N), .(
    NUTS0 = NUTS0[1],
    GLOB_country = GLOB_country[1]
  ),
  by = NUTS3
]

# CAPRI keys by the SAME dominant rule -- but computed SEPARATELY and merged, NOT by adding
# CAPRI_NUTS to the `by` above. Adding it there would split the pixel counts finer and could flip
# which country a border NUTS3 is assigned to, silently changing the existing GLOB_country grouping.
# This is a real assignment rather than a lookup: 411 of 1432 NUTS3 units touch more than one CAPRI
# country at 1km, because the two geometries disagree at borders.
if ("CAPRI_NUTS" %in% names(grid_map_raw)) {
  .capri_lookup <- grid_map_raw[!is.na(NUTS3) & !is.na(CAPRI_NUTS) & nzchar(as.character(CAPRI_NUTS)),
                                .N, by = .(NUTS3, CAPRI_NUTS)][
    order(-N), .(CAPRI_NUTS = as.character(CAPRI_NUTS)[1]), by = NUTS3]
  nuts3_lookup <- merge(nuts3_lookup, .capri_lookup, by = "NUTS3", all.x = TRUE)
  nuts3_lookup[, CAPRI_country := substr(CAPRI_NUTS, 1, 2)]
  cat(sprintf(">>> CAPRI keys: %d NUTS3 -> %d CAPRI regions -> %d CAPRI countries (%d NUTS3 with no CAPRI match)\n",
              nrow(nuts3_lookup), uniqueN(nuts3_lookup$CAPRI_NUTS), uniqueN(nuts3_lookup$CAPRI_country),
              sum(is.na(nuts3_lookup$CAPRI_NUTS))))
} else if (identical(RE_GROUP_COL, "CAPRI_country")) {
  stop("DRIVER_RE_GROUP_COL=CAPRI_NUTS, but the grid mapping carries no CAPRI_NUTS column: ", mapping_file)
}

lsu_data <- merge(lsu_data, grid_map_raw[, .(
  INSPIRE_Europe_buffer_1kmID = as.integer(INSPIRE_Europe_buffer_1kmID), 
  EEA_1kmID = as.integer(EEA_1kmID),
  EEA_10kmID,
  NUTS3,
  pixel_weight = if ("GLOB_5arcminID_area_km2" %in% names(grid_map_raw)) as.numeric(GLOB_5arcminID_area_km2) else 1.0
)], by = "INSPIRE_Europe_buffer_1kmID", all.x = TRUE)

# Drop missing mappings if any
lsu_data <- lsu_data[!is.na(NUTS3)]
lsu_data[is.na(pixel_weight) | pixel_weight == 0, pixel_weight := 1.0]

# Merge with the strictly nested country lookup
lsu_data <- merge(lsu_data, nuts3_lookup, by = "NUTS3", all.x = TRUE)

# Grouping key: NUTS3 normally; EEA_10kmID x NUTS3 for the downscaling prediction grid.
if (nzchar(DOWNSCALE_GRID)) {
  .nb <- nrow(lsu_data); lsu_data <- lsu_data[!is.na(EEA_10kmID)]    # drop the ~3% pixels not in the 10km mapping
  lsu_data[, RESOLUTION := paste0(EEA_10kmID, "_", as.character(NUTS3))]
  cat(sprintf(">>> DOWNSCALE_GRID=%s: %d/%d 1km pixels mapped -> %d units (EEA_10kmID x NUTS3)\n",
              DOWNSCALE_GRID, nrow(lsu_data), .nb, uniqueN(lsu_data$RESOLUTION)))
} else {
  lsu_data[, RESOLUTION := as.character(NUTS3)]
}

# Observation Unit mapping (1km pixel -> Admin Boundary)
# CAPRI_country rides along so it can serve as the RE key; it is nested under NUTS3 by construction
# (assigned from it above), so carrying it changes neither frame's row count.
.capri_cols <- intersect(c("CAPRI_NUTS", "CAPRI_country"), names(lsu_data))
grid_map_pixel <- unique(lsu_data[, c("INSPIRE_Europe_buffer_1kmID", "EEA_1kmID", "RESOLUTION", "NUTS0",
                                      "GLOB_country", "NUTS3", "pixel_weight", .capri_cols), with = FALSE])

# Base Admin unit frame (for joining later)
.af_before <- uniqueN(lsu_data$RESOLUTION)
admin_frame <- unique(lsu_data[, c("RESOLUTION", "NUTS0", "GLOB_country", "NUTS3", .capri_cols), with = FALSE])
if (nrow(admin_frame) != .af_before)
  stop(sprintf("admin_frame is no longer one row per RESOLUTION (%d rows for %d units) -- the CAPRI keys are not nested under NUTS3 as assumed.",
               nrow(admin_frame), .af_before))

# =========================================================================
# 3. STATIC COVARIATES
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
  dcast_cols <- setdiff(colnames(res), "INSPIRE_Europe_buffer_1kmID")
  res <- res[grid_map_pixel, on = "INSPIRE_Europe_buffer_1kmID", nomatch = 0L]

  # Weight by pixel_weight for area-accurate aggregation
  res[, (dcast_cols) := lapply(.SD, function(x) x * pixel_weight), .SDcols = dcast_cols]

  # Sum the pixel counts to get Admin Boundary totals
  res <- res[, lapply(.SD, sum, na.rm = TRUE), by = .(RESOLUTION), .SDcols = dcast_cols]
  cls_cols <- setdiff(colnames(res), "RESOLUTION")
  res[, (cls_cols) := .SD / pmax(1e-9, rowSums(.SD, na.rm = TRUE)), .SDcols = cls_cols]
  setnames(res, cls_cols, paste0(v, "_s", cls_cols))
  return(res)
})
x_pixel_soil_list <- lapply(x_pixel_soil_list, function(dt) setkey(as.data.table(dt), RESOLUTION))
x_pixel_soil <- Reduce(function(x, y) merge(x, y, all = TRUE), x_pixel_soil_list)

soil_dev_list <- lapply(spatial_cat_vars, function(v) {
  p <- paste0(v, "_s")
  cols <- sort(grep(paste0("^", p), colnames(x_pixel_soil), value = TRUE))
  if (length(cols) < 2) {
    return(NULL)
  }
  ref_col <- cols[1]
  other_cols <- cols[-1]
  devs <- x_pixel_soil[, other_cols, with = FALSE] - x_pixel_soil[[ref_col]]
  setnames(devs, other_cols, gsub("_s", "_dev", other_cols))
  return(devs)
})
x_pixel_soil <- cbind(x_pixel_soil, do.call(cbind, soil_dev_list))

# =========================================================================
# 4. MULTI-YEAR PANEL ASSEMBLY
# =========================================================================
cat("\nBuilding Admin Boundary Level Panel...\n")

# LU composition as covariates, ALIGNED to the pixel model's classification.
# We will use the granular LUM map for the appropriate year, falling back to 2000 for years <= 2000 and 2018 for years > 2000.

# Static CHELSA climate covariates from prior_1km. (Yield indices removed per request; spei48
# excluded -- only spei48_2000/2018 in the parquet, not panel-consistent with 2000/2010/2020.)
climate_yield_cols <- c("Growing_Degree_Days_gdd5", "Precipitation_Seasonality_bio15", "Annual_Precipitation_bio12")

T_PAIRS <- lapply(seq_along(MODEL_YEARS), function(i) {
  list(out_year = MODEL_YEARS[i], cov_year = COV_YEARS[i])
})

dat_admin_list <- lapply(T_PAIRS, function(tp) {
  cat("  Processing Tier: Outcome =", tp$out_year, "| Covariates =", tp$cov_year, "\n")

  # 1. OUTCOME (Y) - sum from LSU dataset
  y_cols <- c(
    paste0("cattle_LSU_", tp$out_year),
    paste0("sheep_LSU_", tp$out_year),
    paste0("goat_LSU_", tp$out_year)
  )

  if (all(y_cols %in% colnames(lsu_data))) {
    y_admin <- lsu_data[, .(
      BOV = sum(get(y_cols[1]), na.rm = TRUE), # bovines (cattle)
      SGT = sum(get(y_cols[2]), na.rm = TRUE) + sum(get(y_cols[3]), na.rm = TRUE) # small grazers (sheep + goats)
    ), by = .(RESOLUTION)]
  } else {
    stop(sprintf("LSU counts not found for year %d", tp$out_year))
  }

  # 2. COVARIATES (X-Temporal)
  gdp_col <- paste0("GDP_", tp$cov_year)
  pop_col <- paste0("Pop_", tp$cov_year)

  x_1km_yr <- prior_1km[, .(INSPIRE_Europe_buffer_1kmID, Slope_rad, Elevation, Aspect_cos_mean, Aspect_sin_mean, allPA_share, CISI,
    Pop = get(pop_col), GDP = get(gdp_col)
  )]
  # year-varying GHM threat groups (GHM_HI always, GHM_TI once the gridwork pull lands)
  for (.g in GHM_VARS) x_1km_yr[[.g]] <- prior_1km[[paste0(.g, "_", tp$cov_year)]]
  # Append the pixel model's STATIC climate (CHELSA) + yield covariates, if present. Intensive ->
  # weighted mean below. spei48 is skipped: only spei48_2000/2018 exist (not the 2000/2010/2020 panel).
  static_extra <- intersect(climate_yield_cols, names(prior_1km))
  if (length(static_extra)) x_1km_yr <- cbind(x_1km_yr, prior_1km[, ..static_extra])

  cont_vars <- setdiff(colnames(x_1km_yr), "INSPIRE_Europe_buffer_1kmID")
  # allPA_share is a per-pixel protected fraction -> SUM it (x pixel_weight) to the total protected
  # AREA (km2) of the unit, like the LU areas, rather than the unit-mean share.
  mean_vars <- setdiff(cont_vars, c("Pop", "GDP", "allPA_share")) # intensive -> weighted mean
  sum_vars <- intersect(cont_vars, c("Pop", "GDP", "allPA_share")) # extensive (areas/counts) -> weighted sum

  # Intensive variables -> Weighted Mean
  x_mean <- x_1km_yr[grid_map_pixel, on = "INSPIRE_Europe_buffer_1kmID", nomatch = 0L][
    , lapply(.SD, weighted.mean, w = pixel_weight, na.rm = TRUE),
    by = .(RESOLUTION), .SDcols = mean_vars
  ]

  # Extensive variables -> Weighted Sum (adjusting for split pixels)
  x_sum <- x_1km_yr[grid_map_pixel, on = "INSPIRE_Europe_buffer_1kmID", nomatch = 0L][
    , lapply(.SD, function(x) sum(x * pixel_weight, na.rm = TRUE)),
    by = .(RESOLUTION), .SDcols = sum_vars
  ]
  if ("allPA_share" %in% names(x_sum)) setnames(x_sum, "allPA_share", "allPA_area")  # summed share = area (km2)

  x_admin_yr_cont <- merge(x_mean, x_sum, by = "RESOLUTION")

  # Within-unit weighted SD (heterogeneity) for the configured indicators -> <var>_sd
  if (INCLUDE_SD_DRIVERS) {
    .wsd <- function(x, w) {
      ok <- is.finite(x) & is.finite(w) & w > 0
      if (sum(ok) < 2) {
        return(0)
      }
      wm <- sum(w[ok] * x[ok]) / sum(w[ok])
      sqrt(pmax(sum(w[ok] * (x[ok] - wm)^2) / sum(w[ok]), 0))
    }
    sd_src <- intersect(SD_VARS, mean_vars)
    if (length(sd_src)) {
      x_sd <- x_1km_yr[grid_map_pixel, on = "INSPIRE_Europe_buffer_1kmID", nomatch = 0L][
        , lapply(.SD, .wsd, w = pixel_weight),
        by = .(RESOLUTION), .SDcols = sd_src
      ]
      setnames(x_sd, sd_src, paste0(sd_src, "_sd"))
      x_admin_yr_cont <- merge(x_admin_yr_cont, x_sd, by = "RESOLUTION")
    }
  }

  # Derived terrain SHARE features (flat / steep / lowland / upland fraction of the unit)
  if (INCLUDE_TERRAIN_FEATURES && all(c("Slope_rad", "Elevation") %in% names(x_1km_yr))) {
    .tf <- x_1km_yr[grid_map_pixel, on = "INSPIRE_Europe_buffer_1kmID", nomatch = 0L][
      , .(
        flat_share = weighted.mean(Slope_rad < SLOPE_FLAT_RAD, pixel_weight, na.rm = TRUE),
        steep_share = weighted.mean(Slope_rad > SLOPE_STEEP_RAD, pixel_weight, na.rm = TRUE),
        lowland_share = weighted.mean(Elevation < ELEV_LOW_M, pixel_weight, na.rm = TRUE),
        upland_share = weighted.mean(Elevation > ELEV_HIGH_M, pixel_weight, na.rm = TRUE)
      ),
      by = .(RESOLUTION)
    ]
    x_admin_yr_cont <- merge(x_admin_yr_cont, .tf, by = "RESOLUTION")
  }

  # Total Area for offset (Sum of pixel weights in the admin boundary)
  area_df <- grid_map_pixel[, .(total_area_km2 = sum(pixel_weight, na.rm = TRUE)), by = .(RESOLUTION)]

  # 3. LU COVARIATES — ABSOLUTE AREAS per LU class (km^2), NOT shares.
  # Admin units are non-uniform in size, so absolute LU areas (with the log_area_offset handling
  # exposure) preserve the size signal that shares would discard. Classes = model_class (pixel-aligned).
  # Use the tier's OWN year when a map exists for it, falling back only when it genuinely does not.
  .lum_dir <- file.path(CASCADE_DATA, "02_intermediate")
  .lum_for <- function(y) file.path(.lum_dir, sprintf("LUM_fit_with_energy_levels_and_new_FM_%d_EEA_1kmID.rds", y))
  lum_year <- if (file.exists(.lum_for(tp$cov_year))) tp$cov_year else if (tp$cov_year <= 2000) 2000 else 2018
  if (lum_year != tp$cov_year)
    cat(sprintf("  NOTE: no LUM map for %d -> using %d (covariates lead/lag the outcome by %d years)\n",
                tp$cov_year, lum_year, lum_year - tp$cov_year))
  lum_path <- .lum_for(lum_year)
  lu_areas <- NULL

  if (file.exists(lum_path)) {
    cat(sprintf("  Loading granular %d LUM map for covariates...\n", lum_year))
    temp_map_1km_raw <- as.data.table(readRDS(lum_path))
    temp_map_1km_raw[, EEA_1kmID := as.integer(EEA_1kmID)]
    temp_map_1km_raw <- temp_map_1km_raw[area_km2 > 0 & EEA_1kmID %in% unique(grid_map_pixel$EEA_1kmID)]

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
    temp_lum_wide_1km <- dcast(temp_map_1km_raw, EEA_1kmID ~ get(paste0("LUM_fit_with_energy_levels_and_new_FM_", lum_year)), value.var = "area_km2", fill = 0)
    setnames(temp_lum_wide_1km, old = names(temp_lum_wide_1km)[-1], new = paste0("LUM", names(temp_lum_wide_1km)[-1]))

    # --- Irrigation intersection
    if (INCLUDE_IRRIGATION) {
      temp_lum_wide_1km <- merge(temp_lum_wide_1km, temp_binary_irrigation_map, by = "EEA_1kmID", all.x = TRUE)
      temp_lum_wide_1km[is.na(irrigation_binary), irrigation_binary := 0]
      temp_lum_wide_1km[is.na(area_km2), area_km2 := 0]
      temp_non_irrigated_codes <- mapping_thematic[GLOBIOM_UNFCCC == "Cropland" & GLOB_CSYS != "IR" & !grepl("O", GLOB_CSYS), LUM_Code]
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

    # --- Organic intersection
    if (INCLUDE_ORGANIC) {
      temp_final_dt_1km <- merge(temp_lum_wide_1km, temp_organic_map, by = "EEA_1kmID", all.x = TRUE)
      for (j in org_cols) set(temp_final_dt_1km, which(is.na(temp_final_dt_1km[[j]])), j, 0)

      arable_cols <- intersect(c("LUM2001", "LUM2002", "LUM2003", "LUM2004", "LUM2005", if (INCLUDE_IRRIGATION) "LUM2006"), names(temp_final_dt_1km))
      permanent_cols <- intersect(c("LUM3001", "LUM3002", "LUM3003", "LUM3004", "LUM3005"), names(temp_final_dt_1km))
      pasture_cols <- intersect(c(paste0("LUM", 14:26), "LUM4012", "LUM4013"), names(temp_final_dt_1km))

      temp_final_dt_1km[, total_arable_lum := rowSums(.SD, na.rm = TRUE), .SDcols = arable_cols]
      temp_final_dt_1km[, total_permanent_lum := rowSums(.SD, na.rm = TRUE), .SDcols = permanent_cols]
      temp_final_dt_1km[, total_pasture_lum := rowSums(.SD, na.rm = TRUE), .SDcols = pasture_cols]
      temp_final_dt_1km[, total_agri_lum := total_arable_lum + total_permanent_lum + total_pasture_lum]

      temp_final_dt_1km[, Cropland_organic := ifelse(total_arable_lum + total_permanent_lum > 0, pmin(Cropland_organic / (total_arable_lum + total_permanent_lum), 1.0), 0)]
      temp_final_dt_1km[, Livestock_organic := ifelse(total_pasture_lum > 0, pmin(Livestock_organic / total_pasture_lum, 1.0), 0)]
      temp_final_dt_1km[, Mixed_organic := ifelse(total_agri_lum > 0, pmin(Mixed_organic / total_agri_lum, 1.0), 0)]
      temp_final_dt_1km[, All_organic := ifelse(total_agri_lum > 0, pmin(All_organic / total_agri_lum, 1.0), 0)]

      # Historical organic downweight: scale the (present-day) organic fractions to the year's
      # organic extent. The un-organic remainder stays in the conventional base class -> area-conserving.
      .ow <- organic_area_weight(tp$cov_year)
      if (.ow < 1) {
        for (.j in c("Cropland_organic", "Livestock_organic", "Mixed_organic", "All_organic")) {
          set(temp_final_dt_1km, NULL, .j, temp_final_dt_1km[[.j]] * .ow)
        }
        cat(sprintf(
          "    Organic-area downweight for %d: x%.2f (EU organic share %.1f%% vs %d ref %.1f%%)\n",
          tp$cov_year, .ow, .eu_org_share[as.character(tp$cov_year)], ORGANIC_AREA_REF_YEAR,
          .eu_org_share[as.character(ORGANIC_AREA_REF_YEAR)]
        ))
      }

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
      for (col in neg_names) temp_final_dt_1km[get(col) < 0, (col) := 0]
    }

    # E. Aggregate to Admin Resolution (Area-Weighted by pixel_weight)
    temp_final_dt_1km <- merge(temp_final_dt_1km, grid_map_pixel[, .(EEA_1kmID, RESOLUTION, pixel_weight)], by = "EEA_1kmID", all.x = TRUE)
    temp_final_dt_1km[is.na(pixel_weight), pixel_weight := 1.0]
    lum_cols <- grep("^LUM", colnames(temp_final_dt_1km), value = TRUE)
    temp_final_dt_1km[, (lum_cols) := lapply(.SD, function(x) x * pixel_weight), .SDcols = lum_cols]

    temp_agg_wide <- temp_final_dt_1km[!is.na(RESOLUTION), lapply(.SD, sum, na.rm = TRUE), by = .(RESOLUTION), .SDcols = patterns("^LUM")]
    temp_compiled_long <- melt(temp_agg_wide, id.vars = c("RESOLUTION"), variable.name = "LUM_Code", value.name = "area_km2")
    temp_compiled_long[, LUM_Code := as.integer(sub("LUM", "", LUM_Code))]

    # Map raw LUM codes -> model_class via class_lookup
    temp_compiled_long <- merge(temp_compiled_long, class_lookup, by = "LUM_Code", all.x = TRUE)
    temp_compiled_long[is.na(model_class) | model_class == "NODATA", model_class := "no_choice"]

    lu_admin_wide <- dcast(temp_compiled_long, RESOLUTION ~ model_class, value.var = "area_km2", fill = 0, fun.aggregate = sum)

    # Tag with lu_area_ prefix
    lu_cols <- setdiff(colnames(lu_admin_wide), "RESOLUTION")
    setnames(lu_admin_wide, lu_cols, paste0("lu_area_", lu_cols))
    lu_areas <- lu_admin_wide
  }

  # Final Join for this Tier
  tier_dat <- merge(y_admin, x_admin_yr_cont, by = "RESOLUTION", all.x = TRUE)
  tier_dat <- merge(tier_dat, x_pixel_soil, by = "RESOLUTION", all.x = TRUE)
  tier_dat <- merge(tier_dat, area_df, by = "RESOLUTION", all.x = TRUE)
  if (!is.null(lu_areas)) {
    tier_dat <- merge(tier_dat, lu_areas, by = "RESOLUTION", all.x = TRUE)
  }
  tier_dat <- merge(tier_dat, admin_frame, by = "RESOLUTION", all.x = TRUE)

  tier_dat[, out_year := tp$out_year]
  tier_dat[, cov_year := tp$cov_year]

  # Eagerly free massive 1km intermediate maps for this year to prevent RAM bloat
  rm(list = intersect(c("temp_map_1km_raw", "temp_binary_irrigation_map", "temp_organic_map", 
                        "temp_lum_wide_1km", "temp_final_dt_1km", "temp_agg_wide", "temp_compiled_long", "lu_admin_wide"), ls()))
  gc(verbose = FALSE)

  return(tier_dat)
})

dat_admin <- rbindlist(dat_admin_list, fill = TRUE)

# =========================================================================
# 5. FILTERING AND MATRIX PREPARATION
# =========================================================================
skewed_vars <- c("GDP", "Pop", "allPA_area")   # right-skewed -> log1p (allPA_area = summed protected area; RAI removed 2026-08-04)
# climate/yield columns actually present after the admin aggregation (defensive intersect)
climate_yield_present <- intersect(climate_yield_cols, colnames(dat_admin))
# heterogeneity (SD) drivers actually built (present in dat_admin); empty if INCLUDE_SD_DRIVERS=FALSE
sd_cols <- if (INCLUDE_SD_DRIVERS) intersect(paste0(SD_VARS, "_sd"), colnames(dat_admin)) else character(0)
# derived terrain share features (flat/steep/lowland/upland)
terrain_cols <- if (INCLUDE_TERRAIN_FEATURES) intersect(c("flat_share", "steep_share", "lowland_share", "upland_share"), colnames(dat_admin)) else character(0)
spatial_cont_cols <- c("Slope_rad", "Elevation", "Aspect_cos_mean", "Aspect_sin_mean", GHM_VARS, "CISI", skewed_vars, climate_yield_present, sd_cols, terrain_cols)
soil_dev_cols <- character(0) # soil dummies removed from the covariate set per request (the _dev cols stay in dat_admin but are not modelled)
# absolute LU areas per class (pixel-aligned). Drop no_choice: it's the catch-all for the
# NO_CHOICE_LU classes (water/wetland/natural_other) + unmapped/NODATA land -- not a livestock
# driver, and lumping them dilutes the signal. The differentiated classes are kept.
lu_area_cols <- setdiff(grep("^lu_area_", colnames(dat_admin), value = TRUE), "lu_area_no_choice")
essential_cols <- c(spatial_cont_cols, soil_dev_cols, lu_area_cols, "total_area_km2")

area_before_cov_filter <- sum(dat_admin$total_area_km2, na.rm = TRUE)
rows_before <- nrow(dat_admin)

dat_admin <- dat_admin[complete.cases(dat_admin[, ..essential_cols])]
dat_admin[is.na(dat_admin)] <- 0

area_after_cov_filter <- sum(dat_admin$total_area_km2, na.rm = TRUE)
rows_after <- nrow(dat_admin)

cat(sprintf(
  "\n>>> Covariate Filtering Check:\n    Dropped %d rows and %.1f area due to missing essential covariates.\n",
  rows_before - rows_after, area_before_cov_filter - area_after_cov_filter
))

# =========================================================================
# GROUP COVERAGE FILTER
# =========================================================================
# Controls for border cases where large grouping aggregates (e.g. countries)
# have very low spatial coverage in the modeling dataset (due to border clipping, missing data, etc.),
# which can result in highly unstable group-specific random effects.
# If coverage is below 15% or total modeled area is below 100 km2, we filter the group out.

MIN_GROUP_COVERAGE <- 0.15 # Minimum 15% coverage of original physical extent
MIN_GROUP_AREA_KM2 <- 100.0 # Minimum 100 km2 of total modeled physical area

cat("\n>>> Applying Group Coverage Filter...\n")
# 1. Calculate original physical area per group (GLOB_country) from grid_map_pixel
orig_group_area <- grid_map_pixel[, .(orig_area = sum(pixel_weight, na.rm = TRUE)), by = .(country_key = GLOB_country)]

# We no longer need the large 1km level raw datasets now that dat_admin and orig_group_area are assembled.
# Eagerly free them from RAM before proceeding to matrix preparation and model estimation.
rm(list = intersect(c("prior_1km", "lsu_data", "grid_map_raw", "grid_map_pixel", "x_pixel_soil", "x_pixel_soil_list", "soil_dev_list"), ls()))
gc(verbose = FALSE)

# 2. Calculate modeled physical area per group in dat_admin (latest year to avoid inflation)
modeled_group_area <- dat_admin[cov_year == max(cov_year), .(modeled_area = sum(total_area_km2, na.rm = TRUE)), by = .(country_key = get(RE_GROUP_COL))]

# 3. Merge and compute coverage percent
group_coverage <- merge(modeled_group_area, orig_group_area, by = "country_key", all.y = TRUE)
group_coverage[is.na(modeled_area), modeled_area := 0.0]
group_coverage[, coverage_pct := modeled_area / orig_area]

# Identify low-coverage groups
groups_to_drop <- group_coverage[coverage_pct < MIN_GROUP_COVERAGE | modeled_area < MIN_GROUP_AREA_KM2]

if (nrow(groups_to_drop) > 0) {
  cat("    Dropped low-coverage or small-area grouping aggregates:\n")
  print(groups_to_drop[, .(country_key,
    modeled_area_km2 = round(modeled_area, 1),
    orig_area_km2 = round(orig_area, 1),
    coverage_pct = round(coverage_pct * 100, 1)
  )])

  # Filter dat_admin to drop these groups
  rows_before_group_filter <- nrow(dat_admin)
  dat_admin <- dat_admin[!get(RE_GROUP_COL) %in% groups_to_drop$country_key]
  cat(sprintf(
    "    Group Coverage Filter: Dropped %d rows belonging to low-coverage groups.\n",
    rows_before_group_filter - nrow(dat_admin)
  ))
} else {
  cat("    All grouping aggregates passed coverage and area thresholds.\n")
}

# Save FULL dataset for inspection
timestamp_str <- Sys.Date()
dat_admin_file <- paste0("output/dat_admin_FULL_", MODEL_LABEL, "_", timestamp_str, ".rds") # MODEL_LABEL-tagged, matching dat_pixel_file in the pixel driver
saveRDS(dat_admin, dat_admin_file)
cat(sprintf("\n>>> Saved FULL dat_admin to %s\n", dat_admin_file))

# Transform skewed and calculate offset
# The GRAZING-LAND offset. Matched by pattern, not by GLOBIOM's name: the same land is "Pasture_*"
# under GLOBIOM_UNFCCC+mngmt and "Grassland_*" under BMLEH_Los1_label. Hardcoding "Pasture" made the
# grep return NOTHING under BMLEH, so total_pasture_area became 0 and the offset collapsed to a
# CONSTANT log(1e-4) for every region. That is not a degraded model, it is a different one -- the
# log-area offset is what makes these absolute counts comparable across admin units of wildly
# different size, and losing it silently is the exact failure the count-offset bug already cost us
# once. Fail loudly instead.
GRAZE_PAT <- Sys.getenv("DRIVER_GRAZE_PATTERN", "^lu_area_(Pasture|Grassland)")
pasture_cols <- grep(GRAZE_PAT, lu_area_cols, value = TRUE)
if (!length(pasture_cols))
  stop(sprintf("no grazing-land columns match '%s' among: %s\n  The count model's offset IS grazing area; without it every region would get a constant offset.",
               GRAZE_PAT, paste(sub("lu_area_", "", lu_area_cols), collapse = ", ")))
cat(sprintf(">>> grazing-land offset built from %d column(s): %s\n",
            length(pasture_cols), paste(sub("lu_area_", "", pasture_cols), collapse = ", ")))
non_pasture_cols <- setdiff(lu_area_cols, pasture_cols)

dat_admin[, total_pasture_area := rowSums(.SD, na.rm = TRUE), .SDcols = pasture_cols]
dat_admin[, log_area_offset := log(pmax(total_pasture_area, 1e-4))]
if (dat_admin[, sd(log_area_offset, na.rm = TRUE)] < 1e-8)
  stop("the grazing-area offset is constant across regions -- it carries no information. Check the grazing columns above.")

for (.v in skewed_vars) {
  dat_admin[[.v]] <- log1p(dat_admin[[.v]])
  setnames(dat_admin, .v, paste0("log1p_", .v))
}

# Proportion method for Pasture
prop_pasture_cols <- character(0)
for (.v in pasture_cols) {
  prop_name <- sub("lu_area_", "prop_", .v)
  dat_admin[[prop_name]] <- dat_admin[[.v]] / pmax(dat_admin$total_pasture_area, 1e-6)
  prop_pasture_cols <- c(prop_pasture_cols, prop_name)
}

# Drop a reference category. Named by env or, failing that, the LARGEST grazing class by area --
# rather than a hardcoded "prop_Pasture_HI" that simply does not exist under another classification,
# in which case nothing was dropped and the proportions stayed collinear with the intercept.
ref_cat <- Sys.getenv("DRIVER_GRAZE_REF", "")
if (!nzchar(ref_cat) || !ref_cat %in% prop_pasture_cols) {
  .tot <- vapply(pasture_cols, function(v) sum(dat_admin[[v]], na.rm = TRUE), 0)
  ref_cat <- sub("lu_area_", "prop_", names(which.max(.tot)))
}
if (ref_cat %in% prop_pasture_cols) {
  prop_pasture_cols <- setdiff(prop_pasture_cols, ref_cat)
  cat(sprintf("\n>>> Proportion Method: Dropped %s to serve as the reference baseline.\n", ref_cat))
}

# log1p the absolute LU areas for non-pasture (Cropland, Forest, Urban, etc.)
for (.v in non_pasture_cols) {
  dat_admin[[.v]] <- log1p(dat_admin[[.v]])
  setnames(dat_admin, .v, paste0("log1p_", .v))
}

non_pasture_cols <- if (length(non_pasture_cols)) paste0("log1p_", non_pasture_cols) else character(0)

# Recombine for the final predictor matrix
lu_area_cols <- c(prop_pasture_cols, non_pasture_cols)

spatial_cont_cols_trans <- c(
  "Slope_rad", "Elevation", "Aspect_cos_mean", "Aspect_sin_mean", GHM_VARS, "CISI",
  paste0("log1p_", skewed_vars),   # incl. log1p_allPA_area (summed protected area)
  climate_yield_present, # static CHELSA climate + yields, untransformed (matches the pixel model)
  sd_cols, # within-unit heterogeneity (weighted SD), untransformed -> standardized in the sampler
  terrain_cols # derived terrain shares (flat/steep/lowland/upland), in [0,1]
)

# --- Year fixed effects (multi-year panel only) ---
# Reference (effects) coding: 2 dummies for 3 years (baseline = earliest year). They absorb the
# COMMON temporal level shifts between years so those shifts don't leak into the driver slopes
# (e.g. 2020 having both more livestock AND higher GDP would otherwise inflate the GDP slope and
# wrongly steer the spatial allocation). They are EXCLUDED from the horseshoe and from the country
# REs (global fixed time controls), and are 0/1 so the sampler's >2-unique rule leaves them
# un-standardized. PREDICT/DOWNSCALE: set them to their sample means (time-integrated, average-year
# allocation); a common year effect cancels anyway under renormalization to the GLOBIOM coarse totals.
year_dummy_cols <- character(0)
if (length(MODEL_YEARS) > 1L) {
  .yrs <- sort(unique(dat_admin$out_year))
  for (.yy in .yrs[-1L]) {
    .cn <- paste0("year_", .yy)
    dat_admin[[.cn]] <- as.numeric(dat_admin$out_year == .yy)
    year_dummy_cols <- c(year_dummy_cols, .cn)
  }
  cat(sprintf(
    ">>> Panel: added %d year dummy/-ies (ref=%d): %s\n",
    length(year_dummy_cols), .yrs[1L], paste(year_dummy_cols, collapse = ", ")
  ))
}

# NO continuous `time` column (removed 2026-08-06). It used to sit here alongside the year dummies,
# which made the design EXACTLY RANK-DEFICIENT: with 3 years the intercept + 2 dummies already
# saturate the time dimension, so time = 2000 + 10*year_2010 + 20*year_2020 EXACTLY (max|diff| = 0;
# X_mat rank 37 of 38, condition number 2e15, VIF = Inf for time/year_2010/year_2020). Worse, `time`
# was in re_idx, so that singularity was replicated across all 34 country RE blocks -- a perfectly
# flat likelihood direction. It is why the covariate effects never converged (theta_w median ESS 23
# at 12000x4, unchanged by more iterations, by re_slab_c2, or by re_center) while the PREDICTIONS
# were fine (predictions are invariant along a flat direction). The sampler's rank guard never fired.
# The year dummies retain FREE year effects (no functional form imposed). If a linear, extrapolable
# trend is wanted instead, drop the dummies and restore `time` -- but never both.
X_mat <- cbind(
  intercept = 1,
  as.matrix(dat_admin[, c(spatial_cont_cols_trans, soil_dev_cols, lu_area_cols, year_dummy_cols), with = FALSE])
)
# Guard: never ship a rank-deficient design again (the old failure was silent).
.rk <- qr(X_mat)$rank
if (.rk < ncol(X_mat)) warning(sprintf(
  "X_mat is RANK DEFICIENT: rank %d of %d columns -- flat likelihood direction(s), effects will not converge.", .rk, ncol(X_mat)))
X_mat[!is.finite(X_mat)] <- 0
rownames(X_mat) <- dat_admin$RESOLUTION

# --- Covariate SELECTION (env): pick exactly which covariates enter the model ---
# DRIVER_COVARIATES = comma-separated whitelist of FINAL X_mat column names (use ONLY these); else
# DRIVER_EXCLUDE_COVARIATES = drop these from the default set. Intercept + year dummies always kept.
# (Discover the available names with DRIVER_DIAG_EXIT=TRUE, which prints all X_mat columns.)
.always_keep <- c("intercept", year_dummy_cols)
.cov_wl <- trimws(strsplit(Sys.getenv("DRIVER_COVARIATES", ""), ",")[[1]])
.cov_wl <- .cov_wl[nzchar(.cov_wl)]
# 2026-08-04 spec change: RAI dropped from the data entirely; GHM_HI (+ GHM_TI once available) are
# now WANTED as drivers, so they leave this exclusion list. The rest of the earlier
# "accessibility/human-footprint + raw terrain mean" exclusion (CISI, Slope/Elevation/Aspect means,
# and the _sd heterogeneity twins) is UNCHANGED -- terrain still enters via the flat/steep/lowland/
# upland shares. Set DRIVER_EXCLUDE_COVARIATES explicitly to override.
.cov_ex <- trimws(strsplit(Sys.getenv("DRIVER_EXCLUDE_COVARIATES",
  "CISI,GHM_HI_sd,GHM_TI_sd,CISI_sd,Slope_rad,Elevation,Aspect_cos_mean,Aspect_sin_mean"), ",")[[1]])
.cov_ex <- .cov_ex[nzchar(.cov_ex)]
if (length(.cov_wl)) {
  .unknown <- setdiff(.cov_wl, colnames(X_mat))
  if (length(.unknown)) cat(sprintf(">>> DRIVER_COVARIATES: %d name(s) not in X_mat (ignored): %s\n", length(.unknown), paste(.unknown, collapse = ", ")))
  X_mat <- X_mat[, colnames(X_mat) %in% union(.always_keep, .cov_wl), drop = FALSE]
  cat(sprintf(">>> DRIVER_COVARIATES whitelist applied -> %d covariate columns: %s\n", ncol(X_mat), paste(colnames(X_mat), collapse = ", ")))
} else if (length(.cov_ex)) {
  .drop <- setdiff(intersect(.cov_ex, colnames(X_mat)), .always_keep)
  X_mat <- X_mat[, !(colnames(X_mat) %in% .drop), drop = FALSE]
  cat(sprintf(">>> DRIVER_EXCLUDE_COVARIATES: dropped %d -> %d columns remain\n", length(.drop), ncol(X_mat)))
}

# Column index bookkeeping for the year FE: keep them LINEAR but out of the horseshoe and the REs.
year_idx <- which(colnames(X_mat) %in% year_dummy_cols) # integer(0) for a single timestep
# --- RE STRUCTURE ---------------------------------------------------------------
# Subset of random slopes based on country-level heterogeneity (e.g. GDP, Population).
# DRIVER_RE_VARS overrides; "all" restores full slopes across all covariates.
.re_spec_count <- Sys.getenv("DRIVER_RE_VARS",
  paste0("log1p_GDP,log1p_Pop,GHM_HI,GHM_TI,log1p_allPA_area,flat_share,Slope_rad_sd,Growing_Degree_Days_gdd5,",
         # the largest grazing class, by name under whichever classification is in use
         "log1p_", names(which.max(vapply(pasture_cols, function(v) sum(dat_admin[[v]], na.rm = TRUE), 0)))))
if (identical(tolower(trimws(.re_spec_count)), "all")) {
  re_col_idx <- setdiff(seq_len(ncol(X_mat)), year_idx)      # legacy: every cov except year FE
} else {
  .re_want <- trimws(strsplit(.re_spec_count, ",")[[1]]); .re_want <- .re_want[nzchar(.re_want)]
  .re_miss <- setdiff(.re_want, colnames(X_mat))
  if (length(.re_miss)) cat(sprintf(">>> DRIVER_RE_VARS: %d name(s) not in X_mat (ignored): %s\n",
                                    length(.re_miss), paste(.re_miss, collapse = ", ")))
  re_col_idx <- sort(setdiff(intersect(c(which(colnames(X_mat) == "intercept"),
                                         which(colnames(X_mat) %in% .re_want)),
                                       seq_len(ncol(X_mat))), year_idx))
}
cat(sprintf(">>> RE covariates: %d of %d (%s)\n", length(re_col_idx), ncol(X_mat),
            paste(colnames(X_mat)[re_col_idx], collapse = ", ")))
hs_col_idx <- setdiff(2:ncol(X_mat), year_idx) # horseshoe on continuous cov (not intercept, not year FE)

# --- BART covariate partition (only used when use_bart) ---------------------------------------
# Topo/climate/socioecon drivers where a smooth nonlinearity is plausible (validated set from
# count-bart-solution). Everything else stays LINEAR. Overridable via DRIVER_BART_VARS (comma list).
BART_VARS <- strsplit(Sys.getenv("DRIVER_BART_VARS",
  "Slope_rad,Elevation,Growing_Degree_Days_gdd5,Annual_Precipitation_bio12,log1p_GDP,log1p_Pop"), ",")[[1]]
bart_col_idx <- if (use_bart) which(colnames(X_mat) %in% trimws(BART_VARS)) else integer(0)
linear_col_idx <- setdiff(seq_len(ncol(X_mat)), bart_col_idx)  # BART cols are NOT in the linear block
if (use_bart) cat(sprintf(">>> BART on %d cols: %s\n", length(bart_col_idx),
                          paste(colnames(X_mat)[bart_col_idx], collapse = ", ")))

offset_vec <- as.numeric(dat_admin$log_area_offset)

# Downscaling prediction grid: save the covariate space at EEA_10kmID x NUTS3 and stop (no outcome /
# no MCMC). X_mat here is built with the SAME covariate construction as the fitted model, so the
# fitted totals-beta and composition-delta apply directly; predict + renormalize to the NUTS3 totals.
if (nzchar(DOWNSCALE_GRID)) {
  .dg <- list(
    X_mat = X_mat, offset_vec = offset_vec, col_names = colnames(X_mat), years = MODEL_YEARS,
    unit = dat_admin$RESOLUTION, time = dat_admin$out_year, NUTS3 = dat_admin$NUTS3, NUTS0 = dat_admin$NUTS0,
    GLOB_country = dat_admin[[RE_GROUP_COL]]
  )
  .dgf <- file.path("output", sprintf("downscale_grid_X_%s_%s.rds", DOWNSCALE_GRID, paste(MODEL_YEARS, collapse = "_")))
  saveRDS(.dg, .dgf)
  cat(sprintf(
    "\n>>> DOWNSCALE_GRID: saved prediction-grid covariates -> %s (X %d units x %d cols)\n",
    .dgf, nrow(X_mat), ncol(X_mat)
  ))
  stop("DOWNSCALE_GRID: prediction grid built, stopping before MCMC.")
}

# Prepare Grouping for Random Effects
if (use_re) {
  # droplevels(): GLOB_country can arrive as a factor with UNUSED levels (e.g. a country whose rows
  # were filtered out). Without dropping them, re_group_names (levels) > the sampler's actual groups
  # AND group_idx has a gap -> misaligned/over-counted REs and the post-processing "1485 vs 1530"
  # mismatch. droplevels makes group_idx contiguous and re_group_names == #groups in postb_total.
  re_groups <- droplevels(as.factor(dat_admin[[RE_GROUP_COL]]))
  group_idx_vec <- as.integer(re_groups)
  re_group_names <- levels(re_groups)
  cat(sprintf(">>> RE Model: Detected %d groups using key: %s\n", length(re_group_names), RE_GROUP_COL))
} else {
  group_idx_vec <- NULL
}

# Early exit for covariate-build diagnostics (matches DRIVER_DIAG_EXIT in the pixel driver).
if (isTRUE(as.logical(Sys.getenv("DRIVER_DIAG_EXIT", "FALSE")))) {
  cat(sprintf("\n>>> DIAG: X_mat = %d rows x %d cols.\n", nrow(X_mat), ncol(X_mat)))
  cat(sprintf("    LU area cols (%d): %s\n", length(lu_area_cols), paste(head(lu_area_cols, 12), collapse = ", ")))
  cat(sprintf("    climate/yield cols (%d): %s\n", length(climate_yield_present), paste(climate_yield_present, collapse = ", ")))
  cat(sprintf(
    "    soil_dev cols: %d | spatial_cont_trans: %d | all finite: %s | any all-zero col: %s\n",
    length(soil_dev_cols), length(spatial_cont_cols_trans), all(is.finite(X_mat)),
    any(colSums(abs(X_mat)) == 0)
  ))
  cat(sprintf("    all X_mat cols: %s\n", paste(colnames(X_mat), collapse = ", ")))
  cat(sprintf(
    "    year FE cols (%d): %s | re_idx %d cols, hs_idx %d cols\n",
    length(year_dummy_cols), paste(year_dummy_cols, collapse = ", "), length(re_col_idx), length(hs_col_idx)
  ))
  stop("DRIVER_DIAG_EXIT: covariate build complete, stopping before MCMC.")
}

# Dump the assembled MODEL INPUTS (X_mat, per-target Y, offset, group index, names) so the fit +
# convergence diagnostics can be iterated standalone WITHOUT re-reading the 1.3GB prior_1km inputs.
if (isTRUE(as.logical(Sys.getenv("DRIVER_DUMP_INPUTS", "FALSE")))) {
  .dump <- list(
    X_mat = X_mat, offset_vec = offset_vec,
    group_idx_vec = group_idx_vec, re_group_names = if (use_re) re_group_names else NULL,
    Y = setNames(lapply(target_livestock, function(t) round(as.matrix(dat_admin[[t]]))), target_livestock),
    years = MODEL_YEARS, col_names = colnames(X_mat),
    year_dummy_cols = year_dummy_cols, year_idx = year_idx, re_col_idx = re_col_idx, hs_col_idx = hs_col_idx
  )
  .dp <- file.path("output", sprintf("count_model_inputs_%s.rds", paste(MODEL_YEARS, collapse = "_")))
  saveRDS(.dump, .dp)
  cat(sprintf(
    "\n>>> DRIVER_DUMP_INPUTS: saved model inputs to %s (X %dx%d, %d groups)\n",
    .dp, nrow(X_mat), ncol(X_mat), length(unique(group_idx_vec))
  ))

  # F (count): build + self-check the RECIPE so fit & predict share one assembly (predict_count_prior from
  # raw admin data). NO focal -> just transforms + [intercept, time, features, year-dummies] + log-area offset
  # + group keying. Self-diagnosing (tryCatch): reports max|dX|/|d offset| vs this X_mat; saves the recipe.
  suppressWarnings(tryCatch({
    source("codes/prior_model_predict.R")
    .recipe <- list(
      transforms = list(list(fn = "log1p", cols = c(skewed_vars, sub("^log1p_", "", lu_area_cols)), prefix = "log1p_")),
      time_col = "out_year", add_time = "time" %in% colnames(X_mat),
      year_map = setNames(as.list(as.integer(sub("^year_", "", year_dummy_cols))), year_dummy_cols),
      offset_col = "total_area_km2", offset_floor = 1e-4,
      feature_cols = setdiff(colnames(X_mat), c("intercept", "time")), col_order = colnames(X_mat),
      group_col = RE_GROUP_COL, group_levels_factor = if (use_re) re_group_names else NULL,
      group_levels_appear = if (use_re) unique(group_idx_vec) else NULL,
      train_group_idx = if (use_re) group_idx_vec else NULL)
    .sc <- count_recipe_selfcheck(.recipe, dat_admin, X_mat, offset_vec)
    cat(sprintf(">>> COUNT RECIPE self-check: ncol match=%s | max|dX|=%.2e | max|d offset|=%.2e | group match=%s\n",
                .sc$ncol_match, .sc$max_dX, .sc$max_doffset, .sc$group_match))
    saveRDS(.recipe, file.path("output", sprintf("count_model_recipe_%s.rds", paste(MODEL_YEARS, collapse = "_"))))
    cat(">>> saved count recipe (build_count_prior_model(fit, recipe, ...) after fitting ships a predictor)\n")
  }, error = function(e) cat(sprintf(">>> COUNT RECIPE build skipped (%s)\n", conditionMessage(e)))))
  stop("DRIVER_DUMP_INPUTS: model inputs dumped, stopping before MCMC.")
}

# =========================================================================
# 6. MODEL ESTIMATION
# =========================================================================

# We fit models separately for cattle and sheep_goat (target_livestock defined in the config).
beta_medians_list <- list()

for (tgt in target_livestock) {
  # Per-target saved-model dir under output/saved_model_outputs (matches the pixel driver's layout).
  # RUN_MODE namespace replaces the old date-stamped flat file: test clears + re-fits; production
  # persists -> a re-run resumes from the checkpoint.
  model_disk_path <- file.path(
    "output/saved_model_outputs",
    paste0(MODEL_LABEL, "_", tgt, if (use_re) "_RE_" else "_pooled_", RE_GROUP_COL)
  )
  dir.create(model_disk_path, recursive = TRUE, showWarnings = FALSE)
  if (RUN_MODE == "test") {
    .stale <- list.files(model_disk_path, full.names = TRUE)
    if (length(.stale) > 0) {
      cat(sprintf(">>> RUN_MODE=test: clearing %d stale file(s) in %s\n", length(.stale), model_disk_path))
      unlink(.stale, recursive = TRUE)
    }
  }
  checkpoint_path <- file.path(model_disk_path, "fit.rds")

  # ---------------------------------------------------------------------------------------------
  # CHECKPOINT FINGERPRINT (2026-08-21). The path is keyed only on MODEL_LABEL / target / RE flag /
  # RE_GROUP_COL, so a production re-run with DIFFERENT SAMPLER SETTINGS used to reload the old fit
  # unconditionally and report it as the new one -- the note below ("assumes the checkpoint was
  # produced with the current N_CHAINS / niter / nburn") documented the hazard without guarding it.
  # Same failure class as the nested_cut store hash, which silently returned a symmetric fit for a
  # diagonal arm. Anything that changes what the sampler DOES belongs in this fingerprint.
  .cfg_path <- file.path(model_disk_path, "fit_config.rds")
  .cfg_now  <- list(niter = niter, nburn = nburn, thin_keep = thin_keep, n_chains = N_CHAINS,
                    use_re = use_re, re_group = RE_GROUP_COL, use_bart = use_bart,
                    re_prec_pooled = COUNT_RE_PREC_POOLED,
                    re_slab_c2 = 16,
                    use_country_shrinkage = COUNT_USE_COUNTRY_SHRINKAGE,
                    slab_c2_country = COUNT_SLAB_C2_COUNTRY,
                    n = nrow(X_mat), p = ncol(X_mat), cols = colnames(X_mat))
  if (file.exists(checkpoint_path) && file.exists(.cfg_path)) {
    .cfg_old <- tryCatch(readRDS(.cfg_path), error = function(e) NULL)
    if (!identical(.cfg_old, .cfg_now)) {
      cat(">>> Checkpoint config MISMATCH -- settings changed since it was written. Discarding it\n")
      cat("    and re-fitting (delete-and-refit, not resume). Changed: ",
          paste(names(which(vapply(names(.cfg_now), function(z)
            !identical(.cfg_old[[z]], .cfg_now[[z]]), logical(1)))), collapse = ", "), "\n", sep = "")
      unlink(checkpoint_path); unlink(.cfg_path)
    }
  } else if (file.exists(checkpoint_path) && !file.exists(.cfg_path)) {
    cat(">>> Checkpoint has NO config fingerprint (written before 2026-08-21) -- cannot verify it\n")
    cat("    matches the current settings. Discarding and re-fitting.\n")
    unlink(checkpoint_path)
  }

  if (file.exists(checkpoint_path)) {
    cat(sprintf("\n>>> Checkpoint found for %s: %s. Loading existing estimation...\n", tgt, checkpoint_path))
    res_full <- readRDS(checkpoint_path)

    # Reconstruct res_list to enable convergence assessment.
    # NOTE: assumes the checkpoint was produced with the current N_CHAINS / niter / nburn
    # (env-driven above) -- keep DRIVER_NCHAINS/NITER/NBURN consistent across run + reload.
    n_ret <- (niter - nburn) %/% thin_keep   # thinned retained draws per chain

    res_list <- lapply(seq_len(N_CHAINS), function(ch) {
      chain_res <- list(
        var_names = dimnames(res_full$postb_pooled)[[1]],
        nuts0_names = res_full$nuts0_names
      )

      # Detect if chains are kept separate (4D for postb_pooled, 3D for post_sigma_re)
      is_4d_pooled <- (length(dim(res_full$postb_pooled)) == 4)

      if (!is.null(res_full$postb_pooled)) {
        if (is_4d_pooled) {
          tmp <- res_full$postb_pooled[, , , ch, drop = FALSE]
          dim(tmp) <- dim(tmp)[1:3]
          dimnames(tmp) <- dimnames(res_full$postb_pooled)[1:3]
          chain_res$postb_pooled <- tmp
        } else {
          idx <- ((ch - 1) * n_ret + 1):(ch * n_ret)
          chain_res$postb_pooled <- res_full$postb_pooled[, , idx, drop = FALSE]
        }
      }
      if (!is.null(res_full$post_sigma_re)) {
        if (is_4d_pooled) {
          tmp <- res_full$post_sigma_re[, , ch, drop = FALSE]
          dim(tmp) <- dim(tmp)[1:2]
          dimnames(tmp) <- dimnames(res_full$post_sigma_re)[1:2]
          chain_res$post_sigma_re <- tmp
        } else {
          idx <- ((ch - 1) * n_ret + 1):(ch * n_ret)
          if (length(dim(res_full$post_sigma_re)) == 3) {
            chain_res$post_sigma_re <- res_full$post_sigma_re[, , idx, drop = FALSE]
          } else {
            chain_res$post_sigma_re <- res_full$post_sigma_re[, idx, drop = FALSE]
          }
        }
      }
      if (!is.null(res_full$postb_total)) {
        if (is_4d_pooled) {
          tmp <- res_full$postb_total[, , , , ch, drop = FALSE]
          dim(tmp) <- dim(tmp)[1:4]
          dimnames(tmp) <- dimnames(res_full$postb_total)[1:4]
          chain_res$postb_total <- tmp
        } else {
          idx <- ((ch - 1) * n_ret + 1):(ch * n_ret)
          chain_res$postb_total <- res_full$postb_total[, , , idx, drop = FALSE]
        }
      }
      if (!is.null(res_full$post_log_lik)) {
        if (is_4d_pooled) {
          chain_res$post_log_lik <- res_full$post_log_lik[, ch]
        } else {
          idx <- ((ch - 1) * n_ret + 1):(ch * n_ret)
          chain_res$post_log_lik <- res_full$post_log_lik[idx]
        }
      }
      if (!is.null(res_full$post_log_lik_pointwise)) {
        if (is_4d_pooled) {
          tmp <- res_full$post_log_lik_pointwise[, , ch, drop = FALSE]
          dim(tmp) <- dim(tmp)[1:2]
          dimnames(tmp) <- dimnames(res_full$post_log_lik_pointwise)[1:2]
          chain_res$post_log_lik_pointwise <- tmp
        } else {
          idx <- ((ch - 1) * n_ret + 1):(ch * n_ret)
          chain_res$post_log_lik_pointwise <- res_full$post_log_lik_pointwise[idx, , drop = FALSE]
        }
      }
      if (!is.null(res_full$horseshoe)) {
        chain_res$horseshoe <- list()
        if (!is.null(res_full$horseshoe$post_kappa_mean)) {
          if (is_4d_pooled) {
            tmp <- res_full$horseshoe$post_kappa_mean[, , ch, drop = FALSE]
            dim(tmp) <- dim(tmp)[1:2]
            dimnames(tmp) <- dimnames(res_full$horseshoe$post_kappa_mean)[1:2]
            chain_res$horseshoe$post_kappa_pooled <- tmp
          } else {
            idx <- ((ch - 1) * n_ret + 1):(ch * n_ret)
            chain_res$horseshoe$post_kappa_pooled <- res_full$horseshoe$post_kappa_mean[, idx, drop = FALSE]
          }
        }
      }
      chain_res
    })
  } else {
    cat(sprintf("\n>>> Fitting Count Model for %s <<<\n", tgt))

    # LSU totals are continuous (fractional) and huge (up to ~1.3M). Convert to integer counts in
    # units of LSU_PER_COUNT (default 100): round(LSU / unit). This makes the NB/CRT integer
    # assumption meaningful AND brings the counts into the range where the dispersion r is
    # identified (raw-LSU magnitude leaves r unidentified -> r<->intercept ridge -> poor mixing).
    Y_target <- round(as.matrix(dat_admin[[tgt]]) / LSU_PER_COUNT)
    colnames(Y_target) <- tgt
    cat(sprintf(
      ">>> %s modeled in units of %g LSU/count: mean %.1f, max %.0f, zeros %.1f%%\n",
      tgt, LSU_PER_COUNT, mean(Y_target), max(Y_target), 100 * mean(Y_target == 0)
    ))

    # --- Empirical-Bayes dispersion prefit for BART (count-bart-solution) ---
    # NB-BART is metastable if r is estimated JOINTLY with BART (the r-collapse loop: a transient
    # BART mu-inflation collapses r -> -log(r) feedback inflates U -> bad mode). Fix r via a fast
    # NB-LINEAR prefit: mean_param removes the r<->intercept ridge, so r IS well identified there
    # (validated r=1.84 for BOV, over-pred 1.15x). Then run BART with r held fixed (r_warmup=1e9) =
    # the stable Poisson mechanism with the correct dispersion. Linear runs need no prefit.
    r_eb <- NULL
    if (use_bart) {
      cat(sprintf(">>> BART: empirical-Bayes r prefit (NB-linear) for %s ...\n", tgt))
      .pf <- mncount_rcpp(
        X = X_mat, Y = Y_target, family = "negbin", offset = offset_vec, intercept = FALSE,
        use_re = use_re, group_idx = group_idx_vec, re_idx = re_col_idx, use_bart = FALSE,
        linear_idx = 1:ncol(X_mat), r_method = "crt", mean_param = TRUE, re_regularize = TRUE,
        re_slab_c2 = 100, support_prior_strength = 1, niter = max(600, floor(niter / 2)),
        nburn = max(400, floor(nburn / 2)), standardize = TRUE, method = c("center", "scale"),
        use_horseshoe = TRUE, horseshoe_idx = hs_col_idx, equation_specific_hs = TRUE,
        tau0_mu = 0.5, tau0_dev = 0.2, chain_id = 1
      )
      r_eb <- as.numeric(.pf$r_disp_final)
      cat(sprintf(">>> EB r for %s = %.3f (held fixed in the BART chains)\n", tgt, r_eb))
    }

    res_list <- with_progress({
      p_bar <- progressor(steps = niter * N_CHAINS)

      future_lapply(seq_len(N_CHAINS), function(i) {
        source("codes/count_rcpp.R")

        mncount_rcpp(
          X = X_mat,
          Y = Y_target,
          family = "negbin",
          offset = offset_vec,
          intercept = FALSE,
          use_re = use_re,
          group_idx = group_idx_vec,
          re_idx = re_col_idx, # country REs on all covariates EXCEPT the year FE
          use_bart = use_bart,
          linear_idx = linear_col_idx,          # BART cols excluded from the linear block
          # --- count-BART stable config (count-bart-solution memory): mean_param removes the
          # r<->intercept ridge; bart_absorb_mean=FALSE avoids the exp-link level leak; fixed EB r
          # (r_init + r_warmup=1e9) avoids the r-collapse loop; slim trees stored for grid predict. ---
          mean_param = TRUE,
          bart_idx = if (use_bart) bart_col_idx else NULL,
          n_trees_bart = n_trees_bart,
          bart_base = 0.95, bart_power = 3.0, bart_k = 3.0, bart_f_cap = 3.0,
          bart_absorb_mean = FALSE,
          store_bart_trees = use_bart, do_slim_trees = use_bart,
          r_init = if (use_bart) r_eb else NULL,
          r_warmup = if (use_bart) 1e9 else NULL,
          # --- ported from the MNL work (audit-driven), all validated on simulated NB panels ---
          r_method = "crt", # exact CRT-Gibbs dispersion (replaces the slow/fragile log-RW Metropolis)
          re_prec_pooled = COUNT_RE_PREC_POOLED,   # one RE scale per predictor across outcome columns
          use_country_shrinkage = COUNT_USE_COUNTRY_SHRINKAGE, # Country-gatekeeper hierarchical shrinkage tau_m
          tau0_country = 1.0,
          slab_c2_country = COUNT_SLAB_C2_COUNTRY, # regularized Finnish slab cap on country gatekeeper: tau_m <= sqrt(4.0)=2.0
          re_regularize = TRUE, # regularised-horseshoe RE variance (slice): tames the heavy half-Cauchy tail
          re_slab_c2 = 16,  # effective SD bounded by sqrt(16)=4.0 (harmonized with composition model)
          support_prior_strength = 1, # participation-ratio RE shrinkage (gentler than 2: at strength 2 the sparsest of ~34 groups got ~70x shrinkage)
          init_jitter = 0.1, # per-chain overdispersed init -> honest Rhat across the chains
          niter = niter,
          nburn = nburn,
          thin = thin_keep,                    # RAM: store every k-th draw
          standardize = TRUE,
          method = c("center", "scale"),
          calc_loo = calc_loo_flag,            # RAM: LOO/WAIC off by default (heavy pointwise matrix + PSIS)
          use_horseshoe = TRUE,
          horseshoe_idx = hs_col_idx, # continuous covariates only (excludes intercept + year FE)
          equation_specific_hs = TRUE,
          estimate_c2 = TRUE, slab_df = 4, slab_s2 = 15,
          tau0_mu = 0.5,
          tau0_dev = 0.2,
          chain_id = i,
          progress_cb = function(...) p_bar()
        )
      }, future.seed = TRUE, future.stdout = NA)
    })

    cat("\n>>> Saving Checkpoint...\n")
    res_full <- combine_chains(res_list, keep_chains = TRUE)
    saveRDS(res_full, checkpoint_path)
    saveRDS(.cfg_now, .cfg_path)      # fingerprint alongside, so a settings change invalidates it
    Sys.sleep(5) # Let the filesystem flush before writing convergence diagnostics
  }

  # =========================================================================
  # 7. CONVERGENCE ASSESSMENT
  # =========================================================================
  cat("\n>>> Assessing Convergence...\n")
  assess_convergence(
    res_list = res_list,
    name = paste0(MODEL_LABEL, "_", tgt, "_", if (use_re) "RE" else "pooled"),
    categories = c(tgt, "baseline"),
    baseline_idx = 2,
    type = if (use_re) "re" else "pooled"
  )

  # =========================================================================
  # 8. PARAMETER ANALYTICS & SUMMARY
  # =========================================================================
  cat("\n=========================================================================\n")
  cat(sprintf("Analyzing Parameter Estimates for: %s\n", tgt))
  cat("=========================================================================\n")

  # --- Integrate out year fixed effects into the intercept ---
  if (length(year_dummy_cols) > 0) {
    cat(sprintf(">>> Integrating out %d year fixed effects into the intercept (mean over %d years)\n", 
                length(year_dummy_cols), length(MODEL_YEARS)))
    
    integrate_years <- function(arr) {
      if (is.null(arr)) return(arr)
      dn <- dimnames(arr)[[1]]
      if (is.null(dn)) return(arr)
      
      idx_int <- which(dn == "intercept")
      idx_years <- which(dn %in% year_dummy_cols)
      if (length(idx_int) != 1 || length(idx_years) == 0) return(arr)
      
      .dims <- dim(arr)
      n_dims <- length(.dims)
      
      # Safe slice extraction for the intercept
      args_int <- c(list(arr, idx_int), rep(list(quote(expr=)), n_dims - 1), list(drop = FALSE))
      int_slice <- do.call("[", args_int)
      
      # Average the time shocks into the intercept (time-integrated baseline)
      for (y_idx in idx_years) {
        args_y <- c(list(arr, y_idx), rep(list(quote(expr=)), n_dims - 1), list(drop = FALSE))
        y_slice <- do.call("[", args_y)
        int_slice <- int_slice + (y_slice / length(MODEL_YEARS))
      }
      
      # Assign integrated intercept back
      arr <- do.call("[<-", c(args_int[1:(length(args_int)-1)], list(value = int_slice)))
      
      # Drop year dummy rows safely across all dimensions
      args_drop <- c(list(arr, -idx_years), rep(list(quote(expr=)), n_dims - 1), list(drop = FALSE))
      arr <- do.call("[", args_drop)
      return(arr)
    }
    
    res_full$postb_pooled <- integrate_years(res_full$postb_pooled)
    res_full$postb_total <- integrate_years(res_full$postb_total)
    if (!is.null(res_full$horseshoe$post_kappa_mean)) {
      res_full$horseshoe$post_kappa_mean <- integrate_years(res_full$horseshoe$post_kappa_mean)
    }
  }

  # 1. Coefficient Draws
  # res_full$postb_pooled is either 4D [covariates, target, iterations, chains] or 3D [covariates, target, iterations]
  raw_pooled <- res_full$postb_pooled
  k_vars <- dim(raw_pooled)[1]

  if (length(dim(raw_pooled)) == 4) {
    draws_beta <- matrix(raw_pooled[, 1, , ], nrow = k_vars)
  } else {
    draws_beta <- matrix(raw_pooled[, 1, ], nrow = k_vars)
  }
  rownames(draws_beta) <- dimnames(raw_pooled)[[1]]

  # Posterior Statistics
  post_mean <- rowMeans(draws_beta)
  post_sd <- apply(draws_beta, 1, sd)
  ci_lower <- apply(draws_beta, 1, quantile, probs = 0.025)
  ci_upper <- apply(draws_beta, 1, quantile, probs = 0.975)

  # 2. Horseshoe Shrinkage (if present)
  hs_shrink <- rep(NA_real_, nrow(draws_beta))
  if (!is.null(res_full$horseshoe) && !is.null(res_full$horseshoe$post_kappa_mean)) {
    raw_kappa <- res_full$horseshoe$post_kappa_mean
    if (length(dim(raw_kappa)) == 3) {
      hs_shrink <- rowMeans(matrix(raw_kappa, nrow = k_vars))
    } else {
      hs_shrink <- rowMeans(raw_kappa)
    }
  }

  rhat_val <- rep(NA_real_, nrow(draws_beta))
  ess_val  <- rep(NA_real_, nrow(draws_beta))
  if (requireNamespace("posterior", quietly = TRUE) && length(dim(raw_pooled)) == 4) {
    for (i in seq_len(nrow(draws_beta))) {
      arr <- raw_pooled[i, 1, , , drop=TRUE] # [iterations, chains]
      if (is.matrix(arr)) {
        rhat_val[i] <- posterior::rhat(arr)
        ess_val[i]  <- posterior::ess_bulk(arr)
      }
    }
  }

  # 3. Create Summary Table
  sum_tab <- data.table(
    Covariate = rownames(draws_beta) %||% colnames(X_mat) %||% paste0("V", 1:nrow(draws_beta)),
    Post_Mean = post_mean,
    Post_SD = post_sd,
    CI_2.5 = ci_lower,
    CI_97.5 = ci_upper,
    HS_Shrinkage = hs_shrink,
    Rhat = rhat_val,
    ESS = ess_val
  )

  # Format numeric outputs for display
  print_tab <- copy(sum_tab)
  num_cols <- c("Post_Mean", "Post_SD", "CI_2.5", "CI_97.5", "HS_Shrinkage")
  print_tab[, (num_cols) := lapply(.SD, function(x) round(x, 4)), .SDcols = num_cols]

  print(print_tab)

  # Save to output CSV
  csv_file <- paste0("output/parameter_summary_", MODEL_LABEL, "_", tgt, "_", Sys.Date(), ".csv")
  write.csv(sum_tab, csv_file, row.names = FALSE)
  cat(sprintf("\n>>> Saved parameter summary to %s\n", csv_file))

  # 4. Dispersion Parameter (r) for Negative-Binomial
  if (res_full$family == "negbin" && !is.null(res_full$post_r)) {
    # res_full$post_r is either 3D [target, iterations, chains] or 2D [target, iterations]
    raw_r <- res_full$post_r
    if (length(dim(raw_r)) == 3) {
      draws_r <- as.vector(raw_r[1, , ])
    } else {
      draws_r <- as.vector(raw_r[1, ])
    }
    r_mean <- mean(draws_r)
    r_sd <- sd(draws_r)
    r_ci <- quantile(draws_r, probs = c(0.025, 0.975))

    cat("\n-------------------------------------------------------------------------\n")
    cat(sprintf("Dispersion Parameter (r) for %s:\n", tgt))
    cat(sprintf(
      "  Mean: %.4f | SD: %.4f | 95%% CI: [%.4f, %.4f]\n",
      r_mean, r_sd, r_ci[1], r_ci[2]
    ))
    cat("-------------------------------------------------------------------------\n")
  }

  # 5. WAIC and LOO-CV
  if (!is.null(res_full$waic) || !is.null(res_full$loo)) {
    cat("\n-------------------------------------------------------------------------\n")
    cat("Model Fit Diagnostics:\n")
    if (!is.null(res_full$waic)) {
      cat("\nWAIC Diagnostic:\n")
      print(res_full$waic)
    }
    if (!is.null(res_full$loo)) {
      cat("\nLOO-CV Diagnostic:\n")
      print(res_full$loo)
    }
    cat("-------------------------------------------------------------------------\n")
  }

  # 6. EXTRACT POSTERIOR BETA MEDIANS & RHAT (mirroring pixel-level model)
  # -------------------------------------------------------------------------
  cov_names_count <- dimnames(res_full$postb_total)[[1]]

  # combine_chains(keep_chains=TRUE) returns a 5D postb_total [cov,target,groups,draws,chains] for the
  # RE model. Collapse draws x chains -> 4D so the RE branch below handles it (the per-chain structure
  # for Rhat is recovered from postb_pooled's chain dimension). Fixes the "data length != k x G" error.
  if (length(dim(res_full$postb_total)) == 5) {
    .dm <- dim(res_full$postb_total)
    dim(res_full$postb_total) <- c(.dm[1], .dm[2], .dm[3], .dm[4] * .dm[5])
  }

  if (length(dim(res_full$postb_total)) == 4) {
    # RE model: postb_total is [covariates, target, groups, draws*chains]
    # GROUP LABELLING (fixed 2026-08-05). postb_total[, , k] is the k-th group in APPEARANCE order
    # -- the sampler builds its per-group arrays from `groups <- unique(group_idx)` (count_rcpp.R),
    # NOT from the sorted factor levels. group_idx_vec = as.integer(factor) is in SORTED-level ids,
    # so slice k belongs to re_group_names[unique(group_idx_vec)[k]], not re_group_names[k].
    # Using the sorted levels directly misattributed 28 of 34 countries (e.g. the slice labelled
    # "Belgium" was BosniaHerzg, "Croatia" was Switzerland) in the exported per-country betas and
    # the per-group convergence table. predict_count was never affected (it is passed
    # group_levels = unique(group_idx) explicitly; see recipe$group_levels_appear).
    group_names_vec <- if (exists("re_group_names") && !is.null(re_group_names)) {
      .appear <- if (exists("group_idx_vec") && !is.null(group_idx_vec)) unique(group_idx_vec) else seq_along(re_group_names)
      if (length(.appear) != length(re_group_names))
        stop(sprintf("group labelling: %d appearance-order groups vs %d level names -- cannot key postb_total safely.",
                     length(.appear), length(re_group_names)))
      re_group_names[.appear]
    } else if (!is.null(res_full$nuts0_names)) {
      res_full$nuts0_names
    } else {
      paste0("group_", seq_len(dim(res_full$postb_total)[3]))
    }

    d_ks <- dim(res_full$postb_total)[1]
    d_pop <- dim(res_full$postb_total)[2] # always 1 for count model
    d_Ns <- dim(res_full$postb_total)[3]
    d_draws <- dim(res_full$postb_total)[4]
    N_params <- d_ks * d_pop * d_Ns

    # Coordinate frame matching column-major layout [ks, lu, Ns]
    check_point_rhat_df_tgt <- data.table::as.data.table(expand.grid(
      ks = cov_names_count,
      LSTYP = tgt,
      Ns = group_names_vec,
      stringsAsFactors = FALSE
    ))
    setnames(check_point_rhat_df_tgt, "Ns", RE_GROUP_COL)

    flat_mat <- matrix(res_full$postb_total, nrow = N_params, ncol = d_draws)

    cat(sprintf("  Calculating posterior medians for %s (RE, %d groups)...\n", tgt, d_Ns))
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      check_point_rhat_df_tgt$value_median <- matrixStats::rowMedians(flat_mat)
    } else {
      check_point_rhat_df_tgt$value_median <- apply(flat_mat, 1, median)
    }

    # Rhat requires [draws, chains] — detect chain dimension from res_full
    # postb_pooled: [ks, lu, draws, chains] when keep_chains=TRUE
    if (length(dim(res_full$postb_pooled)) == 4) {
      d_chains <- dim(res_full$postb_pooled)[4]
      d_draws_per_chain <- d_draws / d_chains
    } else {
      d_chains <- 1L
      d_draws_per_chain <- d_draws
    }

    cat(sprintf("  Calculating Rhat for %s (RE, %d params)...\n", tgt, N_params))
    rhat_vals <- numeric(N_params)
    for (ii in seq_len(N_params)) {
      chain_mat <- matrix(flat_mat[ii, ], nrow = d_draws_per_chain, ncol = d_chains)
      rhat_vals[ii] <- tryCatch(posterior::rhat(chain_mat), error = function(e) NA_real_)
    }
    check_point_rhat_df_tgt$value_rhat <- rhat_vals

    rm(flat_mat)
    gc()
  } else {
    # Pooled model: postb_total is [covariates, target, draws, chains] (4D) or
    # [covariates, target, draws] (3D)
    d_ks <- dim(res_full$postb_total)[1]
    d_pop <- dim(res_full$postb_total)[2]
    d_draws <- if (length(dim(res_full$postb_total)) >= 3) dim(res_full$postb_total)[3] else 1L
    d_chains <- if (length(dim(res_full$postb_total)) == 4) dim(res_full$postb_total)[4] else 1L
    N_params <- d_ks * d_pop

    check_point_rhat_df_tgt <- data.table::as.data.table(expand.grid(
      ks = cov_names_count,
      LSTYP = tgt,
      stringsAsFactors = FALSE
    ))

    flat_mat <- matrix(res_full$postb_total, nrow = N_params, ncol = d_draws * d_chains)

    cat(sprintf("  Calculating posterior medians for %s (pooled)...\n", tgt))
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      check_point_rhat_df_tgt$value_median <- matrixStats::rowMedians(flat_mat)
    } else {
      check_point_rhat_df_tgt$value_median <- apply(flat_mat, 1, median)
    }

    cat(sprintf("  Calculating Rhat for %s (pooled, %d params)...\n", tgt, N_params))
    rhat_vals <- numeric(N_params)
    for (ii in seq_len(N_params)) {
      chain_mat <- matrix(flat_mat[ii, ], nrow = d_draws, ncol = d_chains)
      rhat_vals[ii] <- tryCatch(posterior::rhat(chain_mat), error = function(e) NA_real_)
    }
    check_point_rhat_df_tgt$value_rhat <- rhat_vals

    rm(flat_mat)
    gc()
  }

  # Save per-target convergence + median table
  conv_file <- paste0(
    "output/", SAMPLER, "_admin_convergence_check_", MODEL_LABEL, "_", tgt, "_", timestamp_str, ".rds"
  )  # mirrors pixel's <sampler>_pixel_convergence_check_<label>_<ts>.rds
  saveRDS(as.data.frame(check_point_rhat_df_tgt), conv_file)
  cat(sprintf("  >>> Saved convergence check to %s\n", conv_file))

  # Build the LUM-style beta_point_df for this target
  if (RE_GROUP_COL %in% colnames(check_point_rhat_df_tgt)) {
    beta_df <- as.data.frame(check_point_rhat_df_tgt)[
      ,
      c(RE_GROUP_COL, "LSTYP", "ks", "value_median")
    ]
  } else {
    beta_df <- as.data.frame(check_point_rhat_df_tgt)[
      ,
      c("LSTYP", "ks", "value_median")
    ]
  }
  beta_df$value <- beta_df$value_median
  beta_df$value_median <- NULL
  beta_df$LSTYP <- as.character(beta_df$LSTYP)
  beta_df$ks <- as.character(beta_df$ks)

  beta_medians_list[[tgt]] <- beta_df
}

# =========================================================================
# UNIFIED BETA MEDIAN EXPORT (mirroring pixel-level model)
# =========================================================================
cat("\n>>> Combining Posterior Beta Medians (LUM Case Style)...\n")
beta_means_df <- do.call(rbind, beta_medians_list)
rownames(beta_means_df) <- NULL

# Column ordering: [RE_GROUP_COL (if present),] LSTYP, ks, value
has_re_col <- RE_GROUP_COL %in% colnames(beta_means_df)
col_order <- if (has_re_col) c(RE_GROUP_COL, "LSTYP", "ks", "value") else c("LSTYP", "ks", "value")
beta_means_df <- beta_means_df[, col_order]

timestamp_str <- Sys.Date() # consistent with the first definition + the pixel driver (was format(Sys.time(),...))
out_file <- paste0(
  "output/", SAMPLER, "_admin_beta_median_",   # mirrors pixel's <sampler>_pixel_beta_median_<label>_<ts>.csv
  MODEL_LABEL,
  if (use_re) paste0("_RE_", RE_GROUP_COL) else "_pooled",
  "_", timestamp_str, ".csv"
)
write.csv(beta_means_df, out_file, row.names = FALSE)
cat(sprintf(">>> Saved unified posterior beta medians to %s\n", out_file))

# =========================================================================
# 9. SUBTYPE COMPOSITION (δ): drivers -> D/O/F split per species (the 6 subclasses)
# Run the composition pipeline as SEPARATE clean R processes: it uses mnlogit_rcpp_sym, a different
# C++ core than the count sampler loaded in THIS process, so process isolation avoids a symbol clash.
# Steps read the dat_admin_FULL this run just saved + the Eurostat extracts -> subclass parameters.
# =========================================================================
if (FIT_COMPOSITION) {
  cat("\n=========================================================================\n")
  cat(">>> Section 9: subtype composition (delta) for BOV/SGT x D/O/F x {conv,organic} ...\n")
  .rs <- file.path(R.home("bin"), "Rscript")
  .logf <- "output/composition_run.log"
  if (file.exists(.logf)) unlink(.logf)
  .steps <- list(
    c(s = "prep/prepare_composition_training.R", e = ""),               # D/O/F (NUTS2) + organic (NUTS3) training tables
    c(s = "composition/fit_composition.R", e = "D_SPECIES=bov"),
    c(s = "composition/fit_composition.R", e = "D_SPECIES=sgt"),
    c(s = "composition/fit_organic.R", e = ""),                                # organic-vs-conventional spatial pattern (NUTS3, off the Livestock_organic map)
    c(s = "composition/consolidate_composition.R", e = ""),
    c(s = "composition/build_subclass_parameters.R", e = "")
  ) # merge totals-beta + D/O/F delta + organic delta (Eurostat-anchored) -> gamma per country/subclass/system/driver
  .ok <- TRUE
  for (.st in .steps) {
    if (!file.exists(.st["s"])) {
      cat(sprintf("    SKIP (missing script): %s\n", .st["s"]))
      .ok <- FALSE
      break
    }
    cat(sprintf("    > %s %s\n", .st["s"], .st["e"]))
    .out <- tryCatch(system2(.rs, .st["s"],
      env = if (nzchar(.st["e"])) .st["e"] else character(0),
      stdout = TRUE, stderr = TRUE
    ), error = function(err) structure("", status = 1L))
    cat(.out, file = .logf, sep = "\n", append = TRUE)
    .rc <- attr(.out, "status")
    if (is.null(.rc)) .rc <- 0L
    if (!identical(as.integer(.rc), 0L)) {
      cat(sprintf("    !! failed (rc=%s) — see %s\n", .rc, .logf))
      .ok <- FALSE
      break
    }
  }
  if (.ok) {
    cat(">>> Composition complete -> output/composition/subclass_allocation_parameters.csv\n")
    if (file.exists("output/composition/subclass_allocation_parameters.csv")) {
      .p <- data.table::fread("output/composition/subclass_allocation_parameters.csv")
      cat(">>> Visible drivers per subclass:\n")
      print(.p[, .(n_drivers = .N, n_visible = sum(visible)), by = subclass])
    }
  } else {
    cat(">>> Composition step failed — the TOTALS above are complete and unaffected.\n")
    cat("    Debug by running the composition scripts standalone (see ", .logf, ").\n")
  }
}

# Promote this finished test run to the production bucket if requested (matches the pixel driver).
if (RUN_MODE == "test" && PROMOTE_TO_PRODUCTION) promote_test_to_production()

cat("\nEstimation and Assessment Complete.\n")

if (FIT_COMPOSITION) {
  cat("\nRunning Composition Sub-Models...\n")
  source("composition/fit_composition.R")
}

cat("\nBuilding Diagnostic HTML Report...\n")
source("postprocess/build_count_html.R")

