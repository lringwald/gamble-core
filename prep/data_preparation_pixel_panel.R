# -------------------------------------------------------------------------
# Panel Data Preparation: Pixel MNLogit (CLC 2000 + 2018)
# -------------------------------------------------------------------------
# Builds a two-period panel at pixel resolution for level-share MNLogit.
#
# Key design decisions:
#   - LU from CLC Annual TS serves BOTH as response Y AND as covariate X
#     (same-period level specification, not lagged)
#   - Time-varying covariates matched to closest forward vintage:
#       t=2000 → Pop_2000, GDP_2000, GHM_Ovr_2000, GHM_HI_2000
#       t=2018 → Pop_2020, GDP_2020, GHM_Ovr_2020, GHM_HI_2020
#   - Time-invariant covariates (topo, soil, PAs, CISI, RAI) shared
#   - Two cross-sections stacked into one panel
#
# Outputs (saved to temp/):
#   x_input_CLC_panel_{reso}_{y0}_{y1}_{date}.rds  — panel covariate matrix
#   lu_levels_CLC_panel_{reso}_{y0}_{y1}_{date}.rds — stacked LU shares (Y)
# -------------------------------------------------------------------------

rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(data.table)
library(arrow)
library(assertthat)
library(rlang)

source("codes/spatial_utils.R")


# =========================================================================
# 0. CONFIGURATION
# =========================================================================

PIXEL_RES     <- 5L
reso_km       <- PIXEL_RES
reso          <- paste0(PIXEL_RES, "km")
cell_area_km2 <- as.numeric(reso_km^2)
years         <- c(2000L, 2018L)        # CLC observation years
date_suffix   <- gsub("-", "", Sys.Date())

# Grid ID column names
cell_id_col  <- paste0("EEA_", reso, "ID")
x_col        <- paste0("x_", reso, "ID")
y_col        <- paste0("y_", reso, "ID")
cell_id_cols <- c(cell_id_col, x_col, y_col)

# Input paths
gridwork_dir <- "../LAMASUS_gridwork/output"

path_prior_inputs    <- file.path(gridwork_dir, "prior_model_1km_master_inputs.parquet")
path_grid_mapping    <- file.path(gridwork_dir, "one_kmID_master_mapping_2026-03-30.parquet")
# DATA ROOT. gamble-core holds the MODEL CODE; the bulk inputs live in cascadinggamble. GAMBLE_INPUT_DIR
# points at wherever they are, defaulting to the in-repo `input/` so an existing self-contained
# checkout keeps working unchanged.
INPUT_DIR            <- Sys.getenv("GAMBLE_INPUT_DIR", "input")
path_clc_class_map   <- file.path(INPUT_DIR, "CLC_Code1_LEVEL123_ETL2_mapping.csv")

# CLC Annual TS paths (one per year)
path_clc_ts <- setNames(
  file.path(gridwork_dir, paste0("CLC_Annual_TS_", years, ".rds")),
  as.character(years)
)

# Time-varying covariate vintage mapping (CLC year → closest forward vintage)
vintage_map <- list(
  "2000" = list(Pop = "Pop_2000", GDP = "GDP_2000",
                GHM_Ovr = "GHM_Ovr_2000", GHM_HI = "GHM_HI_2000"),
  "2018" = list(Pop = "Pop_2020", GDP = "GDP_2020",
                GHM_Ovr = "GHM_Ovr_2020", GHM_HI = "GHM_HI_2020")
)

# All 10 modelled LU classes
all_classes <- c("ACRP", "GRSL", "HCRP", "HEAS", "PAST",
                 "PCRP", "SPVA", "URBN", "WOFO", "WTLN")


# =========================================================================
# 1. LOAD GRID MAPPING & PRIOR INPUTS
# =========================================================================

message("Loading 1km → 5km grid mapping...")
grid_mapping <- read_parquet(path_grid_mapping) |>
  mutate(across(contains("km"), as.integer)) |>
  select(INSPIRE_Europe_buffer_1kmID, EEA_1kmID, all_of(cell_id_cols)) |>
  distinct(INSPIRE_Europe_buffer_1kmID, .keep_all = TRUE) |>
  # Enforce one consistent (x, y) centroid per 5km cell
  group_by(.data[[cell_id_col]]) |>
  mutate(
    !!x_col := first(.data[[x_col]]),
    !!y_col := first(.data[[y_col]])
  ) |>
  ungroup() |>
  distinct()

message("Loading 1km prior model inputs...")
prior_inputs_1km <- read_parquet(path_prior_inputs) |>
  mutate(
    EEA_1kmID = as.integer(EEA_1kmID),
    across(where(is.double), ~ ifelse(is.nan(.), NA_real_, .))
  )

# Attach pixel IDs to every 1km row
prior_inputs_joined <- prior_inputs_1km |>
  inner_join(grid_mapping, by = "EEA_1kmID")

message("  ", nrow(prior_inputs_joined), " 1km cells matched to ",
        n_distinct(prior_inputs_joined[[cell_id_col]]), " ", reso, " cells.")


# =========================================================================
# 2. AGGREGATE PRIOR COVARIATES TO PIXEL
# =========================================================================

message("Aggregating prior inputs to ", reso, "...")

# --- Helper functions for modal aggregation ---
dominant_int <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_integer_)
  as.integer(names(sort(table(x), decreasing = TRUE))[1L])
}

dominant_chr <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1L]
}

# --- 2a. Time-varying continuous covariates (keep ALL vintages; subset later) ---
prior_timevar_cols <- c(
  paste0("Pop_",     c(2020, 2015, 2010, 2005, 2000, 1995)),
  paste0("GHM_Ovr_", c(2020, 2015, 2010, 2005, 2000, 1995)),
  paste0("GHM_HI_",  c(2020, 2015, 2010, 2005, 2000, 1995)),
  paste0("GDP_",     c(2020, 2015, 2010, 2005, 2000, 1995))
)
# Keep only columns that actually exist in the parquet
prior_timevar_cols <- intersect(prior_timevar_cols, colnames(prior_inputs_joined))

prior_timevar_5km <- prior_inputs_joined |>
  group_by(across(all_of(cell_id_cols))) |>
  summarise(across(all_of(prior_timevar_cols), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

# --- 2b. Time-invariant continuous covariates ---
prior_invariant_cols <- c(
  "CISI", "RAI",
  "allPA_share", "CvPA_16_expansion_share", "CvPA_22_expansion_share",
  "CvPA_23_expansion_share", "N2K_share", "RstA_16_share", "RstA_18_share",
  "RstA_23_share", "StPA_share", "StPA_16_expansion_share",
  "StPA_22_expansion_share", "StPA_23_expansion_share"
)
prior_invariant_cols <- intersect(prior_invariant_cols, colnames(prior_inputs_joined))

prior_invariant_5km <- prior_inputs_joined |>
  group_by(across(all_of(cell_id_cols))) |>
  summarise(across(all_of(prior_invariant_cols), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

# --- 2c. Topography: mean + cross-cell variance ---
prior_topo_5km <- prior_inputs_joined |>
  mutate(Slope_deg = Slope_rad * (180 / pi)) |>
  group_by(across(all_of(cell_id_cols))) |>
  summarise(
    area_w_mean_slope_rad = mean(Slope_rad,       na.rm = TRUE),
    area_w_var_slope_rad  = var(Slope_rad,        na.rm = TRUE),
    area_w_mean_slope_deg = mean(Slope_deg,       na.rm = TRUE),
    area_w_var_slope_deg  = var(Slope_deg,        na.rm = TRUE),
    area_w_mean_alti      = mean(Elevation,       na.rm = TRUE),
    area_w_var_alti       = var(Elevation,        na.rm = TRUE),
    area_w_mean_hills     = mean(Hillshade,       na.rm = TRUE),
    area_w_var_hills      = var(Hillshade,        na.rm = TRUE),
    sin_mean              = mean(Aspect_sin_mean, na.rm = TRUE),
    cos_mean              = mean(Aspect_cos_mean, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    mean_aspect_deg = (atan2(sin_mean, cos_mean) * 180 / pi) %% 360,
    circ_variance   = 1 - sqrt(sin_mean^2 + cos_mean^2)
  )

# --- 2d. ESDAC soil classes: dominant (modal) value ---
prior_soil_cols <- c("WRBLV1", "WR", "WM2", "WM1", "VS", "ROO", "OC_TOP",
                     "IL", "ERODI", "DR", "DIMP", "CEC_TOP", "BS_TOP",
                     "AWC_TOP", "AWC_SUB")
prior_soil_cols <- intersect(prior_soil_cols, colnames(prior_inputs_joined))

prior_soil_5km <- prior_inputs_joined |>
  group_by(across(all_of(cell_id_cols))) |>
  summarise(across(all_of(prior_soil_cols), dominant_int), .groups = "drop")

# --- 2e. Administrative regions: dominant value ---
prior_region_cols <- c("NUTS0", "NUTS1", "NUTS2", "NUTS3", "EEA_biogeo")
prior_region_cols <- intersect(prior_region_cols, colnames(prior_inputs_joined))

prior_regions_5km <- prior_inputs_joined |>
  mutate(across(all_of(prior_region_cols), as.character)) |>
  group_by(.data[[cell_id_col]]) |>
  summarise(across(all_of(prior_region_cols), dominant_chr), .groups = "drop")

# --- Combine all aggregated prior covariates ---
prior_covariates_pixel <- prior_timevar_5km |>
  left_join(prior_invariant_5km, by = cell_id_cols) |>
  left_join(prior_topo_5km,      by = cell_id_cols) |>
  left_join(prior_soil_5km,      by = cell_id_cols) |>
  left_join(prior_regions_5km,   by = cell_id_col)

message("  Aggregated prior covariates: ", nrow(prior_covariates_pixel), " ", reso, " cells × ",
        ncol(prior_covariates_pixel), " columns.")

# Free large objects
rm(prior_inputs_1km, prior_inputs_joined,
   prior_timevar_5km, prior_invariant_5km, prior_topo_5km,
   prior_soil_5km, prior_regions_5km)
gc()


# =========================================================================
# 3. CLC ANNUAL TS → PIXEL LU LEVELS
# =========================================================================

message("Loading CLC class mapping...")
clc_class_map <- read.csv(path_clc_class_map)
clc_to_lu     <- clc_class_map |> select(Code1, ETL2_aug_abbr)

# Mapping from LAMASUS buffer ID → pixel cell
grid_map_clc <- grid_mapping |> select(INSPIRE_Europe_buffer_1kmID, all_of(cell_id_cols))

# Process each CLC year
lu_levels_list <- list()

for (yr in years) {
  yr_chr <- as.character(yr)
  message("Processing CLC Annual TS for ", yr, "...")
  
  clc_raw <- readRDS(path_clc_ts[yr_chr]) |>
    mutate(INSPIRE_Europe_buffer_1kmID = as.integer(INSPIRE_Europe_buffer_1kmID)) |>
    setDT()
  
  # Map Code1 → broad LU class
  clc_mapped <- clc_raw |>
    left_join(clc_to_lu, by = "Code1") |>
    rename(lu_class = ETL2_aug_abbr) |>
    filter(!is.na(lu_class), lu_class %in% all_classes) |>
    # Map to 5km cell
    inner_join(grid_map_clc, by = "INSPIRE_Europe_buffer_1kmID") |>
    # Aggregate to 5km
    group_by(across(all_of(c(cell_id_cols, "lu_class")))) |>
    summarise(area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop")
  
  # Pivot to wide: one row per cell, one column per LU class
  lu_wide <- clc_mapped |>
    pivot_wider(
      names_from  = lu_class,
      values_from = area_km2
    )
  
  # Ensure all expected classes are present (as NA if missing)
  for (cls in all_classes) {
    if (!cls %in% colnames(lu_wide)) lu_wide[[cls]] <- NA_real_
  }
  
  # Compute shares (normalise so rows sum to 1)
  lu_wide <- lu_wide |>
    mutate(
      total_area = rowSums(pick(all_of(all_classes))),
      across(all_of(all_classes), ~ .x / total_area)
    )
  
  lu_wide$times <- yr
  lu_levels_list[[yr_chr]] <- lu_wide
  
  message("  ", yr, ": ", nrow(lu_wide), " ", reso, " cells with LU data.")
  rm(clc_raw, clc_mapped)
  gc()
}


# =========================================================================
# 4. LU COMPLEXITY & FOCAL STATISTICS (per period)
# =========================================================================

message("Computing LU complexity and focal statistics...")

process_lu_features <- function(lu_wide, yr) {
  
  # --- 4a. Complexity metrics ---
  lu_shares_mat <- lu_wide |> select(all_of(all_classes)) |> as.matrix()
  
  complexity <- lu_wide |>
    select(all_of(cell_id_cols)) |>
    mutate(
      richness  = rowSums(lu_shares_mat > 0),
      shannon   = -rowSums(ifelse(lu_shares_mat > 0, lu_shares_mat * log(lu_shares_mat), 0)),
      simpson   = 1 - rowSums(lu_shares_mat^2),
      evenness  = ifelse(richness > 1, shannon / log(richness), 0),
      normalized_complexity = ifelse(richness > 1, shannon / log(richness), 0)
    )
  
  # --- 4b. Focal/neighborhood statistics ---
  lu_focal <- focal_df(
    df             = lu_wide,
    x_id_col       = x_col,
    y_id_col       = y_col,
    value_cols     = all_classes,
    w              = 3,
    fun            = "mean",
    include_center = FALSE,
    decay_type     = "inverse",
    decay_param    = 1
  )
  
  # Combine
  lu_features <- lu_wide |>
    left_join(complexity, by = cell_id_cols) |>
    left_join(lu_focal,   by = c(x_col, y_col))
  
  return(lu_features)
}

lu_features_list <- list()
for (yr in years) {
  yr_chr <- as.character(yr)
  lu_features_list[[yr_chr]] <- process_lu_features(lu_levels_list[[yr_chr]], yr)
  message("  ", yr, ": LU features computed (",
          ncol(lu_features_list[[yr_chr]]), " columns).")
}

rm(lu_levels_list)
gc()


# =========================================================================
# 5. CONSOLIDATE PANEL
# =========================================================================

message("Consolidating panel dataset...")

# Time-invariant covariate columns from prior_covariates_pixel
topo_cols <- c("area_w_mean_slope_rad", "area_w_var_slope_rad",
               "area_w_mean_slope_deg", "area_w_var_slope_deg",
               "area_w_mean_alti", "area_w_var_alti",
               "area_w_mean_hills", "area_w_var_hills",
               "sin_mean", "cos_mean", "mean_aspect_deg", "circ_variance")

# Build each cross-section with time-matched vintages
panel_list <- list()

for (yr in years) {
  yr_chr <- as.character(yr)
  vmap   <- vintage_map[[yr_chr]]
  
  # Select time-matched vintage columns and rename to generic names
  prior_timematched <- prior_covariates_pixel |>
    select(all_of(cell_id_cols),
           # Time-varying (rename to generic)
           Pop     = all_of(vmap$Pop),
           GDP     = all_of(vmap$GDP),
           GHM_Ovr = all_of(vmap$GHM_Ovr),
           GHM_HI  = all_of(vmap$GHM_HI),
           # Time-invariant continuous
           all_of(prior_invariant_cols),
           # Topography
           all_of(topo_cols)
    )
  
  # Combine LU features + prior covariates + soil + admin
  panel_yr <- lu_features_list[[yr_chr]] |>
    left_join(prior_timematched, by = cell_id_cols) |>
    left_join(
      prior_covariates_pixel |> select(all_of(cell_id_col), all_of(prior_soil_cols)),
      by = cell_id_col
    ) |>
    left_join(
      prior_covariates_pixel |> select(all_of(cell_id_col), all_of(prior_region_cols)),
      by = cell_id_col
    )
  
  panel_list[[yr_chr]] <- panel_yr
  message("  t=", yr, ": ", nrow(panel_yr), " rows × ", ncol(panel_yr), " cols")
}

# Stack the two cross-sections
panel_full <- bind_rows(panel_list)

# Replace NA with 0 for numeric covariates
panel_full <- panel_full |>
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

message("  Panel: ", nrow(panel_full), " rows × ", ncol(panel_full), " cols (",
        n_distinct(panel_full[[cell_id_col]]), " unique cells × ",
        length(years), " periods)")


# =========================================================================
# 5b. VALIDATION
# =========================================================================

message("Running validation checks...")

# No duplicate ns × times
assert_that(
  nrow(panel_full) == nrow(distinct(panel_full, .data[[cell_id_col]], times)),
  msg = "Duplicate cell × period rows detected!"
)

# LU shares sum to ~1
share_sums <- rowSums(panel_full[, all_classes])
assert_that(

  all(abs(share_sums - 1) < 0.05, na.rm = TRUE),
  msg = paste0("LU shares do not sum to ~1. Range: [",
               round(min(share_sums, na.rm = TRUE), 4), ", ",
               round(max(share_sums, na.rm = TRUE), 4), "]")
)

# Consistent column count across periods
n_per_year <- panel_full |>
  group_by(times) |>
  summarise(n_rows = n(), .groups = "drop")
message("  Rows per period:")
print(n_per_year)

# Covariate summary
cov_cols <- setdiff(
  colnames(panel_full),
  c(cell_id_cols, "times", "total_area", all_classes, prior_region_cols)
)
message("  Covariates (", length(cov_cols), "): ", paste(cov_cols, collapse = ", "))


# =========================================================================
# 6. EXPORT
# =========================================================================

message("Exporting...")
dir.create("temp", showWarnings = FALSE)

file_tag <- paste0("CLC_panel_", reso, "_", paste0(years, collapse = "_"), "_", date_suffix)

# --- X input (full panel with covariates) ---
# Rename cell_id_col to 'ns' for compatibility with estimation scripts
x_out <- panel_full |>
  rename(ns = all_of(cell_id_col)) |>
  select(
    where(is.character), times, ns,
    starts_with("x_"), starts_with("y_"),
    where(is.numeric)
  )

saveRDS(x_out, file = paste0("temp/x_input_", file_tag, ".rds"))

# --- Y response (LU shares only, for convenience) ---
y_out <- panel_full |>
  rename(ns = all_of(cell_id_col)) |>
  select(ns, times, all_of(all_classes))

saveRDS(y_out, file = paste0("temp/lu_levels_", file_tag, ".rds"))

message("\nDone. Outputs saved to temp/:")
message("  x_input_",    file_tag, ".rds  (", nrow(x_out), " × ", ncol(x_out), ")")
message("  lu_levels_",  file_tag, ".rds  (", nrow(y_out), " × ", ncol(y_out), ")")
