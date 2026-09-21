rm(list=ls())
require(tidyverse)
require(data.table)
require(arrow)

# Helper function to find the most recently modified file matching a pattern
get_latest_file <- function(dir_path, pattern) {
  files <- list.files(path = dir_path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  info <- file.info(files)
  latest_file <- rownames(info[order(info$mtime, decreasing = TRUE), ])[1]
  return(basename(latest_file))
}

DS_DIR <- "../LAMASUS_downscaling/"
GRIDWORK_DIR <- "../LAMASUS_gridwork/output"
PIXEL_RES <- 10

# Dynamically pick the latest mapping file
mapping_file <- get_latest_file(GRIDWORK_DIR, "^one_kmID_master_mapping_.*\\.parquet$")
if (is.null(mapping_file)) stop("No mapping file found in ", GRIDWORK_DIR)
cat(sprintf("Using latest mapping file: %s\n", mapping_file))
message(">>> Loading thematic mapping and base year maps...")
mapping_thematic <- read.csv("../LAMASUS_downscaling/aux_files/LAMASUS_LUM_thematic_mapping.csv") %>% setDT()


mapping_grid <- arrow::read_parquet(file.path(GRIDWORK_DIR, mapping_file)) %>% as.data.table()


temp_map_1km <- readRDS(paste0(DS_DIR,"input/LUM_fit_with_energy_levels_and_new_FM_2018_EEA_1kmID.rds")) %>%
  mutate(EEA_1kmID = as.integer(EEA_1kmID)) %>%
  setDT()

# Filter to current mapping grid scope
temp_map_1km <- temp_map_1km[EEA_1kmID %in% unique(mapping_grid$EEA_1kmID)]

# Store total expected area for validation
map_tot_area_check   <- sum(temp_map_1km$area_km2, na.rm = TRUE)
map_tot_pixels_check <- map_tot_area_check 

temp_binary_irrigation_map <- readRDS(paste0(DS_DIR,"/input/irrigation_binary_EEA_1kmID.rds")) %>%
  mutate(EEA_1kmID = as.integer(EEA_1kmID)) %>%
  filter(irrigation_binary == 1) %>%
  setDT()

temp_binary_irrigation_map <- temp_binary_irrigation_map[EEA_1kmID %in% unique(mapping_grid$EEA_1kmID)]


# --- Load and Pre-process Calibrated Organic Areas ---
message(">>> Loading calibrated organic certificate areas...")
temp_organic_map <- readRDS(paste0(GRIDWORK_DIR,"/organic_certificaties_final_calibrated_master.rds"))
# Columns: INSPIRE_Europe_buffer_1kmID, Cropland_organic, Livestock_organic, Mixed_organic, All_organic
# These are in km2 per 1km grid cell.

# Join with mapping_grid to get EEA_1kmID
temp_organic_map <- temp_organic_map[unique(mapping_grid[,.(EEA_1kmID,INSPIRE_Europe_buffer_1kmID)]), on="INSPIRE_Europe_buffer_1kmID"]
org_cols <- c("Cropland_organic", "Livestock_organic", "Mixed_organic", "All_organic")
temp_organic_map <- temp_organic_map[, .SD, .SDcols = c("EEA_1kmID", org_cols)]

# Fill NAs with 0
for (j in org_cols) set(temp_organic_map, which(is.na(temp_organic_map[[j]])), j, 0)

# ---------------------------------------------------------------
# 2. Refine Mapping Grid
# ---------------------------------------------------------------
message(">>> Refining mapping grid and calculating spatial weights...")
message(">>> Collapsing mapping grid to 1km resolution...")
mapping_grid <- unique(mapping_grid[, .(
  EEA_1kmID, 
  ns = get(paste0("EEA_", PIXEL_RES, "kmID")), 
  NUTS2
)])

# 2. Check for duplicates: if this > 0, your 1km cells span multiple NUTS2/ns areas
dups <- mapping_grid[, .N, by = EEA_1kmID][N > 1]
if(nrow(dups) > 0) {
  message("!!! Warning: Some 1km cells mapped to multiple regions. Taking the first match.")
  mapping_grid <- mapping_grid[, .SD[1], by = EEA_1kmID]
}

setnames(mapping_thematic, 
         old = c("GLOBIOM_UNFCCC", "GLOBIOM_mngmt"), 
         new = c("GLOB_lu_class", "GLOB_CSYS"))

# ---------------------------------------------------------------
# 3. Pivot & Irrigation Integration
# ---------------------------------------------------------------
message(">>> Integrating irrigation data...")

# Pivot wider to handle thematic categories as columns
temp_wide <- dcast(
  temp_map_1km,
  EEA_1kmID ~ LUM_fit_with_energy_levels_and_new_FM_2018,
  value.var = "area_km2",
  fill = 0
)

# Add LUM_ prefix to category columns
setnames(temp_wide, old = names(temp_wide)[-1], new = paste0("LUM", names(temp_wide)[-1]))

# Join irrigation status
temp_wide <- merge(temp_wide, temp_binary_irrigation_map, by = "EEA_1kmID", all.x = TRUE)
temp_wide[is.na(irrigation_binary), irrigation_binary := 0]
temp_wide[is.na(area_km2),          area_km2 := 0]

# Identify cropland and non-irrigated categories
temp_cropland_codes      <- mapping_thematic[GLOB_lu_class == "Cropland" & !grepl("O",GLOB_CSYS), LUM_Code]
temp_non_irrigated_codes <- mapping_thematic[GLOB_lu_class == "Cropland" & GLOB_CSYS != "IR" & !grepl("O",GLOB_CSYS), LUM_Code]

temp_cropland_cols      <- paste0("LUM", temp_cropland_codes)
temp_non_irrigated_cols <- paste0("LUM", temp_non_irrigated_codes)

# Calculate irrigated cropland area (LUM2006)
temp_wide[, LUM2006 := irrigation_binary * area_km2 * rowSums(.SD, na.rm = TRUE), .SDcols = temp_non_irrigated_cols]
temp_wide[, area_km2 := NULL]

# Subtract irrigated area from non-irrigated cropland pools
temp_wide[, (temp_non_irrigated_cols) := {
  dt_sub <- .SD
  row_sums <- rowSums(dt_sub, na.rm = TRUE)
  weights <- dt_sub / row_sums
  weights[is.na(weights)] <- 0 
  dt_sub - (LUM2006 * weights)
}, .SDcols = temp_non_irrigated_cols]

temp_wide[, irrigation_binary := NULL]

final_dt <- merge(temp_wide, temp_organic_map, by = "EEA_1kmID", all.x = TRUE)

# 3. Fill NA indicators with 0 (cells with no organic data remain 100% conventional)
org_cols <- c("Cropland_organic", "Livestock_organic", "Mixed_organic", "All_organic")
for (j in org_cols) set(final_dt, which(is.na(final_dt[[j]])), j, 0)
# --- Calculate Pixel-Level Constraints & Shares ---
message(">>> Calculating pixel-level organic shares constrained by LUM area...")
arable_cols  <- c("LUM2001", "LUM2002", "LUM2003","LUM2004","LUM2005","LUM2006") 
pasture_cols <- c(paste0("LUM", 14:26), c("LUM4012", "LUM4013"))

# Calculate total available area per pool
final_dt[, total_arable_lum  := rowSums(.SD, na.rm=TRUE), .SDcols = arable_cols]
final_dt[, total_pasture_lum := rowSums(.SD, na.rm=TRUE), .SDcols = pasture_cols]
final_dt[, total_agri_lum    := total_arable_lum + total_pasture_lum]

# Convert absolute km2 to shares [0, 1] relative to LUM availability
final_dt[, Cropland_organic  := ifelse(total_arable_lum > 0,  pmin(Cropland_organic / total_arable_lum, 1.0), 0)]
final_dt[, Livestock_organic := ifelse(total_pasture_lum > 0, pmin(Livestock_organic / total_pasture_lum, 1.0), 0)]
final_dt[, Mixed_organic     := ifelse(total_agri_lum > 0,    pmin(Mixed_organic / total_agri_lum, 1.0), 0)]
final_dt[, All_organic       := ifelse(total_agri_lum > 0,    pmin(All_organic / total_agri_lum, 1.0), 0)]

all_lum_cols <- names(final_dt)[grepl("^LUM", names(final_dt))]


# 2. Sequential allocation loop
for (col in all_lum_cols) {
  
  # We use temporary columns to track the math inside this specific loop
  final_dt[, primary_org_area := 0]
  
  # STEP 1: Cropland/Livestock First
  if (col %in% arable_cols) {
    final_dt[, primary_org_area := get(col) * Cropland_organic]
    final_dt[, (col) := get(col) - primary_org_area]
  } else if (col %in% pasture_cols) {
    final_dt[, primary_org_area := get(col) * Livestock_organic]
    final_dt[, (col) := get(col) - primary_org_area]
  }
  
  # STEP 2: Mixed on the remaining area
  final_dt[, mixed_org_area := get(col) * Mixed_organic]
  final_dt[, (col) := get(col) - mixed_org_area] # LUM column is now ONLY Conventional
  # STEP 3: All on the remaining area
  
  final_dt[, all_org_area := get(col) * All_organic]
  final_dt[, (col) := get(col) - all_org_area] # LUM column is now ONLY Conventional
  
  
  # STEP 3: Combine them into one "Organic" column
  final_dt[, paste0(col, "999") := primary_org_area + mixed_org_area + all_org_area]
  
  # Clean up temp columns for next loop iteration
  final_dt[, c("primary_org_area", "mixed_org_area", "all_org_area") := NULL]
}

# Remove agricultural pool columns
final_dt[, c("total_arable_lum", "total_pasture_lum", "total_agri_lum") := NULL]


# 3. Final safety: remove tiny floating point negatives
for (j in names(final_dt)[grepl("^LUM", names(final_dt))]) {
  set(final_dt, which(final_dt[[j]] < 0), j, 0)
}

message(">>> Aggregating to {PIXEL_RES} resolution...")
final_dt <- merge(
  final_dt, 
  mapping_grid, 
  by = "EEA_1kmID"
)

# 2. Aggregate all LUM columns at once while wide
# lapply(.SD, sum) is vectorized and very fast in data.table
message(">>> Summarizing wide table...")
agg_wide <- final_dt[, lapply(.SD, sum, na.rm = TRUE), 
                     by = .(ns, NUTS2), 
                     .SDcols = patterns("^LUM")]

# 3. NOW melt the summarized table
# The number of rows here will be (Unique ns/NUTS2 combinations * Number of LUM codes)
# This is likely 10-50x smaller than 531 million.
temp_long <- melt(
  agg_wide,
  id.vars = c("ns", "NUTS2"),
  variable.name = "LUM_Code",
  value.name = "LUM_area_km2"
)

# 4. Clean up the LUM_Code string to integer
temp_long[, LUM_Code := as.integer(sub("LUM", "", LUM_Code))]
mapping_thematic[,full_lum_class:=ifelse(GLOB_CSYS!="",paste0(GLOB_lu_class,"_",GLOB_CSYS),GLOB_lu_class)]
mapping_thematic <- unique(mapping_thematic[,.(LUM_Code,full_lum_class)])


temp_long <- temp_long[mapping_thematic, on = "LUM_Code"]

temp_long <- temp_long[,.(area_km2=sum(LUM_area_km2)), by=.(ns,NUTS2,full_lum_class)]

# --- Final Wide Transformation ---
message(">>> Creating final wide table for prior modeling...")
lum_wide_10km <- dcast(
  temp_long,
  ns + NUTS2 ~ full_lum_class,
  value.var = "area_km2",
  fill = 0
)
# Standardize column names (no spaces)
setnames(lum_wide_10km, old = names(lum_wide_10km), new = make.names(names(lum_wide_10km)))

# Automatically drop zero-sum columns
cat_cols <- setdiff(names(lum_wide_10km), c("ns", "NUTS2"))
zero_cols <- cat_cols[lum_wide_10km[, lapply(.SD, sum), .SDcols = cat_cols] == 0]

if (length(zero_cols) > 0) {
  message(sprintf(">>> Dropping %d zero-sum categories: %s", length(zero_cols), paste(zero_cols, collapse = ", ")))
  lum_wide_10km[, (zero_cols) := NULL]
}

# Save output for model runner
output_file <- file.path(GRIDWORK_DIR, sprintf("LUM_complete_%dkm_wide_2018.parquet", PIXEL_RES))
arrow::write_parquet(lum_wide_10km[, NODATA := NULL], output_file)
message(sprintf(">>> Final wide table saved to %s", output_file))

# Summary check
print(temp_long[,.(area_km2=sum(area_km2)), by=.(full_lum_class)])
