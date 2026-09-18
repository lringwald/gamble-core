#!/usr/bin/env Rscript
# =============================================================================
# Area audit: LUM -> HRL crop split -> Eurostat cascade, checked PER 1km CELL
# =============================================================================
# The pipeline already asserts area conservation, but only on the GLOBAL total -- a fault that moved
# area between cells would pass that untouched. This re-runs the real steps on the real data and
# checks conservation cell by cell, then quantifies the two places where area is legitimately NOT
# preserved, so they are known quantities rather than surprises:
#   * HRL clipping   where HRL crop area exceeds LUM cropland in a cell, crops are scaled to fit
#   * grid/pixel     1km cells outside the model grid are dropped (not this script's concern)
# =============================================================================
suppressMessages({library(data.table)})
GW <- Sys.getenv("GAMBLE_GRIDWORK_DIR", "../LAMASUS_gridwork/output")
YEAR <- Sys.getenv("AUDIT_YEAR", "2018")
source("codes/target_class_split.R")
ok <- function(cond, msg) cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg))

# ---- 1. LUM ------------------------------------------------------------------------------------
lum <- as.data.table(readRDS(file.path(GW, sprintf("LUM_fit_with_energy_levels_and_new_FM_%s.rds", YEAR))))
setnames(lum, grep("bufferID|1kmID", names(lum), value = TRUE)[1], "join_id")
setnames(lum, grep("code|Code", names(lum), value = TRUE)[1], "LUM_Code")
setnames(lum, grep("area", names(lum), value = TRUE)[1], "area")
lum <- lum[area > 0, .(join_id = as.integer(join_id), LUM_Code = as.integer(LUM_Code), area)]
A0 <- lum[, .(a0 = sum(area)), by = join_id]
cat(sprintf("\nLUM %s: %s cells | %.1f km2 total | per-cell mean %.3f\n",
            YEAR, format(nrow(A0), big.mark=","), sum(A0$a0), mean(A0$a0)))

# ---- 2. class mapping ----------------------------------------------------------------------------
map <- fread("aux_files/LUM_Code_to_macro_model_mapping.csv")
CLS <- Sys.getenv("AUDIT_CLASS_COL", "BMLEH_Los1_label")
m <- unique(map[, .(LUM_Code = as.integer(LUM_Code), cls = trimws(as.character(get(CLS))))])
m <- m[!is.na(LUM_Code)]
d <- merge(lum, m, by = "LUM_Code", all.x = TRUE)
unm <- d[is.na(cls) | !nzchar(cls), sum(area)]
cat(sprintf("\nclass mapping on %s: %.1f km2 unmapped (%.4f%%)\n", CLS, unm, 100*unm/sum(lum$area)))
d[is.na(cls) | !nzchar(cls), cls := "NODATA"]
d <- d[, .(area = sum(area)), by = .(join_id, model_class = cls)][, focal_class := model_class]
A1 <- d[, .(a1 = sum(area)), by = join_id]
k <- merge(A0, A1, by = "join_id", all = TRUE); k[is.na(k)] <- 0
ok(max(abs(k$a0 - k$a1)) < 1e-8, sprintf("mapping conserves area PER CELL (max dev %.2g)", max(abs(k$a0-k$a1))))

# ---- 3. HRL crop split ---------------------------------------------------------------------------
LEG <- data.table(
  Crop_Type_Code = c(1110,1120,1130,1140,1150,1210,1220,1310,1320,1410,1420,1430,1440,2100,2200,2310,2320),
  crop_class = c("Wheat","Barley","Maize","Rice","Other_cereals","Fresh_vegetables","Dry_pulses",
                 "Potatoes","Sugar_beet","Sunflower","Soybeans","Rapeseed","Flax_cotton_hemp",
                 "Grapes","Olives","Fruits","Nuts"),
  crop_group = c(rep("arable",13), rep("permanent",4)))
cr <- as.data.table(readRDS(file.path(GW, sprintf("Crop_Types/Crop_Types_%s_1km.rds", YEAR))))
setnames(cr, grep("bufferID", names(cr), value=TRUE)[1], "join_id")
setnames(cr, grep("area", names(cr), value=TRUE)[1], "carea")
cr[, join_id := as.integer(join_id)]
cr <- merge(cr, LEG, by = "Crop_Type_Code")[carea > 0]
clip_tot <- 0; new_rows <- list()
for (grp in c("arable","permanent")) {
  rc <- if (grp=="arable") "Cropland_arable_other" else "Cropland_permanent_other"
  A <- d[model_class == rc, .(A = sum(area)), by = join_id]; if (!nrow(A)) next
  Hc <- cr[crop_group == grp, .(H_c = sum(carea)), by = .(join_id, crop_class)]
  Ht <- Hc[, .(H = sum(H_c)), by = join_id]
  j <- merge(A, Ht, by = "join_id", all.x = TRUE); j[is.na(H), H := 0]
  j[, `:=`(scale = fifelse(H > 0, pmin(1, A/H), 0), resid = pmax(A - H, 0))]
  clip_tot <- clip_tot + j[H > A, sum(H - A)]
  new_rows[[grp]] <- rbind(
    merge(Hc, j[, .(join_id, scale)], by="join_id")[scale > 0, .(join_id, model_class=crop_class, area=H_c*scale)],
    j[resid > 1e-9, .(join_id, model_class = rc, area = resid)])
}
d2 <- rbind(d[!model_class %in% c("Cropland_arable_other","Cropland_permanent_other"), .(join_id, model_class, area)],
            rbindlist(new_rows), fill = TRUE)
d2[, focal_class := model_class]
A2 <- d2[, .(a2 = sum(area)), by = join_id]
k <- merge(A1, A2, by = "join_id", all = TRUE); k[is.na(k)] <- 0
cat(sprintf("\nHRL crop split:\n"))
ok(max(abs(k$a1 - k$a2)) < 1e-6, sprintf("conserves area PER CELL (max dev %.2g)", max(abs(k$a1-k$a2))))
cat(sprintf("  HRL area CLIPPED because it exceeded LUM cropland in the cell: %.0f km2\n", clip_tot))
cat(sprintf("  (that area is HRL's, never entered LUM, and is correctly not added)\n"))

# ---- 4. Eurostat cascade -------------------------------------------------------------------------
e <- new.env(); sys.source("projects/BMLEH_Los1_CAPRI/gamble_model/target_rules.R", e)
sh <- fread("output/gamble_model/BMLEH_Los1_CAPRI/crop_shares_nuts2.csv")
suppressMessages(library(arrow))
gf <- sort(list.files(GW, "^one_kmID_master_mapping_.*\\.parquet$", full.names=TRUE))
g <- as.data.table(read_parquet(gf[length(gf)]))
geo <- unique(g[, .(join_id = as.integer(INSPIRE_Europe_buffer_1kmID), geo = as.character(NUTS2))])
geo <- geo[!is.na(join_id) & !duplicated(join_id)]
d3 <- apply_target_classification(copy(d2), get("BMLEH_TARGET_RULES", e), sh, geo, verbose = FALSE)
A3 <- d3[, .(a3 = sum(area)), by = join_id]
k <- merge(A2, A3, by = "join_id", all = TRUE); k[is.na(k)] <- 0
cat(sprintf("\nEurostat cascade:\n"))
ok(max(abs(k$a2 - k$a3)) < 1e-6, sprintf("conserves area PER CELL (max dev %.2g)", max(abs(k$a2-k$a3))))
ok(abs(sum(A0$a0) - sum(A3$a3)) < 1, sprintf("END-TO-END total preserved: %.1f -> %.1f km2", sum(A0$a0), sum(A3$a3)))
miss <- setdiff(A0$join_id, A3$join_id)
cat(sprintf("  cells lost end-to-end: %d (%.4f%% of area)\n", length(miss), 100*A0[join_id %in% miss, sum(a0)]/sum(A0$a0)))
cat(sprintf("  geo coverage: %.2f%% of cascade-split area had an exact NUTS2 share\n",
    100 * mean(d2[model_class %in% get("BMLEH_TARGET_RULES", e)[action=="split", from_class], join_id] %in% geo$join_id)))
cat(sprintf("\nfinal: %d classes | cropland %.0f km2 of %.0f km2 total\n",
    uniqueN(d3$model_class), d3[grepl("^Cropland", model_class), sum(area)], sum(A3$a3)))

# ---- 5. does the produced class list match the target? -----------------------------------------
tm <- tryCatch(fread("../cascadinggamble-core/data/BMLEH_Los1_thematic_mapping.csv"), error = function(e) NULL)
if (!is.null(tm)) {
  tgt <- sort(unique(tm$BMLEH_Los1_DS_reporting_fine)); got <- sort(unique(d3$model_class))
  cat(sprintf("\nclass list vs target (%d target / %d produced):\n", length(tgt), length(got)))
  cat(sprintf("  in target, not produced: %s\n", paste(setdiff(tgt, got), collapse = ", ")))
  cat(sprintf("  produced, not in target: %s\n", paste(setdiff(got, tgt), collapse = ", ")))
  ar <- d3[, .(km2 = sum(area)), by = model_class][order(-km2)]
  cat("\n  smallest produced classes (candidates for being dropped by the zero-area filter):\n")
  print(head(ar[order(km2)][, .(model_class, km2 = round(km2, 1))], 6))
}

fwrite(d3[, .(km2 = sum(area)), by = model_class][order(-km2)],
       "output/gamble_model/BMLEH_Los1_CAPRI/class_areas_1km.csv")
cat("\nwrote output/gamble_model/BMLEH_Los1_CAPRI/class_areas_1km.csv\n")
