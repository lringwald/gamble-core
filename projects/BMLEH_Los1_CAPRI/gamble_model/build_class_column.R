#!/usr/bin/env Rscript
# =============================================================================
# BMLEH_Los1_CAPRI — generate the BMLEH_Los1_label column in the LUM mapping
# =============================================================================
# The project needs classes that AgMIP_label does not expose, even though the LUM map itself knows
# them. Rather than patch this downstream, give the project its own class column and run the driver
# with DRIVER_CLASS_COLS=BMLEH_Los1_label.
#
# Division of labour, deliberately:
#   THIS COLUMN  handles what the LUM map already distinguishes (short rotation, grassland
#                intensity). No external statistic is involved, so it belongs at mapping level.
#   target_rules.R handles what needs crop statistics to resolve (the HRL crop classes and the
#                arable/permanent residuals). Those cannot be decided per LUM code.
#
# Base is AgMIP_label, changed in exactly two places:
#
#   1. LUM_Code 7 "Short rotation" -> Cropland_permanent_energy
#      AgMIP calls this Forests_managed, so short-rotation coppice was invisible: it sat inside
#      managed forest with long rotation (702k km2) and the energy class never materialised. It is
#      3,358 km2. GLOBIOM_subclass calls it Forests_SR, but switching CLASS_COLS to that column
#      wholesale is not an option -- the crop cascade keys on the AgMIP cropland names.
#      The target name comes from the authoritative table (data/BMLEH_Los1_thematic_mapping.csv):
#      CAPRI_CODE NECR "New energy crops", LUM_Code 7 -> Cropland_PERMANENT_energy. Short-rotation
#      coppice is a permanent crop in CAPRI, not an arable one; this was Cropland_arable_energy here
#      until 2026-09-11, which is a class that does not exist in the target.
#
#   2. Pasture/grassland codes -> Grassland_intensive | Grassland_extensive
#      CONFIRMED against data/BMLEH_Los1_thematic_mapping.csv, which lists the codes explicitly:
#        GRAI intensive  15|16|19|20|4012              (+ organic twins)
#        GRAE extensive  17|18|21|22|23|24|25|4013     (+ organic twins)
#      That matches the assignment below exactly, including Rough grazing (23) as EXTENSIVE -- which
#      GLOBIOM_mngmt marks HI, so deriving intensity from that column would have got it wrong.
#      LUM_Code 26 "Unmanaged semi-natural and natural grassland" is NOT grassland here -- it stays
#      Natural_unmanaged, exactly as AgMIP has it.
#
# NOT handled here, by decision: Natural_unmanaged_fallow / _set_aside. Natural_unmanaged is the
# single class in the prior and the downscaling engine allocates it to fallow/set-aside.
#
#   Rscript projects/BMLEH_Los1_CAPRI/gamble_model/build_class_column.R
# =============================================================================
suppressMessages(library(data.table))
MAP <- Sys.getenv("BMLEH_MAPPING", "aux_files/LUM_Code_to_macro_model_mapping.csv")
if (!file.exists(MAP)) stop("mapping not found: ", MAP, " (run from the gamble-core root)")

# Pure renames onto the target's own vocabulary (data/BMLEH_Los1_thematic_mapping.csv).
# Forests is deliberately NOT among them: the target carries a single "Forests" class, but we keep
# Forests_primary and Forests_managed separate -- the managed/primary distinction is real and the
# downscaler can collapse it, whereas we could not recover it after collapsing here.
RENAME_TO_TARGET <- c(Built_up_area = "Artificial")
GRASS_INTENSIVE <- c(15, 16, 19, 20, 4012)          # "very high density" / "high density" / "intensive heterogeneous"
GRASS_EXTENSIVE <- c(17, 18, 21, 22, 23, 24, 25, 4013)  # moderate/low density, rough grazing, silvo-pastoral,
                                                        # managed semi-natural, "extensive heterogeneous"
ENERGY_CODE <- 7   # -> Cropland_permanent_energy (CAPRI NECR)

m <- fread(MAP, colClasses = c(LUM_Code = "character"))
if (!"AgMIP_label" %in% names(m)) stop("mapping has no AgMIP_label column to base this on")
m[, .code := suppressWarnings(as.integer(LUM_Code))]
# organic twins carry code <base>999 and are dropped unless INCLUDE_ORGANIC; give them their base
# code's class anyway so the column is complete and the _O suffix logic can work if it is switched on
m[, .base := fifelse(!is.na(.code) & .code > 999 & .code %% 1000 == 999, .code %/% 1000, .code)]

m[, BMLEH_Los1_label := as.character(AgMIP_label)]
# carry the base class onto the organic rows, where AgMIP_label is blank
blank <- m[, !nzchar(trimws(BMLEH_Los1_label)) & !is.na(.base)]
if (any(blank)) {
  lk <- m[!is.na(.code) & .code == .base, .(.base, cls = BMLEH_Los1_label)]
  m[blank, BMLEH_Los1_label := lk$cls[match(.base, lk$.base)]]
  m[is.na(BMLEH_Los1_label), BMLEH_Los1_label := ""]
}
for (.from in names(RENAME_TO_TARGET))
  m[BMLEH_Los1_label == .from, BMLEH_Los1_label := RENAME_TO_TARGET[[.from]]]
m[.base %in% ENERGY_CODE,     BMLEH_Los1_label := "Cropland_permanent_energy"]
m[.base %in% GRASS_INTENSIVE, BMLEH_Los1_label := "Grassland_intensive"]
m[.base %in% GRASS_EXTENSIVE, BMLEH_Los1_label := "Grassland_extensive"]

chg <- m[trimws(as.character(AgMIP_label)) != trimws(BMLEH_Los1_label) & nzchar(trimws(AgMIP_label))]
cat(sprintf("\n%d LUM codes differ from AgMIP_label:\n", nrow(chg)))
print(chg[, .(LUM_Code, LUM_label = substr(LUM_label, 1, 46), AgMIP_label, BMLEH_Los1_label)])
cat(sprintf("\n%d organic rows filled from their base code (blank in AgMIP_label)\n", sum(blank)))

m[, c(".code", ".base") := NULL]
setcolorder(m, c(setdiff(names(m), "BMLEH_Los1_label"), "BMLEH_Los1_label"))
fwrite(m, MAP)
cat(sprintf("\nwrote %s\n  distinct BMLEH_Los1_label: %s\n", MAP,
            paste(sort(setdiff(unique(m$BMLEH_Los1_label), "")), collapse = ", ")))
