#!/usr/bin/env Rscript
# =============================================================================
# Eurostat crop SHARE preparation for the BMLEH_Los1 crop-type carve-out.
#
# WHY: HRL Crop Types gives the SPATIAL pattern but only a 17-class legend, while BMLEH_Los1
# asks for 33 crop classes. Several HRL classes therefore cover more than one target class
# (Wheat -> soft + durum, Maize -> grain + fodder, Other_cereals -> rye + oats + other). This
# script produces the REGIONAL SHARES that split them, from Eurostat NUTS2 crop areas.
#
# Same cascade as prep/prepare_eurostat_livestock_shares.R: a spatial layer carries the pattern,
# Eurostat carries the composition, and every fallback is reported rather than silently applied.
#
#   Rscript prep/prepare_eurostat_crop_shares.R
#
# Outputs (output/eurostat/):
#   crop_shares_nuts2.csv   group x geo x target_class -> share (sums to 1 within group x geo)
#   crop_coverage.csv       which resolution level each group x geo was answered at
#
# Years: 2017-2019 averaged, matching the HRL layer (Crop_Types_Avg_2017_2019_1km.rds). Averaging
# smooths crop rotation, which is the point -- a single year misstates the arable mix.
# =============================================================================
suppressMessages({library(eurostat); library(data.table)})
YEARS   <- as.integer(strsplit(Sys.getenv("EU_CROP_YEARS", "2017,2018,2019"), ",")[[1]])
OUT_DIR <- Sys.getenv("EU_CROP_OUT", "output/eurostat")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- the splits ------------------------------------------------------------------------------
# Each entry: the HRL class (or "<residual>") that needs splitting -> BMLEH target classes and the
# Eurostat codes that measure them. A target may take SEVERAL codes (summed) where Eurostat is
# finer than BMLEH. Codes verified present in apro_cpshr; NUTS2 coverage is in crop_coverage.csv.
SPLITS <- list(
  Wheat = list(
    Cropland_arable_softwheat  = "C1110",   # Common wheat and spelt
    Cropland_arable_durumwheat = "C1120"),  # Durum wheat
  Maize = list(
    Cropland_arable_maize        = "C1500", # Grain maize and corn-cob-mix
    Cropland_arable_fodder_maize = "G3000"),# Green maize
  Other_cereals = list(
    Cropland_arable_rye          = "C1200", # Rye and winter cereal mixtures (maslin)
    Cropland_arable_oats         = "C1410", # Oats
    Cropland_arable_other_cereal = c("C1420","G9100")),
  Fruits = list(                            # HRL "Fruits" also carries citrus; T0000 separates it
    Cropland_permanent_citrus      = "T0000",
    Cropland_permanent_other_fruit = "F0000")
)
# Classes HRL cannot see at all. These are carved out of the arable RESIDUAL in proportion to
# their Eurostat area, so the residual is decomposed rather than dumped into one bucket.
RESIDUAL_SPLIT <- list(
  Cropland_arable_fodder_other     = c("G2100","G2900","G9900"),
  Cropland_arable_fodder_rootcrops = "R9000",
  Cropland_arable_other_oil        = c("I1140","I1150","I1190"),
  Cropland_arable_tobacco          = "I3000",
  Cropland_permanent_flowers       = "N0000",
  Cropland_permanent_nurseries     = "L0000"
)
# NOT sourced here, by decision:
#   energy        <- Forests_SR in the LUM map (short-rotation coppice)
#   fallow,       <- LUM map directly; Natural_unmanaged is the downscaling fallback
#   set_aside
#   grassland ext/int <- LUM management classes
#   other_industrial  <- the leftover Cropland_arable_other residual
# STILL OPEN (Eurostat aggregate only, needs a finer dataset):
#   tomatoes vs other_veg (V0000_S0000), apples vs other_fruit (F0000), grapes vs wine (W1000)

cat(sprintf("\nfetching apro_cpshr for %s ...\n", paste(YEARS, collapse=", ")))
raw <- as.data.table(get_eurostat("apro_cpshr", time_format = "num"))
if (!"strucpro" %in% names(raw)) stop("apro_cpshr: no strucpro dimension; the API layout changed.")
# the year column is `time` on a filtered fetch and `TIME_PERIOD` on a bulk one -- normalise
tcol <- intersect(c("time","TIME_PERIOD","period"), names(raw))[1]
if (is.na(tcol)) stop("apro_cpshr: no recognisable time column (", paste(names(raw), collapse=", "), ")")
setnames(raw, tcol, "yr")
raw[, yr := as.integer(substr(as.character(yr), 1, 4))]
# AREA in 1000 ha. AR_THS_HA ("area under cultivation") is the right measure for rotational
# crops, but NON-ROTATIONAL ones are only published as MAR_THS_HA ("main area") -- flowers
# (N0000) and nurseries (L0000) have ZERO AR_THS_HA rows, so an AR-only filter silently drops
# them and they come out as 0% of the residual everywhere. Prefer AR, fall back to MAR per code.
d <- raw[strucpro %in% c("AR_THS_HA","MAR_THS_HA") & yr %in% YEARS & !is.na(values)]
has_ar <- d[strucpro == "AR_THS_HA", unique(crops)]
d <- d[strucpro == fifelse(crops %in% has_ar, "AR_THS_HA", "MAR_THS_HA")]
cat(sprintf("  codes taken from MAR_THS_HA (no AR available): %s\n",
            paste(setdiff(unique(d$crops), has_ar), collapse=", ")))
if (!nrow(d)) stop("no AR_THS_HA rows for ", paste(YEARS, collapse=","))
d[, crops := as.character(crops)][, geo := as.character(geo)]
# average over the window: a single year misstates the mix under rotation
a <- d[, .(area = mean(values, na.rm = TRUE), nyr = .N), by = .(geo, crops)]
a[, lev := fifelse(nchar(geo) == 2L, "NUTS0", fifelse(nchar(geo) == 4L, "NUTS2", "NUTS1"))]
cat(sprintf("  %s rows | %d geos | %d codes | years present %s\n", format(nrow(a), big.mark=","),
            uniqueN(a$geo), uniqueN(a$crops), paste(sort(unique(d$yr)), collapse="/")))

# ---- share builder with an explicit fallback chain ---------------------------------------------
# NUTS2 -> its NUTS0 -> EU. Recorded per (group, geo) so a share that came from an EU average is
# never mistaken for a regional measurement.
shares_for <- function(gname, targets) {
  codes <- unlist(targets); tgt <- rep(names(targets), lengths(targets))
  sub <- a[crops %in% codes]
  if (!nrow(sub)) { warning("no data for group ", gname); return(NULL) }
  sub[, target := tgt[match(crops, codes)]]
  agg <- sub[, .(area = sum(area)), by = .(geo, lev, target)]
  wide <- dcast(agg, geo + lev ~ target, value.var = "area", fill = 0)
  tn <- names(targets); for (t in setdiff(tn, names(wide))) wide[, (t) := 0]
  wide[, tot := rowSums(.SD), .SDcols = tn]
  eu <- wide[geo == "EU" | geo == "EU27_2020"][1]
  if (is.null(eu) || !nrow(eu) || is.na(eu$tot) || eu$tot <= 0) {
    s <- wide[lev == "NUTS0", lapply(.SD, sum, na.rm = TRUE), .SDcols = tn]
    eu <- cbind(data.table(geo = "EU_derived", lev = "EU"), s)[, tot := rowSums(.SD), .SDcols = tn]
  }
  n2 <- wide[lev == "NUTS2"]; n0 <- wide[lev == "NUTS0"]
  # A region is only trusted if it reports on a comparable SET of targets to its country. Reporting
  # exactly one target does not mean that target is 100% of the mix -- it means the others were not
  # published. Without this, 91 regions with a single residual code gave it the entire residual
  # (tobacco came out at 19% EU-wide). Require >=60% of the parent's reported targets, and >=2.
  nrep <- function(r) sum(unlist(r[, ..tn]) > 0)
  MINFRAC <- as.numeric(Sys.getenv("EU_CROP_MIN_TARGET_FRAC", "0.6"))
  out <- rbindlist(lapply(seq_len(nrow(n2)), function(i) {
    r <- n2[i]; src <- "NUTS2"
    p <- n0[geo == substr(n2$geo[i], 1, 2)]
    need <- if (nrow(p)) max(2L, ceiling(MINFRAC * nrep(p))) else 2L
    if (r$tot <= 0 || nrep(r) < need) { if (nrow(p) && p$tot > 0) { r <- p; src <- "NUTS0" }
                                        else { r <- eu; src <- "EU" } }
    v <- unlist(r[, ..tn]); v <- v / sum(v)
    data.table(group = gname, geo = n2$geo[i], source = src, target = tn, share = as.numeric(v))
  }))
  # countries with no NUTS2 rows at all still need an answer at NUTS0
  miss <- setdiff(n0$geo, substr(n2$geo, 1, 2))
  if (length(miss)) out <- rbind(out, rbindlist(lapply(miss, function(g) {
    r <- n0[geo == g]; src <- "NUTS0"
    if (!nrow(r) || r$tot <= 0) { r <- eu; src <- "EU" }
    v <- unlist(r[, ..tn]); v <- v / sum(v)
    data.table(group = gname, geo = g, source = src, target = tn, share = as.numeric(v)) })))
  out
}
res <- rbindlist(c(lapply(names(SPLITS), function(g) shares_for(g, SPLITS[[g]])),
                   list(shares_for("<arable residual>", RESIDUAL_SPLIT))), use.names = TRUE)
res <- res[is.finite(share)]
fwrite(res, file.path(OUT_DIR, "crop_shares_nuts2.csv"))
cov <- res[, .(n_geo = uniqueN(geo)), by = .(group, source)]
cov <- dcast(cov, group ~ source, value.var = "n_geo", fill = 0)
fwrite(cov, file.path(OUT_DIR, "crop_coverage.csv"))
cat(sprintf("\nwrote %s (%s rows) and %s\n", file.path(OUT_DIR,"crop_shares_nuts2.csv"),
            format(nrow(res), big.mark=","), file.path(OUT_DIR,"crop_coverage.csv")))
cat("\n--- resolution each split was answered at (number of geos) ---\n"); print(cov)
cat("\n--- share per target: EU AREA-WEIGHTED vs mean-of-regions ---\n")
# the area-weighted figure is the one to sanity-check against published EU totals; the mean of
# regional shares is pulled toward zero by every region where the crop is simply absent
euw <- rbindlist(c(lapply(names(SPLITS), function(g) {
    tg <- SPLITS[[g]]; ar <- sapply(tg, function(cd) a[crops %in% cd & lev=="NUTS0", sum(area)])
    data.table(group=g, target=names(tg), eu_weighted=as.numeric(ar/sum(ar))) }),
  list({ tg <- RESIDUAL_SPLIT; ar <- sapply(tg, function(cd) a[crops %in% cd & lev=="NUTS0", sum(area)])
    data.table(group="<arable residual>", target=names(tg), eu_weighted=as.numeric(ar/sum(ar))) })))
mm <- res[, .(mean_of_regions = round(mean(share),4)), by=.(group,target)]
print(merge(euw, mm, by=c("group","target"))[order(group, -eu_weighted)][
      , .(group, target, eu_weighted=round(eu_weighted,4), mean_of_regions)])
bad <- res[, .(s = sum(share)), by = .(group, geo)][abs(s - 1) > 1e-8]
cat(sprintf("\nshares summing to 1 within group x geo: %s\n", if (!nrow(bad)) "OK" else sprintf("FAILED for %d", nrow(bad))))
