#!/usr/bin/env Rscript
# =============================================================================
# Eurostat crop SHARE preparation for the BMLEH_Los1 crop-type carve-out.
#
# LIVES HERE, in the project's gamble_model/ folder: the split map below is entirely
# BMLEH_Los1's target classification (data/BMLEH_Los1_thematic_mapping.csv). gamble-core hosts the
# samplers and generic scripts; everything project-specific about producing THIS project's
# econometric prior -- custom processing, the model spec, the run wrapper -- sits here together.
#
# WHY: HRL Crop Types gives the SPATIAL pattern but only a 17-class legend, while BMLEH_Los1
# asks for 33 crop classes. Several HRL classes therefore cover more than one target class
# (Wheat -> soft + durum, Maize -> grain + fodder, Other_cereals -> rye + oats + other). This
# script produces the REGIONAL SHARES that split them, from Eurostat NUTS2 crop areas.
#
# Same cascade as prep/prepare_eurostat_livestock_shares.R: a spatial layer carries the pattern,
# Eurostat carries the composition, and every fallback is reported rather than silently applied.
#
#   Rscript scripts/BMLEH_Los1_CAPRI/gamble_model/prepare_eurostat_crop_shares.R   (from repo root)
#
# Outputs (results/gamble_model/BMLEH_Los1_CAPRI/): these are PRIOR-PRODUCTION artefacts, so they
# follow the gamble output convention (alongside NCUT_OUT_DIR = results/prior_nested_cut) rather
# than the downscaling side's data/02_intermediate.
#   crop_shares_nuts2.csv   group x geo x target_class -> share (sums to 1 within group x geo)
#   crop_coverage.csv       which resolution level each group x geo was answered at
#
# Years: 2017-2019 averaged, matching the HRL layer (Crop_Types_Avg_2017_2019_1km.rds). Averaging
# smooths crop rotation, which is the point -- a single year misstates the arable mix.
#
# ---- THE RESOLUTION TRADE-OFF, AND WHY IT IS NOT A CHOICE --------------------------------------
# Eurostat forces a trade: apro_cpshr has 79 codes at NUTS2 (spatial detail, thematic coarse);
# apro_cpsh1 has 214 codes at NUTS0 (thematic detail, spatially coarse). Neither alone reaches the
# BMLEH target.
#
# They COMPOSE, because the aggregates a fine code sits inside are themselves published at NUTS2
# (F0000 in 252 regions, W1000 in 242, V0000_S0000 in 198). So each level contributes what it has:
#
#     area(apples, pixel) = HRL "Fruits" area          <- 1 km spatial pattern
#                         x share(F0000 | NUTS2)       <- regional fruit composition
#                         x share(F1110 | F0000, NUTS0)<- national apple fraction
#
# The NUTS0 factor is constant within a country, so apples do not vary regionally BEYOND what the
# fruit aggregate and the HRL pattern already say. That is a real limitation, not a hidden one:
# every such share is tagged source = "NUTS0-only" in the output so it cannot be mistaken for a
# measured regional value.
# =============================================================================
suppressMessages({library(eurostat); library(data.table)})
YEARS   <- as.integer(strsplit(Sys.getenv("EU_CROP_YEARS", "2017,2018,2019"), ",")[[1]])
# This script lives in the project folder and runs from either root: the cascadinggamble repo
# (where the project folder finally belongs) or gamble-core (where it is developed). Anchor the
# output on whichever we are in rather than hard-failing on one of them.
.IN_CASCADE <- file.exists("scripts/gamble_model_run.R")
.IN_CORE    <- file.exists("codes/mnlogit_rcpp_sym.R")
if (!.IN_CASCADE && !.IN_CORE)
  stop("run from the cascadinggamble repo root or the gamble-core root")
# PREP PRODUCTS LIVE UNDER prep/, NOT gamble_model/. The root here switches by which REPO you are
# in, while fitted runs are written to results/gamble_model/<PROJECT>/ by estimate_prior.R -- so the
# old "gamble_model/<PROJECT>" subpath collided with runs whenever both landed under the same root,
# with nothing in the path saying which was which. `prep/` separates them in both repos.
OUT_DIR <- Sys.getenv("EU_CROP_OUT",
  file.path(if (.IN_CASCADE) "results" else "output", "prep/BMLEH_Los1_CAPRI"))
# Read-side compatibility: a table written before this move is still where it was. Consumers glob
# both (see run.sh), and this keeps the OLD directory usable if it is the only one populated.
.OUT_LEGACY <- file.path(if (.IN_CASCADE) "results" else "output", "gamble_model/BMLEH_Los1_CAPRI")
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
  Cropland_arable_other_industrial = c("I4000","I5000","I9000"),  # hops, aromatics, other industrial
  Cropland_permanent_flowers       = "N0000",
  Cropland_permanent_nurseries     = "L0000",
  # CAPRI OCRO "Other crops production activity" -> the arable area no identified crop accounts for.
  # There is no Eurostat code for it: everything left uncovered is either an aggregate (C1000, PECR,
  # I1100) or a sub-split of a code already in use. So it is DERIVED below as
  # ARA - (all identified arable crops), clamped at zero, and injected under this name. Without it the
  # residual split hands 100% of unattributed arable to seven minor crops, which is measurably too
  # much -- green fodder came out at 1.26x its Eurostat area.
  Cropland_arable_other_crop       = "<OCRO derived>"
)

# ---- the PERMANENT residual --------------------------------------------------------------------
# LUM permanent cropland that HRL could not attribute to a specific permanent crop. Without this it
# was a bare rename into other_fruit, which buried ~79k km2 -- three quarters of all permanent land
# -- in one class and left olives at 0.07 of their Eurostat area. HRL sees only 44.6k km2 of the
# 118.2k km2 Eurostat reports, so the residual is where most real olive and vine area actually is.
# These four codes partition permanent crops exactly: Eurostat defines F0000 as "fruits, berries and
# nuts EXCLUDING citrus fruits, grapes and strawberries", so citrus (T0000) and grapes (W1000) sit
# beside it while nuts (F4000, 11.8k km2) sit INSIDE it -- no double counting, and they sum to
# 118,237 km2, the permanent total. All four are published at NUTS2, so this splits at the same
# resolution as the arable residual rather than falling back to national shares.
PERM_RESIDUAL_SPLIT <- list(
  Cropland_permanent_olives      = "O1000",
  Grapes                         = "W1000",   # HRL's own class name: feeds the later Grapes_detail split
  Cropland_permanent_citrus      = "T0000",
  Cropland_permanent_other_fruit = "F0000"    # feeds the later Fruits_detail (apples) split
)

# ---- FINE splits, NUTS0 only -------------------------------------------------------------------
# apro_cpshr (NUTS2) carries only aggregates for these; the detail is in apro_cpsh1, which has 214
# codes but is published at COUNTRY level. So these shares are constant within a country -- a real
# limitation, recorded as source = "NUTS0-only" so it is never mistaken for regional variation.
# Each: the aggregate `total`, the measured `parts`, and the target taking total - sum(parts).
SPLITS_FINE <- list(
  Fruits_detail = list(total = "F0000",              # excludes citrus, grapes, strawberries
                       parts = list(Cropland_permanent_apples = "F1110"),
                       remainder = "Cropland_permanent_other_fruit"),
  Veg_detail    = list(total = "V0000",              # fresh vegetables incl. melons
                       parts = list(Cropland_arable_tomatoes = "V3100"),
                       remainder = "Cropland_arable_other_veg"),
  Grapes_detail = list(total = "W1000",
                       parts = list(Cropland_permanent_wine = "W1100"),   # grapes for wines
                       remainder = "Cropland_permanent_grapes")           # table/other grapes
)
# NOT sourced here, by decision:
#   energy        <- Forests_SR in the LUM map (short-rotation coppice)
#   fallow,       <- LUM map directly; Natural_unmanaged is the downscaling fallback
#   set_aside
#   grassland ext/int <- LUM management classes
#   other_industrial  <- the leftover Cropland_arable_other residual
# STILL OPEN (Eurostat aggregate only, needs a finer dataset):
#   tomatoes vs other_veg (V0000_S0000), apples vs other_fruit (F0000), grapes vs wine (W1000)

# Arable crop codes that, summed, should account for arable land. Used only to derive OCRO.
OCRO_ACCOUNTED <- c("C1000","C1500","R1000","R2000","I1110","I1120","I1130","I1140","I1150","I1190",
                    "I3000","I4000","I5000","I9000","N0000","L0000","G1000","G2000","G3000",
                    "G2100","G2900","G9900","R9000","P0000","V0000_S0000","Q0000")
OCRO_TOTAL <- "ARA"   # arable land
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
d[, `:=`(geo = as.character(geo), crops = as.character(crops))]
# Resolve AR-vs-MAR PER (geo, crop), not per crop globally. Availability varies by COUNTRY, not just
# by code: Germany publishes wine grapes (W1000) only as MAR_THS_HA while other countries publish AR.
# A global rule pins W1000 to AR and then silently drops Germany -- a real producer -- from the split.
d[, has_ar := any(strucpro == "AR_THS_HA"), by = .(geo, crops)]
.mar <- unique(d[has_ar == FALSE, .(geo, crops)])
d <- d[strucpro == fifelse(has_ar, "AR_THS_HA", "MAR_THS_HA")][, has_ar := NULL]
cat(sprintf("  geo x crop pairs taken from MAR_THS_HA (no AR published there): %d%s\n", nrow(.mar),
            if (nrow(.mar)) sprintf(" e.g. %s", paste(head(paste0(.mar$geo, ":", .mar$crops), 6), collapse=", ")) else ""))
if (!nrow(d)) stop("no AR_THS_HA rows for ", paste(YEARS, collapse=","))
d[, crops := as.character(crops)][, geo := as.character(geo)]
# average over the window: a single year misstates the mix under rotation
a <- d[, .(area = mean(values, na.rm = TRUE), nyr = .N), by = .(geo, crops)]
a[, lev := fifelse(nchar(geo) == 2L, "NUTS0", fifelse(nchar(geo) == 4L, "NUTS2", "NUTS1"))]
# ---- derive the OCRO residual: arable land minus everything identified on it -------------------
.tot <- a[crops == OCRO_TOTAL, .(geo, lev, tot = area)]
.acc <- a[crops %in% OCRO_ACCOUNTED, .(acc = sum(area)), by = .(geo, lev)]
.oc  <- merge(.tot, .acc, by = c("geo", "lev"), all.x = TRUE)
.oc[is.na(acc), acc := 0][, area := pmax(tot - acc, 0)]
cat(sprintf("  OCRO derived for %d geos | EU total %.0f kha (%.1f%% of arable) | %d geos clamped to 0\n",
            nrow(.oc), .oc[lev == "NUTS0", sum(area)],
            100 * .oc[lev == "NUTS0", sum(area)] / pmax(.oc[lev == "NUTS0", sum(tot)], 1e-9),
            .oc[area <= 0, .N]))
a <- rbind(a, .oc[area > 0, .(geo, crops = "<OCRO derived>", area, nyr = 3L, lev)], fill = TRUE)
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
  # The count test alone is not enough. It asks how MANY targets a region reports, not whether it
  # reports the ones that matter: the arable residual has 7 targets, and a region publishing 5 minor
  # ones while omitting fodder -- 75% of the EU mix -- passes the count test and then reads as 0%
  # fodder. That hit 55 of 80 regions. So also require the reported targets to COVER most of the
  # parent's area mass, which is the quantity actually at stake.
  nrep <- function(r) sum(unlist(r[, ..tn]) > 0)
  MINFRAC <- as.numeric(Sys.getenv("EU_CROP_MIN_TARGET_FRAC", "0.6"))
  MINMASS <- as.numeric(Sys.getenv("EU_CROP_MIN_PARENT_MASS", "0.8"))
  .rej <- c(count = 0L, mass = 0L)
  out <- rbindlist(lapply(seq_len(nrow(n2)), function(i) {
    r <- n2[i]; src <- "NUTS2"
    p <- n0[geo == substr(n2$geo[i], 1, 2)]
    need <- if (nrow(p)) max(2L, ceiling(MINFRAC * nrep(p))) else 2L
    # share of the parent's mix carried by the targets this region actually reports
    covered <- if (nrow(p) && p$tot > 0) sum((unlist(p[, ..tn]) / p$tot)[unlist(r[, ..tn]) > 0]) else 1
    if (r$tot <= 0 || nrep(r) < need || covered < MINMASS) {
      if (r$tot > 0 && nrep(r) >= need) .rej["mass"]  <<- .rej["mass"] + 1L
      else                              .rej["count"] <<- .rej["count"] + 1L
      if (nrow(p) && p$tot > 0) { r <- p; src <- "NUTS0" }
                                        else { r <- eu; src <- "EU" } }
    v <- unlist(r[, ..tn]); v <- v / sum(v)
    data.table(group = gname, geo = n2$geo[i], source = src, target = tn, share = as.numeric(v))
  }))
  cat(sprintf("  %-18s %3d NUTS2 regions | rejected: %d too few targets, %d covering <%.0f%% of parent mass\n",
              gname, nrow(n2), .rej["count"], .rej["mass"], 100 * MINMASS))
  # EVERY country gets a 2-char row. This is the fallback the pixel cascade relies on, and it has to
  # exist even for countries already covered at NUTS2, because the grid's region codes and Eurostat's
  # do not agree region-for-region:
  #   * Eurostat files single-region countries under EEZZ / LVZZ ("not regionalised"). Those are
  #     4 chars, so they look like NUTS2 codes, but they match NO region on the grid (which has
  #     EE00 / LV00) -- and with no country row, Estonia and Latvia fell through to the EU-mean mix.
  #   * NUTS vintages drift: the grid has UKN1, Eurostat publishes UKN0.
  # Emitting the country row unconditionally makes both cases resolve nationally instead of EU-wide.
  ctry <- setdiff(unique(c(n0$geo[nchar(n0$geo) == 2L], substr(n2$geo, 1, 2))), c("EU", "XX"))
  if (length(ctry)) out <- rbind(out, rbindlist(lapply(ctry, function(g) {
    r <- n0[geo == g]; src <- "NUTS0"
    if (!nrow(r) || r$tot <= 0) {            # no national row -> aggregate the country's own regions
      k <- n2[substr(geo, 1, 2) == g]
      if (nrow(k) && sum(k$tot) > 0) {
        r <- cbind(data.table(geo = g, lev = "NUTS0"),
                   k[, lapply(.SD, sum, na.rm = TRUE), .SDcols = tn])[, tot := rowSums(.SD), .SDcols = tn]
        src <- "NUTS0-agg"
      } else { r <- eu; src <- "EU" }
    }
    v <- unlist(r[, ..tn]); v <- v / sum(v)
    data.table(group = gname, geo = g, source = src, target = tn, share = as.numeric(v)) })))
  out
}
res <- rbindlist(c(lapply(names(SPLITS), function(g) shares_for(g, SPLITS[[g]])),
                   list(shares_for("<arable residual>", RESIDUAL_SPLIT)),
                   list(shares_for("<permanent residual>", PERM_RESIDUAL_SPLIT))), use.names = TRUE)

# ---- fine splits from apro_cpsh1 (country level) ------------------------------------------------
cat("\nfetching apro_cpsh1 (country-level detail: apples, tomatoes, wine grapes) ...\n")
r1 <- as.data.table(get_eurostat("apro_cpsh1", time_format = "num"))
t1 <- intersect(c("time","TIME_PERIOD","period"), names(r1))[1]; setnames(r1, t1, "yr")
r1[, yr := as.integer(substr(as.character(yr), 1, 4))]
d1 <- r1[strucpro %in% c("AR_THS_HA","MAR_THS_HA") & yr %in% YEARS & !is.na(values)]
d1[, `:=`(geo = as.character(geo), crops = as.character(crops))]
d1[, has_ar := any(strucpro == "AR_THS_HA"), by = .(geo, crops)]   # per country, see above
.mar1 <- unique(d1[has_ar == FALSE, .(geo, crops)])
d1 <- d1[strucpro == fifelse(has_ar, "AR_THS_HA", "MAR_THS_HA")][, has_ar := NULL]
cat(sprintf("  geo x crop pairs taken from MAR_THS_HA: %d%s\n", nrow(.mar1),
            if (nrow(.mar1)) sprintf(" e.g. %s", paste(head(paste0(.mar1$geo, ":", .mar1$crops), 6), collapse=", ")) else ""))
a1 <- d1[nchar(as.character(geo)) == 2L, .(area = mean(values, na.rm = TRUE)), by = .(geo, crops)]
fine <- rbindlist(lapply(names(SPLITS_FINE), function(g) {
  sp <- SPLITS_FINE[[g]]; pn <- names(sp$parts); pc <- unlist(sp$parts)
  tot <- a1[crops == sp$total, .(geo, tot = area)]
  prt <- a1[crops %in% pc, .(p = sum(area)), by = geo]
  j <- merge(tot, prt, by = "geo", all.x = TRUE); j[is.na(p), p := 0]
  j <- j[tot > 0]
  if (!nrow(j)) return(NULL)
  # A part exceeding its own aggregate means total and part did not come from the same measure: the
  # AR/MAR choice is made per (geo, crop), so a country publishing the total only as MAR and the part
  # only as AR ends up comparing the two. Clamp and report rather than emit a negative remainder.
  # As of 2026-09 this fires for DE/Grapes only, where the two measures are the same quantity
  # (W1000 MAR ~100.1 vs W1100 AR ~100.2 kha -- German grapes are essentially all wine grapes), so
  # the clamp lands on the right answer. Re-read this warning if it ever names another country.
  bad <- j[p > tot]
  if (nrow(bad)) cat(sprintf("   [warn] %s: part > total for %d country(ies), clamped to 100%%: %s\n",
                             g, nrow(bad), paste(sprintf("%s (%.3g vs %.3g)", bad$geo, bad$p, bad$tot), collapse=", ")))
  j[, p := pmin(p, tot)]
  rbind(data.table(group = g, geo = j$geo, source = "NUTS0-only", target = pn, share = j$p / j$tot),
        data.table(group = g, geo = j$geo, source = "NUTS0-only", target = sp$remainder,
                   share = (j$tot - j$p) / j$tot))
}), use.names = TRUE)
res <- rbind(res, fine, use.names = TRUE)
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
euw_fine <- rbindlist(lapply(names(SPLITS_FINE), function(g) {
  sp <- SPLITS_FINE[[g]]
  tt <- a1[crops == sp$total, sum(area)]; pp_ <- a1[crops %in% unlist(sp$parts), sum(area)]
  if (!length(tt) || tt <= 0) return(NULL)
  data.table(group = g, target = c(names(sp$parts), sp$remainder),
             eu_weighted = c(pp_/tt, (tt-pp_)/tt)) }))
euw <- rbindlist(c(lapply(names(SPLITS), function(g) {
    tg <- SPLITS[[g]]; ar <- sapply(tg, function(cd) a[crops %in% cd & lev=="NUTS0", sum(area)])
    data.table(group=g, target=names(tg), eu_weighted=as.numeric(ar/sum(ar))) }),
  list({ tg <- RESIDUAL_SPLIT; ar <- sapply(tg, function(cd) a[crops %in% cd & lev=="NUTS0", sum(area)])
    data.table(group="<arable residual>", target=names(tg), eu_weighted=as.numeric(ar/sum(ar))) }),
  list({ tg <- PERM_RESIDUAL_SPLIT; ar <- sapply(tg, function(cd) a[crops %in% cd & lev=="NUTS0", sum(area)])
    data.table(group="<permanent residual>", target=names(tg), eu_weighted=as.numeric(ar/sum(ar))) })))
euw <- rbind(euw, euw_fine, use.names = TRUE)
mm <- res[, .(mean_of_regions = round(mean(share),4)), by=.(group,target)]
print(merge(euw, mm, by=c("group","target"))[order(group, -eu_weighted)][
      , .(group, target, eu_weighted=round(eu_weighted,4), mean_of_regions)])
bad <- res[, .(s = sum(share)), by = .(group, geo)][abs(s - 1) > 1e-8]
cat(sprintf("\nshares summing to 1 within group x geo: %s\n", if (!nrow(bad)) "OK" else sprintf("FAILED for %d", nrow(bad))))
