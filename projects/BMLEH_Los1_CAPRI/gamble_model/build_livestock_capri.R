#!/usr/bin/env Rscript
# =============================================================================
# BMLEH_Los1_CAPRI — livestock parameters, CAPRI-keyed, for the downscaler
# =============================================================================
#   Rscript projects/BMLEH_Los1_CAPRI/gamble_model/build_livestock_capri.R
#
# Reads the count-model plug-in table (species x subclass x D/O/F x system x country x driver, with
# gamma) and re-expresses it on CAPRI activity codes and CAPRI country codes.
#
# TWO THINGS IT DOES NOT DO, deliberately:
#  * It does not split one of our subclasses across the CAPRI activities beneath it. CAPRI is finer
#    on cattle by yield (DCOH/DCOL), weight (HEIH/HEIL, BULH/BULL) and sex (CAMR/CAFR/CAMF/CAFF);
#    the model carries none of those. Every CAPRI activity under a subclass therefore receives that
#    subclass's gamma UNCHANGED, with `shared_with` saying how many activities share it. That is the
#    correct statement -- "these activities have the same driver response as far as this prior knows"
#    -- rather than an invented allocation.
#  * It does not invent pigs. PIGF/SOWS/PKPL are reported as uncovered.
#
# Output: livestock_capri_parameters.csv
#   capri | capri_label | country | system | driver | gamma | subclass | shared_with | country_purity
# =============================================================================
suppressMessages({library(data.table); library(arrow)})
SRC  <- Sys.getenv("LS_PARAMS", "output/composition/subclass_country_parameters.csv")
OUT  <- Sys.getenv("LS_OUT", "output/gamble_model/BMLEH_Los1_CAPRI")
if (!file.exists(SRC)) stop("count-model parameters not found: ", SRC)
source("projects/BMLEH_Los1_CAPRI/gamble_model/livestock_capri_mapping.R")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

p <- fread(SRC)
cat(sprintf("source: %s rows | %d countries | %d subclasses | %d drivers\n",
            format(nrow(p), big.mark = ","), uniqueN(p$country), uniqueN(p$subclass), uniqueN(p$driver)))

# ---- countries: GLOB names -> CAPRI 2-char codes ------------------------------------------------
# The two do not agree cell-for-cell, so the relabel is by DOMINANT overlap and its purity is carried
# through to the output. Germany is 0.99; Switzerland maps at only 0.41, which a consumer must be able
# to see rather than discover later.
gw <- Sys.getenv("GAMBLE_GRIDWORK_DIR", "../LAMASUS_gridwork/output")
gf <- sort(list.files(gw, "^one_kmID_master_mapping_.*\\.parquet$", full.names = TRUE))
g  <- as.data.table(read_parquet(gf[length(gf)]))
g[, `:=`(cap2 = substr(as.character(CAPRI_NUTS), 1, 2), glob = as.character(GLOB_country))]
g <- g[!is.na(cap2) & nzchar(cap2) & !is.na(glob) & nzchar(glob)]
# IS THE SOURCE ALREADY CAPRI-KEYED? Since 2026-09-14 the count model can be fitted with
# DRIVER_RE_GROUP_COL=CAPRI_NUTS, in which case `country` already holds CAPRI 2-char codes and the
# crosswalk above is not merely unnecessary but WRONG: it inner-joins codes against GLOB names, so
# only the handful that coincide survive and the table silently collapses to a country or two (896
# rows over ONE country, when it should span 32). Detect it and pass the native keys straight
# through, at purity 1 -- a key that was fitted in CAPRI space is not a relabel of anything.
.src_countries <- unique(as.character(p$country))
.is_capri_keyed <- length(.src_countries) > 0 &&
  all(nchar(.src_countries) == 2) &&
  mean(.src_countries %in% unique(g$cap2)) > 0.8
if (.is_capri_keyed) {
  cat(sprintf("source is ALREADY CAPRI-keyed (%d codes: %s...) -- skipping the GLOB->CAPRI relabel\n",
              length(.src_countries), paste(head(sort(.src_countries), 6), collapse = ", ")))
  m <- copy(p)
  # cells = 1, NOT NA: the aggregation below is a CELL-WEIGHTED mean
  # (sum(gamma*cells)/sum(cells)), so an NA weight turns every gamma into NA. With native CAPRI keys
  # each country is its own code, nothing merges, and an equal weight reduces the mean to the value
  # itself -- which is the intended no-op.
  m[, `:=`(capri_country = country, country_purity = 1, cells = 1L)]
} else {
  # Only NOW build the GLOB -> CAPRI crosswalk. It used to be computed (and reported) unconditionally,
  # so a CAPRI-keyed run printed "33 GLOB names -> 32 CAPRI codes" and "DROPPED Switzerland->IT 0.41"
  # describing a relabel it then did not perform -- a message that reads as a real caveat about the
  # output when it is a caveat about nothing.
  MINPUR <- as.numeric(Sys.getenv("LS_MIN_PURITY", "0.70"))
  xw <- g[, .N, by = .(glob, cap2)][order(-N)]
  xw <- xw[, .(capri_country = cap2[1], country_purity = N[1] / sum(N), cells = sum(N)), by = glob]
  # A country with no CAPRI home of its own gets absorbed into its largest neighbour, which is not a
  # relabel but a mislabel: Switzerland has no CAPRI code and mapped to IT at 0.41 purity, which would
  # have put Swiss parameters inside Italy. Refuse below a purity floor rather than emit that.
  drop <- xw[country_purity < MINPUR]
  if (nrow(drop)) cat(sprintf("  DROPPED (purity < %.2f, no CAPRI home of their own): %s\n", MINPUR,
                              paste(sprintf("%s->%s %.2f", drop$glob, drop$capri_country, drop$country_purity), collapse = ", ")))
  xw <- xw[country_purity >= MINPUR]
  cat(sprintf("country crosswalk: %d GLOB names -> %d CAPRI codes | purity min %.2f median %.2f\n",
              nrow(xw), uniqueN(xw$capri_country), min(xw$country_purity), median(xw$country_purity)))

  m <- merge(p, xw, by.x = "country", by.y = "glob", all.x = FALSE)
}

# ---- subclasses -> CAPRI activities -------------------------------------------------------------
mm <- merge(m, LIVESTOCK_CAPRI_MAP, by = "subclass", allow.cartesian = TRUE)
mm[, shared_with := uniqueN(capri), by = .(subclass, country, system, driver)]
# Several of our countries can legitimately share one CAPRI code -- BL is Belgium+Luxembourg by
# definition. Combine them by a cell-weighted mean rather than emitting duplicate rows for the same
# key, and say how many sources went in.
# Rename the SOURCE country column first: naming the by-group `country` as well shadows it, and
# merged_from silently read 1 everywhere while the merge itself was working correctly.
setnames(mm, "country", "src_country")
out <- mm[, .(gamma = sum(gamma * cells) / sum(cells),
              subclass = paste(sort(unique(subclass)), collapse = "+"),
              shared_with = shared_with[1],
              merged_from = uniqueN(src_country),
              merged_names = paste(sort(unique(src_country)), collapse = "+"),
              country_purity = round(min(country_purity), 4)),
          by = .(capri, capri_label, country = capri_country, system, driver)]
setcolorder(out, c("capri","capri_label","country","system","driver","gamma",
                   "subclass","shared_with","merged_from","merged_names","country_purity"))
setorder(out, capri, country, system, driver)
fwrite(out, file.path(OUT, "livestock_capri_parameters.csv"))

cov <- rbind(
  out[, .(status = "mapped", n_country = uniqueN(country), n_driver = uniqueN(driver),
          shared_with = shared_with[1]), by = .(capri, capri_label, subclass)],
  LIVESTOCK_CAPRI_UNCOVERED[, .(capri, capri_label = label, subclass = NA_character_,
          status = "NOT COVERED", n_country = 0L, n_driver = 0L, shared_with = NA_integer_)],
  fill = TRUE)
fwrite(cov, file.path(OUT, "livestock_capri_coverage.csv"))

cat(sprintf("\nwrote %s  (%s rows)\n", file.path(OUT, "livestock_capri_parameters.csv"),
            format(nrow(out), big.mark = ",")))
cat(sprintf("      %s\n\n", file.path(OUT, "livestock_capri_coverage.csv")))
cat("coverage:\n"); print(cov[, .(capri, subclass, status, shared_with)], nrows = 30)
cat(sprintf("\n%d of the 16 CAPRI leaf activities are covered; %d are not (%s).\n",
            uniqueN(out$capri), nrow(LIVESTOCK_CAPRI_UNCOVERED),
            paste(LIVESTOCK_CAPRI_UNCOVERED$capri, collapse = ", ")))
