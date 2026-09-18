# =====================================================================================
# build_globiom_subclass.R -- add a COMPLETE class column to the LUM thematic mapping.
#
# WHY. `GLOBIOM_mngmt` carries management only for MANAGED land (HI/LI/IRO/O ...) and is
# BLANK on all 38 unmanaged rows. So the composed class `GLOBIOM_UNFCCC + GLOBIOM_mngmt`
# resolves Forests/Cropland/Pasture into intensities but collapses every natural cover
# into one bucket: `Natural_unmanaged` (12.7% of area) lumps natural grassland, moors,
# sclerophyllous scrub, transitional woodland-shrub AND three LUM residual codes.
# Being the only Natural_* class it also sits as a ROOT SINGLETON in the nest tree, so it
# gets no within-nest structure at all.
#
# `GLOBIOM_subclass` closes that gap: ONE self-contained class column (same role as
# `AgMIP_label`) that is the composed managed name where management exists, and a real
# unmanaged subclass where it does not. Use via DRIVER_CLASS_COLS="GLOBIOM_subclass".
#
# THE TAXONOMY IS NOT INVENTED HERE -- it is read off `BIOCLIMA_DS_reporting`, which
# already distinguishes exactly these covers (and, usefully, already decides where the
# ambiguous rows go: LUM 324 transitional woodland-shrub -> shrubland, and the residual
# "Leftover forest" -> shrubland, "Leftover cropland"/"Leftover grassland" -> grassland).
#
# BACKWARD COMPATIBILITY IS CHECKED, NOT ASSUMED: for every row with non-blank mngmt the
# script asserts GLOBIOM_subclass == paste(GLOBIOM_UNFCCC, GLOBIOM_mngmt, sep="_"), so
# every currently-existing managed class name is reproduced EXACTLY. Idempotent.
#
# NOT DONE HERE (deliberate): `Cropland_other` has the same blank-mngmt gap and
# BIOCLIMA_DS_reporting would split it into permanent low/medium/high. Left alone --
# it is outside the nature/no_choice scope and would interact with the `_O` organic
# suffix (Cropland_other_LI_O ...). The rule below is the place to add it.
# =====================================================================================
suppressMessages(library(data.table))

MAP_FILE <- Sys.getenv("LUM_MAPPING_FILE",
                       "../LAMASUS_downscaling/aux_files/LAMASUS_LUM_thematic_mapping.csv")

build_globiom_subclass <- function(m) {
  m <- as.data.table(copy(m))
  req <- c("GLOBIOM_UNFCCC", "GLOBIOM_mngmt", "BIOCLIMA_DS_reporting")
  miss <- setdiff(req, names(m))
  if (length(miss)) stop("mapping is missing required column(s): ", paste(miss, collapse = ", "))

  .blank <- function(x) is.na(x) | trimws(as.character(x)) == ""
  unf <- trimws(as.character(m$GLOBIOM_UNFCCC))
  mng <- trimws(as.character(m$GLOBIOM_mngmt))
  rep <- trimws(as.character(m$BIOCLIMA_DS_reporting))

  # 1) MANAGED: exactly the historical composed name.
  sub <- ifelse(!.blank(mng) & !.blank(unf), paste(unf, mng, sep = "_"), NA_character_)

  # 2) UNMANAGED NATURAL: refine from the reporting taxonomy.
  #    Natural_other is kept WHOLE -- 331-335 individually are each well under a percent,
  #    and the smallest existing class (Forests_SR, 0.084%) is already the worst-behaved.
  nat <- is.na(sub) & unf == "Natural_unmanaged"
  sub[nat & grepl("Grassland_1_unmanaged", rep, fixed = TRUE)] <- "Natural_grassland"
  sub[nat & grepl("Natural_1_shrubland",   rep, fixed = TRUE)] <- "Natural_shrubland"

  # 3) EVERYTHING ELSE (Urban, Waterbodies, Wetlands, NODATA, Natural_protected,
  #    Natural_other, Cropland_other): keep the UNFCCC class unchanged.
  sub[is.na(sub)] <- unf[is.na(sub)]

  # --- assertions -------------------------------------------------------------------
  keep <- !.blank(mng) & !.blank(unf)
  bad  <- which(keep & sub != paste(unf, mng, sep = "_"))
  if (length(bad)) stop("managed rows would be RENAMED (backward-compat broken): ",
                        paste(utils::head(m$LUM_Code[bad], 5), collapse = ", "))
  still <- which(sub == "Natural_unmanaged")
  if (length(still)) warning(sprintf(
    "%d row(s) stayed Natural_unmanaged (no reporting taxonomy): LUM_Code %s",
    length(still), paste(m$LUM_Code[still], collapse = ", ")))
  if (any(.blank(sub))) stop("GLOBIOM_subclass came out blank on some rows")

  m[, GLOBIOM_subclass := sub]
  m[]
}

if (sys.nframe() == 0L) {
  stopifnot(file.exists(MAP_FILE))
  m0 <- fread(MAP_FILE)
  bk <- sub("\\.csv$", sprintf("_backup_%s.csv", format(Sys.Date())), MAP_FILE)
  if (!file.exists(bk)) { fwrite(m0, bk); message("backup written: ", bk) }

  m1 <- build_globiom_subclass(m0)
  cat("\n--- rows where GLOBIOM_subclass DIFFERS from the old composed class ---\n")
  old <- ifelse(trimws(as.character(m0$GLOBIOM_mngmt)) %in% c("", "NA"),
                trimws(as.character(m0$GLOBIOM_UNFCCC)),
                paste(trimws(as.character(m0$GLOBIOM_UNFCCC)),
                      trimws(as.character(m0$GLOBIOM_mngmt)), sep = "_"))
  # NB: build the frame FIRST, then filter. Subsetting `old` inside the data.table `j`
  # recycles it against the full table and silently mispairs was/now.
  d <- data.table(LUM_Code = m1$LUM_Code, LUM_label = substr(m1$LUM_label, 1, 34),
                  was = old, now = m1$GLOBIOM_subclass)[was != now]
  print(d)
  cat("\nchanged rows:", nrow(d), "of", nrow(m1), "\n")

  fwrite(m1, MAP_FILE)
  message("GLOBIOM_subclass written to ", MAP_FILE)
}
