# =============================================================================
# Eurostat livestock SHARE preparation for the subtype-allocation cascade.
# See docs/livestock_subtype_allocation_spec.md.
#
# Produces, for the model's NUTS0 countries x timesteps:
#   - production D/O/F shares per species (BOV clean; SGT breeders/followers + D/O hook)
#   - NUTS2 dairy-cow share (the one sub-national refinement)
#   - organic shares (2000=0, 2010=2012*backcast, 2018/2020 direct)
#   - a COVERAGE report: which country x year x attribute is observed vs fallback.
# Inputs: input/eurostat/*.tsv (downloaded via the dissemination API).
# Outputs: output/eurostat/{production_shares_nuts0,dairy_share_nuts2,organic_shares_nuts0,
#          coverage_report}.csv
# =============================================================================
suppressMessages(library(data.table))

# DATA ROOT. gamble-core holds the MODEL CODE; the bulk inputs live in cascadinggamble. GAMBLE_INPUT_DIR
# points at wherever they are, defaulting to the in-repo `input/` so an existing self-contained
# checkout keeps working unchanged.
INPUT_DIR    <- Sys.getenv("GAMBLE_INPUT_DIR", "input")
EUROSTAT_DIR <- Sys.getenv("EUROSTAT_DIR", file.path(INPUT_DIR, "eurostat"))
OUT_DIR      <- "output/eurostat"; dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
MODEL_YEARS  <- c(2000, 2010, 2020)             # outcome timesteps (2018/2020 -> 2020 for population)
ORG_BACKCAST_2010 <- 0.85                       # 2010 organic = 0.85 * 2012 share (≈2yr @ ~8%/yr)
                                                #   override with the EU organic index if available.
DAT_ADMIN    <- Sys.getenv("DRIVER_DAT_ADMIN",  # model geo (NUTS0 + NUTS3) for the coverage check
                 "output/dat_admin_FULL_Admin_Livestock_test_GLOBIOM_2026-06-25.rds")

# --- Eurostat animal-code -> production-category map -------------------------
# BOV (apro_mt_lscatl): dairy = A2300F, suckler/meat = A2300G, total = A2000, followers = residual.
BOV_TOTAL <- "A2000"; BOV_D <- "A2300F"; BOV_O <- "A2300G"
# SGT dairy/meat split IS in the population tables:
#   sheep: A4110KC = MILK (dairy) ewes, A4110KD = non-milk (meat) ewes, A4110K = all ewes (fallback);
#          A4120 = other sheep (followers); A4100 = sheep total.
#   goat:  A4210 = she-goats total, else A4210KA + A4210KB; EU she-goats are ~all dairy -> SGTD.
#          A4220 = other goats (followers); A4200 = goat total.
SGT_EWE_MILK <- "A4110KC"; SGT_EWE_MEAT <- "A4110KD"; SGT_EWES_ALL <- "A4110K"
SGT_OTHER_SHEEP <- "A4120"; SGT_SHEEP_TOT <- "A4100"
SGT_SHEGOATS <- "A4210"; SGT_SHEGOATS_A <- "A4210KA"; SGT_SHEGOATS_B <- "A4210KB"
SGT_OTHER_GOATS <- "A4220"; SGT_GOAT_TOT <- "A4200"
# organic (org_lstspec): bovines A2000, sheep A4100, goats A4200.
ORG_BOV <- "A2000"; ORG_SHEEP <- "A4100"; ORG_GOAT <- "A4200"

# ---------------------------------------------------------------------------
# 1. Eurostat TSV reader -> tidy long (dims..., year, value)
# ---------------------------------------------------------------------------
read_estat <- function(name) {
  path <- file.path(EUROSTAT_DIR, paste0(name, ".tsv"))
  raw  <- fread(path, sep = "\t", header = TRUE, colClasses = "character")
  key  <- names(raw)[1]                                   # e.g. "freq,month,animals,unit,geo\\TIME_PERIOD"
  dims <- strsplit(sub("\\\\.*", "", key), ",")[[1]]
  parts <- tstrsplit(raw[[1]], ",", fixed = TRUE)
  for (i in seq_along(dims)) raw[[dims[i]]] <- trimws(parts[[i]])
  raw[[key]] <- NULL
  yrs  <- grep("^[12][0-9]{3}$", names(raw), value = TRUE)
  long <- melt(raw, id.vars = dims, measure.vars = yrs, variable.name = "year", value.name = "v")
  long[, year := as.integer(as.character(year))]
  long[, value := suppressWarnings(as.numeric(sub("[[:space:]].*$", "", trimws(v))))]  # strip ":" + flag letters
  long[v %in% c(":", ""), value := NA_real_]
  long[, v := NULL]
  long[]
}

# collapse the month dim PER CELL: prefer the Nov/Dec annual survey, fall back to May/June
# where Nov/Dec is absent (small producers report only one round, and it varies by year/country --
# picking a single month for the whole dataset silently drops them).
collapse_month <- function(dt) {
  if (!"month" %in% names(dt)) return(dt)
  pref <- c("M11_M12", "M12", "M00", "M11", "M05_M06", "M06")
  dt <- dt[month %in% pref]
  dt[, mrank := match(month, pref)]
  setorder(dt, mrank)
  idv <- setdiff(names(dt), c("month", "value", "mrank"))
  dt[, .(value = { vv <- value[!is.na(value)]; if (length(vv)) vv[1] else NA_real_ }), by = idv]
}

# ---------------------------------------------------------------------------
# 2. Load the data
# ---------------------------------------------------------------------------
cat(">>> Reading Eurostat extracts ...\n")
cattle <- collapse_month(read_estat("apro_mt_lscatl"))
sheep  <- collapse_month(read_estat("apro_mt_lssheep"))
goat   <- collapse_month(read_estat("apro_mt_lsgoat"))
org    <- read_estat("org_lstspec")
reg    <- read_estat("agr_r_animal")                       # NUTS2 regional

# unit normalisation -> HEADS. apro_mt = THS_HD (thousand head); org = HD/NR.
to_heads <- function(dt) { dt <- copy(dt)
  dt[unit == "THS_HD", value := value * 1000]; dt[unit == "NR", value := value]; dt }
cattle <- to_heads(cattle); sheep <- to_heads(sheep); goat <- to_heads(goat); org <- to_heads(org); reg <- to_heads(reg)

getv <- function(dt, an, yr) dt[animals == an & year == yr, .(geo, value)]   # geo->value at year yr

# ---------------------------------------------------------------------------
# 3. Production D/O/F shares (NUTS0) per timestep
# ---------------------------------------------------------------------------
cat(">>> Building production D/O/F shares ...\n")
prod_rows <- list()
for (yr in MODEL_YEARS) {
  # BOV
  tot <- getv(cattle, BOV_TOTAL, yr); dd <- getv(cattle, BOV_D, yr); oo <- getv(cattle, BOV_O, yr)
  b <- Reduce(function(a, b) merge(a, b, by = "geo", all = TRUE),
              list(setnames(copy(tot), "value", "tot"), setnames(copy(dd), "value", "d"), setnames(copy(oo), "value", "o")))
  b[, f := pmax(tot - d - o, 0)]
  b[, `:=`(BOVD = d / tot, BOVO = o / tot, BOVF = f / tot)]
  prod_rows[[paste0("BOV", yr)]] <- b[, .(species = "BOV", geo, year = yr, D = BOVD, O = BOVO, F = BOVF,
                                          src = ifelse(is.na(d) | is.na(tot), "MISSING", "apro_mt"))]
  # SGT dairy/meat/followers from the milk-ewe + she-goat breakdown:
  #   SGTD = milk ewes (A4110KC) + she-goats (dairy);  SGTO = non-milk ewes;
  #   SGTF = other sheep (A4120) + other goats (A4220);  total = sheep + goats.
  em <- getv(sheep, SGT_EWE_MILK, yr); ed <- getv(sheep, SGT_EWE_MEAT, yr); ek <- getv(sheep, SGT_EWES_ALL, yr)
  osh <- getv(sheep, SGT_OTHER_SHEEP, yr); ts <- getv(sheep, SGT_SHEEP_TOT, yr)
  sga <- getv(goat, SGT_SHEGOATS_A, yr); sgb <- getv(goat, SGT_SHEGOATS_B, yr); sg0 <- getv(goat, SGT_SHEGOATS, yr)
  og <- getv(goat, SGT_OTHER_GOATS, yr); tg <- getv(goat, SGT_GOAT_TOT, yr)
  .lst <- Map(function(x, nm) setnames(copy(x), "value", nm),
              list(em, ed, ek, osh, ts, sga, sgb, sg0, og, tg),
              c("em", "ed", "ek", "osh", "ts", "sga", "sgb", "sg0", "og", "tg"))
  s <- Reduce(function(a, b) merge(a, b, by = "geo", all = TRUE), .lst)
  for (.c in c("em", "ed", "ek", "osh", "sga", "sgb", "sg0", "og")) s[is.na(get(.c)), (.c) := 0]
  s[, shegoat := fifelse(sg0 > 0, sg0, sga + sgb)]                       # she-goats total (dairy proxy)
  s[, milk_ewes := em]                                                  # dairy ewes (0 if unreported)
  s[, meat_ewes := fifelse(ed > 0, ed, pmax(ek - em, 0))]               # fallback: total ewes - milk ewes
  s[, tot := rowSums(cbind(ts, tg), na.rm = TRUE)]
  s[, `:=`(SGTD = (milk_ewes + shegoat) / tot,
           SGTO = meat_ewes / tot,
           SGTF = pmax(tot - (milk_ewes + shegoat) - meat_ewes, 0) / tot)]   # followers = residual
  prod_rows[[paste0("SGT", yr)]] <- s[, .(species = "SGT", geo, year = yr, D = SGTD, O = SGTO, F = SGTF,
                                          src = ifelse(tot > 0, "apro_mt", "MISSING"))]
}
production <- rbindlist(prod_rows)

# ---------------------------------------------------------------------------
# 4. NUTS2 dairy-cow share (refinement) from agr_r_animal
# ---------------------------------------------------------------------------
cat(">>> Building NUTS2 dairy share ...\n")
dairy_nuts2 <- rbindlist(lapply(MODEL_YEARS, function(yr) {
  tt <- reg[animals == BOV_TOTAL & year == yr, .(geo, tot = value)]
  dd <- reg[animals == BOV_D     & year == yr, .(geo, d = value)]
  m  <- merge(tt, dd, by = "geo", all = TRUE)[nchar(geo) == 4]          # NUTS2 = 4-char codes
  m[, dairy_share := d / tot][, year := yr][, .(nuts2 = geo, year, dairy_share)]
}))

# ---------------------------------------------------------------------------
# 5. Organic shares (NUTS0) with the temporal rule
# ---------------------------------------------------------------------------
cat(">>> Building organic shares (2000=0, 2010=0.85*2012, 2018/2020 direct) ...\n")
org_share_year <- function(yr_org, yr_tot) {                # organic count / total count, same heads
  ob <- org[animals == ORG_BOV   & year == yr_org, .(geo, ob = value)]
  os <- org[animals == ORG_SHEEP & year == yr_org, .(geo, os = value)]
  og <- org[animals == ORG_GOAT  & year == yr_org, .(geo, og = value)]
  tb <- cattle[animals == BOV_TOTAL & year == yr_tot, .(geo, tb = value)]
  ts <- sheep [animals == SGT_SHEEP_TOT & year == yr_tot, .(geo, ts = value)]
  tg <- goat  [animals == SGT_GOAT_TOT  & year == yr_tot, .(geo, tg = value)]
  m <- Reduce(function(a, b) merge(a, b, by = "geo", all = TRUE), list(ob, os, og, tb, ts, tg))
  m[, BOV_org := ob / tb]
  m[, SGT_org := (os + og) / (ts + tg)]
  m[, .(geo, BOV_org, SGT_org)]
}
org_rows <- list()
for (yr in MODEL_YEARS) {
  if (yr <= 2000) {                                          # 2000 -> fully conventional
    geos <- unique(cattle[year == yr, geo]); m <- data.table(geo = geos, BOV_org = 0, SGT_org = 0, src = "assumed_0")
  } else if (yr <= 2010) {                                   # 2010 -> 2012 backcast
    m <- org_share_year(2012, yr); m[, `:=`(BOV_org = BOV_org * ORG_BACKCAST_2010,
                                            SGT_org = SGT_org * ORG_BACKCAST_2010, src = "2012_backcast")]
  } else {                                                   # direct (nearest available)
    yo <- if (yr %in% org$year) yr else max(org$year[org$year <= yr])
    m <- org_share_year(yo, yr); m[, src := paste0("org_", yo)]
  }
  org_rows[[as.character(yr)]] <- m[, year := yr][]
}
organic <- rbindlist(org_rows, fill = TRUE)
organic[is.na(BOV_org), BOV_org := 0][is.na(SGT_org), SGT_org := 0]   # missing organic -> 0 (conventional)

# ---------------------------------------------------------------------------
# 6. Coverage check against the MODEL domain (NUTS0 x year x attribute)
# ---------------------------------------------------------------------------
cat(">>> Coverage check against model domain ...\n")
stopifnot(file.exists(DAT_ADMIN))
da <- as.data.table(readRDS(DAT_ADMIN))
model_nuts0 <- sort(unique(da$NUTS0)); model_nuts3 <- sort(unique(da$RESOLUTION))
cat(sprintf("    model domain: %d NUTS0 countries, %d NUTS3 units\n", length(model_nuts0), length(model_nuts3)))

cov <- CJ(geo = model_nuts0, year = MODEL_YEARS, species = c("BOV", "SGT"))
cov <- merge(cov, production[, .(geo, year, species, prod_src = src,
                                 prod_ok = is.finite(D) & is.finite(O) & is.finite(F))],
             by = c("geo", "year", "species"), all.x = TRUE)
org_long <- melt(organic, id.vars = c("geo", "year", "src"), measure.vars = c("BOV_org", "SGT_org"),
                 variable.name = "species", value.name = "org_val")
org_long[, species := ifelse(species == "BOV_org", "BOV", "SGT")]
cov <- merge(cov, org_long[, .(geo, year, species, org_src = src, org_ok = is.finite(org_val))],
             by = c("geo", "year", "species"), all.x = TRUE)
cov[is.na(prod_ok), prod_ok := FALSE][is.na(org_ok), org_ok := FALSE]
cov[is.na(prod_src), prod_src := "MISSING"][is.na(org_src), org_src := "MISSING"]
# NUTS2 dairy coverage: does each country have ANY NUTS2 dairy rows?
n2_have <- unique(substr(dairy_nuts2[is.finite(dairy_share), nuts2], 1, 2))
cov[, nuts2_dairy := substr(geo, 1, 2) %in% n2_have]

# ---------------------------------------------------------------------------
# 7. Write outputs + console summary
# ---------------------------------------------------------------------------
fwrite(production, file.path(OUT_DIR, "production_shares_nuts0.csv"))
fwrite(dairy_nuts2, file.path(OUT_DIR, "dairy_share_nuts2.csv"))
fwrite(organic, file.path(OUT_DIR, "organic_shares_nuts0.csv"))
fwrite(cov, file.path(OUT_DIR, "coverage_report.csv"))

cat("\n==================== COVERAGE SUMMARY ====================\n")
cat(sprintf("production shares present: %d / %d country*year*species cells\n", sum(cov$prod_ok), nrow(cov)))
cat(sprintf("organic shares present:    %d / %d (post temporal-rule; 0 = treated conventional)\n",
            sum(cov$org_ok), nrow(cov)))
miss_prod <- cov[prod_ok == FALSE, .(geo, year, species)]
if (nrow(miss_prod)) { cat("\nPRODUCTION MISSING (need fallback ladder -> GLOBIOM D/O/F):\n"); print(miss_prod) }
no_n2 <- sort(unique(cov[nuts2_dairy == FALSE, geo]))
cat(sprintf("\ncountries WITHOUT NUTS2 dairy (fall back to national dairy share): %s\n",
            if (length(no_n2)) paste(no_n2, collapse = ", ") else "none"))
cat(sprintf("\nWrote 4 files to %s/\n", OUT_DIR))
