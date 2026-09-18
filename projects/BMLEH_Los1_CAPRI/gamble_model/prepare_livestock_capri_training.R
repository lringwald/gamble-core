#!/usr/bin/env Rscript
# =============================================================================
# BMLEH_Los1_CAPRI — NUTS2 cattle training table in the CAPRI classification
# =============================================================================
# Replaces the 3-category D/O/F composition target with 8 CAPRI-aligned categories. Eurostat's
# regional table (agr_r_animal) carries the full cattle detail at NUTS2, so this is estimated rather
# than apportioned:
#
#   DCOW  A2300F                 dairy cows            -> CAPRI DCOH + DCOL   (no yield split)
#   SCOW  A2300G                 non-dairy cows        -> CAPRI SCOW          EXACT
#   HEIR  A2230C + A2220C        heifers, not for slaughter -> CAPRI HEIR     EXACT
#   HEIF  A2230B + A2220B        heifers, for slaughter -> CAPRI HEIH + HEIL  (no weight split)
#   BULL  A2120 + A2130          male bovines >= 1y    -> CAPRI BULH + BULL   (no weight split)
#   CAMR  A2110C                 male calves, not for slaughter -> CAPRI CAMR EXACT
#   CAFR  A2210C                 female calves, not for slaughter -> CAPRI CAFR EXACT
#   CAFF  A2010B                 bovines < 1y for slaughter -> CAPRI CAMF + CAFF (no sex split)
#
# The four unsplit pairs are left unsplit by decision: yield, weight and the fattening-calf sex split
# are not in any regional source, and inventing them would put structure in the prior that no data
# supports. Four of the twelve CAPRI cattle activities come out exactly; the rest come out as pairs.
#
# The partition CLOSES against the A2000 total -- median ratio 1.0000, 10th-90th 0.9992-1.0006,
# 96.9% of regions within 2% -- which is what makes the multinomial well-posed.
#
#   Rscript projects/BMLEH_Los1_CAPRI/gamble_model/prepare_livestock_capri_training.R
# =============================================================================
suppressMessages(library(data.table))
EU   <- Sys.getenv("EUROSTAT_DIR", file.path(Sys.getenv("GAMBLE_INPUT_DIR", "input"), "eurostat"))
YEAR <- as.integer(Sys.getenv("LS_YEAR", "2020"))
DRV  <- Sys.getenv("LS_DRIVERS", "output/composition/bov_training_nuts2.csv")
OUT  <- Sys.getenv("LS_OUT", "output/composition/bov_training_nuts2_capri.csv")
MINTOT <- as.numeric(Sys.getenv("LS_MIN_TOTAL", "1"))

read_estat <- function(name) {
  raw <- fread(file.path(EU, paste0(name, ".tsv")), sep = "\t", header = TRUE, colClasses = "character")
  key <- names(raw)[1]; dims <- strsplit(sub("\\\\.*", "", key), ",")[[1]]
  parts <- tstrsplit(raw[[1]], ",", fixed = TRUE)
  for (i in seq_along(dims)) raw[[dims[i]]] <- trimws(parts[[i]])
  raw[[key]] <- NULL
  yrs <- grep("^[12][0-9]{3}$", names(raw), value = TRUE)
  long <- melt(raw, id.vars = dims, measure.vars = yrs, variable.name = "year", value.name = "v")
  long[, year := as.integer(as.character(year))]
  long[, value := suppressWarnings(as.numeric(sub("[[:space:]].*$", "", trimws(v))))]
  long[!is.na(value)]
}
d <- read_estat("agr_r_animal")[year == YEAR & nchar(geo) == 4L, .(geo, animals, value)]
w <- dcast(d, geo ~ animals, value.var = "value", fill = NA_real_)
g  <- function(cd) if (cd %in% names(w)) w[[cd]] else rep(NA_real_, nrow(w))
s2 <- function(a, b) { x <- rowSums(cbind(g(a), g(b)), na.rm = TRUE); x[is.na(g(a)) & is.na(g(b))] <- NA_real_; x }
w[, `:=`(nDCOW = g("A2300F"), nSCOW = g("A2300G"),
         nHEIR = s2("A2230C","A2220C"), nHEIF = s2("A2230B","A2220B"),
         nBULL = s2("A2120","A2130"),
         nCAMR = g("A2110C"), nCAFR = g("A2210C"), nCAFF = g("A2010B"))]
CATS <- c("nDCOW","nSCOW","nHEIR","nHEIF","nBULL","nCAMR","nCAFR","nCAFF")
w[, tot_parts := rowSums(.SD, na.rm = TRUE), .SDcols = CATS]
w[, tot_eurostat := g("A2000")]
cat(sprintf("agr_r_animal %d: %d NUTS2 regions\n", YEAR, nrow(w)))

keep <- w[complete.cases(w[, ..CATS]) & tot_parts >= MINTOT]
cat(sprintf("  complete on all 8 categories and >= %g head: %d regions\n", MINTOT, nrow(keep)))
chk <- keep[!is.na(tot_eurostat) & tot_eurostat > 0]
cat(sprintf("  partition closes: median %.4f | within 2%%: %.1f%% (%d checked)\n",
            median(chk$tot_parts/chk$tot_eurostat),
            100*mean(abs(chk$tot_parts/chk$tot_eurostat - 1) < 0.02), nrow(chk)))

# Reuse the driver columns already assembled for the D/O/F training table, so the only thing changing
# is the OUTCOME. Anything the driver table does not cover cannot be fitted and is reported.
drv <- fread(DRV)
dcols <- setdiff(names(drv), c("nD","nO","nF"))
out <- merge(keep[, c("geo", CATS), with = FALSE], drv[, ..dcols], by.x = "geo", by.y = "nuts2")
setnames(out, "geo", "nuts2")
lost <- setdiff(keep$geo, drv$nuts2)
cat(sprintf("  joined to drivers: %d regions (%d dropped for having no driver row: %s)\n",
            nrow(out), length(lost), paste(head(lost, 8), collapse = ", ")))
fwrite(out, OUT)
cat(sprintf("\nwrote %s\n  %d regions x %d categories (was %d regions x 3)\n",
            OUT, nrow(out), length(CATS), nrow(drv)))
cat("\ncategory shares, EU-wide:\n")
tt <- colSums(out[, ..CATS]); print(round(100*tt/sum(tt), 2))
