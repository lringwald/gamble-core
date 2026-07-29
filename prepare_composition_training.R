# =============================================================================
# Step 1 of the driver-conditioned composition (δ) fit.
# Builds NUTS2 training tables: multinomial D/O/F head counts (agr_r_animal) + the
# count-model drivers aggregated to NUTS2. For BOV and SGT. Year 2010 (the saved dat_admin).
# Output: output/composition/{bov,sgt}_training_nuts2.csv
# =============================================================================
suppressMessages(library(data.table))
OUT <- "output/composition"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
YR <- 2010

read_estat <- function(name){ raw<-fread(file.path("input/eurostat",paste0(name,".tsv")),sep="\t",header=TRUE,colClasses="character")
 key<-names(raw)[1]; dims<-strsplit(sub("\\\\.*","",key),",")[[1]]; parts<-tstrsplit(raw[[1]],",",fixed=TRUE)
 for(i in seq_along(dims)) raw[[dims[i]]]<-trimws(parts[[i]]); raw[[key]]<-NULL
 yrs<-grep("^[12][0-9]{3}$",names(raw),value=TRUE); long<-melt(raw,id.vars=dims,measure.vars=yrs,variable.name="year",value.name="v")
 long[,year:=as.integer(as.character(year))]; long[,value:=suppressWarnings(as.numeric(sub("[[:space:]].*$","",trimws(v))))]
 long[v %in% c(":",""),value:=NA_real_]; long[,v:=NULL]; long[] }

reg <- read_estat("agr_r_animal")
gv2 <- function(an) { x <- reg[animals==an & year==YR & nchar(geo)==4, .(nuts2=geo, value)]; x }  # NUTS2 = 4-char

# ---- responses (NUTS2 head counts, THS_HD; ratios so unit cancels) ----
# BOV: D=dairy cows A2300F, O=suckler A2300G, F=total A2000 - cows
bov <- Reduce(function(a,b) merge(a,b,by="nuts2",all=TRUE), list(
  setnames(gv2("A2000"),"value","tot"), setnames(gv2("A2300F"),"value","nD"), setnames(gv2("A2300G"),"value","nO")))
bov[, nF := pmax(tot - nD - nO, 0)]
bov <- bov[is.finite(nD) & is.finite(nO) & is.finite(tot) & tot > 0, .(nuts2, nD, nO, nF)]

# SGT: D=milk ewes A4110KC + she-goats(A4210K|KA+KB, dairy), O=non-milk ewes A4110KD,
#      F = (sheep A4100 + goat A4200) - breeding females
sgt0 <- Reduce(function(a,b) merge(a,b,by="nuts2",all=TRUE), list(
  setnames(gv2("A4100"),"value","tsh"), setnames(gv2("A4200"),"value","tgo"),
  setnames(gv2("A4110KC"),"value","emilk"), setnames(gv2("A4110KD"),"value","emeat"),
  setnames(gv2("A4210K"),"value","sg"), setnames(gv2("A4210KA"),"value","sga"), setnames(gv2("A4210KB"),"value","sgb")))
for (c in c("tsh","tgo","emilk","emeat","sg","sga","sgb")) sgt0[is.na(get(c)), (c):=0]
sgt0[, shegoat := fifelse(sg>0, sg, sga+sgb)]
sgt0[, tot := tsh + tgo]
sgt0[, `:=`(nD = emilk + shegoat, nO = emeat, nF = pmax(tot - emilk - shegoat - emeat, 0))]
sgt <- sgt0[tot > 0, .(nuts2, nD, nO, nF)]

# ---- drivers: aggregate dat_admin (NUTS3) -> NUTS2 ----
# auto-detect the most recent dat_admin_FULL the count model saved (timestamped with the run date),
# or override with DRIVER_DAT_ADMIN. Step 1 (run_prior_module_count_model.R) writes this file.
.daf <- Sys.getenv("DRIVER_DAT_ADMIN", "")
if (.daf == "") {
  .cands <- list.files("output", pattern = "^dat_admin_FULL_.*\\.rds$", full.names = TRUE)
  if (!length(.cands)) stop("No dat_admin_FULL_*.rds in output/ — run run_prior_module_count_model.R first.")
  .daf <- .cands[which.max(file.mtime(.cands))]
}
cat(sprintf(">>> drivers from: %s\n", .daf))
d <- as.data.table(readRDS(.daf))
d[, nuts2 := substr(RESOLUTION, 1, 4)]
# EXTENSIVE (sum -> log1p, EXACTLY like the totals/grid X): LU areas + allPA_area + GDP + Pop.
area_cols <- c(setdiff(grep("^lu_area_", names(d), value=TRUE), "lu_area_no_choice"),
               intersect(c("allPA_area", "GDP", "Pop"), names(d)))
# INTENSIVE (weighted mean): terrain + climate + heterogeneity + shares.
meanv <- intersect(c("Slope_rad","Elevation","Aspect_cos_mean","Aspect_sin_mean","GHM_HI","CISI",
                     "Growing_Degree_Days_gdd5","Precipitation_Seasonality_bio15","Annual_Precipitation_bio12"), names(d))
meanv <- c(meanv, grep("_sd$", names(d), value=TRUE),
           intersect(c("flat_share","steep_share","lowland_share","upland_share"), names(d)))
d[is.na(total_area_km2), total_area_km2 := 0]
drv <- d[, c(
  lapply(.SD[, area_cols, with=FALSE], function(x) log1p(sum(x, na.rm=TRUE))),                            # areas/GDP/Pop: sum->log1p
  lapply(.SD[, meanv, with=FALSE], function(x) weighted.mean(x, total_area_km2, na.rm=TRUE)),             # rest: area-wmean
  country = .(NUTS0[1])), by = nuts2]
lu <- area_cols   # for the write_train driver count below

write_train <- function(resp, tag) {
  m <- merge(resp, drv, by="nuts2")
  m <- m[is.finite(nD+nO+nF) & (nD+nO+nF) > 0]
  fwrite(m, file.path(OUT, paste0(tag, "_training_nuts2.csv")))
  cat(sprintf("%s: %d NUTS2 regions x %d drivers | mean shares D=%.2f O=%.2f F=%.2f\n",
              toupper(tag), nrow(m), length(c(area_cols,meanv)),
              mean(m$nD/(m$nD+m$nO+m$nF)), mean(m$nO/(m$nD+m$nO+m$nF)), mean(m$nF/(m$nD+m$nO+m$nF))))
}
write_train(bov, "bov")
write_train(sgt, "sgt")
cat(sprintf("drivers (%d): %s\n", length(c(area_cols,meanv)), paste(c(area_cols,meanv), collapse=", ")))
