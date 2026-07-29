# =============================================================================
# Build the FINAL plug-in parameter table: per (country, subclass, driver) coefficient gamma.
# gamma = beta_total_country (the totals' per-country slope, the LEVEL) + delta_subclass (the
# composition contrast, the SPLIT). Predict: subclass weight_i ∝ exp(sum_driver gamma*X_i + offset_i),
# renormalize within NUTS3 (intercept/r/country-baseline cancel there -> slopes are what matter).
# Output: output/composition/subclass_country_parameters.csv
# =============================================================================
suppressMessages(library(data.table))
OUT <- "output/composition/subclass_country_parameters.csv"

# country order = droplevels(GLOB_country) levels, as the count model used for the REs
daf <- list.files("output","^dat_admin_FULL_.*rds$",full.names=TRUE)
d <- as.data.table(readRDS(daf[which.max(file.mtime(daf))]))
countries <- levels(droplevels(as.factor(d[["GLOB_country"]])))

# map composition driver name -> totals (X_mat) name: lu_area_X -> log1p_lu_area_X; GDP/Pop -> log1p_*
to_tot <- function(nm) { nm <- ifelse(grepl("^lu_area_", nm), paste0("log1p_", nm), nm)
  ifelse(nm %in% c("GDP","Pop","allPA_area"), paste0("log1p_", nm), nm) }   # composition raw name -> totals X name

rows <- list()
for (sp in c("BOV","SGT")) {
  # --- totals: per-country per-driver beta (mean over draws+chains) ---
  # auto-detect the saved totals fit dir (MODEL_LABEL is resolution-based now, e.g. NUTS3_Livestock_*)
  .dirs <- list.dirs("output/saved_model_outputs", recursive = FALSE)
  .d <- .dirs[grepl(sprintf("_%s_RE_", sp), basename(.dirs)) & file.exists(file.path(.dirs, "fit.rds"))]
  if (!length(.d)) stop(sprintf("No saved totals fit dir for %s under output/saved_model_outputs/", sp))
  tf <- file.path(.d[which.max(file.mtime(file.path(.d, "fit.rds")))], "fit.rds")
  r <- readRDS(tf); pt <- r$postb_total
  dn <- dimnames(pt)[[1]]
  bt <- if (length(dim(pt))==5) apply(pt[,1,,,], c(1,2), mean) else apply(pt[,1,,], c(1,2), mean)  # [driver x country]
  rownames(bt) <- dn; ng <- ncol(bt)
  cn <- if (ng == length(countries)) countries else as.character(seq_len(ng))

  # --- composition: delta per subclass (fixed slopes), mapped to totals names ---
  cf <- readRDS(sprintf("output/composition/%s_composition_fit.rds", tolower(sp)))
  cd <- cf$drivers; cd_tot <- to_tot(cd)
  delta <- list(D = rowMeans(cf$dD), O = rowMeans(cf$dO), F = rowMeans(cf$dF))
  for (lv in names(delta)) names(delta[[lv]]) <- cd_tot

  for (ci in seq_len(ng)) for (lv in c("D","O","F")) {
    dvec <- setNames(numeric(length(dn)), dn)               # delta aligned to totals drivers (0 if absent)
    dvec[intersect(names(delta[[lv]]), dn)] <- delta[[lv]][intersect(names(delta[[lv]]), dn)]
    rows[[paste(sp,ci,lv)]] <- data.table(
      species = sp, country = cn[ci], subclass = paste0(sp, lv), driver = dn,
      beta_total = bt[, ci], delta = dvec, gamma = bt[, ci] + dvec)
  }
}
tab <- rbindlist(rows)
fwrite(tab, OUT)
cat(sprintf("Wrote %d rows -> %s\n", nrow(tab), OUT))
cat(sprintf("  %d countries x 6 subclasses x %d drivers\n", uniqueN(tab$country), uniqueN(tab$driver)))
cat("\nexample (BOVD, country 1, top |gamma| drivers):\n")
print(tab[subclass=="BOVD" & country==tab$country[1]][order(-abs(gamma))][1:6, .(driver, beta_total=round(beta_total,3), delta=round(delta,3), gamma=round(gamma,3))])
