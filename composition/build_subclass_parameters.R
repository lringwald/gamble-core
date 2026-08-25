# =============================================================================
# Build the FINAL plug-in parameter table: per (country, subclass, system, driver) coefficient gamma.
#   gamma = beta_total_country (totals per-country slope = the LEVEL)
#         + delta_DOF          (D/O/F composition contrast, NUTS2 Eurostat  -> the subtype SPLIT)
#         + delta_org_signed   (organic vs conventional contrast, NUTS3 map -> the SYSTEM SPLIT)
# 12 cells per species = {D,O,F} x {conventional, organic}.  Predict: cell weight_i proportional to
# exp(sum_driver gamma*X_i), renormalize within NUTS3 (species/D-O-F baseline cancels; the ORGANIC
# level does NOT cancel -> it is anchored per country to the Eurostat national organic share).
# Output: output/composition/subclass_country_parameters.csv
#   columns: species, country, subclass, dof, system, driver, beta_total, delta_dof, delta_org, gamma
# =============================================================================
suppressMessages(library(data.table))
OUT <- "output/composition/subclass_country_parameters.csv"
ORG_YEAR <- as.integer(Sys.getenv("ORG_ANCHOR_YEAR", "2010"))   # Eurostat organic-share year to anchor to

# country order = droplevels(GLOB_country) levels, as the count model used for the REs
daf <- list.files("output","^dat_admin_FULL_.*rds$",full.names=TRUE)
d <- as.data.table(readRDS(daf[which.max(file.mtime(daf))]))
countries <- levels(droplevels(as.factor(d[["GLOB_country"]])))
iso2name  <- unique(d[, .(NUTS0 = as.character(NUTS0), name = as.character(GLOB_country))])   # NUTS0 -> country name
name_by_iso <- setNames(iso2name$name, iso2name$NUTS0)

# map composition/organic driver name -> totals (X_mat) name: lu_area_X -> log1p_lu_area_X; GDP/Pop/allPA_area -> log1p_*
to_tot <- function(nm) { nm <- ifelse(grepl("^lu_area_", nm), paste0("log1p_", nm), nm)
  ifelse(nm %in% c("GDP","Pop","allPA_area"), paste0("log1p_", nm), nm) }

# ---- ORGANIC (optional): fitted pattern (dOrg) + Eurostat national levels + per-country anchor ----
of <- "output/composition/organic_fit.rds"; otr <- "output/composition/organic_training_nuts3.csv"
osh <- "output/eurostat/organic_shares_nuts0.csv"
HAVE_ORG <- all(file.exists(of, otr, osh))
if (HAVE_ORG) {
  ofit  <- readRDS(of); dorg_raw <- setNames(rowMeans(ofit$dOrg), ofit$drivers)   # over c("intercept", sel), raw scale
  otrain <- fread(otr); otrain[, name := name_by_iso[country]]
  esh <- fread(osh)[year == ORG_YEAR, .(NUTS0 = geo, BOV = as.numeric(BOV_org), SGT = as.numeric(SGT_org))]
  esh[, name := name_by_iso[NUTS0]]
  eur_target <- function(sp, nm) { v <- esh[name == nm][[sp]]; if (!length(v) || all(is.na(v))) NA_real_ else mean(v, na.rm=TRUE) }
  # per-country anchor c: weighted-mean plogis(2*dOrg.X + c) over the country's NUTS3 == Eurostat share.
  #   pattern eta_i = 2 * [1, X_i[sel]] %*% dOrg ; weight w_i = grassland area (n_org+n_conv).
  Xcols <- intersect(names(dorg_raw), c("intercept", names(otrain)))               # intercept + sel drivers present
  anchor_c <- function(sp, nm) {
    sub <- otrain[name == nm]; tgt <- eur_target(sp, nm)
    if (!nrow(sub) || is.na(tgt)) return(NA_real_)
    tgt <- min(max(tgt, 1e-4), 1 - 1e-4)
    Xi <- cbind(intercept = 1, as.matrix(sub[, setdiff(Xcols, "intercept"), with = FALSE]))
    eta <- 2 * as.numeric(Xi %*% dorg_raw[colnames(Xi)]); w <- sub$n_org + sub$n_conv
    f <- function(c) sum(w * plogis(eta + c)) / sum(w) - tgt
    if (f(-40) > 0) return(-40); if (f(40) < 0) return(40)                         # target outside reachable range
    uniroot(f, c(-40, 40))$root
  }
  dorg_tot <- setNames(numeric(0), character(0))                                   # dOrg mapped to totals driver names
  nm_tot <- to_tot(names(dorg_raw)); dorg_tot <- setNames(as.numeric(dorg_raw), nm_tot)
} else {
  cat(sprintf(">>> organic pieces missing (%s) -> emitting 6 conventional cells only.\n",
              paste(c(of, otr, osh)[!file.exists(c(of, otr, osh))], collapse = ", ")))
}

systems <- if (HAVE_ORG) c("conventional","organic") else "conventional"
rows <- list()
for (sp in c("BOV","SGT")) {
  # --- totals: per-country per-driver beta (mean over draws+chains) ---
  .dirs <- list.dirs("output/saved_model_outputs", recursive = FALSE)
  .d <- .dirs[grepl(sprintf("_%s_RE_", sp), basename(.dirs)) & file.exists(file.path(.dirs, "fit.rds"))]
  if (!length(.d)) stop(sprintf("No saved totals fit dir for %s under output/saved_model_outputs/", sp))
  tf <- file.path(.d[which.max(file.mtime(file.path(.d, "fit.rds")))], "fit.rds")
  r <- readRDS(tf); pt <- r$postb_total; dn <- dimnames(pt)[[1]]
  bt <- if (length(dim(pt))==5) apply(pt[,1,,,], c(1,2), mean) else apply(pt[,1,,], c(1,2), mean)  # [driver x country]
  rownames(bt) <- dn; ng <- ncol(bt)
  cn <- if (ng == length(countries)) countries else as.character(seq_len(ng))

  # --- D/O/F composition delta (fixed slopes), mapped to totals names ---
  cf <- readRDS(sprintf("output/composition/%s_composition_fit.rds", tolower(sp)))
  delta <- list(D = rowMeans(cf$dD), O = rowMeans(cf$dO), F = rowMeans(cf$dF))
  for (lv in names(delta)) names(delta[[lv]]) <- to_tot(cf$drivers)

  # --- organic delta mapped to totals names (0-vec if no organic) ---
  dorg_vec <- setNames(numeric(length(dn)), dn)
  if (HAVE_ORG) dorg_vec[intersect(names(dorg_tot), dn)] <- dorg_tot[intersect(names(dorg_tot), dn)]

  for (ci in seq_len(ng)) {
    cc <- anchor <- NULL
    for (lv in c("D","O","F")) {
      dvec <- setNames(numeric(length(dn)), dn)                            # D/O/F delta aligned to totals drivers
      dvec[intersect(names(delta[[lv]]), dn)] <- delta[[lv]][intersect(names(delta[[lv]]), dn)]
      for (sys in systems) {
        dorg_signed <- setNames(numeric(length(dn)), dn)
        if (sys != "conventional") {
          if (is.null(anchor)) anchor <- anchor_c(sp, cn[ci])              # solve once per (sp, country)
          c_anch <- if (is.na(anchor)) 0 else anchor
          dorg_signed <- dorg_vec; dorg_signed["intercept"] <- dorg_signed["intercept"] + c_anch/2   # organic:  +dOrg, +c/2
        } else if (HAVE_ORG) {
          if (is.null(anchor)) anchor <- anchor_c(sp, cn[ci])
          c_anch <- if (is.na(anchor)) 0 else anchor
          dorg_signed <- -dorg_vec; dorg_signed["intercept"] <- dorg_signed["intercept"] - c_anch/2   # conventional: -dOrg, -c/2
        }
        rows[[paste(sp,ci,lv,sys)]] <- data.table(
          species = sp, country = cn[ci], subclass = paste0(sp, lv), dof = lv, system = sys, driver = dn,
          beta_total = bt[, ci], delta_dof = dvec, delta_org = dorg_signed,
          gamma = bt[, ci] + dvec + dorg_signed)
      }
    }
  }
}
tab <- rbindlist(rows)
fwrite(tab, OUT)
cat(sprintf("Wrote %d rows -> %s\n", nrow(tab), OUT))
cat(sprintf("  %d countries x %d subclasses (%d D/O/F x %d system) x %d drivers\n",
            uniqueN(tab$country), uniqueN(tab[, .(subclass, system)]), 3L, length(systems), uniqueN(tab$driver)))
if (HAVE_ORG) {
  cat("\norganic level check — modelled vs Eurostat national share (should match after anchor):\n")
  chk <- rbindlist(lapply(c("BOV","SGT"), function(sp) rbindlist(lapply(unique(esh$name), function(nm) {
    sub <- otrain[name == nm]; if (!nrow(sub)) return(NULL)
    Xi <- cbind(intercept = 1, as.matrix(sub[, setdiff(Xcols, "intercept"), with = FALSE]))
    c_a <- anchor_c(sp, nm); if (is.na(c_a)) return(NULL)
    eta <- 2 * as.numeric(Xi %*% dorg_raw[colnames(Xi)]); w <- sub$n_org + sub$n_conv
    data.table(species = sp, country = nm, modelled = round(sum(w*plogis(eta+c_a))/sum(w), 3),
               eurostat = round(eur_target(sp, nm), 3)) }))))
  print(chk[!is.na(eurostat)][order(-eurostat)][1:10])
}
