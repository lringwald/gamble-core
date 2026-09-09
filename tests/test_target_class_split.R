# The BMLEH cascade on REAL Eurostat shares: composition, area conservation, geo fallback.
suppressMessages(library(data.table))
source("codes/target_class_split.R")
source("projects/BMLEH_Los1_CAPRI/gamble_model/target_rules.R")
np <- 0L; nf <- 0L
ok <- function(c_, m) { if (isTRUE(c_)) { np <<- np+1L; cat(sprintf("[PASS] %s\n", m)) }
                        else { nf <<- nf+1L; cat(sprintf("[FAIL] %s\n", m)) } }
sh <- fread("output/gamble_model/BMLEH_Los1_CAPRI/crop_shares_nuts2.csv")

# 4 pixels: an Italian region (durum country), an Austrian one, a German one, and an UNKNOWN geo
geo <- data.table(join_id = 1:4, geo = c("ITF4","AT11","DE21","ZZ99"))
dt <- rbindlist(lapply(1:4, function(i) data.table(join_id = i, focal_class = NA_character_,
  model_class = c("Wheat","Maize","Other_cereals","Fruits","Fresh_vegetables","Grapes",
                  "Cropland_arable_other","Barley","Forests_SR","Nuts","Cropland_permanent_other"),
  area = c(10,10,10,10,10,10,10,10,10,10,10))))
A0 <- dt[, sum(area)]
out <- apply_target_classification(dt, BMLEH_TARGET_RULES, sh, geo, verbose = FALSE)

ok(abs(out[, sum(area)] - A0) < 1e-9, sprintf("area conserved (%.1f -> %.1f)", A0, out[, sum(area)]))
ok(!any(grepl("^(Wheat|Maize|Other_cereals|Fruits|Grapes|Fresh_vegetables|Barley|Forests_SR|Nuts)$",
              out$model_class)), "no HRL/LUM source class survives the cascade")
ok(all(grepl("^Cropland_", unique(out$model_class))), "every output class is a BMLEH target class")
ok(abs(out[model_class == "Cropland_arable_energy", sum(area)] - 40) < 1e-9,   # 4 pixels x 10
   "Forests_SR -> energy, area intact")

w <- out[grepl("wheat", model_class), .(a = sum(area)), by = .(join_id, model_class)]
it <- w[join_id == 1]; at <- w[join_id == 2]
d_it <- it[model_class == "Cropland_arable_durumwheat", a] / it[, sum(a)]
d_at <- at[model_class == "Cropland_arable_durumwheat", a] / at[, sum(a)]
cat(sprintf("\n  durum share: ITF4 %.3f | AT11 %.3f\n", d_it, d_at))
ok(d_it > 0.9, "Puglia (ITF4) is nearly all durum -- regional share applied")
ok(d_at < 0.2, "Austria (AT11) is nearly all soft wheat")
ok(d_it > d_at, "the split VARIES by region (NUTS2 shares are live)")

# the two-stage fruit cascade: NUTS2 citrus split, then NUTS0 apple split of the remainder
fr <- out[grepl("permanent_(apples|citrus|other_fruit)", model_class), .(a = sum(area)), by = model_class]
cat(sprintf("  fruit cascade: %s\n", paste(sprintf("%s %.2f", sub("Cropland_permanent_","",fr$model_class), fr$a), collapse=" | ")))
ok(nrow(fr) == 3, "Fruits -> citrus + apples + other_fruit (both stages ran)")
ok(out[model_class == "Cropland_permanent_apples", sum(area)] > 0, "apples produced by the NUTS0 stage")
ok(out[model_class == "Cropland_arable_tomatoes", sum(area)] > 0, "tomatoes produced")
ok(out[model_class == "Cropland_permanent_wine", sum(area)] > 0, "wine grapes produced")

# unknown geo must still be allocated, not dropped
u <- out[join_id == 4, sum(area)]
ok(abs(u - 110) < 1e-9, sprintf("unknown geo ZZ99 falls back rather than dropping (%.1f of 110)", u))
ok(out[model_class == "Cropland_arable_other_industrial", sum(area)] > 0, "residual split reaches other_industrial")

cat(sprintf("\n%d distinct target classes produced\n", uniqueN(out$model_class)))
cat(sprintf("RESULT: %d/%d checks passed\n", np, np + nf))
if (nf > 0L) quit(status = 1)
