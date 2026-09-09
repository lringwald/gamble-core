# =============================================================================
# target_class_split.R — map model classes onto a PROJECT's target classification
# =============================================================================
# Generic mechanism; the project supplies the rules and the share table. Two operations, applied
# as an ORDERED CASCADE so a split can feed a later split (HRL "Fruits" -> citrus + other_fruit at
# NUTS2, then other_fruit -> apples + other_fruit at NUTS0):
#
#   rename : from_class -> to_class                       (1:1, no data needed)
#   split  : from_class -> group, using shares[group, geo] (area allocated by share)
#
# `shares` is long: group | geo | source | target | share, with share summing to 1 within
# (group, geo). `geo_lookup` maps join_id -> the geo code the shares are keyed on; a pixel whose
# geo is missing falls back to its country, then to the share table's own fallback rows.
#
# AREA IS CONSERVED at every step: sum(area) after == before, asserted.
# =============================================================================
suppressMessages(library(data.table))

# rules: data.table(order, action ["rename"|"split"], from_class, to)
#   rename -> `to` is the new class name
#   split  -> `to` is the group name in `shares`
apply_target_classification <- function(dt, rules, shares = NULL, geo_lookup = NULL,
                                        verbose = TRUE, tol = 1e-6) {
  stopifnot(is.data.table(dt), all(c("join_id","model_class","area") %in% names(dt)))
  rules <- as.data.table(rules)[order(order)]
  if (any(rules$action == "split")) {
    if (is.null(shares) || is.null(geo_lookup))
      stop("apply_target_classification: split rules need both `shares` and `geo_lookup`")
    shares <- as.data.table(shares); geo_lookup <- as.data.table(geo_lookup)
    stopifnot(all(c("group","geo","target","share") %in% names(shares)),
              all(c("join_id","geo") %in% names(geo_lookup)))
    setkey(geo_lookup, join_id)
  }
  A0 <- dt[, sum(area, na.rm = TRUE)]
  for (i in seq_len(nrow(rules))) {
    r <- rules[i]
    if (!r$from_class %in% dt$model_class) {
      if (verbose) message(sprintf("  [skip] %-28s absent from the design", r$from_class)); next }
    if (r$action == "rename") {
      n <- dt[model_class == r$from_class, .N]
      dt[model_class == r$from_class, `:=`(model_class = r$to, focal_class = r$to)]
      if (verbose) message(sprintf("  rename %-26s -> %-32s (%d rows)", r$from_class, r$to, n))
      next
    }
    sh <- shares[group == r$to]
    if (!nrow(sh)) stop("no shares for group '", r$to, "'")
    src <- dt[model_class == r$from_class]
    keep <- dt[model_class != r$from_class]
    g <- geo_lookup[.(src$join_id), geo, on = "join_id"]
    src[, geo := g]
    # geo fallback: exact -> country (first 2 chars) -> the group's area-weighted mean
    have <- unique(sh$geo)
    src[, geo_use := fifelse(geo %in% have, geo,
                      fifelse(substr(geo, 1, 2) %in% have, substr(geo, 1, 2), NA_character_))]
    gmean <- sh[, .(share = mean(share)), by = target][, share := share / sum(share)][]
    ok  <- src[!is.na(geo_use)]; bad <- src[is.na(geo_use)]
    out <- list()
    if (nrow(ok))
      out[[1]] <- merge(ok[, .(join_id, geo_use, area)], sh[, .(geo, target, share)],
                        by.x = "geo_use", by.y = "geo", allow.cartesian = TRUE)[
                        , .(join_id, model_class = target, focal_class = target, area = area * share)]
    if (nrow(bad)) {                       # no geo match at all -> the group's mean composition
      b <- bad[rep(seq_len(.N), each = nrow(gmean)), .(join_id, area)]
      b[, `:=`(target = rep(gmean$target, nrow(bad)), share = rep(gmean$share, nrow(bad)))]
      out[[2]] <- b[, .(join_id, model_class = target, focal_class = target, area = area * share)]
    }
    add <- rbindlist(out, fill = TRUE)
    if (verbose) message(sprintf("  split  %-26s -> %-2d class(es) via '%s'  [%d px exact, %d fallback]",
        r$from_class, uniqueN(sh$target), r$to, nrow(ok), nrow(bad)))
    dt <- rbind(keep, add, fill = TRUE)
  }
  A1 <- dt[, sum(area, na.rm = TRUE)]
  if (abs(A1 - A0) > tol * max(1, A0))
    stop(sprintf("apply_target_classification: area not conserved (%.6g -> %.6g)", A0, A1))
  if (verbose) message(sprintf("  area conserved: %.6g", A1))
  dt[]
}
