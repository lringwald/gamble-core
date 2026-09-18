# =============================================================================
# harmonize_to_targets.R — rake a pixel classification onto regional target areas
# =============================================================================
# WHY, measured: HRL types only ~62% of LUM arable area, and its detection rate is strongly spatially
# structured -- R^2 of HRL coverage on the pixel model's own covariates is 0.713, with focal_softwheat
# r=+0.57, yield_index_High +0.48, Precipitation_Seasonality -0.44, focal_olives -0.37. Everything HRL
# misses becomes residual and ~69% of the residual is fodder, so an unharmonised prior learns "fodder
# belongs in Mediterranean low-yield terrain" when it is really seeing where the sensor product fails.
# Green fodder comes out at 29.5% of arable against Eurostat's 19.7%. Country random effects cannot
# absorb this: at R^2 0.71 the detection geography is nearly a function of the drivers themselves.
#
# WHAT: per zone, move area from over-represented to under-represented classes so the zone totals hit
# the targets, then allocate those flows to pixels. Pixel TOTALS are untouched -- this redistributes
# composition within a pixel's existing area, it does not move land.
#
# The spatial step reuses cascadecore::allocate_transitions_spatial (gravity-seeded IPF), the same
# allocator the downstream start-map harmonisation uses, so the prior and the downscaler agree on how
# a zone-level flow becomes a pixel-level one. Without terra/cascade-core it falls back to
# proportional allocation, which conserves the same margins but ignores neighbourhood structure.
# =============================================================================
suppressMessages({library(data.table)})

# obs:     data.table(pixel, zone, class, area)   -- the classification to correct
# targets: data.table(zone, class, area)          -- what the zone totals should be
harmonize_to_targets <- function(obs, targets, classes = NULL, allocator = TRUE,
                                 coords = NULL, pixel_size = 10, tol = 1e-6, verbose = TRUE) {
  obs <- as.data.table(copy(obs)); targets <- as.data.table(copy(targets))
  stopifnot(all(c("pixel","zone","class","area") %in% names(obs)),
            all(c("zone","class","area") %in% names(targets)))
  if (is.null(classes)) classes <- sort(unique(targets$class))
  # Only the named block is raked; everything else passes through untouched, so a caller can correct
  # arable composition without disturbing forest or grassland.
  keep <- obs[!class %in% classes]; work <- obs[class %in% classes]
  if (!nrow(work)) { warning("no observed area in the named classes"); return(obs) }

  P0 <- work[, .(p0 = sum(area)), by = pixel]
  zo <- work[, .(obs = sum(area)), by = .(zone, class)]
  zt <- targets[class %in% classes, .(tgt = sum(area)), by = .(zone, class)]
  z  <- merge(zo, zt, by = c("zone","class"), all = TRUE)
  z[is.na(obs), obs := 0][is.na(tgt), tgt := 0]
  # Targets are SCALED to the observed zone total: the block's extent comes from LUM and is not up for
  # revision here, only its composition. Raking to raw target levels would silently change land area.
  z[, `:=`(tot_o = sum(obs), tot_t = sum(tgt)), by = zone]
  z <- z[tot_o > 0]
  z[, tgt_s := fifelse(tot_t > 0, tgt * tot_o / tot_t, obs)]
  z[, `:=`(surplus = pmax(obs - tgt_s, 0), deficit = pmax(tgt_s - obs, 0))]

  # zone-level from->to flows: each surplus class sends to deficit classes in proportion to deficit
  zz <- z[, {
    sp <- .SD[surplus > 0]; df <- .SD[deficit > 0]; D <- sum(df$deficit)
    stay <- .SD[, .(from_class = class, to_class = class, value = obs - surplus)]
    mv <- if (nrow(sp) && nrow(df) && D > 0)
      CJ(i = seq_len(nrow(sp)), j = seq_len(nrow(df)))[
        , .(from_class = sp$class[i], to_class = df$class[j], value = sp$surplus[i] * df$deficit[j] / D)]
      else NULL
    rbind(stay[value > 0], mv)[value > 0]
  }, by = zone]
  if (verbose) {
    mv <- zz[from_class != to_class, sum(value)]
    cat(sprintf("  harmonise: %d zones | %.0f km2 moved between classes (%.1f%% of the block)\n",
                uniqueN(zz$zone), mv, 100 * mv / work[, sum(area)]))
  }
  zz[, value := value / sum(value), by = .(zone, from_class)]     # -> shares

  pix <- work[, .(area = sum(area)), by = .(pixel, zone, from_class = class)]
  done <- NULL
  if (allocator && !is.null(coords) && requireNamespace("terra", quietly = TRUE)) {
    af <- tryCatch({
      src <- "../cascadinggamble-core/external/cascade-core/R/gravity_allocator.R"
      if (!exists("allocate_transitions_spatial") && file.exists(src)) source(src)
      get("allocate_transitions_spatial")
    }, error = function(e) NULL)
    if (!is.null(af)) {
      pd <- merge(pix, coords, by = "pixel")
      done <- tryCatch(as.data.table(af(pixel_dt = pd, shares_dt = zz, zone_cols = "zone",
          pixel_id_col = "pixel", from_col = "from_class", to_col = "to_class",
          share_col = "value", area_col = "area", coord_cols = c("x","y"), pixel_size = pixel_size)),
          error = function(e) { if (verbose) cat("  (gravity allocator failed: ", conditionMessage(e),
                                                 " -- falling back to proportional)\n", sep=""); NULL })
    }
  }
  if (is.null(done)) {   # proportional fallback: same margins, no neighbourhood structure
    done <- merge(pix, zz, by = c("zone","from_class"), allow.cartesian = TRUE)
    done[, area := area * value]
  }
  out <- done[, .(area = sum(area)), by = .(pixel, zone, class = to_class)]

  P1 <- out[, .(p1 = sum(area)), by = pixel]
  chk <- merge(P0, P1, by = "pixel", all = TRUE); chk[is.na(chk)] <- 0
  dev <- max(abs(chk$p0 - chk$p1))
  if (dev > tol * max(1, max(chk$p0)))
    stop(sprintf("harmonise changed pixel totals (max dev %.3g) -- it must only move composition", dev))
  if (verbose) {
    zf <- merge(out[, .(got = sum(area)), by = .(zone, class)], z[, .(zone, class, tgt_s)],
                by = c("zone","class"), all = TRUE)
    zf[is.na(zf)] <- 0
    cat(sprintf("  pixel totals preserved (max dev %.2g) | zone targets hit to %.3g km2\n",
                dev, max(abs(zf$got - zf$tgt_s))))
  }
  rbind(keep, out[, .(pixel, zone, class, area)], fill = TRUE)
}
