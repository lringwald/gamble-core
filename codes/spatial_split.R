# =============================================================================
# spatial_split.R — train / held-out splits that account for spatial dependence
# =============================================================================
# A random pixel split is OPTIMISTIC on this grid: neighbouring pixels are strongly
# autocorrelated, so nearly every held-out pixel has a training neighbour. Worse, the model's
# `focal_*` covariates ARE the mean of a pixel's 8 neighbours (compute_focal_coord, queen-8,
# centre excluded) -- so under a random split a held-out pixel's own predictors are built from
# training pixels' RESPONSES. That is direct leakage of the thing being predicted.
#
# This routine builds four kinds of split, and the default ("hybrid") builds two held-out sets at
# once so a single fit yields BOTH an interpolation score and a spatial-transfer score. The gap
# between them is the optimism, measured on the design at hand rather than assumed.
#
#   random   uniform over pixels -- the legacy behaviour, kept for comparability
#   country  uniform WITHIN each RE group, so every country appears in training
#   block    whole contiguous square tiles
#   hybrid   half country-stratified random + half blocks, scored separately   [default]
#
# GEOMETRY IS DERIVED, NOT CHOSEN. The block side is a multiple of the focal window:
#
#     res_m = min(diff(sort(unique(coord_X))))       <- same rule compute_focal_coord uses
#     side  = block_mult * (2 * focal_reach + 1) cells
#
# `focal_reach` is 1 because the focal is queen-8 and that is HARDCODED in compute_focal_coord
# (`CJ(dx = c(-r, 0, r), dy = c(-r, 0, r))` minus the centre). Deriving res_m the same way the
# focal does means the block geometry cannot drift out of step with the covariate it exists to
# defeat, and it scales automatically with resolution: block_mult = 5 gives 150 km blocks on a
# 10 km design and 75 km on a 5 km design.
#
# THE BUFFER IS THE PART THAT ACTUALLY WORKS. Block size alone does not remove the leakage: a
# held-out pixel on a block EDGE still has training pixels among its 8 neighbours, and the focal
# is built from exactly those. So a ring of `buffer_cells` around every test block is dropped from
# training entirely -- it is neither trained on nor scored. Without it, a "spatial" split still
# measures something close to interpolation while claiming otherwise.
#
# Returns a list with `train`, `test`, `test_random`, `test_block`, `buffer` (row indices into the
# supplied vectors) and `meta` (the realised geometry and fractions -- realised, because block
# sampling lands NEAR a target fraction, not on it, and a fit whose "20% hold-out" is really 14%
# must not be silently mislabelled).
# =============================================================================

make_spatial_split <- function(coord_X, coord_Y, group,
                               test_frac    = 0.2,
                               mode         = c("hybrid", "random", "country", "block"),
                               block_mult   = 5L,
                               focal_reach  = 1L,
                               buffer_cells = 1L,
                               block_share  = 0.5,   # of test_frac, the part taken as blocks
                               min_train_per_group = 30L,
                               seed         = 20260916L) {
  mode <- match.arg(mode)
  n <- length(group)
  stopifnot(n > 0, test_frac > 0, test_frac < 1)

  # ---- random / country need no geometry -----------------------------------------------------
  if (mode %in% c("random", "country") || is.null(coord_X) || is.null(coord_Y)) {
    if (mode == "block" || mode == "hybrid")
      stop("mode '", mode, "' needs coord_X / coord_Y; the design dump carries them as coord_X/coord_Y ",
           "(rebuild the design if they are NULL).")
    set.seed(seed)
    te <- if (mode == "random") {
      sort(sample.int(n, round(n * test_frac)))
    } else {
      # stratify by group: every country present on BOTH sides, so "unseen group" never confounds
      # an arm comparison. The scorer drops test rows whose group is absent from training.
      sort(unlist(lapply(split(seq_len(n), group), function(ix)
        sample(ix, min(length(ix) - 1L, max(0L, round(test_frac * length(ix)))))), use.names = FALSE))
    }
    tr <- setdiff(seq_len(n), te)
    return(list(train = tr, test = te, test_random = te, test_block = integer(0),
                buffer = integer(0),
                meta = list(mode = mode, test_frac_requested = test_frac,
                            test_frac_realised = length(te) / n,
                            n_train = length(tr), n_test = length(te),
                            n_test_random = length(te), n_test_block = 0L, n_buffer = 0L)))
  }

  stopifnot(length(coord_X) == n, length(coord_Y) == n)
  if (any(!is.finite(coord_X)) || any(!is.finite(coord_Y)))
    stop("coord_X / coord_Y contain non-finite values; cannot tile the grid.")

  # ---- grid geometry, inferred exactly as compute_focal_coord infers it -----------------------
  ux <- sort(unique(coord_X)); uy <- sort(unique(coord_Y))
  if (length(ux) < 2L || length(uy) < 2L) stop("need at least two distinct X and Y coordinates to infer the grid.")
  res_m <- min(c(diff(ux), diff(uy)))
  if (!is.finite(res_m) || res_m <= 0) stop("could not infer a positive grid spacing from the coordinates.")

  side_cells <- as.integer(block_mult) * (2L * as.integer(focal_reach) + 1L)
  side_m     <- side_cells * res_m

  ix <- as.integer(round((coord_X - min(coord_X)) / res_m))
  iy <- as.integer(round((coord_Y - min(coord_Y)) / res_m))
  cellkey <- as.numeric(ix) * 1e6 + as.numeric(iy)          # unique per occupied cell

  bx <- ix %/% side_cells; by <- iy %/% side_cells
  tile <- paste0(bx, "_", by)
  tiles <- unique(tile)

  set.seed(seed)
  target_block_n  <- if (mode == "block") round(n * test_frac) else round(n * test_frac * block_share)
  target_random_n <- if (mode == "block") 0L else round(n * test_frac * (1 - block_share))

  # ---- choose whole tiles, refusing any that would strand a country --------------------------
  # Online guard: a tile is skipped if taking it would leave some group with fewer than
  # `min_train_per_group` pixels. Dropping a country out of training is not a small error -- the
  # scorer silently discards test rows whose group never appeared in training, so those pixels
  # would vanish from the score with no warning at all.
  grp_total <- table(group)
  grp_left  <- grp_total
  ord <- sample(tiles)
  chosen <- character(0); n_block <- 0L; n_skipped_guard <- 0L
  tile_rows <- split(seq_len(n), tile)
  for (tl in ord) {
    if (n_block >= target_block_n) break
    rows <- tile_rows[[tl]]
    tg <- table(group[rows])
    left_after <- grp_left
    left_after[names(tg)] <- left_after[names(tg)] - tg
    if (any(left_after[names(tg)] < min_train_per_group)) { n_skipped_guard <- n_skipped_guard + 1L; next }
    chosen <- c(chosen, tl); n_block <- n_block + length(rows); grp_left <- left_after
  }
  test_block <- if (length(chosen)) sort(unlist(tile_rows[chosen], use.names = FALSE)) else integer(0)

  # ---- buffer ring: cells within `buffer_cells` of a test-block cell, excluded from training ---
  buffer <- integer(0)
  if (length(test_block) && buffer_cells > 0L) {
    b <- as.integer(buffer_cells)
    tkey <- numeric(0)
    for (dx in -b:b) for (dy in -b:b)
      tkey <- c(tkey, (as.numeric(ix[test_block]) + dx) * 1e6 + (as.numeric(iy[test_block]) + dy))
    near <- which(cellkey %in% unique(tkey))
    buffer <- sort(setdiff(near, test_block))
  }

  # ---- random half, drawn from what is left, stratified by group -----------------------------
  eligible <- setdiff(seq_len(n), c(test_block, buffer))
  test_random <- integer(0)
  if (target_random_n > 0L && length(eligible)) {
    frac <- min(0.9, target_random_n / length(eligible))
    test_random <- sort(unlist(lapply(split(eligible, group[eligible]), function(ix2) {
      k <- min(length(ix2) - 1L, max(0L, round(frac * length(ix2))))
      if (k <= 0L) integer(0) else sample(ix2, k)
    }), use.names = FALSE))
  }

  test  <- sort(c(test_block, test_random))
  train <- setdiff(seq_len(n), c(test, buffer))

  # ---- final verification: no group may be missing from training ------------------------------
  missing_grp <- setdiff(as.character(unique(group)), as.character(unique(group[train])))
  if (length(missing_grp))
    warning("spatial split leaves no training pixels for group(s): ", paste(missing_grp, collapse = ", "),
            ". Their held-out rows will be DROPPED by the scorer. Lower block_mult or raise test_frac.",
            call. = FALSE)

  list(train = train, test = test, test_random = test_random, test_block = test_block, buffer = buffer,
       meta = list(
         mode = mode, res_m = res_m, res_km = res_m / 1000,
         focal_reach_cells = focal_reach, block_mult = block_mult,
         block_side_cells = side_cells, block_side_km = side_m / 1000,
         buffer_cells = buffer_cells,
         n_tiles_total = length(tiles), n_tiles_test = length(chosen),
         n_tiles_skipped_by_group_guard = n_skipped_guard,
         test_frac_requested = test_frac, block_share_requested = block_share,
         test_frac_realised = length(test) / n,
         block_frac_realised = length(test_block) / n,
         random_frac_realised = length(test_random) / n,
         buffer_frac = length(buffer) / n,
         n_train = length(train), n_test = length(test),
         n_test_block = length(test_block), n_test_random = length(test_random),
         n_buffer = length(buffer),
         groups_missing_from_train = missing_grp))
}

# Pretty-print what a split actually did. Realised numbers, not requested ones.
print_spatial_split <- function(sp) {
  m <- sp$meta
  if (identical(m$mode, "random") || identical(m$mode, "country")) {
    cat(sprintf("  split: %s | train %d / test %d (%.1f%%)\n",
                m$mode, m$n_train, m$n_test, 100 * m$test_frac_realised))
    return(invisible(NULL))
  }
  cat(sprintf("  split: %s | grid %.0f km | blocks %d cells = %.0f km (mult %d x queen-%d window) | buffer %d cell(s)\n",
              m$mode, m$res_km, m$block_side_cells, m$block_side_km, m$block_mult,
              2 * m$focal_reach_cells + 1, m$buffer_cells))
  cat(sprintf("         tiles %d of %d held out%s\n", m$n_tiles_test, m$n_tiles_total,
              if (m$n_tiles_skipped_by_group_guard > 0)
                sprintf(" (%d skipped: would have stranded a country)", m$n_tiles_skipped_by_group_guard) else ""))
  cat(sprintf("         train %d | test %d = %d block (%.1f%%) + %d random (%.1f%%) | buffer dropped %d (%.1f%%)\n",
              m$n_train, m$n_test, m$n_test_block, 100 * m$block_frac_realised,
              m$n_test_random, 100 * m$random_frac_realised, m$n_buffer, 100 * m$buffer_frac))
  invisible(NULL)
}
