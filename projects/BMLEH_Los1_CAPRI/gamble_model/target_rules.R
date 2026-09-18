# =============================================================================
# BMLEH_Los1_CAPRI — model class -> target classification rules
# =============================================================================
# Consumed by codes/target_class_split.R (generic cascade). ORDER MATTERS: a split may feed a
# later split, which is how the NUTS2/NUTS0 resolution trade-off is composed --
# HRL "Fruits" -> citrus + other_fruit at NUTS2, then other_fruit -> apples + other_fruit at NUTS0.
#
# Target: data/BMLEH_Los1_thematic_mapping.csv (43 fine classes, 33 of them crops).
# Shares:  results/gamble_model/BMLEH_Los1_CAPRI/crop_shares_nuts2.csv
# =============================================================================
suppressMessages(library(data.table))

BMLEH_TARGET_RULES <- rbindlist(list(
  # ---- 1:1 renames: HRL identifies these uniquely, no statistics needed -------------------------
  # Energy (short-rotation coppice) and the grassland intensity split are NOT here: the LUM map
  # already distinguishes them per code, so they are resolved at mapping level by the
  # BMLEH_Los1_label column (build_class_column.R). This file only handles what genuinely needs a
  # crop statistic to resolve. Run the driver with DRIVER_CLASS_COLS=BMLEH_Los1_label.
  data.table(order = 1:10, action = "rename",
    from_class = c("Barley","Rice","Dry_pulses","Potatoes","Sugar_beet","Sunflower","Soybeans",
                   "Rapeseed","Flax_cotton_hemp","Olives"),
    to         = c("Cropland_arable_barley","Cropland_arable_rice","Cropland_arable_pulses",
                   "Cropland_arable_potatoes","Cropland_arable_sugar_beet","Cropland_arable_sunflower",
                   "Cropland_arable_soya","Cropland_arable_rape","Cropland_arable_textile",
                   "Cropland_permanent_olives")),

  # ---- NUTS2 splits (apro_cpshr): regional composition -----------------------------------------
  data.table(order = 12:15, action = "split",
    from_class = c("Wheat","Maize","Other_cereals","Fruits"),
    to         = c("Wheat","Maize","Other_cereals","Fruits")),

  # ---- the PERMANENT residual: LUM permanent land HRL could not attribute ----------------------
  # Must run BEFORE the detail splits below, because its outputs feed them: it emits "Grapes" (HRL's
  # own class name) so rule 19 refines it into wine/table, and Cropland_permanent_other_fruit so
  # rule 17 refines it into apples. This was a bare rename into other_fruit until 2026-09-10, which
  # buried ~79k km2 -- three quarters of all permanent land -- in one class and left olives at 0.07
  # of their Eurostat area. HRL sees only 44.6k of the 118.2k km2 Eurostat reports, so most real
  # olive and vine area is in here, not in the HRL-identified classes.
  data.table(order = 16, action = "split",
    from_class = "Cropland_permanent_other", to = "<permanent residual>"),

  # ---- NUTS0 splits (apro_cpsh1): national composition, applied to the NUTS2 output -------------
  # Fruits_detail deliberately runs AFTER Fruits: it refines the other_fruit remainder that the
  # NUTS2 citrus split produced, so both resolutions contribute.
  data.table(order = 17:19, action = "split",
    from_class = c("Cropland_permanent_other_fruit","Fresh_vegetables","Grapes"),
    to         = c("Fruits_detail","Veg_detail","Grapes_detail")),

  # ---- the arable residual: everything HRL could not see ---------------------------------------
  data.table(order = 20, action = "split",
    from_class = "Cropland_arable_other", to = "<arable residual>"),

  # ---- remaining 1:1, after the splits ---------------------------------------------------------
  # Nuts has no BMLEH class of its own, so it folds into other_fruit. Deliberately AFTER
  # Fruits_detail, unlike the permanent residual: Eurostat's F0000 genuinely contains apples (so the
  # residual should be split), whereas area HRL positively identified as nuts is not apples and must
  # not be.
  data.table(order = 21, action = "rename",
    from_class = "Nuts", to = "Cropland_permanent_other_fruit")
), use.names = TRUE)

# NOT handled here (they are not crop classes and come straight from the LUM map):
#   Grassland_extensive / Grassland_intensive  <- BMLEH_Los1_label, from the LUM density classes
#   Cropland_arable_energy                     <- BMLEH_Los1_label, from LUM code 7 short rotation
#   Natural_unmanaged_fallow / _set_aside      <- LUM map; Natural_unmanaged is the
#                                                 downscaling fallback
#   Artificial, Forests, Natural_other, Waterbodies_inland, Wetlands_natural <- LUM directly
