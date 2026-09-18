# =============================================================================
# BMLEH_Los1_CAPRI — livestock: our subclasses -> CAPRI activity codes
# =============================================================================
# Our model resolves species x production subtype x system:
#   BOVD/BOVO/BOVF, SGTD/SGTO/SGTF  x  {conventional, organic}          (D/O/F composition)
#   BOVDCOW/BOVSCOW/BOVHEIR/BOVHEIF/BOVBULL/BOVCAMR/BOVCAFR/BOVCAFF      (8-category CAPRI cattle)
# BOTH cattle vocabularies are carried below. The 8-category one resolves FOUR CAPRI activities
# exactly, where D/O/F resolves none exactly.
# CAPRI resolves 16 leaf activities, which are FINER on cattle in two directions our model does not
# carry -- yield class (DCOH/DCOL), and weight/sex class (HEIH/HEIL, BULH/BULL, CAMR/CAFR/CAMF/CAFF).
#
# So the mapping is honest about being one-to-MANY. Each of our subclasses maps to a SET of CAPRI
# activities; the split WITHIN that set needs herd statistics this model never saw, and is left to the
# downscaler rather than invented here. Where CAPRI's own aggregates line up with our subtypes they
# are given too, because those are exact:
#   BEFM (beef meat)  = SCOW + HEIH + HEIL + BULH + BULL   == our BOVO
#   CACR (all dairy)  = DCOH + DCOL + HEIR + CAMR + CAFR + CAMF + CAFF == our BOVD + BOVF
#   CATA (all cattle) = every cattle activity              == BOVD + BOVO + BOVF
#
# NOT COVERED: PIGF (pig fattening) and SOWS (pig breeding). The count model fits BOV and SGT only --
# pigs are not a grazing species and were never in the herd totals -- so there is no prior for them
# here and none is fabricated. PKPL ("other animals") is likewise outside the model.
# =============================================================================
suppressMessages(library(data.table))

LIVESTOCK_CAPRI_MAP <- rbindlist(list(
  # ---- bovine dairy: the milking herd -----------------------------------------------------------
  data.table(subclass = "BOVD", capri = c("DCOH","DCOL"),
             capri_label = c("Dairy cows high yield","Dairy cows low yield"),
             split_needed = "yield class"),
  # ---- bovine other/suckler: the beef herd ------------------------------------------------------
  data.table(subclass = "BOVO", capri = c("SCOW","HEIH","HEIL","BULH","BULL"),
             capri_label = c("Other (suckler) cows","Heifers fattening high weight",
                             "Heifers fattening low weight","Male adult cattle high weight",
                             "Male adult cattle low weight"),
             split_needed = "weight class"),
  # ---- bovine followers: replacements and calves ------------------------------------------------
  data.table(subclass = "BOVF", capri = c("HEIR","CAMR","CAFR","CAMF","CAFF"),
             capri_label = c("Heifers breeding","Raising male calves","Raising female calves",
                             "Fattening male calves","Fattening female calves"),
             split_needed = "sex / raising vs fattening"),
  # ---- bovine, CAPRI-ALIGNED 8-category target -------------------------------------------------
  # When the composition is fitted on the 8 CAPRI cattle categories (LIVESTOCK_SPEC.md) the mapping
  # gets much tighter than the D/O/F block above: FOUR activities come out EXACTLY, and the other
  # four are pairs CAPRI splits by a dimension no regional source carries (yield, weight, calf sex).
  # These subclass names (BOV + category) coexist with the BOVD/BOVO/BOVF rows: the join takes
  # whichever vocabulary the composition fit actually produced, so both remain runnable.
  data.table(subclass = "BOVDCOW", capri = c("DCOH","DCOL"),
             capri_label = c("Dairy cows high yield","Dairy cows low yield"),
             split_needed = "yield class"),
  data.table(subclass = "BOVSCOW", capri = "SCOW",
             capri_label = "Other (suckler) cows", split_needed = NA_character_),
  data.table(subclass = "BOVHEIR", capri = "HEIR",
             capri_label = "Heifers breeding", split_needed = NA_character_),
  data.table(subclass = "BOVHEIF", capri = c("HEIH","HEIL"),
             capri_label = c("Heifers fattening high weight","Heifers fattening low weight"),
             split_needed = "weight class"),
  data.table(subclass = "BOVBULL", capri = c("BULH","BULL"),
             capri_label = c("Male adult cattle high weight","Male adult cattle low weight"),
             split_needed = "weight class"),
  data.table(subclass = "BOVCAMR", capri = "CAMR",
             capri_label = "Raising male calves", split_needed = NA_character_),
  data.table(subclass = "BOVCAFR", capri = "CAFR",
             capri_label = "Raising female calves", split_needed = NA_character_),
  data.table(subclass = "BOVCAFF", capri = c("CAMF","CAFF"),
             capri_label = c("Fattening male calves","Fattening female calves"),
             split_needed = "calf sex"),
  # ---- small grazers: CAPRI is NOT finer here, so these are the clean ones ----------------------
  data.table(subclass = "SGTD", capri = "SHGM", capri_label = "Milk ewes and goats", split_needed = NA_character_),
  data.table(subclass = c("SGTO","SGTF"), capri = "SHGF",
             capri_label = "Sheep and goat fattening", split_needed = "two of ours -> one of theirs")
), use.names = TRUE)

# CAPRI aggregates that our subtypes reproduce exactly (no within-set split needed)
# Given twice, once per cattle vocabulary. Which one applies depends on the composition fit; the
# D/O/F rows are kept so a three-category fit still resolves its aggregates.
LIVESTOCK_CAPRI_AGG <- rbindlist(list(
  data.table(capri    = c("BEFM","CACR","CATA"),
             label    = c("Beef meat activities","All dairy","All cattle activities"),
             subclass = c("BOVO","BOVD+BOVF","BOVD+BOVO+BOVF"),
             vocab    = "DOF"),
  data.table(capri    = c("BEFM","CACR","CATA"),
             label    = c("Beef meat activities","All dairy","All cattle activities"),
             subclass = c("BOVSCOW+BOVHEIF+BOVBULL",
                          "BOVDCOW+BOVHEIR+BOVCAMR+BOVCAFR+BOVCAFF",
                          "BOVDCOW+BOVSCOW+BOVHEIR+BOVHEIF+BOVBULL+BOVCAMR+BOVCAFR+BOVCAFF"),
             vocab    = "CAPRI8")))

# CAPRI activities with no counterpart in this model
LIVESTOCK_CAPRI_UNCOVERED <- data.table(
  capri = c("PIGF","SOWS","PKPL"),
  label = c("Pig fattening","Pig breeding","Other animals"),
  reason = "the count model fits grazing species (BOV, SGT) only; pigs were never in the herd totals")
