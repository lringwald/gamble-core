# BMLEH_Los1_CAPRI — livestock prior: settled specification

## Classification: plain `BMLEH_Los1_label`, no intensity

    DRIVER_CLASS_COLS=BMLEH_Los1_label
    DRIVER_RE_GROUP_COL=CAPRI_NUTS

LU covariates come out as `Cropland_arable_other`, `Cropland_permanent_other`,
`Cropland_permanent_energy`, `Grassland_extensive`, `Grassland_intensive`, `Forests_*`, etc.

DECIDED, and not to be revisited on the grounds below: this loses the GLOBIOM cropland INTENSITY
split (Cropland_HI / LI / IR -> Cropland_arable_other + Cropland_permanent_other). That is accepted
because the LU downscaling in BMLEH_Los1_CAPRI produces exactly these classes and no intensity, and
`gamma` is keyed by DRIVER NAME -- a prior trained on covariates the downscaler cannot construct is
not usable, however well it fits. Concatenating `BMLEH_Los1_label,GLOBIOM_mngmt` to keep intensity
was considered and rejected for the same reason.

Grassland is unaffected: `Pasture_HI`/`Pasture_LI` and `Grassland_intensive`/`Grassland_extensive`
are the same land, so the GRAZING-AREA OFFSET -- which is the model's exposure term -- carries over
exactly. That offset is now matched by pattern over `(Pasture|Grassland)` and fails loudly rather
than collapsing to a constant, which is what it silently did under the hardcoded GLOBIOM name.

## Composition: 8 CAPRI-aligned cattle categories, not D/O/F

Estimated from `agr_r_animal` at NUTS2 (not apportioned). Partition closes against the A2000 total:
median ratio 1.0000, 98.6% of regions within 2%, 219 of 227 regions complete.

    DCOW A2300F -> DCOH+DCOL      SCOW A2300G -> SCOW (exact)
    HEIR A2230C+A2220C -> HEIR (exact)   HEIF A2230B+A2220B -> HEIH+HEIL
    BULL A2120+A2130 -> BULH+BULL        CAMR A2110C -> CAMR (exact)
    CAFR A2210C -> CAFR (exact)          CAFF A2010B -> CAMF+CAFF

The four remaining pairs (yield, weight x2, fattening-calf sex) are NOT split: no regional source
carries them and inventing them would put structure in the prior that no data supports.
PIGF / SOWS / PKPL are not covered at all -- the count model fits grazing species only.

## Known limitation, accepted

Only 2000 and 2018 granular LUM maps exist, so `lum_year <- if (cov_year <= 2000) 2000 else 2018`.
The 2010 tier is therefore explained by land use from EIGHT YEARS LATER. 2020 uses 2018, a two-year
lag, which is fine. Worth stating in any write-up rather than leaving implicit.
