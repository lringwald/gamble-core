#!/usr/bin/env bash
# =============================================================================
# run/fit_flat_bart.sh — GLOBIOM flat MNL with BART on the terrain surface
# =============================================================================
# Usage: run/fit_flat_bart.sh [smoke|prod]      Default: prod
#
#   smoke  300 sweeps (burn 100) x 1 chain, "test" namespace   — checks the partition, not the fit
#   prod  8000 sweeps (burn 2000) x 4 chains, "production"     — ~18-24 h
#
# WHAT IS NON-PARAMETRIC HERE: the terrain block ONLY --
#   Slope_rad, Elevation, Aspect_cos_mean, Aspect_sin_mean, lon, lat
# Everything else (climate, soil, socio-economic, accessibility, yields) stays LINEAR and
# interpretable, and the focal LU-lags stay linear because they are autoregressive. Measured:
# pure-topo-in-BART beats the linear+BART hybrid, and symmetric trees + prevalence-scaled leaf
# shrinkage are super-additive (+75.4 held-out together vs +38.4 for symmetric alone).
#
# BART ONLY EXISTS ON THE FLAT PATH. codes/nested_cut.R has no BART support, so this does NOT
# produce a nested fit -- it is the single-level MNL over all 27 GLOBIOM classes at once.
#
# DRIVER_ADD_COORDS is what puts lon/lat in the design. Without it DRIVER_BART_COLS=topo resolves
# to the 4 terrain columns and says so in the log instead of failing, which is easy to miss.
#
# The design is rebuilt from the parquet every run (~30 min): the flat driver has no dump-loading
# path, so a design dump cannot be reused here the way run/nested.R reuses one.
#
# KILL-SAFE: posteriors stream to disk in batches and "production" hot-starts each chain from its
# final state on disk, so re-running this script with the SAME arguments RESUMES an interrupted
# run. "smoke" uses the test namespace, which is cleared at the start and never hot-starts.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."

MODE="${1:-prod}"
case "$MODE" in
  smoke) NITER=300;  NBURN=100;  NCHAINS=1; RUN_MODE=test ;;
  prod)  NITER=8000; NBURN=2000; NCHAINS=4; RUN_MODE=production ;;
  *) echo "FATAL: mode must be smoke|prod"; exit 1 ;;
esac
THIN=4

[ -f drivers/run_lu_pixel_model.R ] || { echo "FATAL: run from the repo root"; exit 1; }
# lon/lat come from an sf reprojection of the EPSG:3035 coordinates; dbarts carries the ensembles.
Rscript -e 'q(status = as.integer(!all(sapply(c("sf","dbarts","qs2"), requireNamespace, quietly = TRUE))))' \
  || { echo "FATAL: need the sf, dbarts and qs2 packages"; exit 1; }
# Same guard run/fit_nested.sh and entrypoint.sh use: a checkout without the 2026-08-19 symmetric-HS
# fix misallocates probability across classes while still reporting a symmetric prior.
grep -q "kronecker(Msym" codes/mnlogit_rcpp_sym.R \
  || { echo "FATAL: codes/mnlogit_rcpp_sym.R lacks the symmetric-HS fix"; exit 1; }

# Reading the 1 km parquet straight off Google Drive intermittently dies with
# `IOError ... [errno 60] Operation timed out` under arrow's parallel reads; prefer the local copy.
LOCAL_PARQUET=/Users/leopoldringwald/gamble_local_data/prior_model_1km_master_inputs.parquet
[ -f "$LOCAL_PARQUET" ] && export GAMBLE_MASTER_PARQUET="${GAMBLE_MASTER_PARQUET:-$LOCAL_PARQUET}"

TS=$(date +%Y%m%d_%H%M)
LOG="output/flat_bart_topo_${MODE}_${TS}.log"
mkdir -p output

cat <<EOF
>>> GLOBIOM flat MNL, BART on terrain + coordinates
    mode        : $MODE ($RUN_MODE namespace$([ "$RUN_MODE" = production ] && echo ", resumable"))
    sweeps      : $NITER (burn $NBURN, thin $THIN) x $NCHAINS chain(s)
    BART on     : Slope_rad, Elevation, Aspect_cos_mean, Aspect_sin_mean, lon, lat
    linear      : intercept + 28 focal LU-lags + climate/soil/socio/yields
    parquet     : ${GAMBLE_MASTER_PARQUET:-<driver resolves>}
    log         : $LOG
EOF

DRIVER_CLASS_COLS=GLOBIOM_subclass \
DRIVER_PROMOTE_NATURAL_OTHER=TRUE \
DRIVER_MODEL_YEARS=2018 DRIVER_FOCAL_YEARS=2010 DRIVER_COV_YEARS=2020 \
DRIVER_USE_BART=TRUE DRIVER_BART_COLS=topo \
DRIVER_ADD_COORDS=TRUE DRIVER_COORD_TYPE=lonlat \
DRIVER_NITER="$NITER" DRIVER_NBURN="$NBURN" DRIVER_THIN="$THIN" DRIVER_NCHAINS="$NCHAINS" \
DRIVER_RUN_MODE="$RUN_MODE" \
nohup Rscript drivers/run_lu_pixel_model.R > "$LOG" 2>&1 &
echo ">>> pid $!   tail -f \"$LOG\""
