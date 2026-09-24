#!/usr/bin/env bash
# =============================================================================
# GLOBIOM_subclass — build the design, then fit the prior in resumable batches
# =============================================================================
#   ./projects/GLOBIOM_subclass/gamble_model/run.sh design   # ~30 min, writes the dump
#   ./projects/GLOBIOM_subclass/gamble_model/run.sh fit      # segment 1: 4000 burn + 1000 kept
#   ./projects/GLOBIOM_subclass/gamble_model/run.sh more     # +500 sweeps, continues the LATEST
#   ./projects/GLOBIOM_subclass/gamble_model/run.sh more <dir>   # continue a specific segment
#   ./projects/GLOBIOM_subclass/gamble_model/run.sh score    # re-score what is on disk, no sampling
#   ./projects/GLOBIOM_subclass/gamble_model/run.sh smoke    # design + 5k-pixel 1-chain fit
#   ./projects/GLOBIOM_subclass/gamble_model/run.sh all
#
# Run from the gamble-core root. Design switches that are load-bearing and non-default:
#   DRIVER_CLASS_COLS=GLOBIOM_subclass
#     27 modelled classes plus the CURATED GLOBIOM_subclass_nest column (Cropland / Forests /
#     Natural / Pasture / Urban / Waterbodies).
#   DRIVER_PROMOTE_NATURAL_OTHER=TRUE
#     Moves Natural_other OUT of the non-choosable residual so it is modelled.
#   NO DRIVER_NO_CHOICE IS SET, DELIBERATELY.
#     A curated nest column SUPERSEDES the no_choice collapse (run_lu_pixel_model.R ~447-458):
#     NO_CHOICE_LU is emptied and the TREE does the grouping instead. Copying BMLEH's
#     DRIVER_NO_CHOICE line here would fight the curated tree, not help it.
#   DRIVER_FOCAL_YEARS=2010 with DRIVER_MODEL_YEARS=2018
#     The focal (neighbourhood composition) covariate is LAGGED by 8 years. That lag is what makes
#     the spatial Y-lag defensible: a t-1 neighbourhood cannot be simultaneously determined with
#     the outcome, which is the objection a contemporaneous focal invites.
#
# Sampler settings are the STANDARD (docs/sampler_model_specification.md sec. 8.0) and live in
# estimate_prior.R, not here. Measurement arms are env overrides, e.g.
#   GB_RE_ASIS=FALSE ./run.sh fit          # the ASIS arm
#   GB_RE_SUPPORT=0  ./run.sh fit          # the knob sec. 8.1b leaves unreconciled
#   GB_BART=TRUE     ./run.sh fit          # the BART treatment arm
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/../../.."
[ -f codes/mnlogit_rcpp_sym.R ] || { echo "run from the gamble-core root"; exit 1; }

P=projects/GLOBIOM_subclass/gamble_model
DUMP=output/designs/pixel_model_inputs_GLOBIOM_subclass.rds
RES=results/gamble_model/GLOBIOM_subclass
mkdir -p "$(dirname "$DUMP")"

# The local parquet copy exists because reading straight off Google Drive intermittently dies with
# `IOError ... [errno 60] Operation timed out` under arrow's parallel reads. Use it when present;
# otherwise let the driver resolve the canonical copy itself.
# NOTE: this copy must be kept in step with the canonical one. It was an Aug-17 vintage still
# carrying the pre-rename LAMASUS_1km_bufferID key until 2026-09-16; a stale copy makes the driver
# ask for a column the file does not have.
LOCAL_PARQUET=/Users/leopoldringwald/gamble_local_data/prior_model_1km_master_inputs.parquet
if [ -f "$LOCAL_PARQUET" ]; then
  export GAMBLE_MASTER_PARQUET="${GAMBLE_MASTER_PARQUET:-$LOCAL_PARQUET}"
else
  echo "    (no local parquet at $LOCAL_PARQUET -- letting the driver resolve the canonical copy)"
fi

design() {
  echo ">>> design: GLOBIOM_subclass (27 classes + curated 6-nest tree), 10 km, focal lagged to 2010"
  # BOTH dump flags are required: DRIVER_DUMP_EXIT is nested INSIDE the DUMP_INPUTS block, so on its
  # own it neither dumps nor exits -- it silently runs a full fit instead.
  DRIVER_DUMP_INPUTS=TRUE DRIVER_DUMP_EXIT=TRUE \
  DRIVER_DUMP_PATH="$DUMP" \
  DRIVER_CLASS_COLS=GLOBIOM_subclass \
  DRIVER_PROMOTE_NATURAL_OTHER=TRUE \
  DRIVER_LUM_LINEAGE=globiom \
  DRIVER_PIXEL_RES="${DRIVER_PIXEL_RES:-10}" \
  DRIVER_MODEL_YEARS=2018 \
  DRIVER_FOCAL_YEARS=2010 \
  DRIVER_COV_YEARS=2020 \
  DRIVER_PIXEL_INTERSECT="${DRIVER_PIXEL_INTERSECT:-NUTS3}" \
  Rscript drivers/run_lu_pixel_model.R
  # FAIL LOUDLY: callers pipe this through `tail`, and a pipeline reports the LAST command's status,
  # so a build that died on a data error still exited 0 and looked finished. Assert the dump exists.
  [ -f "$DUMP" ] || { echo ">>> DESIGN BUILD FAILED: no dump at $DUMP (see the error above)"; exit 1; }
  echo ">>> design written to $DUMP"
}

# Budget: 4000 burn-in + 1000 kept, then extend with `more` (500 each, NBURN=0 since the chain is
# already burnt in). Every run streams its draws AND a per-chain final_state, so the budget is not a
# one-shot decision: take 1000, read the diagnostics, extend if anything is ragged.
# BURN-IN RAISED 2000 -> 4000 (2026-09-17). The 2026-09-16 fit did not converge at 2000 burn-in
# (21% of mu above Rhat 1.05, min ESS 6). Two things changed since: the collinear yield block is
# gone (VIF 312/108/83 -> 9.1/8.5, which should help) and the RE block grew from intercept+4 to
# intercept+6 (which cuts the other way -- more per-group parameters mix slower). A longer burn-in
# is the cheap hedge against the second.
# The first segment carries the burn-in, so it is the expensive one -- 5000 sweeps x 4 chains on
# ~64k pixels x 27 classes. Budget hours, not minutes, and run it backgrounded (see README).
fit() { echo ">>> fit: flat MNL, standard config, 4 chains (4000 burn + 1000 kept)"
        GB_INPUT="${GB_INPUT:-${DESIGN_PATH:-$DUMP}}" \
        GB_NITER="${GB_NITER:-5000}" GB_NBURN="${GB_NBURN:-4000}" GB_THIN="${GB_THIN:-1}" \
        Rscript "$P/estimate_prior.R" "$@"
        local latest
        latest=$(ls -1d "$RES"/prior_* 2>/dev/null | sort -r | head -1)
        if [ -n "$latest" ] && [ -d "$latest/posterior" ]; then
          echo ">>> Posterior sweep diagnostics on $latest"
          Rscript postprocess/diagnose_posterior.R "$latest" || true
        fi
      }

# Continue sampling from a finished segment's saved chain state. Default target is the newest
# RESUMABLE segment, so repeated `more` calls chain forward one segment at a time.
more() {
  local from="${1:-}"
  if [ -z "$from" ]; then
    # Pick the newest RESUMABLE segment, not merely the newest directory. An aborted run leaves an
    # empty folder behind, and by name that husk sorts newest.
    for d in $(ls -1d "$RES"/prior_* 2>/dev/null | sort -r); do
      if ls "$d"/posterior/final_state_chain_*.qs >/dev/null 2>&1; then from="$d"; break; fi
      echo "    (skipping $(basename "$d") -- no chain state; aborted or still running)"
    done
  fi
  if [ -z "$from" ] || [ ! -d "$from" ]; then
    echo "no RESUMABLE segment under $RES."
    echo "a resumable segment is one that finished and wrote posterior/final_state_chain_*.qs:"
    ls -1d "$RES"/prior_* 2>/dev/null | while read -r d; do
      printf "    %-28s states %s  batches %s\n" "$(basename "$d")" \
        "$(ls "$d"/posterior/final_state_chain_*.qs 2>/dev/null | wc -l | tr -d ' ')" \
        "$(ls "$d"/posterior/posterior_batch_*.qs 2>/dev/null | wc -l | tr -d ' ')"
    done
    exit 1
  fi
  echo ">>> continuing from $from"
  GB_INPUT="${GB_INPUT:-${DESIGN_PATH:-$DUMP}}" GB_RESUME="$from" \
    GB_NITER="${GB_NITER:-500}" GB_THIN="${GB_THIN:-1}" Rscript "$P/estimate_prior.R"
}

# Re-score what is already on disk (pools every chain and every ancestor segment). Useful after a
# run was killed mid-flight: batches without a final_state are unusable for RESUME but perfectly
# scorable.
score() {
  local from="${1:-}"
  if [ -z "$from" ]; then from=$(ls -1d "$RES"/prior_* 2>/dev/null | sort -r | head -1); fi
  [ -n "$from" ] && [ -d "$from" ] || { echo "no segment to score under $RES"; exit 1; }
  echo ">>> scoring $from (no sampling)"
  GB_INPUT="${GB_INPUT:-${DESIGN_PATH:-$DUMP}}" GB_RESUME="$from" GB_SCORE_ONLY=TRUE Rscript "$P/estimate_prior.R"
}

case "${1:-all}" in
  design) design ;;
  fit)    fit ;;
  more)   more "${2:-}" ;;
  score)  score "${2:-}" ;;
  smoke)  design; fit --smoke ;;
  all)    design; fit ;;
  *) echo "usage: $0 [design|fit|more [dir]|score [dir]|smoke|all]"; exit 1 ;;
esac
