#!/usr/bin/env bash
# =============================================================================
# BMLEH_Los1_CAPRI — build the design, then fit the prior
# =============================================================================
#   ./projects/BMLEH_Los1_CAPRI/gamble_model/run.sh design   # ~10 min, writes the dump + GeoTIFF
#   ./projects/BMLEH_Los1_CAPRI/gamble_model/run.sh fit      # moderate segment, 4 chains
#   ./projects/BMLEH_Los1_CAPRI/gamble_model/run.sh more     # continue the LATEST segment
#   ./projects/BMLEH_Los1_CAPRI/gamble_model/run.sh more <dir>   # continue a specific one
#   ./projects/BMLEH_Los1_CAPRI/gamble_model/run.sh all
#   ./projects/BMLEH_Los1_CAPRI/gamble_model/run.sh smoke    # design + a 5k-pixel 1-chain fit
#
# Run from the gamble-core root. Two switches are load-bearing and non-default:
#   DRIVER_CROP_SPLIT=TRUE     CLASS_SCHEME derives to BMLEH and the HRL split only auto-enables for AGMIP.
#   DRIVER_RE_GROUP_COL=CAPRI_NUTS
#     The random-effect groups become CAPRI country codes (the driver slices CAPRI_NUTS to 2 chars for
#     the RE key), which is what the downscaling engine keys on. The default GLOB_country gives country
#     NAMES ("CzechRep", "UK"), and those are not translatable after the fact: CAPRI country codes are
#     not ISO (BL, IR, EL, CS, KO, MO) and the pixel-level mapping between the two is many-to-many --
#     the geometries disagree at borders, giving 150 distinct (CAPRI, GLOB) pairs. A fit grouped on
#     GLOB_country therefore cannot be relabelled into CAPRI keys; it has to be re-fitted.
#     geo_region stays on NUTS2 independently, for the Eurostat crop shares.
#   DRIVER_NO_CHOICE="Waterbodies_marine"
#     Only marine water is non-choosable. Waterbodies_inland (INLW), Wetlands_natural (TWET) and
#     Natural_other (OLND) are BMLEH target classes with CAPRI codes, so they must be MODELLED
#     classes -- inside no_choice the prior emits no transitions for them at all.
#   DRIVER_FOCAL_YEARS=2018
#     Focal must be CONTEMPORANEOUS. HRL crop types exist only for 2018, so a t-1 lag focal carries
#     no crop-typed neighbourhood at all -- every focal_<crop> column disappears and only the
#     residual-derived ones survive.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/../../.."
[ -f codes/mnlogit_rcpp_sym.R ] || { echo "run from the gamble-core root"; exit 1; }

P=projects/BMLEH_Los1_CAPRI/gamble_model
DUMP=output/designs/pixel_model_inputs_BMLEH_Los1_CAPRI.rds
# The share table is written by prepare_eurostat_crop_shares.R, whose OUT_DIR resolves to
# <root>/gamble_model/BMLEH_Los1_CAPRI -- "output" when run from the gamble-core root, "results" when
# run from the cascadinggamble root. This is that path.
# DO NOT repoint this at output/eurostat/crop_shares_nuts2.csv: that file is older, unrelated output
# produced by the SUPERSEDED generator (pre 2026-09-09), and it lacks the "<permanent residual>" and
# derived-OCRO groups that target_rules.R now requires -- the design build fails on
# "share table has no rows for split group(s): <permanent residual>".
SHARES=output/gamble_model/BMLEH_Los1_CAPRI/crop_shares_nuts2.csv
mkdir -p "$(dirname "$DUMP")"

# The local parquet copy exists because reading straight off Google Drive intermittently dies with
# `IOError ... [errno 60] Operation timed out` under arrow's parallel reads. Keep it in step with the
# canonical copy: a stale one makes the driver ask for a column the file does not have.
LOCAL_PARQUET=/Users/leopoldringwald/gamble_local_data/prior_model_1km_master_inputs.parquet
if [ -f "$LOCAL_PARQUET" ]; then
  export GAMBLE_MASTER_PARQUET="${GAMBLE_MASTER_PARQUET:-$LOCAL_PARQUET}"
fi

design() {
  [ -f "$SHARES" ] || { echo "no share table; run: Rscript $P/prepare_eurostat_crop_shares.R"; exit 1; }
  echo ">>> design: BMLEH_Los1_label + HRL crop split + target cascade"
  DRIVER_DUMP_INPUTS=TRUE DRIVER_DUMP_EXIT=TRUE DRIVER_DUMP_TIF=TRUE \
  DRIVER_DUMP_PATH="$DUMP" \
  DRIVER_CLASS_COLS=BMLEH_Los1_label \
  DRIVER_CROP_SPLIT=TRUE \
  DRIVER_RE_GROUP_COL=CAPRI_NUTS \
  DRIVER_LUM_LINEAGE=globiom \
  DRIVER_NO_CHOICE="Waterbodies_marine" \
  DRIVER_FOCAL_YEARS=2018 \
  DRIVER_TARGET_RULES="$P/target_rules.R" \
  DRIVER_TARGET_SHARES="$SHARES" \
  DRIVER_PIXEL_RES="${DRIVER_PIXEL_RES:-10}" \
  DRIVER_PIXEL_INTERSECT="${DRIVER_PIXEL_INTERSECT:-NUTS3}" \
  Rscript drivers/run_lu_pixel_model.R
  # FAIL LOUDLY. Callers pipe this through `tail`, and a pipeline reports the LAST command's status:
  # a design build that died on a data error still exited 0 and was reported as finished. Assert the
  # dump exists rather than trusting the exit code.
  [ -f "$DUMP" ] || { echo ">>> DESIGN BUILD FAILED: no dump at $DUMP (see the error above)"; exit 1; }
  echo ">>> design written to $DUMP"
}
RES=results/gamble_model/BMLEH_Los1_CAPRI
# Production sampler settings, all three established by measurement rather than default:
#   BM_RE_ASIS=TRUE   RE-scale ASIS. Without it the RE variance does not converge -- sigma_re Rhat
#                     1.618 and joint log_lik Rhat 1.878, against 1.007 and 1.037 with it. This is
#                     the single setting that decides whether the posterior is quotable.
#   BM_SLAB_C2=FALSE  c2 fixed rather than sampled. It never binds (cap sigma<=3.8 vs sigma_re 0.108)
#                     and sampling it adds the worst-mixing scalar in the model for nothing.
#   BM_THIN=1         thinning never improves ESS, it only discards draws; storage is not the
#                     constraint here (~60 MB per 500 draws x 4 chains).
# Budget: 2000 burn-in + 500 kept, then extend with `more` (500 more each time, NBURN=0 since the
# chain is already burnt in). Every run streams its draws AND a per-chain final_state, and a resumed
# chain restarts within one draw's noise of where it stopped -- measured 91 nats against an 83-nat
# draw-to-draw floor, with per-chain identity preserved 4/4 on the parameter state. So there is no
# reason to guess the budget up front: take 500, look at the diagnostics, extend if kappa is ragged.
fit()  { echo ">>> fit: flat MNL, no BART, RE-ASIS on"
         BM_INPUT="$DUMP" BM_RE_ASIS="${BM_RE_ASIS:-TRUE}" \
         BM_SLAB_C2="${BM_SLAB_C2:-FALSE}" BM_SLAB_C2_VAL="${BM_SLAB_C2_VAL:-14.69}" \
         BM_NITER="${BM_NITER:-5000}" BM_NBURN="${BM_NBURN:-4000}" BM_THIN="${BM_THIN:-1}" \
         Rscript "$P/estimate_prior.R" "$@"; }
# Continue sampling from a finished segment's saved chain state. Default target is the newest
# segment, so repeated `more` calls chain forward one segment at a time.
more() {
  local from="${1:-}"
  if [ -z "$from" ]; then
    # Pick the newest RESUMABLE segment, not merely the newest directory. An aborted run leaves an
    # empty folder behind, and by name that husk sorts newest -- which is how `more` came to pick a
    # directory with no chain state at all while a finished segment sat right beside it.
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
  BM_INPUT="$DUMP" BM_RESUME="$from" BM_RE_ASIS="${BM_RE_ASIS:-TRUE}" \
    BM_SLAB_C2="${BM_SLAB_C2:-FALSE}" BM_SLAB_C2_VAL="${BM_SLAB_C2_VAL:-14.69}" \
    BM_NITER="${BM_NITER:-500}" BM_THIN="${BM_THIN:-1}" Rscript "$P/estimate_prior.R"
}

case "${1:-all}" in
  design) design ;;
  fit)    fit ;;
  more)   more "${2:-}" ;;
  smoke)  design; fit --smoke ;;
  all)    design; fit ;;
  *) echo "usage: $0 [design|fit|more [dir]|smoke|all]"; exit 1 ;;
esac
