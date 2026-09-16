#!/usr/bin/env bash
# =============================================================================
# run/fit_nested.sh — batch runner for the nested land-use model
# =============================================================================
# Usage: run/fit_nested.sh [smoke|prod] [factorized|iv] [intercept|intercept+socio] [symhs|diag]
# Default:                  prod         factorized      intercept                   symhs
#
# symhs is the STANDARD (BMLEH_Los1, 2026-09-16); diag is the measurement arm.
# See docs/sampler_model_specification.md sec. 8.0.
#
#   smoke  6000 px / 200 iter  (~10-15 min)      prod  64,173 px / 1000 iter (~24 h)
#
# Same code path as run/nested.R -- this is the non-interactive form. Backgrounds itself with nohup
# and prints the log path; the node store makes it kill-safe and resumable.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."

MODE="${1:-prod}"; VARIANT="${2:-factorized}"; RE_COLS="${3:-intercept}"; HS="${4:-symhs}"
DESIGN="${NCUT_DESIGN:-output/pixel_model_inputs.rds}"

[ -f "$DESIGN" ] || { echo "FATAL: no design dump at $DESIGN (build one with run/flat.R)"; exit 1; }
case "$MODE"    in smoke) SUB=6000; NITER=200 ;; prod) SUB=0; NITER=1000 ;;
                   *) echo "FATAL: mode must be smoke|prod"; exit 1 ;; esac
case "$VARIANT" in factorized) USE_IV=FALSE ;; iv) USE_IV=TRUE ;;
                   *) echo "FATAL: variant must be factorized|iv"; exit 1 ;; esac
case "$RE_COLS" in intercept|intercept+socio) ;;
                   *) echo "FATAL: RE block must be intercept|intercept+socio"; exit 1 ;; esac
case "$HS" in
  diag)  SYM=FALSE ;;
  symhs) SYM=TRUE
         # symmetric_hs had two defects (baseline-dependent block kernel; stale mu_R in the complement
         # redraw), both fixed 2026-08-19. Without the fix it misallocates probability across nests.
         grep -q "kronecker(Msym" codes/mnlogit_rcpp_sym.R || {
           echo "FATAL: codes/mnlogit_rcpp_sym.R lacks the 2026-08-19 symmetric-HS fix"; exit 1; } ;;
  *) echo "FATAL: last argument must be diag|symhs"; exit 1 ;;
esac

TS=$(date +%Y%m%d_%H%M)
TAGS="${VARIANT}_${RE_COLS//+/_}_${HS}"
# The store fingerprints data/niter/M but NOT the sampler version, so the tag carries every knob that
# changes the fit -- never let two different configurations share one store path.
# NCUT_STORE_DIR overrides the derived path -- needed to RESUME a run that was started under a
# different naming convention (e.g. the 2026-08-20 symmetric run). Only ever point it at a store
# built by the SAME sampler version and the same knobs.
STORE="${NCUT_STORE_DIR:-output/ncut_store_${MODE}_${TAGS}}"
LOG="output/ncut_${MODE}_${TAGS}_${TS}.log"

cat <<EOF
>>> nested land-use model
    mode/variant : $MODE / $VARIANT (use_iv=$USE_IV)
    RE block     : $RE_COLS      shrinkage: $HS (symmetric_hs=$SYM)
    design       : $DESIGN ($(date -r "$DESIGN" '+%Y-%m-%d %H:%M'))
    store / log  : $STORE
                   $LOG
EOF

NCUT_USE_IV="$USE_IV" NCUT_RE_COLS="$RE_COLS" NCUT_SUB="$SUB" NCUT_SYM_HS="$SYM" \
NCUT_STORE_DIR="$STORE" \
nohup Rscript run_nested_cut.R GLOBIOM "$DESIGN" 25 "$NITER" TRUE auto TRUE 1 1 > "$LOG" 2>&1 &
echo ">>> pid $!   tail -f \"$LOG\""
