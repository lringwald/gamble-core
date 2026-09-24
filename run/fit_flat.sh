#!/usr/bin/env bash
# =============================================================================
# run/fit_flat.sh — batch runner for the flat pixel land-use model
# =============================================================================
# Usage: run/fit_flat.sh [smoke|prod] [linear|bart] [design.rds]
# Default:               prod         linear        (auto-detected)
#
# Examples:
#   run/fit_flat.sh smoke
#   run/fit_flat.sh prod linear output/designs/pixel_model_inputs.rds
#   run/fit_flat.sh prod bart
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."

MODE="${1:-prod}"
BART_MODE="${2:-linear}"
DESIGN="${3:-${DESIGN_PATH:-${GB_INPUT:-${BM_INPUT:-}}}}"

if [ -z "$DESIGN" ] || [ ! -f "$DESIGN" ]; then
  # Auto-detect latest staged design dump
  cand=$(ls -t output/designs/pixel_model_inputs*.rds output/pixel_model_inputs*.rds 2>/dev/null | head -1 || true)
  if [ -n "$cand" ] && [ -f "$cand" ]; then
    DESIGN="$cand"
    echo ">>> Auto-detected design dump: $DESIGN"
  else
    echo "FATAL: No design dump found. Provide a path or build one with run/flat.R (or TASK=flat_design)."
    exit 1
  fi
fi

case "$MODE" in
  smoke) NITER=600;  NBURN=200;  NCHAINS=1; NPIX=5000; RUN_MODE=smoke ;;
  prod)  NITER=6000; NBURN=2000; NCHAINS=4; NPIX=0;    RUN_MODE=production ;;
  *) echo "FATAL: mode must be smoke|prod"; exit 1 ;;
esac

case "$BART_MODE" in
  linear) USE_BART=FALSE ;;
  bart)   USE_BART=TRUE ;;
  *) echo "FATAL: 2nd argument must be linear|bart"; exit 1 ;;
esac

TS=$(date +%Y%m%d_%H%M)
TAGS="flat_${MODE}_${BART_MODE}"
LOG="output/${TAGS}_${TS}.log"
mkdir -p output results

cat <<EOF
>>> flat pixel land-use model
    mode      : $MODE (niter=$NITER, nburn=$NBURN, chains=$NCHAINS, subsample=$NPIX)
    terrain   : $BART_MODE (use_bart=$USE_BART)
    design    : $DESIGN ($(date -r "$DESIGN" '+%Y-%m-%d %H:%M'))
    log       : $LOG
EOF

DESIGN_PATH="$DESIGN" NITER="$NITER" NBURN="$NBURN" N_CHAINS="$NCHAINS" \
SUBSAMPLE="$NPIX" USE_BART="$USE_BART" RUN_MODE="$RUN_MODE" \
nohup Rscript drivers/run_flat_fit.R > "$LOG" 2>&1 &
echo ">>> pid $!   tail -f \"$LOG\""
