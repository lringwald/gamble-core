#!/usr/bin/env bash
# =============================================================================
# entrypoint.sh — routine entrypoint for gamble-core (DeepOrigin / WKUBE / Celery)
# =============================================================================
# Dispatches containerized tasks in the foreground using environment variables
# injected by the Routine Configuration Schema.
# =============================================================================
set -eo pipefail

# Ensure we are in the repository root
cd "$(dirname "$0")"

# Allow arbitrary command passthrough if specified (e.g. bash, Rscript custom.R)
if [ $# -gt 0 ] && command -v "$1" >/dev/null 2>&1 && [ "$1" != "nested" ] && [ "$1" != "flat_design" ] && [ "$1" != "flat_fit" ] && [ "$1" != "count" ] && [ "$1" != "report" ] && [ "$1" != "test" ]; then
  exec "$@"
fi

# -----------------------------------------------------------------------------
# Configuration mapping (Environment Variables -> Driver Knobs)
# -----------------------------------------------------------------------------
TASK="${1:-${TASK:-nested}}"
CLASSIFICATION="${CLASSIFICATION:-GLOBIOM_subclass}"
VARIANT="${VARIANT:-factorized}"
RE_BLOCK="${RE_BLOCK:-intercept}"
SYMMETRIC_HS="${SYMMETRIC_HS:-TRUE}"
NITER="${NITER:-1000}"
SUBSAMPLE="${SUBSAMPLE:-0}"
M="${M:-25}"
N_CHAINS="${N_CHAINS:-1}"
N_CORES="${N_CORES:-${N_CHAINS}}"
DESIGN_PATH="${DESIGN_PATH:-output/designs/pixel_model_inputs.rds}"
OUTPUT_DIR="${OUTPUT_DIR:-output}"

# Ensure output directories exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/designs"
mkdir -p "$OUTPUT_DIR/report"

echo "======================================================================"
echo " Starting gamble-core Routine"
echo " Time:            $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
echo " Task:            $TASK"
echo " Output Dir:      $OUTPUT_DIR"
echo "======================================================================"

# If GAMBLE_MASTER_PARQUET is passed, export it
if [ -n "${GAMBLE_MASTER_PARQUET:-}" ]; then
  echo ">>> Master Parquet: $GAMBLE_MASTER_PARQUET"
  export GAMBLE_MASTER_PARQUET
fi

case "$TASK" in
  nested)
    echo ">>> Task: NESTED land-use model"
    echo "    Variant:        $VARIANT"
    echo "    RE block:       $RE_BLOCK"
    echo "    Symmetric HS:   $SYMMETRIC_HS"
    echo "    Iterations:     $NITER"
    echo "    Subsample:      $SUBSAMPLE (0 = full)"
    echo "    Chains/Cores:   $N_CHAINS / $N_CORES"
    echo "    Design file:    $DESIGN_PATH"

    # Verify design file exists
    if [ ! -f "$DESIGN_PATH" ]; then
      echo "------------------------------------------------------------------"
      echo "ERROR: Design dump not found at: $DESIGN_PATH"
      echo "To build a design dump, run this routine with TASK=flat_design first,"
      echo "or mount a precomputed pixel_model_inputs.rds into the container."
      echo "------------------------------------------------------------------"
      exit 1
    fi

    # Sanity check for symmetric horseshoe fix
    if [ "$SYMMETRIC_HS" = "TRUE" ] || [ "$SYMMETRIC_HS" = "true" ]; then
      grep -q "kronecker(Msym" codes/mnlogit_rcpp_sym.R || {
        echo "ERROR: codes/mnlogit_rcpp_sym.R lacks the symmetric-HS fix."
        exit 1
      }
      SYM_VAL="TRUE"
    else
      SYM_VAL="FALSE"
    fi

    USE_IV_VAL=$([ "$VARIANT" = "iv" ] && echo "TRUE" || echo "FALSE")
    STORE_DIR="${OUTPUT_DIR}/ncut_store_${VARIANT}_${RE_BLOCK//+/_}"

    export NCUT_USE_IV="$USE_IV_VAL"
    export NCUT_RE_COLS="$RE_BLOCK"
    export NCUT_SUB="$SUBSAMPLE"
    export NCUT_SYM_HS="$SYM_VAL"
    export NCUT_STORE_DIR="$STORE_DIR"

    # Run in the FOREGROUND so WKUBE / Celery tracks process completion
    Rscript drivers/run_nested_cut.R \
      GLOBIOM \
      "$DESIGN_PATH" \
      "$M" \
      "$NITER" \
      TRUE \
      auto \
      TRUE \
      "$N_CHAINS" \
      "$N_CORES"
    ;;

  flat_design)
    echo ">>> Task: Assembling DESIGN dump only"
    echo "    Classification: $CLASSIFICATION"
    export DRIVER_CLASS_COLS="$CLASSIFICATION"
    export DRIVER_PROMOTE_NATURAL_OTHER="TRUE"
    export DRIVER_DUMP_INPUTS="TRUE"
    export DRIVER_DUMP_EXIT="TRUE"

    Rscript drivers/run_lu_pixel_model.R
    ;;

  flat_fit)
    echo ">>> Task: Fitting FLAT pixel model (no nesting)"
    echo "    Classification: $CLASSIFICATION"
    echo "    Iterations:     $NITER"
    echo "    Chains:         $N_CHAINS"
    export DRIVER_CLASS_COLS="$CLASSIFICATION"
    export DRIVER_PROMOTE_NATURAL_OTHER="TRUE"
    export DRIVER_NITER="$NITER"
    export DRIVER_NCHAINS="$N_CHAINS"

    Rscript drivers/run_lu_pixel_model.R
    ;;

  count)
    echo ">>> Task: Fitting livestock COUNT model (BOV + SGT)"
    echo "    Iterations:     $NITER"
    export DRIVER_RUN_MODE="production"
    export DRIVER_NITER="$NITER"

    Rscript drivers/run_ls_count_model.R
    ;;

  report)
    echo ">>> Task: Generating fit report"
    Rscript postprocess/nested_report.R
    ;;

  test)
    echo ">>> Task: Running sampler gates and test suites"
    Rscript tests/run_all.R fast
    ;;

  bmleh|bmleh_all)
    echo ">>> Task: BMLEH_Los1_CAPRI project model (design + fit)"
    if [ ! -f "projects/BMLEH_Los1_CAPRI/gamble_model/run.sh" ]; then
      echo "------------------------------------------------------------------"
      echo "ERROR: projects/BMLEH_Los1_CAPRI/gamble_model/run.sh not found!"
      echo "Note: 'projects/' is git-ignored by default. To run in a remote routine,"
      echo "ensure 'projects/BMLEH_Los1_CAPRI' is committed to your branch or mounted."
      echo "------------------------------------------------------------------"
      exit 1
    fi
    chmod +x projects/BMLEH_Los1_CAPRI/gamble_model/run.sh
    export BM_NITER="${BM_NITER:-${NITER:-5000}}"
    ./projects/BMLEH_Los1_CAPRI/gamble_model/run.sh all
    ;;

  bmleh_smoke)
    echo ">>> Task: BMLEH_Los1_CAPRI project model (smoke run: design + 5k px fit)"
    if [ ! -f "projects/BMLEH_Los1_CAPRI/gamble_model/run.sh" ]; then
      echo "ERROR: projects/BMLEH_Los1_CAPRI/gamble_model/run.sh not found!"
      exit 1
    fi
    chmod +x projects/BMLEH_Los1_CAPRI/gamble_model/run.sh
    ./projects/BMLEH_Los1_CAPRI/gamble_model/run.sh smoke
    ;;

  bmleh_design)
    echo ">>> Task: BMLEH_Los1_CAPRI (design only)"
    chmod +x projects/BMLEH_Los1_CAPRI/gamble_model/run.sh
    ./projects/BMLEH_Los1_CAPRI/gamble_model/run.sh design
    ;;

  bmleh_fit)
    echo ">>> Task: BMLEH_Los1_CAPRI (fit only)"
    chmod +x projects/BMLEH_Los1_CAPRI/gamble_model/run.sh
    export BM_NITER="${BM_NITER:-${NITER:-5000}}"
    ./projects/BMLEH_Los1_CAPRI/gamble_model/run.sh fit
    ;;

  *)
    echo "ERROR: Unknown TASK: '$TASK'"
    echo "Supported tasks: nested | flat_design | flat_fit | count | report | test | bmleh | bmleh_smoke | bmleh_design | bmleh_fit"
    exit 1
    ;;
esac

echo "======================================================================"
echo " gamble-core Routine finished successfully at $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
echo "======================================================================"
