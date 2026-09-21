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

# --- results must land on the MOUNTED DRIVE, not in the container ---------------------
# The platform uploads only what appears under the FUSE mount; a run that writes into the
# repo's own output/ finishes with "No files created ... to upload" and looks like a job
# that produced nothing.
#
# OUTPUT_DIR CANNOT DO THIS. No R file reads it -- every driver writes to a RELATIVE
# "output/..." path (output/saved_model_outputs, output/designs, output/intermediate,
# output/plots, output/diagnostics), about 25 of them, hardcoded. So the redirection has to
# happen at the FILESYSTEM: give the name `output` to the drive and the existing paths
# follow, with no change to any R code.
GAMBLE_WORK_DIR="${GAMBLE_WORK_DIR:-/mnt/wdrv}"
if [ -d "$GAMBLE_WORK_DIR" ] && [ ! -L output ]; then
  mkdir -p "$GAMBLE_WORK_DIR/output"
  # carry over whatever the checkout shipped (in practice just .gitkeep) before handing the
  # name over, so nothing tracked is lost
  [ -d output ] && cp -a output/. "$GAMBLE_WORK_DIR/output/" 2>/dev/null || true
  rm -rf output
  ln -s "$GAMBLE_WORK_DIR/output" output
  echo ">>> output/ -> $GAMBLE_WORK_DIR/output  (results land on the mounted drive)"
fi

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

# Master Parquet resolution (env override -> /data auto-detect)
if [ -n "${GAMBLE_MASTER_PARQUET:-}" ]; then
  echo ">>> Master Parquet: $GAMBLE_MASTER_PARQUET"
  export GAMBLE_MASTER_PARQUET
elif [ -f "/data/prior_model_1km_master_inputs.parquet" ]; then
  echo ">>> Auto-detected Master Parquet: /data/prior_model_1km_master_inputs.parquet"
  export GAMBLE_MASTER_PARQUET="/data/prior_model_1km_master_inputs.parquet"
fi

# Cascade data resolution (/data auto-detect)
if [ -z "${GAMBLE_CASCADE_DATA:-}" ] && [ -d "/data/cascadinggamble-core/data" ]; then
  export GAMBLE_CASCADE_DATA="/data/cascadinggamble-core/data"
fi

# Design-dump resolution (explicit DESIGN_PATH -> drive/-/data auto-detect). A cold container
# has no design: output/ ships with nothing but .gitkeep, so TASK=nested hit its preflight and
# exited before doing any work. If a dump is already staged on the drive, use it rather than
# making the caller spell out the path.
if [ ! -f "$DESIGN_PATH" ]; then
  for cand_dir in "$GAMBLE_WORK_DIR" /data; do
    [ -d "$cand_dir" ] || continue
    cand=$(find "$cand_dir" -maxdepth 3 -name 'pixel_model_inputs*.rds' 2>/dev/null | sort | head -1)
    if [ -n "${cand:-}" ]; then
      echo ">>> Auto-detected design dump: $cand"
      DESIGN_PATH="$cand"; break
    fi
  done
fi

# The BMLEH project scripts ARE tracked: they were force-added in d93583e (2026-09-19) and
# .gitignore's `projects/` rule does not apply to files git already tracks. So if run.sh is
# missing here, the CODE is not the problem -- this checkout or image predates that commit,
# or the build context dropped it. The old message here said the opposite and sent a debug
# session after a non-existent commit step; say the checkable thing instead.
BM_RUN=projects/BMLEH_Los1_CAPRI/gamble_model/run.sh
require_bmleh() {
  [ -f "$BM_RUN" ] && { chmod +x "$BM_RUN"; return 0; }
  echo "------------------------------------------------------------------"
  echo "ERROR: $BM_RUN not found in this image."
  echo
  echo "These scripts ARE committed (d93583e, 2026-09-19) even though .gitignore lists"
  echo "'projects/' -- ignore rules do not apply to already-tracked files. Do NOT re-add them."
  echo
  echo "Check instead which commit this image was built from:"
  echo "    git -C . log -1 --oneline    # if .git is present"
  echo "    ls projects/                 # empty => built before d93583e, or context-filtered"
  echo "If it predates d93583e, rebuild the image. If you build locally, confirm"
  echo ".dockerignore does not exclude projects/."
  echo "------------------------------------------------------------------"
  exit 1
}

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
    require_bmleh
    export BM_NITER="${BM_NITER:-${NITER:-5000}}"
    ./projects/BMLEH_Los1_CAPRI/gamble_model/run.sh all
    ;;

  bmleh_smoke)
    echo ">>> Task: BMLEH_Los1_CAPRI project model (smoke run: design + 5k px fit)"
    require_bmleh
    ./projects/BMLEH_Los1_CAPRI/gamble_model/run.sh smoke
    ;;

  bmleh_design)
    echo ">>> Task: BMLEH_Los1_CAPRI (design only)"
    require_bmleh
    ./projects/BMLEH_Los1_CAPRI/gamble_model/run.sh design
    ;;

  bmleh_fit)
    echo ">>> Task: BMLEH_Los1_CAPRI (fit only)"
    require_bmleh
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
