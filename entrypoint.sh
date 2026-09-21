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
# PROJECT names the project-specific build (e.g. BMLEH_Los1_CAPRI). It is the FIRST thing
# consulted when picking a staged design dump, because a project build and a classification
# are not the same axis: pixel_model_inputs_BMLEH_Los1_CAPRI.rds is selected by project, while
# pixel_model_inputs_GLOBIOM_subclass.rds is selected by classification.
PROJECT="${PROJECT:-}"
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

# --- Redirect output/ and results/ to the mounted drive -------------------------
# The platform uploads only what appears under the FUSE mount. Two directories
# must be redirected:
#   output/  — design dumps, intermediate files, reports
#   results/ — posterior batches and chain states (estimate_prior.R writes here)
# Without the results/ symlink every posterior draw lands in the container's
# ephemeral filesystem and is lost at job end ("No files found to upload").
GAMBLE_WORK_DIR="${GAMBLE_WORK_DIR:-/mnt/wdrv}"
if [ -d "$GAMBLE_WORK_DIR" ]; then
  if [ ! -L output ]; then
    mkdir -p "$GAMBLE_WORK_DIR/output"
    [ -d output ] && cp -a output/. "$GAMBLE_WORK_DIR/output/" 2>/dev/null || true
    rm -rf output
    ln -s "$GAMBLE_WORK_DIR/output" output
    echo ">>> output/  -> $GAMBLE_WORK_DIR/output"
  fi
  if [ ! -L results ]; then
    mkdir -p "$GAMBLE_WORK_DIR/results"
    [ -d results ] && cp -a results/. "$GAMBLE_WORK_DIR/results/" 2>/dev/null || true
    rm -rf results
    ln -s "$GAMBLE_WORK_DIR/results" results
    echo ">>> results/ -> $GAMBLE_WORK_DIR/results"
  fi
  echo ">>> Mounted drive: $GAMBLE_WORK_DIR  (all model output lands here)"
fi

# --- Design-dump resolution -----------------------------------------------------
# A cold container has no design: output/ ships with nothing but .gitkeep, and TASK=nested is the
# routine's default, so the preflight fires before any work happens. If a dump is already staged on
# the drive, use it rather than making the caller spell the path out.
# (This was written once and lost in a merge on 2026-09-21; job 8721 failed for exactly that.)
if [ ! -f "$DESIGN_PATH" ]; then
  # The trailing `true` is load-bearing. This script runs under `set -eo pipefail`, and without it
  # the loop's exit status is that of its last test -- `[ -d /data ]`, which is FALSE whenever /data
  # is not mounted. pipefail carries that through the pipe, the assignment fails, and the routine
  # dies right here with no message. `bash -n` accepts it either way.
  cands=$(for cand_dir in "$GAMBLE_WORK_DIR" /data; do
            [ -d "$cand_dir" ] && find "$cand_dir" -maxdepth 4 -name 'pixel_model_inputs*.rds' 2>/dev/null
            true
          done | sort)
  n_cand=$(printf '%s' "$cands" | grep -c . || true)
  if [ "${n_cand:-0}" -eq 1 ]; then
    DESIGN_PATH="$cands"
    echo ">>> Auto-detected design dump: $DESIGN_PATH"
  elif [ "${n_cand:-0}" -gt 1 ]; then
    # NEVER GUESS BETWEEN DESIGNS. Alphabetical order would hand a GLOBIOM run the BMLEH dump and
    # fit it without a word -- a wrong answer is worse than a failed job. Match the classification
    # the task actually asked for; if that is ambiguous too, stop and show the candidates.
    # PROJECT first, then CLASSIFICATION -- the two selectors name different axes, and a project
    # build (pixel_model_inputs_BMLEH_Los1_CAPRI.rds) does not carry a classification in its name.
    match=""
    [ -n "$PROJECT" ] && match=$(printf '%s\n' "$cands" | grep -F "$PROJECT" | head -1 || true)
    sel="PROJECT=$PROJECT"
    if [ -z "$match" ]; then
      match=$(printf '%s\n' "$cands" | grep -F "$CLASSIFICATION" | head -1 || true)
      sel="CLASSIFICATION=$CLASSIFICATION"
    fi
    if [ -n "$match" ]; then
      DESIGN_PATH="$match"
      echo ">>> Auto-detected design dump for $sel: $DESIGN_PATH"
      echo ">>>   ($n_cand dumps on the drive; the others were ignored)"
    else
      echo "------------------------------------------------------------------"
      echo "ERROR: $n_cand design dumps are staged and none names PROJECT='$PROJECT' or"
      echo "       CLASSIFICATION='$CLASSIFICATION':"
      printf '%s\n' "$cands" | sed 's/^/         /'
      echo "Set PROJECT, or DESIGN_PATH explicitly. Refusing to guess: picking the"
      echo "wrong design fits a different model and reports success."
      echo "------------------------------------------------------------------"
      exit 1
    fi
  fi
  # Say what was SEARCHED and what is actually there. Neither the routine author nor anyone reading
  # this log can list the mount by hand, so a bare "not found" ends the investigation; the .rds files
  # that DO exist usually name the problem (wrong directory, or a project-specific filename).
  if [ ! -f "$DESIGN_PATH" ]; then
    echo ">>> No pixel_model_inputs*.rds under $GAMBLE_WORK_DIR or /data (searched 4 levels)."
    for cand_dir in "$GAMBLE_WORK_DIR" /data; do
      [ -d "$cand_dir" ] || continue
      echo ">>>   .rds files present under $cand_dir:"
      find "$cand_dir" -maxdepth 4 -name '*.rds' 2>/dev/null | head -10 | sed 's/^/>>>     /'
      [ -z "$(find "$cand_dir" -maxdepth 4 -name '*.rds' 2>/dev/null | head -1)" ] && echo ">>>     (none)"
    done
  fi
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
