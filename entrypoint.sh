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

# The task list, named ONCE. It used to be spelled out again inside the passthrough test, which
# silently omitted all four bmleh* tasks -- harmless only because `command -v bmleh` happens to fail.
# The knob registry, generated from config/knobs.json (tools/gen_config.py).
[ -f config/knobs.generated.sh ] && . config/knobs.generated.sh

GAMBLE_TASKS="bmleh bmleh_all bmleh_smoke bmleh_design bmleh_fit nested flat_design flat_fit count report test"
is_task() { case " $GAMBLE_TASKS " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }

# =============================================================================
# EVERY SETTING, AND WHERE IT CAME FROM
# =============================================================================
# The recurring failure in this pipeline is not a wrong value, it is a value that exists in one
# layer and never reaches the next: a routine schema default that never entered the container, an
# OUTPUT_DIR nothing read, a design dump chosen by filename. A log that prints only the VALUE
# cannot tell "you asked for this" from "nothing arrived and I fell back" -- which is exactly the
# distinction that took two failed jobs to establish.
#
# So every knob is recorded WITH ITS PROVENANCE, and written out as a manifest before any work
# starts. One artifact then answers "what did this run actually use", instead of reading a log and
# hoping.
KNOB_NAMES=(); KNOB_VALUES=(); KNOB_SOURCES=()
knob() {                       # knob NAME DEFAULT  -- env wins, otherwise the default
  local n="$1" d="$2" v s
  v="${!n-}"
  if [ -n "$v" ]; then s="environment"; else v="$d"; s="default"; printf -v "$n" '%s' "$v"; fi
  KNOB_NAMES+=("$n"); KNOB_VALUES+=("$v"); KNOB_SOURCES+=("$s")
}
knob_record() {                # knob_record NAME VALUE SOURCE -- for knobs resolved specially
  KNOB_NAMES+=("$1"); KNOB_VALUES+=("$2"); KNOB_SOURCES+=("$3")
}
knob_update() {                # a knob decided LATER (the design auto-detect) corrects its record,
  local i                      # so the manifest shows what was used, not what was first guessed
  for i in "${!KNOB_NAMES[@]}"; do
    if [ "${KNOB_NAMES[$i]}" = "$1" ]; then KNOB_VALUES[$i]="$2"; KNOB_SOURCES[$i]="$3"; return 0; fi
  done
  knob_record "$1" "$2" "$3"
}
json_escape() { printf '%s' "$1" | sed 's/\\/\\\\/g; s/"/\\"/g'; }
write_manifest() {
  local out="$1" i n v src sep=""
  { printf '{\n'
    printf '  "generated_at": "%s",\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    printf '  "task": "%s",\n' "$(json_escape "$TASK")"
    printf '  "task_source": "%s",\n' "$(json_escape "$TASK_SRC")"
    printf '  "git_commit": "%s",\n' "$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
    printf '  "settings": {\n'
    for i in "${!KNOB_NAMES[@]}"; do
      n="${KNOB_NAMES[$i]}"; v="$(json_escape "${KNOB_VALUES[$i]}")"; src="${KNOB_SOURCES[$i]}"
      printf '%s    "%s": { "value": "%s", "source": "%s" }' "$sep" "$n" "$v" "$src"; sep=$',\n'
    done
    printf '\n  }\n}\n'
  } > "$out"
}
# ---- registry + profiles -----------------------------------------------------------------------
# NO ASSOCIATIVE ARRAYS: bash 3.2 rejects `declare -A`, and then silently treats t[key] as index 0,
# so a lookup APPEARS to work while returning the wrong entry. Prefixed variables instead.
PROFILE_KEYS=""
profile_set()  { printf -v "PROFILE_VAL_$1" '%s' "$2"; PROFILE_KEYS="$PROFILE_KEYS $1"; }
profile_has()  { case " $PROFILE_KEYS " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }
profile_get()  { eval "printf '%s' \"\${PROFILE_VAL_$1}\""; }
in_registry()  { local x; for x in "${KNOB_REGISTRY[@]}"; do [ "$x" = "$1" ] && return 0; done; return 1; }

load_profile() {
  local pf line k v
  [ -n "$PROFILE" ] || return 0
  pf="config/profiles/${PROFILE}.env"
  if [ ! -f "$pf" ]; then
    echo "ERROR: no profile '$PROFILE' ($pf does not exist)."
    echo "  Available: $(ls config/profiles/*.env 2>/dev/null | xargs -n1 basename 2>/dev/null | sed 's/\.env$//' | tr '\n' ' ')"
    exit 1
  fi
  while IFS= read -r line || [ -n "$line" ]; do
    case "$line" in ''|\#*) continue ;; esac
    k="${line%%=*}"; v="${line#*=}"
    k="$(printf '%s' "$k" | tr -d '[:space:]')"
    if ! in_registry "$k"; then
      echo "ERROR: profile '$PROFILE' sets '$k', which is not declared in config/knobs.json."
      echo "       Rejected rather than ignored -- a setting nobody reads is precisely the failure"
      echo "       this mechanism exists to prevent."
      exit 1
    fi
    profile_set "$k" "$v"
  done < "$pf"
  echo ">>> profile: $PROFILE  ($pf)"
}

# Resolve every registry knob: environment > profile > declared default. Each is exported under the
# name the R drivers actually read, and a model-scope knob moved off its default is announced --
# those defaults are the validated standard, so deviating is a modelling claim, not a preference.
SPEC_DEVIATIONS=()
resolve_registry() {
  local kn def envvar scope cur src
  for kn in "${KNOB_REGISTRY[@]}"; do
    eval "def=\$KNOB_DEFAULT_$kn"; eval "envvar=\$KNOB_ENV_$kn"; eval "scope=\$KNOB_SCOPE_$kn"
    cur="${!kn-}"
    if   [ -n "$cur" ];         then src="environment"
    elif extra_has "$kn";       then cur="$(extra_get "$kn")";   src="EXTRA override"
    elif profile_has "$kn";     then cur="$(profile_get "$kn")"; src="profile:$PROFILE"
    else cur="$def";                 src="default"
    fi
    printf -v "$kn" '%s' "$cur"
    # A knob resting on its BLANK SENTINEL means "not set". The routine form marks every declared
    # property required and refuses an empty string, so an optional knob needs a value the form can
    # hold; not exporting it here is what makes the driver fall back to its own default, exactly as
    # if the field had been left empty.
    eval "blank=\$KNOB_BLANK_$kn"
    if [ -n "$blank" ] && [ "$cur" = "$blank" ]; then
      src="$src (unset: '$blank')"
      KNOB_NAMES+=("$kn"); KNOB_VALUES+=("$cur"); KNOB_SOURCES+=("$src")
      continue
    fi
    KNOB_NAMES+=("$kn"); KNOB_VALUES+=("$cur"); KNOB_SOURCES+=("$src")
    [ -n "$cur" ] && export "$envvar=$cur"
    if [ "$scope" = "model" ] && [ "$cur" != "$def" ]; then
      SPEC_DEVIATIONS+=("$kn: $def -> $cur   ($src)")
    fi
  done
}
print_spec_deviations() {
  [ ${#SPEC_DEVIATIONS[@]} -eq 0 ] && { echo " Model spec: STANDARD (every model knob at its validated default)"; return 0; }
  echo "----------------------------------------------------------------------"
  echo " MODEL SPEC DEVIATIONS -- this run does not estimate the standard model."
  echo " Results are not comparable with runs at the defaults."
  local d; for d in "${SPEC_DEVIATIONS[@]}"; do echo "   $d"; done
  echo "----------------------------------------------------------------------"
}

print_settings() {
  local i n v src
  echo "----------------------------------------------------------------------"
  echo " Settings (value <- where it came from)"
  for i in "${!KNOB_NAMES[@]}"; do
    n="${KNOB_NAMES[$i]}"; v="${KNOB_VALUES[$i]}"; src="${KNOB_SOURCES[$i]}"
    printf "   %-22s %-46s <- %s\n" "$n" "${v:-(empty)}" "$src"
  done
  echo "----------------------------------------------------------------------"
}

# PROJECT TASKS ARE DISCOVERED, NOT ENUMERATED. Any projects/<NAME>/gamble_model/run.sh yields
# <NAME>_design, <NAME>_fit, <NAME>_smoke, <NAME>_all and <NAME>_more, so the task list never has to
# be kept in step with the projects directory by hand. GLOBIOM_subclass has had a runnable run.sh
# all along and no task, purely because nobody added one.
PROJECT_ACTIONS="design fit smoke all more"
gamble_project_tasks() {
  for d in projects/*/; do
    [ -f "${d}gamble_model/run.sh" ] || continue
    n=$(basename "$d")
    for a in $PROJECT_ACTIONS; do printf '%s_%s ' "$n" "$a"; done
  done
}
# Splits <PROJECT>_<action> into PROJ_NAME / PROJ_ACTION, and ONLY when that project really exists
# -- so `flat_design`, which also ends in _design, falls through to the core task list.
resolve_project_task() {
  for a in $PROJECT_ACTIONS; do
    case "$1" in
      *_"$a")
        n="${1%_$a}"
        if [ -f "projects/$n/gamble_model/run.sh" ]; then PROJ_NAME="$n"; PROJ_ACTION="$a"; return 0; fi
        ;;
    esac
  done
  return 1
}

# Allow arbitrary command passthrough if specified (e.g. bash, Rscript custom.R)
if [ $# -gt 0 ] && ! is_task "$1" && command -v "$1" >/dev/null 2>&1; then
  exec "$@"
fi

# -----------------------------------------------------------------------------
# Configuration mapping (Environment Variables -> Driver Knobs)
# -----------------------------------------------------------------------------
# WHERE THE TASK CAME FROM, not just what it is. A schema default only helps if the platform
# actually injects it: job 8721 ran `nested` while the routine's schema said `bmleh`, because
# neither $1 nor $TASK was set and the fallback won without saying so. "Task: nested" on its own
# cannot tell "you asked for nested" apart from "nothing reached me". Now it can.
if   [ $# -gt 0 ] && [ -n "${1:-}" ]; then TASK="$1";           TASK_SRC="command argument"
elif [ -n "${TASK:-}" ];             then                       TASK_SRC="TASK env"
elif [ -n "${task:-}" ];             then TASK="$task";         TASK_SRC="task env (lowercase)"
elif [ -n "${WORKFLOW:-}" ];         then TASK="$WORKFLOW";     TASK_SRC="WORKFLOW env"
elif [ -n "${workflow:-}" ];         then TASK="$workflow";     TASK_SRC="workflow env"
elif [ -n "${GAMBLE_TASK:-}" ];      then TASK="$GAMBLE_TASK";  TASK_SRC="GAMBLE_TASK env"
else
  # NO SILENT FALLBACK. This used to default to `nested`, which meant a routine whose task never
  # reached the container quietly started a different, long-running workflow -- twice, on jobs 8721
  # and 8733, the second of which got as far as fitting a BMLEH design through a generically-coded
  # GLOBIOM nest tree. A misconfiguration should cost a second, not a day of cores.
  echo "------------------------------------------------------------------"
  echo "ERROR: no task given, so there is nothing to run."
  echo
  echo "Nothing arrived as a command argument, or in TASK / task / WORKFLOW /"
  echo "workflow / GAMBLE_TASK. A default set in a routine's config schema is"
  echo "NOT enough -- the platform has to pass the value into the container:"
  echo "    command:  ./entrypoint.sh BMLEH_Los1_CAPRI_smoke"
  echo "    or env:   TASK=BMLEH_Los1_CAPRI_smoke"
  echo
  echo "  Core tasks:    $GAMBLE_TASKS"
  echo "  Project tasks: $(gamble_project_tasks)"
  echo "------------------------------------------------------------------"
  exit 1
fi

# The bmleh* spellings predate discovery. Keep them working, and resolve them to the explicit name
# so a config pinned to the old form neither breaks nor stays vague about which project it means.
case "$TASK" in
  bmleh|bmleh_all) TASK="BMLEH_Los1_CAPRI_all" ;;
  bmleh_smoke)     TASK="BMLEH_Los1_CAPRI_smoke" ;;
  bmleh_design)    TASK="BMLEH_Los1_CAPRI_design" ;;
  bmleh_fit)       TASK="BMLEH_Los1_CAPRI_fit" ;;
esac

if ! is_task "$TASK" && ! resolve_project_task "$TASK"; then
  echo "ERROR: unknown task '$TASK' (from $TASK_SRC)."
  echo "  Core tasks:    $GAMBLE_TASKS"
  echo "  Project tasks: $(gamble_project_tasks)"
  exit 1
fi
knob_record TASK "$TASK" "$TASK_SRC"
knob EXTRA "none"
# EXTRA is the escape hatch that keeps the routine form small: "NITER=8000, USE_BART=TRUE" instead
# of declaring twenty-five fields that would each have to hold a value and would then outrank the
# profile they were meant to inherit from. A key that is not a declared knob is REJECTED -- typing
# NITTER=8000 into a free-text box and having it silently ignored is the worst of both worlds.
EXTRA_KEYS=""
extra_set() { printf -v "EXTRA_VAL_$1" '%s' "$2"; EXTRA_KEYS="$EXTRA_KEYS $1"; }
extra_has() { case " $EXTRA_KEYS " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }
extra_get() { eval "printf '%s' \"\${EXTRA_VAL_$1}\""; }
parse_extra() {
  local item k v
  [ -n "$EXTRA" ] && [ "$EXTRA" != "none" ] || return 0
  # A here-string, not a pipe: the loop must run in THIS shell or extra_set writes into a subshell
  # and every override is silently lost. It also supplies the trailing newline that `read` needs --
  # without it the LAST item of "A=1, B=2" is never read, which is how B=2 went missing.
  while IFS= read -r item || [ -n "$item" ]; do
    item="$(printf '%s' "$item" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
    [ -n "$item" ] || continue
    case "$item" in
      *=*) ;;
      *) echo "ERROR: EXTRA item '$item' is not KEY=VALUE."; exit 1 ;;
    esac
    k="${item%%=*}"; v="${item#*=}"
    k="$(printf '%s' "$k" | tr -d '[:space:]')"
    if ! in_registry "$k"; then
      echo "ERROR: EXTRA sets '$k', which is not a declared knob (config/knobs.json)."
      echo "  Declared: ${KNOB_REGISTRY[*]}"
      exit 1
    fi
    extra_set "$k" "$v"
    echo ">>> override: $k=$v  (EXTRA)"
  done <<< "$(printf '%s' "$EXTRA" | tr ',;' '\n\n')"
}

knob PROFILE "none"
# "none"/"auto" are how the form expresses "not set" -- normalise them back to empty before use.
[ "$PROFILE" = "none" ] && PROFILE=""
load_profile
parse_extra
resolve_registry
# PROJECT names the project-specific build (e.g. BMLEH_Los1_CAPRI). It is the FIRST thing
# consulted when picking a staged design dump, because a project build and a classification
# are not the same axis: pixel_model_inputs_BMLEH_Los1_CAPRI.rds is selected by project, while
# pixel_model_inputs_GLOBIOM_subclass.rds is selected by classification.
# "auto" is how the form says "not set" -- see the blank-sentinel note above.
knob PROJECT "auto"; [ "$PROJECT" = "auto" ] && PROJECT=""
knob DESIGN_PATH "auto"
[ "$DESIGN_PATH" = "auto" ] && DESIGN_PATH="output/designs/pixel_model_inputs.rds"
knob OUTPUT_DIR "output"

# --- Redirect output/ and results/ to the mounted drive -------------------------
# The platform uploads only what appears under the FUSE mount. Two directories
# must be redirected:
#   output/  — design dumps, intermediate files, reports
#   results/ — posterior batches and chain states (estimate_prior.R writes here)
# Without the results/ symlink every posterior draw lands in the container's
# ephemeral filesystem and is lost at job end ("No files found to upload").
knob GAMBLE_WORK_DIR "/mnt/wdrv"
if [ -d "$GAMBLE_WORK_DIR" ]; then
  # A bare `ln` failure under `set -e` kills the job with one cryptic line. Say what is actually
  # wrong: the working directory has to be writable by the RUNTIME user, and a root-owned /app in a
  # container started as uid 1000 is not.
  redirect_dir() {                       # redirect_dir <name>
    local n="$1"
    [ -L "$n" ] && return 0
    mkdir -p "$GAMBLE_WORK_DIR/$n"
    [ -d "$n" ] && cp -a "$n/." "$GAMBLE_WORK_DIR/$n/" 2>/dev/null || true
    rm -rf "$n" 2>/dev/null || true
    if ln -s "$GAMBLE_WORK_DIR/$n" "$n" 2>/dev/null; then
      echo ">>> $n/  -> $GAMBLE_WORK_DIR/$n"
    else
      echo "------------------------------------------------------------------"
      echo "ERROR: cannot redirect $n/ to the mounted drive."
      echo "  $PWD is not writable by this user (uid $(id -u)), so the symlink cannot be created."
      echo "  Everything the model writes would stay in the container and be lost at job end."
      echo "  The image must make its working directory writable: the Dockerfile does"
      echo "  'chmod -R a+rwX /app' after COPY. Rebuild if that is missing."
      echo "------------------------------------------------------------------"
      exit 1
    fi
  }
  redirect_dir output
  redirect_dir results
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
      knob_update DESIGN_PATH "$DESIGN_PATH" "auto-detected, matched $sel"
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

# Written BEFORE any work: a job that dies mid-run still leaves a complete record of how it was
# configured. OUTPUT_DIR is the symlink to the mounted drive by this point, so it is uploaded.
print_settings
print_spec_deviations
write_manifest "$OUTPUT_DIR/run_manifest.json"
echo ">>> settings manifest: $OUTPUT_DIR/run_manifest.json"

# DRY_RUN resolves and reports the configuration, then stops before dispatch. It answers "what
# would this routine actually do" without spending a job to find out -- which, for a workflow whose
# cheapest task is minutes and whose dearest is days, is the difference between checking a config
# and gambling one.
if [ "${DRY_RUN:-FALSE}" = "TRUE" ] || [ "${DRY_RUN:-false}" = "true" ]; then
  echo ">>> DRY_RUN: configuration resolved and written; stopping before the task runs."
  exit 0
fi
mkdir -p "$OUTPUT_DIR/designs"
mkdir -p "$OUTPUT_DIR/report"

echo "======================================================================"
echo " Starting gamble-core Routine"
echo " Time:            $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
echo " Task:            $TASK   (from $TASK_SRC)"
echo " Output Dir:      $OUTPUT_DIR"
echo "======================================================================"

# Master Parquet resolution (env override -> /data auto-detect)
[ "${GAMBLE_MASTER_PARQUET:-}" = "auto" ] && GAMBLE_MASTER_PARQUET=""
[ "${GAMBLE_CASCADE_DATA:-}" = "auto" ] && GAMBLE_CASCADE_DATA=""
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
require_project() {
  BM_RUN="$1"
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

  *)
    if resolve_project_task "$TASK"; then
      PROJ_RUN="projects/$PROJ_NAME/gamble_model/run.sh"
      echo ">>> Task: $PROJ_NAME project model ($PROJ_ACTION)"
      require_project "$PROJ_RUN"
      case "$PROJ_ACTION" in
        fit|all|more) export BM_NITER="${BM_NITER:-${NITER:-5000}}" ;;
      esac
      "./$PROJ_RUN" "$PROJ_ACTION"
    else
      echo "ERROR: Unknown TASK: '$TASK' (from $TASK_SRC)"
      echo "  Core tasks:    $GAMBLE_TASKS"
      echo "  Project tasks: $(gamble_project_tasks)"
      exit 1
    fi
    ;;
esac

echo "======================================================================"
echo " gamble-core Routine finished successfully at $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
echo "======================================================================"
