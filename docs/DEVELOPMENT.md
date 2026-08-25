# Development guide — gamble-core

How to develop and run the GAMBLE prior models from this repo. Run everything **from the repo root**
(all paths are root-relative). Outputs land in `output/` (git-ignored, regenerable). Model theory and
equations for the nested model live in [`nested_cut_model.md`](nested_cut_model.md).

---

## 1. The development loop — two entry points

The models need the pixel **design** (outcomes + covariates). Building it from raw 1 km land-use maps is a
heavy step with an **external data dependency**; model development does not repeat it. Pick your entry point:

### Where the raw data lives (important)
The pixel data-prep does **not** read this repo's `input/` (that holds NUTS geometries + Eurostat only). It
reads the big 1 km land-use maps / grid mapping / crop types from **`GRIDWORK_DIR` (~5 GB)** — the output of
the upstream **`LAMASUS_gridwork`** pipeline — plus `DS_DIR`. Defaults assume sibling folders next to this
repo; override with env vars:
```bash
export GAMBLE_GRIDWORK_DIR=/path/to/LAMASUS_gridwork/output   # ~5 GB, required for data-prep
export GAMBLE_DS_DIR=/path/to/LAMASUS_downscaling/
```

### Entry A — from raw grids (you have the gridwork data)
```bash
# data-prep ONLY: reads GRIDWORK_DIR, writes output/pixel_model_inputs.rds, exits before the fit (minutes)
DRIVER_DUMP_INPUTS=TRUE DRIVER_DUMP_EXIT=TRUE <BRANCH env vars, see §1.1> \
  Rscript run_prior_module_pixel_level_model.R

# fit the nested model off that dump (arg2 = .rds path => FROM_INPUTS, auto-detected)
Rscript run_nested_cut.R <BRANCH> output/pixel_model_inputs.rds <M> <NITER> <USE_RE> <IVMODE>
```

### Entry B — from a pre-built design (no gridwork data — the collaborator path)
`output/pixel_model_inputs.rds` (or a `dat_pixel_FULL_*.rds`) is portable — built once by someone with the
gridwork data, then shared. Model development needs only this, not the 5 GB of rasters:
```bash
Rscript run_nested_cut.R <BRANCH> path/to/pixel_model_inputs.rds <M> <NITER> <USE_RE> <IVMODE>
```

### Entry C — no data at all: synthetic validation
```bash
Rscript test_nested_cut.R      # known-λ DGP: λ recovery + held-out draws vs moments vs flat
```

`pixel_model_inputs.rds` holds `X_mat` (with intercept), `Y_pixel`, `group_idx_vec`, `re_group_names`,
`col_names`, and coords — everything the model needs.

### 1.1 Per-branch data-prep env config
`CLASS_SCHEME` is inferred from `DRIVER_CLASS_COLS`; branch defaults follow from it (see the driver header
around lines 45–260). Typical settings:

| branch   | key `DRIVER_*` env for data-prep |
|----------|----------------------------------|
| AGMIP    | `DRIVER_CLASS_COLS=AgMIP_label` (default) · `DRIVER_CROP_SPLIT=TRUE` (default) · `DRIVER_PIXEL_INTERSECT=NUTS3` |
| GLOBIOM  | `DRIVER_CLASS_COLS="GLOBIOM_UNFCCC,GLOBIOM_mngmt"` · `DRIVER_LUM_LINEAGE=globiom` · `DRIVER_INCLUDE_ORGANIC=TRUE` |
| BIOCLIMA | `DRIVER_LUM_LINEAGE=bioclima` · `DRIVER_PIXEL_INTERSECT=CAPRI_NUTS` (per the saved run) |

Always confirm against the driver header (`DRIVER_*` are read at the top of
`run_prior_module_pixel_level_model.R`); the class-name set the dump produces drives which `NEST_TREES`
builder applies (§3).

### `run_nested_cut.R` arguments
```
Rscript run_nested_cut.R <BRANCH> <arg2> <M> <NITER> <USE_RE> <IVMODE> <STREAM> <N_CHAINS>
  BRANCH   AGMIP | GLOBIOM | BIOCLIMA
  arg2     integer  -> subsample n pixels (0 = all), reconstructs design from saved dat+metadata (CFG)
           *.rds    -> FROM_INPUTS: read that pixel_model_inputs.rds directly (recommended for dev)
  M        inclusive-value imputations           (smoke 15 · production 30)
  NITER    Gibbs sweeps                           (smoke 800 · production 12000+)
  USE_RE   TRUE = country random effects · FALSE = pooled (fast smoke)
  IVMODE   draws | moments | auto                 (default auto)
  STREAM   TRUE = disk-stream sub-fits             (default TRUE)
  N_CHAINS 1 (default) · >1 (e.g. 4) -> per-node Rhat/ESS convergence check in fit$convergence
```
Writes `output/nested_cut_<BRANCH>_<timestamp>.rds` (`fit`, `tree`, coords, obs) and prints λ [95% CI] per nest.

### 1.2 Convergence — `N_CHAINS` (profile "b": full chains in parallel)
Set `N_CHAINS > 1` (the interactive `production` preset defaults to **4**) to get a **flat-model-analogue
convergence diagnostic** on the nested fit. Each **leaf** and each **nest** (fit at its mean inclusive value)
is run as `N_CHAINS` independent chains → split-Rhat + bulk-ESS per parameter via the `posterior` package,
summarised per node in **`fit$convergence`** (a data.frame: `node`, `max_rhat`, `pct_rhat_lt_1.01`,
`pct_rhat_lt_1.1`, `median_ess`, `n_params`) and echoed as a one-line worst-case summary. The **M inclusive-
value imputations stay single-chain** (the cut posterior is pooled across imputations, not a chain), so cost
is ≈ `N_CHAINS ×` the leaf + one per-nest convergence fit — not `N_CHAINS × M`. Inspect via
`nested_cut_convergence(fit)` or read `fit$convergence` directly. Use `N_CHAINS = 1` for fast dev iterations
(no Rhat computed).

**Parallel chains — `N_CORES`.** The chains within a leaf / convergence fit are independent cold starts, so
they fork in parallel via `parallel::mclapply` when `N_CORES > 1` (default `min(N_CHAINS, ncpu-1)`; CLI 9th
arg). On macOS/Linux the workers share `X`/`Y` copy-on-write and (with `stream_disk=TRUE`) each holds ~one
batch of RAM, so e.g. 4 chains on 4 cores run in ≈ 1× wall-clock instead of 4×. **Windows has no fork →
automatic serial fallback** (a startup message notes it). The M inclusive-value imputations stay serial
(hot-start chain). Validated: 3 chains, `N_CORES=3` → 2.9× speedup, identical per-node Rhat. Set
`N_CORES = 1` to force serial.

### 1.3 Persist + resume — `store_dir` (kill-safe runs)
Pass `store_dir = "output/ncut_store_<...>"` (the interactive script sets a stable per-branch default;
`run_nested_cut.R` uses env `NCUT_STORE_DIR`, `none` to disable) and each **completed node's draws are written
to `<store_dir>/nodes/<path>.rds` and kept** (not the throwaway temp dirs `stream_disk` uses — that's a
separate, per-fit RAM guard). If the run is killed, **re-launching with the same `store_dir` and the same
config loads every finished node from disk and only refits what didn't complete** — e.g. a crash during the
root IV fit skips the already-done leaves and resumes there. Validated: leaf-cached resume 15× faster,
full resume ~instant, λ/convergence identical.
- **Config guard:** a `manifest.rds` stores a fingerprint of the design + all fit knobs (data checksums,
  `niter`/`nburn`/`thin`/`M`/`use_re`/`iv_mode`/`n_chains`/…). Re-pointing at a store built with a **different
  config hard-stops** with a message — change `store_dir` to a fresh path or delete the folder to refit.
- **Inspect / clean:** `nested_cut_store_status(store_dir)` lists cached nodes. The store lives under
  `output/` (git-ignored); delete the folder to force a clean rerun.
- Nests cache with `children` stripped (rebuilt from their own caches on load), so there's no draw
  duplication across the tree. `stream_disk` and `store_dir` are independent and compose.

### 1.4 Per-fit progress heartbeat — `fit_progress`
Each sub-fit's own iteration bar is normally silenced (there are many fits). With `fit_progress = TRUE`
(default; interactive knob `FIT_PROGRESS_SEC`) every chain / leaf / imputation prints a **throttled
heartbeat** relaying the sampler's `Iteration i / niter [Burn-in|Sampling]`, labelled by node path and chain,
e.g. `[nested_cut]  root/Forests ch2/4  Iteration 2500 / 6000 [Sampling]`. `fit_progress_sec` sets the
throttle (default 30 s). Passing a `chain_id` to the sampler also **removes the stray `Sampling…0%`
`txtProgressBar`** that used to sit stuck at 0. In **parallel** mode (`n_cores>1`) the forked chains'
heartbeats interleave, but each line is self-labelled (`ch1/4`, `ch2/4`, …) so you can watch all chains
advance at once. During burn-in no disk batches are written, so this heartbeat is the *only* signal that a
fit is progressing then. Set `fit_progress = FALSE` to restore full silence.

---

## 2. `iv_mode` — how leaf uncertainty is carried up

The parent nests consume each sub-nest's **inclusive value** (a generated regressor). How its uncertainty
is propagated is the main speed/accuracy knob (details: `nested_cut_model.md` §6):

| `iv_mode`  | cost per internal node | when to use |
|------------|------------------------|-------------|
| `draws`    | **M** fits (exact cut) | reference / ground truth; any node with non-Gaussian IV |
| `moments`  | **1 + 2·moment_rank** fits (sigma-point, 2nd-order) | fast; IV ~Gaussian |
| `auto`     | diagnoses each node, then moments or draws | **default** — safe + fast |

- **`auto`** tests the skew/kurtosis of the IV field's leading directions and **falls back to exact `draws`
  wherever the IV is skewed** (e.g. sparse, quasi-separated leaves like the 14-class Arable nest), logging
  each decision. Inspect what was used: `nested_cut_modes(fit)`.
- Knobs: `moment_rank` (principal directions kept; 1–2), `moment_skew_tol` / `moment_kurt_tol` (fallback
  thresholds, defaults 1.0 / 2.0). Raise `moment_rank` if a node's IV is multi-directional.

### Memory — `stream_disk`
Each sub-fit's dominant RAM cost is the sampler's big per-group posterior array, which nested_cut does not
use (only the pooled draws). `stream_disk = TRUE` (default in the driver / interactive script) writes each
sub-fit's batches to a temp dir so that array **never accumulates in RAM** — peak drops to ~one batch and
becomes **independent of `niter`** (safer at full n / high iterations / 5 km). It reads only the pooled
draws back, then deletes the temp dir. Verified identical to in-RAM (`λ` and CI match exactly). Set
`stream_disk = FALSE` for slightly faster small/synthetic runs (avoids the disk I/O per sub-fit).

---

## 3. Adding / changing a classification tree

Trees are declared in `codes/nest_trees.R`; the model consumes any named list whose leaves are character
vectors of fine-class names. A new classification = **one `NEST_TREES` entry**.

- **Structured names** (`Type_subtype[_level]`, as GLOBIOM/BIOCLIMA) — auto-group by prefix:
  ```r
  NEST_TREES$MYSCHEME <- function(classes) nest_tree_by_prefix(classes, subsplit = "Cropland")
  ```
  (single-member prefixes become root singletons; `subsplit` groups a nest one level deeper.)
- **Bespoke grouping** (e.g. crops → Arable/Permanent, as AgMIP) — give each class an ancestor path:
  ```r
  paths <- c(setNames(rep(list(c("Cropland","Arable")), length(arable)), arable),
             setNames(rep(list(character(0)), length(singles)), singles), ...)
  NEST_TREES$MYSCHEME <- function(classes) nest_tree_from_paths(paths[intersect(names(paths), classes)])
  ```
- **Trees may be arbitrary / ragged depth** — singletons at the root, some branches 1 level, others 2+.
- Always validate coverage:
  ```r
  build_nest_tree("MYSCHEME", classes)   # prints the tree + "[coverage] N leaves / N classes | OK"
  ```
- To make the driver aware of a new branch, add a `CFG` entry in `run_nested_cut.R` (or just run via
  `FROM_INPUTS`, where the branch name only selects the tree builder).

---

## 4. Inspecting a fit

```r
source("codes/nested_cut.R")
fit <- readRDS("output/nested_cut_GLOBIOM_<ts>.rds")$fit
summary_nested_cut(fit)                      # per-nest zero-sum coef matrices + λ (median, CI)
nested_cut_modes(fit)                        # iv_mode actually used per node (auto fallbacks)
predict_nested_cut(fit, X_new, D = 200)      # posterior draws [n × J × D] of fine-class shares
effective_symmetric_table(fit, X)            # flat-comparable [J × P] marginal-effect table (median + CI)
```
The **native parameters are symmetric per nest**; `effective_symmetric_table` is a derived, flat-comparable
summary (zero-sum over fine classes). See `nested_cut_model.md` §9.

---

## 5. The other models

```bash
# pixel land-use MNL (full fit; env-driven config — see the driver header for DRIVER_* vars)
Rscript run_prior_module_pixel_level_model.R

# livestock count model (NB / Poisson)
Rscript run_prior_module_count_model.R
```

Post-processing (recover → fit metrics → all-class grid maps → HTML report) lives in `postprocess/`
(`engine.R` → `render_maps.R`/`render_heatplot.R` → `build_html.R`).

---

## 6. Conventions

- Run from repo root; Rcpp cores in `codes/*.cpp` compile on first use.
- Never commit `output/` (posterior batches, fits, reports — multi-GB, regenerable); `input/` **is** tracked.
- Exclude `.git/` from OneDrive sync (or keep the working clone outside OneDrive).
- Requirements: see `README.md`.
