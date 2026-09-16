# Running the models — quick reference

Everything runs **from the repo root**. All scripts resolve paths like `codes/…` relative to it, so
`cd gamble-core` first (or open `gamble-core.Rproj`, which sets the working directory for you).

There are three ways to run every model. They are the *same code path* — the `run/*.R` files just set
environment variables and call the drivers, so an interactive fit and a batch fit produce identical
results.

| | how | best for |
|---|---|---|
| **A. RStudio** | open `run/<x>.R`, edit the CONTROL PANEL, source it | exploring, changing knobs, inspecting `fit` afterwards |
| **B. Terminal** | `Rscript run/<x>.R` | reproducing exactly what the CONTROL PANEL says |
| **C. Shell wrapper / env vars** | `run/fit_nested.sh …` or `VAR=… Rscript <driver>.R` | long runs, sweeps, anything scripted |

---

## 1. Build the design (do this first)

Assembles X and Y and writes `output/pixel_model_inputs.rds`, which the nested model consumes.
The **classification** knob decides what the model predicts. ~30 min.

**A.** open `run/flat.R` → set `CLASSIFICATION`, keep `BUILD_DESIGN_ONLY <- TRUE` → source
**B.** `Rscript run/flat.R`
C.
```bash
GAMBLE_MASTER_PARQUET=/Users/leopoldringwald/gamble_local_data/prior_model_1km_master_inputs.parquet \
DRIVER_CLASS_COLS="GLOBIOM_subclass" DRIVER_PROMOTE_NATURAL_OTHER=TRUE \
DRIVER_DUMP_INPUTS=TRUE DRIVER_DUMP_EXIT=TRUE \
Rscript run_lu_pixel_model.R
```

Classifications: `GLOBIOM_subclass` (27 classes + the curated 6-nest tree — production),
`GLOBIOM_UNFCCC,GLOBIOM_mngmt` (legacy 19-class), `AgMIP_label` (crop types), `BIOCLIMA_DS_reporting`.

Check what you built:
```r
inp <- readRDS("output/pixel_model_inputs.rds")
dim(inp$X_mat); colnames(inp$Y_pixel); inp$class_nest   # class_nest = the curated tree
```

## 2. Flat pixel model

Single-level MNL over all classes at once — no nesting.

**A.** `run/flat.R` with `BUILD_DESIGN_ONLY <- FALSE` → source **B.** `Rscript run/flat.R`
**C.** as above without the two `DRIVER_DUMP_*` flags, plus `DRIVER_NITER=1000`.

## 3. Nested model

**A.** open `run/nested.R`, set variant / RE block / symmetric HS → source. Fit lands in `fit`.
**B.** `Rscript run/nested.R`
**C.**
```bash
run/fit_nested.sh smoke factorized intercept      symhs  # ~10-15 min sanity check
run/fit_nested.sh prod  factorized intercept      symhs  # production (~6-24 h), backgrounded
run/fit_nested.sh prod  factorized intercept      diag   # diagonal measurement arm
```
or drive `run_nested_cut.R` directly — args are
`BRANCH design.rds M NITER USE_RE IVMODE STREAM N_CHAINS N_CORES`:
```bash
NCUT_USE_IV=FALSE NCUT_RE_COLS=intercept NCUT_SUB=0 NCUT_SYM_HS=TRUE \
NCUT_STORE_DIR=output/ncut_store_myrun \
Rscript run_nested_cut.R GLOBIOM output/pixel_model_inputs.rds 25 1000 TRUE auto TRUE 1 1
```

Choices, and what the measurements say:

| knob | options | note |
|---|---|---|
| variant | `factorized` \| `iv` | factorized is the validated choice; `iv` estimates λ but on the 27-class tree λ fell outside (0,1] on 4 of 6 nests, and it is far slower (root refit per imputation) |
| RE block | `intercept` \| `intercept+socio` | intercept-only wins on held-out *and* converges (RE Rhat 1.1 vs 1.5); the 7 slopes tie on skill for 8x the per-group parameters |
| shrinkage | `diag` \| `symhs` | symmetric = baseline-invariant shrinkage. The old "~120 nats worse held-out" gate is **VOID** — those arms ran the pre-2026-08-19 bug, at Rhat 1.91. Post-fix the only head-to-head is IN-SAMPLE at one node (McFadden +0.216 sym vs +0.260 diag); no held-out comparison has been run. `diag` is the untested default, not a measured winner |
| `N_CHAINS` | 1 \| 4 | 4 gives per-node Rhat/ESS in `fit$convergence`, at ~4x cost |

## 4. Livestock count model

One run fits **both** BOV and SGT.

**A.** `run/count.R` → source **B.** `Rscript run/count.R`
**C.** `DRIVER_RUN_MODE=production DRIVER_NITER=20000 Rscript run_ls_count_model.R`

The count Gibbs mixes slowly; for publication-grade runs use
`experiments/count/long_run_segmented.R`, which is resumable (long single jobs get killed here).

## 5. Reports

HTML overview of nested fits — comparison table, per-nest and per-class composition, observed vs
predicted maps, provenance for each fit. Single self-contained file:

```bash
Rscript postprocess/nested_report.R                       # 3 newest fits
Rscript postprocess/nested_report.R fitA.rds fitB.rds     # specific ones
REPORT_MAP_CLASSES=20 Rscript postprocess/nested_report.R # more mapped classes (default 12)
```
-> `output/report/nested_fit_report.html`

The FLAT model has its own report pipeline (`postprocess/engine.R` -> `render_maps.R` ->
`build_html.R` -> `output/report/prior_model_fit_report.html`); it reconstructs from on-disk posterior
batches and does not read nested fits. The two share a stylesheet, so they look like one product.

## 6. Scoring and tests

```r
source("run/score.R")                                    # compare fits side by side
```
```bash
Rscript codes/score_nested_cut.R output/nested_cut_A.rds output/nested_cut_B.rds
Rscript tests/run_all.R                  # EVERY test + HTML report & figure -> output/report/
Rscript codes/test_suite_lu_pixel.R      # 37-check sampler gate — run before/after touching codes/
Rscript tests/test_nested_cut.R          # nested-framework synthetic validation
```
Scoring is **in-sample**: it measures reproduction, not generalisation, and favours the richer RE
block. For an honest RE comparison use `experiments/mixing/re_idx_tradeoff.R` (country-stratified
held-out).

---

## Gotchas that have actually cost time here

- **Run from the repo root.** `Rscript prep/foo.R`, never `cd prep && Rscript foo.R`.
- **Design build needs BOTH dump flags.** `DRIVER_DUMP_EXIT` alone is nested inside the
  `DRIVER_DUMP_INPUTS` block, so on its own it neither dumps nor exits — it silently runs a full fit.
- **The node store is resumable but version-blind.** It fingerprints data/niter/M, *not* the sampler
  version. Never point a run at a store built by a different `codes/mnlogit_rcpp_sym.R`; give each
  configuration its own `NCUT_STORE_DIR`. To resume a killed run, reuse its store path exactly.
- **Read the parquet from the local copy.** Straight off Google Drive, arrow intermittently dies with
  `IOError … [errno 60] Operation timed out`. Hence `GAMBLE_MASTER_PARQUET`.
- **Check provenance before trusting a fit** — every saved fit records what produced it:
  ```r
  z <- readRDS("output/nested_cut_…rds")
  z$inputs; z$inputs_mtime; z$tree_source; z$use_iv; z$re_cols; z$symmetric_hs
  ```
  A silently reused stale design has cost a full multi-hour run before.
- **Killing jobs:** kill by PID. `pkill -f run_nested_cut.R` will take down every nested run at once,
  including production ones you meant to keep.
