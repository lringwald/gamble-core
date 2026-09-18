# gamble-core

Processing, estimation, and post-processing environment for the **GAMBLE** high-resolution
prior land-use and livestock models (LAMASUS). Consolidates the production pipelines:

- **Pixel land-use model** — zero-sum Bayesian multinomial logit on the 10 km EU grid, random
  effects by region, Pólya-Gamma Gibbs (Rcpp core).
- **Livestock count model** — Bayesian Negative-Binomial / Poisson regression for livestock heads.
- **Nested land-use model** — cut / inclusive-value nested MNL over an arbitrary classification tree
  (AgMIP / GLOBIOM / BIOCLIMA), with a moment-carrier fast path.

## Layout

```
run/              ENTRY POINTS — start here
  flat.R            build the design dump and/or fit the FLAT pixel model (pick the CLASSIFICATION)
  nested.R          fit the NESTED model (factorized | iv, RE block, symmetric HS)
  count.R           fit the livestock count model (BOV + SGT in one pass)
  score.R           score fitted models side by side
  fit_nested.sh     non-interactive form of nested.R:  run/fit_nested.sh prod factorized intercept symhs
codes/            core samplers + Rcpp cores + shared helpers
  mnlogit_rcpp_sym.R + mnlogit_gibbs_core_sym.cpp   pixel MNL (symmetric zero-sum)
  mnlogit_rcpp.R     + mnlogit_gibbs_core.cpp        base MNL (dependency)
  count_rcpp.R       + count_gibbs_core.cpp          livestock count (NB/Poisson)
  nested_cut.R, nest_trees.R, mnlogit_nested_iv.R    nested / inclusive-value framework
  predict_gamble.R, prior_model_predict.R            prediction engines
  spatial_split.R, spatial_utils.R, mundlak.R        spatial & model specification utilities
  mnl_aux_func.R, rotate_draws.R, target_class_split.R
drivers/          pipeline drivers (env-driven; run/ wraps them)
  run_lu_pixel_model.R   (pixel MNL / design assembly)
  run_ls_count_model.R   (livestock count)
  run_nested_cut.R       (nested cut)
  run_bart_gate.R        (linear vs BART gate)
prep/             data assembly & harmonization:
  data_preparation_pixel_panel.R, compile_complete_LUM_map.R,
  harmonize_to_targets.R, build_globiom_subclass.R,
  prepare_eurostat_livestock_shares.R, prepare_composition_training.R
composition/      subtype composition: fit_composition.R, fit_delta_contrast.R, fit_organic.R,
                  build_subclass_parameters.R, consolidate_composition.R
postprocess/      reporting & diagnostics:
  engine.R, render_maps.R, render_heatplot.R, build_html.R (fit report),
  calculate_fit_metrics.R, plot_results_heatmap.R, count_validation_report.R,
  score_nested_cut.R, check_progress.R, diagnose_posterior.R, MNL_*
experiments/      measurement harnesses (mixing/re_idx_tradeoff.R, focal/, count/, nested/)
tests/            synthetic validation & sampler gates:
  run_all.R (unified runner), test_suite_lu_pixel.R (37-check gate),
  test_suite_ls_count.R (count gate), test_*.R
aux_files/        LUM_Code_to_macro_model_mapping.csv — CURATED class + nest taxonomy (source of truth)
docs/             model specs (docs/nested_cut_model.md = nested model + algorithm)
input/            model input data (grids, NUTS geometries, Eurostat)   [tracked]
output/           all model outputs / batches / reports                 [git-ignored, regenerable]
_local/           scratch + superseded scripts & legacy samplers        [git-ignored]
```

## Running

**Quick reference for every model and every way to launch it: [`run/README.md`](run/README.md).**

Everything is run **from the repo root**. Open `gamble-core.Rproj` in RStudio (which sets the working
directory), then edit the CONTROL PANEL at the top of a `run/` script and source it.

```r
# 1. build the design (choose the classification inside the file), then
source("run/flat.R")      # BUILD_DESIGN_ONLY = TRUE -> output/designs/pixel_model_inputs.rds

# 2. fit the nested model on it (variant / RE block / symmetric HS inside the file)
source("run/nested.R")

# 3. livestock count model, and scoring
source("run/count.R")
source("run/score.R")
```

Non-interactive equivalents:

```bash
run/fit_nested.sh smoke factorized intercept symhs    # ~10-15 min sanity run
run/fit_nested.sh prod  factorized intercept symhs    # full run (~24 h), backgrounded + resumable
Rscript tests/run_all.R                               # EVERY test + HTML report & figure
Rscript tests/run_all.R fast                          # same, minus the two long suites
Rscript tests/test_suite_lu_pixel.R                   # sampler gate (37 checks)
Rscript tests/test_nested_cut.R                       # nested-framework validation
```

The nested runners write a **resumable node store**; a killed run picks up at the last completed node.
A store fingerprints the data/niter/M but **not** the sampler version — never point a run at a store
built by a different version of `codes/mnlogit_rcpp_sym.R`.

Outputs go to `output/` (git-ignored). The Rcpp cores in `codes/*.cpp` compile on first use.

## Requirements

R (>= 4.2) with: `Rcpp`, `RcppArmadillo`, `data.table`, `qs2`, `matrixStats`, `posterior`, `abind`,
`future`, `future.apply`, `progressr`, `ggplot2`, `scales`, `terra`/`sf` (geodata), `base64enc` (report).

## Notes

- **Data policy:** `input/` (grids, NUTS, Eurostat) is tracked; all `output/` (posterior batches, fitted
  objects, reports — multi-GB, regenerable) is git-ignored.
- **OneDrive + git:** if this working copy lives inside OneDrive, exclude `.git/` from OneDrive sync (or
  keep the working clone outside OneDrive) to avoid sync/lock corruption of the git objects.
- Model theory & algorithms: `docs/`. Nested cut model (equations + computation): `docs/nested_cut_model.md`.
