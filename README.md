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
codes/            core samplers + Rcpp cores + shared helpers
  mnlogit_rcpp_sym.R + mnlogit_gibbs_core_sym.cpp   pixel MNL (symmetric zero-sum)
  mnlogit_rcpp.R     + mnlogit_gibbs_core.cpp        base MNL (dependency)
  count_rcpp.R       + count_gibbs_core.cpp          livestock count (NB/Poisson)
  nested_cut.R, nest_trees.R, mnlogit_nested_iv.R    nested / inclusive-value framework
  lnm_gibbs.R, mvclr_gibbs.R (+ cores)               alternative compositional samplers
  mnl_aux_func.R, prior_model_predict.R, spatial_utils.R,
  MNL_parameter_heatplot.R, MNL_reporting_suite.R, MNL_re_viz_utils.R
experiments/focal/  compute_focal_coord.R           focal (t-1) neighbourhood covariate
data prep (root):   data_preparation_pixel_panel.R, prepare_eurostat_livestock_shares.R,
                    prepare_composition_training.R, compile_complete_LUM_map.R
composition (root): fit_composition.R, fit_delta_contrast.R, build_subclass_parameters.R,
                    consolidate_composition.R
drivers (root):     run_prior_module_pixel_level_model.R   (pixel MNL)
                    run_prior_module_count_model.R         (livestock count)
                    run_nested_cut.R                       (nested)
postprocess/        engine.R, render_maps.R, render_heatplot.R, build_html.R   (fit report:
                    recover -> metrics -> grid maps -> HTML), + calculate_fit_metrics.R,
                    plot_results_heatmap.R
tests/ (root):      test_nested_cut.R, test_nested_iv.R, codes/master_test_suite_sym.R
docs/               model specs (see docs/nested_cut_model.md for the nested model + algorithm)
input/              model input data (grids, NUTS geometries, Eurostat)  [tracked]
output/             all model outputs / batches / reports   [git-ignored, regenerable]
```

## Running

All scripts use paths **relative to the repo root** — run from the repo root:

```r
# pixel land-use model (env-driven config; see the driver header for DRIVER_* vars)
Rscript run_prior_module_pixel_level_model.R

# livestock count model
Rscript run_prior_module_count_model.R

# nested model:  <BRANCH> [SUBSAMPLE] [M] [NITER] [USE_RE] [IVMODE]
Rscript run_nested_cut.R AGMIP 0 25 1000 TRUE auto

# nested-framework synthetic validation
Rscript test_nested_cut.R
```

Outputs are written under `output/` (created on demand; git-ignored). The Rcpp cores in `codes/*.cpp`
compile on first use.

## Requirements

R (>= 4.2) with: `Rcpp`, `RcppArmadillo`, `data.table`, `qs2`, `matrixStats`, `posterior`, `abind`,
`future`, `future.apply`, `progressr`, `ggplot2`, `scales`, `terra`/`sf` (geodata), `base64enc` (report).

## Notes

- **Data policy:** `input/` (grids, NUTS, Eurostat) is tracked; all `output/` (posterior batches, fitted
  objects, reports — multi-GB, regenerable) is git-ignored.
- **OneDrive + git:** if this working copy lives inside OneDrive, exclude `.git/` from OneDrive sync (or
  keep the working clone outside OneDrive) to avoid sync/lock corruption of the git objects.
- Model theory & algorithms: `docs/`. Nested cut model (equations + computation): `docs/nested_cut_model.md`.
