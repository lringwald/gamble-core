# Design: self-contained `predict()` / `project()` for the pixel (and count) prior models

## Goal
A fitted model becomes ONE object that reproduces the **exact estimation-time treatment** on new
inputs — transforms, focal/neighbour computation, column assembly, group→RE mapping — so:
1. **No train/predict drift** (the same preprocessing code runs at fit and at predict), and
2. the **focal covariates are computed *inside* the predictor from `Y_t`**, which makes projection
   dynamic for free (predict → new `Y` → focal recomputes → predict → …).

Not literally inside the Rcpp sampler (that stays a generic, aspatial MNL/NB fitter). Instead a
thin **model-wrapper** at the driver level bundling `{fit, recipe, spatial_layout, group_map, class_names}`.

## Principle: fit/transform separation (sklearn-style)
- `build_recipe(raw, ...)` runs at **fit time**, *learns* the parameters (group levels & order, final
  column order, focal spec, baseline, bart/linear partition), returns `(X_mat, Y, recipe)`.
- `apply_recipe(recipe, X_t, Y_t)` runs at **predict time**, *replays* those learned parameters.
- The driver builds its training `X_mat` **via `build_recipe`** — so fitting and prediction share one
  code path. That IS the consistency guarantee.

### Standardization is already baked in — recipe does NOT re-standardize
`postb_total`/`postb_pooled` are returned on the **raw** covariate scale (back-transformed), and the
calibrated slim BART trees split on **raw** inputs (thresholds un-scaled via `bart_scaling`). So predict
does `raw X_mat %*% raw betas + BART(raw bart cols)` — **no re-centering/scaling step needed**. The recipe
only needs the *deterministic* design steps: transforms (log1p/exclusions), focal, column assembly.
(If a future fit uses `standardize=FALSE`, the recipe would own the scaling instead — keep it optional.)

## Objects

### `recipe` (learned at fit; pure metadata, cheap to store)
```
recipe = list(
  feature_cols   = <chr>,        # raw non-focal columns expected in X_t
  transforms     = list( list(cols=<chr>, fn="log1p"), ... ),   # declarative, in order
  focal_spec     = list(
      from        = "Y",         # focal is computed from the LU state Y
      classes     = <chr>,       # LU classes -> focal_<class> columns
      # RESOLUTION-AGNOSTIC, coordinate-based (NOT terra raster focal). Neighbours defined from the
      # (x,y) coords supplied WITH the inputs -> works at any pixel size / irregular layout.
      neigh       = list(type = "radius"|"knn", R = <num>, k = <int>, self = FALSE),
      agg         = "mean",      # focal_<class> = mean of neighbours' Y[,class] (SOFT shares)
      slice_key   = "year",      # focal computed WITHIN each time-slice (panel has 3x duplicate (x,y))
      nodata_flag = "focal_NODATA"       # 1 if a pixel has no neighbours in range
  ),
  col_order      = <chr>,        # final X_mat column names, in order (intercept first)
  baseline       = <chr/int>,    # baseline class
  bart_idx       = <int>,        # column indices into X_mat -> BART
  linear_idx     = <int>,        # column indices into X_mat -> linear (complement)
  offset_fn      = NULL|function # count model: log(area) from input
)
```

### `group_map`
```
group_map = list(
  key        = "GLOB_country",
  levels     = <chr>,            # unique(group_idx) IN THE SAMPLER'S APPEARANCE ORDER (see keying note)
  unseen     = "pooled"          # new/unseen groups -> postb_pooled
)
```
**Keying note (bug we already hit):** the sampler indexes REs by `unique(group_idx)` = order of first
appearance, NOT sorted. `group_map$levels` MUST store that same vector; predict maps a pixel's group via
`match(country, levels)`; `NA` (unseen) → pooled coefficients.

### `spatial_layout`  (resolution-agnostic — coordinate-based)
Just the per-pixel **`(x, y)` coordinates supplied with the inputs** — no raster template, no assumed
resolution. `compute_focal` builds neighbours on the fly from the coords via a spatial index
(kd-tree: `RANN`/`FNN`, or `spdep::dnearneigh`): each pixel's neighbours = others within radius `R`
(or `k` nearest). Works at 10km, 1km, or irregular NUTS3×EEA units identically — the neighbourhood is
defined by geometry, not a grid. `R` (or `k`) lives in `focal_spec$neigh` and is chosen to match the
intended window at the input resolution.
**Consequence:** this REPLACES the driver's terra raster-focal with the coordinate-based routine at
**both** fit and predict (so they're identical). Practically that means refitting with the coord-based
focal (or a one-time check that coord-focal == terra-focal on the current training grid before switching).

**Duplicate (x,y) / panel handling (important):** the training frame stacks years (2000/2010/2020) so every
coordinate appears once per year. `compute_focal` therefore operates **per `slice_key` (year)**: build the
spatial index *within each slice*, and exclude **self by ROW INDEX**, not by coordinate — otherwise a
pixel-year would pull in itself-in-other-years (distance 0) and cross-year neighbours. This matches terra
(which focal's each year's raster separately). If two rows share an identical coordinate *within one slice*
(genuine co-location or a data artifact), they are neighbours at distance 0 to each other (aggregated in),
and each still excludes only its own row — flag such within-slice coordinate duplicates at fit as a warning.
For `predict`/`project`, a call operates on a single time slice (one horizon), so the focal is within-slice
by construction.

### wrapper
```
model = list(fit=<sampler output>, recipe=recipe, group_map=group_map,
             spatial_layout=spatial_layout, class_names=<chr>, family="mnl"|"negbin")
class(model) <- "prior_model"
```

## `predict()`
```
predict(model, X_t, Y_t = NULL, group_t, return = c("shares","eta"))
```
Flow:
1. `focal <- compute_focal(model$spatial_layout, Y_t, model$recipe$focal_spec)`  — SHARED routine (below).
2. `raw <- cbind(X_t, focal)`
3. apply `recipe$transforms` (log1p …) in order.
4. assemble `X_mat` in `recipe$col_order` (+ intercept); assert names/dims match.
5. **linear + RE:** for each group g present, `eta[rows,] += X_mat[rows, linear_idx] %*% Bmean[,,map(g)]`
   where `Bmean = apply(fit$postb_total, 1:3, mean)`; unseen groups → `apply(fit$postb_pooled,1:2,mean)`.
6. **BART:** `f <- reconstruct_bart_f_mean(fit$tree_store, X_mat[, bart_idx], bart_meta)` ;  `eta <- eta + f`.
7. (count) `eta <- eta + offset` ; mean `mu = exp(eta)`.  (mnl) `shares = softmax(eta)`.
8. return `shares` (n × J)  [or `eta`].

Posterior-predictive variant: loop draws (`postb_total[,,,d]` + per-draw `tree_store[[d]]`) and average —
reuse the `assess_pixel_hybrid.R` / `MNL_reporting_suite.R` machinery.

## `project()` (dynamic co-evolution)
```
project(model, X_t, Y_init, n_steps = 20, damp = 0.5, tol = 1e-4)
```
```
Y <- Y_init
for s in 1:n_steps:
    P <- predict(model, X_t, Y_t = Y, group_t)      # focal recomputed from current Y
    Y_new <- (1-damp)*Y + damp*P                     # damped update (mean-field)
    if max|Y_new - Y| < tol: break
    Y <- Y_new
return list(final = Y, trajectory = ...)
```
- **damping is required** — an auto-model iterated hard can hit criticality (field collapses to one class);
  `damp` (mean-field partial step) keeps it stable. Expose it; default ~0.5; monitor per-step change.
- For multi-decade runs: `X_t` (climate/socioecon) varies by horizon; `Y` co-evolves via focal.

## The shared focal routine (the consistency-critical piece)
`compute_focal(spatial_layout, Y, focal_spec)` MUST be the SAME function the driver uses at fit. Refactor
the driver's terra focal block into this function and have `build_recipe` call it — then fit and predict
cannot diverge. Watch: window/kernel, `na.rm`/fill, edge expansion, and the empty-neighbourhood
`focal_NODATA` flag. Cheapest safe path = reuse terra with the stored window; a faster
sparse-adjacency reimplementation is a later optimisation only if it's validated bit-identical to terra.

## What changes where (mostly refactor, not new stats)
- `codes/` new file `prior_model_predict.R`: `build_recipe`, `apply_recipe`, `compute_focal`,
  `predict.prior_model`, `project.prior_model`.
- `run_prior_module_pixel_level_model.R`: build training `X_mat` via `build_recipe` (so the path is shared);
  after fitting, assemble+save the `prior_model` wrapper.
- Sampler: **unchanged** (stays aspatial). Ensure it returns what predict needs (it already does:
  `postb_total/pooled`, `tree_store`, `baseline`, `bart_idx`; confirm `do_slim_trees=TRUE` so trees are
  calibrated/raw-scale).
- Back-end reused as-is: `reconstruct_bart_f_mean`, `predict_slim_bart_cpp`, softmax.

## Edge cases / decisions to lock down
- **Unseen groups** at predict → pooled FE (documented above).
- **Missing/NA covariates** in `X_t` → same imputation/`0`-fill rule as fit (put it in the recipe).
- **`Y_t` at prediction** is *shares in [0,1]* (soft), not hard classes → focal on expected composition;
  fine and smoother than hard assignment (matches training where focal = neighbour area shares).
- **Column mismatch** → hard error (fail loudly, don't silently reorder).
- **BART on new grid** → `predict_slim_bart_cpp` handles any rows; unseen covariate ranges extrapolate
  as BART constants (flag if a downscaling grid goes far outside train support).

## Phasing
0. **DONE — coord focal built + validated ≈ terra.** `experiments/focal/compute_focal_coord.R`:
   resolution-agnostic (infers `res` from coords), per-slice (`out_year`), queen-8, `complete_only`
   (mirrors terra `na.rm=FALSE` → incomplete neighbourhood = 0), defined duplicate rule (aggregate
   sub-units per cell, then focal, then broadcast). Equivalence vs stored terra focal on
   dat_pixel_FULL_10km (160k rows, 16 classes, shares): **mean|d|=0.0004** (clean 0.0011, dup 0.0002);
   max|d|~1.0 on a handful of boundary/dup cells (terra's arbitrary xyz-collapse vs the clean aggregate).
   => coord focal is a drop-in for terra; **existing fit is reusable (no refit)**; duplicate (x,y) handled
   AND the latent terra collapse-ambiguity is fixed. IMPORTANT: normalize own AREAS -> shares per row over
   the focal-class set BEFORE focal (driver L624-627), else scale mismatch.
1. **DONE — static predict.** `codes/prior_model_predict.R::predict_shares(fit, X, bart_cols, linear_cols,
   group_idx, group_levels, type=c("mean","predictive"))` — tree-based BART f (reconstruct_bart_f_mean) +
   linear/RE (postb_total, appearance-order keying, unseen->pooled) + softmax; NO re-standardization.
   Unit test experiments/focal/test_predict_phase1.R: (A) tree-f==stored post_f_mean max|d|=3e-13;
   (B) predict_shares==in-sample fitted shares max|d|=6e-14, cor 1.000000 (EXACT); (C) sorted-vs-correct
   keying max|d|=0.31 (keying handled); (D) unseen-group->pooled OK; (E) predictive path cor 0.999 w/ mean.
2. **DONE — recipe + wrapper.** `codes/prior_model_predict.R`: `build_recipe`/`apply_recipe`/`predict_prior`
   + `model=list(fit,recipe,bart_names,linear_names)`. Test `experiments/focal/test_wrapper_phase2.R`:
   apply_recipe==build design max|dX|=0.0; predict_prior==in-sample fitted max|d|=0.0 cor 1.000000; focal
   live (predict responds to Y change, mean|dP|=0.039). Gotchas: drop pre-existing focal_ before recompute;
   pass Y non-rounded; filter zero-LU rows.
3. **DONE — `project_prior()`.** Iterate predict_prior feeding Y back (focal recomputes) + damping.
   Diagnostic: gap halves each step (damp=0.5 contraction), TRUE fixed point predict(Y*)==Y* to 1e-5,
   unique from any init, stable/no-collapse. Assumes outcome_classes==lu_classes (else needs outcome->focal map).
3b. **CLASS_MAP (user-requested fold-in)** — nested output classes inherit parent's fitted predictor
   (coefficient inheritance, "same hierarchical predictor value"); focal aggregates nested Y children->parent.
   Identity class_map == current behaviour.
4. **Wire build_recipe into the driver** so training design comes through the recipe; save the prior_model wrapper.
4. (optional) sparse-adjacency focal for speed, gated on bit-identical-to-terra validation.

## Decisions (locked 2026-07)
- **Resolution-agnostic**: focal built from supplied `(x,y)` coords via a spatial index (radius/knn), not
  terra. Works at any resolution / irregular layout. Replaces terra focal at fit AND predict (refit or
  validate-equivalence once).
- **`project()` focal feedback = SOFT shares** — dictated by consistency (training `Y` is a composition, so
  training focal = neighbour shares; matching treatment ⇒ soft). Also more stable in the iteration.
- **Both** posterior-mean and posterior-predictive supported: `predict(..., type = c("mean","predictive"))`.
  Default `"mean"` (fast); `"predictive"` averages softmax over draws for uncertainty.

## Still to confirm
- The neighbourhood window: radius `R` vs `k`-nearest, and its value — should reproduce the *intended*
  focal window (e.g. queen-8 at the input resolution ≈ `R ≈ 1.5 × cell_size`). Pick per use-case.
- Whether to keep terra as a validation reference before fully switching the training focal to coord-based.
