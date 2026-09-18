# BMLEH_Los1_CAPRI — livestock (count) pipeline, end to end

Run every command from the **gamble-core root**.

THE ORDER IS LOAD-BEARING. The count model writes `dat_admin_FULL_*.rds`, and every downstream step
takes its land-use vocabulary from that file. Running the composition prep BEFORE the count model is
what produced the GLOBIOM-vs-BMLEH driver mismatch: the fits then carried `lu_area_Pasture_HI` while
the totals beta carried `log1p_lu_area_Grassland_extensive`, so no land-use driver joined and gamma
was assembled from the socioeconomic drivers alone — silently, with a full-looking output file.

## 0. Prerequisites (once)

    # BMLEH_Los1_label into the mapping the COUNT model reads (a different copy from the pixel driver's)
    BMLEH_MAPPING="../cascadinggamble-core/data/aux_files/LUM_Code_to_macro_model_mapping.csv" \
      Rscript projects/BMLEH_Los1_CAPRI/gamble_model/build_class_column.R

    # Stage the 1.5GB master parquet on LOCAL disk. Reading it over Google Drive times out
    # (IOError errno 60) partway into the fit -- arrow does column-projected random-access reads and
    # the Drive FUSE layer handles them badly.
    cp ../cascadinggamble-core/data/02_intermediate/prior_model_1km_master_inputs.parquet /tmp/

## 1. Totals count model (BOV + SGT)          ~14 min, 4 chains

    GAMBLE_MASTER_PARQUET=/tmp/prior_model_1km_master_inputs.parquet \
    DRIVER_CLASS_COLS=BMLEH_Los1_label \
    DRIVER_RE_GROUP_COL=CAPRI_NUTS \
    DRIVER_FIT_COMPOSITION=FALSE \
      Rscript run_ls_count_model.R

`DRIVER_RE_GROUP_COL=CAPRI_NUTS` resolves to the column `CAPRI_country` (the 2-char CAPRI country,
by dominant overlap per NUTS3) -- the downscaler keys on CAPRI codes, and a GLOB_country fit is not
relabellable after the fact. `DRIVER_FIT_COMPOSITION=FALSE` because the composition is run
explicitly below, with the CAPRI target; the built-in sub-run would silently use the D/O/F table.

## 2. Composition training tables            (AFTER step 1)

    Rscript prep/prepare_composition_training.R                                # D/O/F + organic
    Rscript projects/BMLEH_Los1_CAPRI/gamble_model/prepare_livestock_capri_training.R   # 8-category cattle

## 3. Composition fits                        production settings

    D_SPECIES=bov D_TRAINING=output/composition/bov_training_nuts2_capri.csv \
      Rscript composition/fit_composition.R        # 8 CAPRI cattle categories
    D_SPECIES=sgt Rscript composition/fit_composition.R                # 3 categories (D/O/F)
    Rscript composition/fit_organic.R

Defaults are `D_NITER=20000 D_NBURN=10000 D_NCHAINS=4`. Override only for smoke runs.
BOV and SGT may have DIFFERENT category counts -- that is supported; categories are read from each
fit's own `cats` field downstream.

## 4. Assemble gamma

    Rscript composition/consolidate_composition.R
    Rscript composition/build_subclass_parameters.R

`build_subclass_parameters.R` prints two things worth reading:
  * the RE grouping it resolved from the fit (`CAPRI_country`) -- a group-count mismatch is now a
    hard stop, not a silent fallback to integer country labels;
  * `N of M composition drivers match the totals design` -- warns below 50%. A low number means the
    composition tables and the totals fit are on different classifications (see the ordering note).

## 5. CAPRI activity codes + report

    Rscript projects/BMLEH_Los1_CAPRI/gamble_model/build_livestock_capri.R
    Rscript postprocess/build_count_html.R

Outputs: `output/gamble_model/BMLEH_Los1_CAPRI/livestock_capri_parameters.csv`
         `output/report/count_model_fit_report_<date>.html`

## Known limitations

* Only 2000 and 2018 granular LUM maps exist, so the 2010 tier is explained by land use from eight
  years later (`LIVESTOCK_SPEC.md`).
* PIGF / SOWS / PKPL are not covered: the count model fits grazing species only.
* Nothing in the cascadinggamble BMLEH downscaling consumes `livestock_capri_parameters.csv` yet --
  `process_CAPRI_LSU_downscaling.R` builds its own area-based priors (grassland for grazers,
  cropland for pigs/poultry). Wiring that up is a separate step.
