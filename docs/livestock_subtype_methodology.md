# A Bayesian Two-Stage Prior Model for High-Resolution Livestock Composition

*Totals, driver-conditioned subtype allocation, and organic farming systems*

GAMBLE prior-model suite — livestock module. This document specifies, from start to finish, the
econometric model behind the livestock subtype parameters
(`output/composition/subclass_country_parameters.csv`) produced by the count-model pipeline
(`run_prior_module_count_model.R`, Sections 1–9). It is written to be read alongside the code; every
equation is annotated with the script/function that implements it. Companion methodological notes:
[`nested_cut_model.md`](nested_cut_model.md) (the pixel land-use analogue) and
[`livestock_subtype_allocation_spec.md`](livestock_subtype_allocation_spec.md) (data provenance).

---

## Abstract

We estimate a spatially explicit prior for the European livestock herd that is consistent with a coarse
economic land-use model (GLOBIOM) yet resolves detail the coarse model does not carry: the number of
grazing animals per NUTS3 region, their split into production **subtypes** (dairy / meat / followers), and
the **farming system** (organic vs. conventional). The estimator is deliberately **modular** ("cut"): a
Stage-1 Negative-Binomial count model carries all *spatial* signal in the herd totals, while Stage-2
multinomial-logit composition models carry the *driver-conditioned* subtype and system splits. The three
fitted objects are combined by a log-linear plug-in rule and renormalized within region. The organic split
is identified from an organic-grassland proxy derived from a 1 km organic-certification map, with its
national level **anchored to Eurostat** by a per-country moment condition. All coefficients are symmetric
(sum-to-zero) so that no category is a reference and the effects are directly comparable across the
tree. Estimation is fully Bayesian (Pólya-Gamma Gibbs) with country random effects and horseshoe
shrinkage.

---

## 1. Notation and overview

| Symbol | Meaning |
|---|---|
| $r \in \{1,\dots,R\}$ | region (NUTS3 admin unit); country $c(r)$ |
| $t$ | year, $t \in \{2000, 2010, 2020\}$ |
| $s \in \{\mathrm{BOV},\mathrm{SGT}\}$ | species (bovines; small grazers = sheep + goats) |
| $Y_{rts}$ | livestock total of species $s$ in region $r$, year $t$ (LSU/head, model-scaled) |
| $x_{rt}\in\mathbb{R}^{p}$ | driver vector (land use, terrain, climate, socioeconomics) |
| $A_{rt}$ | exposure/offset (agricultural area) |
| $g\in\{\mathrm D,\mathrm O,\mathrm F\}$ | production subtype: **D**airy, **O**ther (meat/suckler), **F**ollowers |
| $q\in\{\text{conv},\text{org}\}$ | farming system: conventional, organic |
| $\beta^{(s)}_{c}$ | Stage-1 totals coefficients (per species, per country) |
| $\delta^{(s)}_{g}$ | Stage-2a subtype composition contrast (per species) |
| $\delta^{(s)}_{\mathrm{org}}$ | Stage-2b organic-system contrast (per species) |
| $\gamma^{(s)}_{c,g,q}$ | final plug-in cell coefficient |

**Target.** For every species we produce a $12$-cell allocation, $\{\mathrm D,\mathrm O,\mathrm F\}\times
\{\text{conv},\text{org}\}$, as a set of per-country, per-driver coefficients $\gamma^{(s)}_{c,g,q}$. At
prediction time a grid cell's share of a subtype-system is
$w \propto \exp(\gamma^{(s)\top}_{c,g,q} x)$, renormalized within region.

**Modularity ("cut").** Information flows **one way**: totals → subtype composition → system composition →
plug-in. Downstream stages condition on upstream point summaries but never feed back. This is a Bayesian
*cut* (Plummer, 2015): it trades a small amount of statistical efficiency for robustness — a
mis-specified composition model cannot contaminate the well-identified totals, and each stage can be
estimated with the data and resolution appropriate to it (Section 3).

---

## 2. Data

### 2.1 Spatial units and resolutions

The three stages are estimated at the finest resolution their response supports:

| Stage | Response | Unit | $n$ | Source |
|---|---|---|---|---|
| 1 — totals | $Y_{rts}$ | NUTS3 × year | ~4,090 rows | GLOBIOM/observed LSU on the 1 km grid, aggregated |
| 2a — D/O/F | head counts | NUTS2 | ~190 (BOV) / ~160 (SGT) | Eurostat `agr_r_animal` |
| 2b — organic | grassland-area proxy | NUTS3 | ~4,090 | organic-certification 1 km map |

The subtype split is capped at NUTS2 because that is the finest resolution at which Eurostat reports
subtype head counts; the organic split, by contrast, is estimated at **native NUTS3** because its proxy
response is available on the grid.

### 2.2 Livestock totals

$Y_{rts}$ is the herd of species $s$ (BOV = cattle; SGT = sheep + goats), obtained from the GLOBIOM/observed
LSU layer on the 1 km grid and summed to NUTS3 for three census years. Totals are optionally scaled by a
constant `y_weight` (default 100) to improve Negative-Binomial mixing without changing the fitted slopes.

### 2.3 Subtype head counts (D/O/F) — Eurostat

Subtypes are defined from Eurostat `agr_r_animal` (annual, NUTS2):

- **BOV**: $\mathrm D=$ dairy cows (`A2300F`), $\mathrm O=$ suckler/meat cows (`A2300G`),
  $\mathrm F=$ total (`A2000`) $-$ cows (followers/young stock).
- **SGT**: $\mathrm D=$ milk ewes (`A4110KC`) + dairy she-goats (`A4210K`), $\mathrm O=$ non-milk ewes
  (`A4110KD`), $\mathrm F=$ (sheep `A4100` + goats `A4200`) $-$ breeding females.

*Note.* The subtype letter **O** denotes **Other/meat**, not organic. Conflating it with the organic
"O" was a defect in an earlier build; the two dimensions are now orthogonal (Sections 5.2–5.3).

### 2.4 Organic system — national levels + spatial proxy

There is **no sub-national organic head count**. The organic split therefore combines two sources:

- **Level (national).** Eurostat `org_lstspec` gives the organic share of each species' herd per country,
  $\pi^{\mathrm{eur}}_{c,s}$ (`output/eurostat/organic_shares_nuts0.csv`: `BOV_org`, `SGT_org`; year 2010
  for the reference run). By construction this carries no within-country detail.
- **Pattern (spatial).** The 1 km organic-certification map yields an *organic-grassland* area per NUTS3.
  We use it as a proxy response: $n^{\mathrm{org}}_{r} = A^{\text{Pasture\_HIO}}_{r}+A^{\text{Pasture\_LIO}}_{r}$
  against conventional $n^{\mathrm{conv}}_{r} = A^{\text{Pasture\_HI}}_{r}+A^{\text{Pasture\_LI}}_{r}$.

**Proxy assumption.** The share of *organic grassland* is taken to be an unbiased spatial signal for the
share of *organic grazing livestock* within a country (they are then re-levelled to Eurostat, §6.3). This
is the model's central organic identifying assumption and is discussed in §9.

### 2.5 Drivers

Drivers $x_{rt}$ (constructed identically across stages so the plug-in rule of §6 is valid):

- **Land-use composition** — area (km²) of each GLOBIOM class (`lu_area_*`), entered as $\log(1+\text{area})$.
- **Terrain** — flat/steep/lowland/upland shares, slope/elevation heterogeneity (sd).
- **Climate** — growing-degree-days, precipitation level and seasonality.
- **Socioeconomic** — GDP, population (as $\log(1+\cdot)$), protected-area share.

Human-footprint/accessibility variables (GHM, CISI, RAI) and raw terrain means are **excluded** from the
composition stages by default (`D_EXCL_ACCESS=TRUE`): they are partly outcome-confounded (developed lowlands
*are* dairy) and induce quasi-separation. Extensive quantities (areas, GDP, Pop) enter as
sums→$\log(1+\cdot)$; intensive quantities (terrain, climate) as area-weighted means.

---

## 3. Stage 1 — livestock totals (Negative-Binomial count model)

Implemented in `mncount_rcpp` (`codes/count_rcpp.R` + `codes/count_gibbs_core.cpp`), driven by
`run_prior_module_count_model.R` Sections 1–8.

### 3.1 Specification

Per species, the herd count is Negative-Binomial with a log-mean that combines fixed driver effects, a
country random intercept, and an area offset:

$$
Y_{rt}\mid \mu_{rt} \;\sim\; \mathrm{NB}\!\left(\mu_{rt},\, \rho\right),
\qquad
\log \mu_{rt} \;=\; \underbrace{x_{rt}^{\top}\beta}_{\text{drivers}} \;+\; \underbrace{u_{c(r)}}_{\text{country RE}} \;+\; \underbrace{\log A_{rt}}_{\text{offset}},
$$

with dispersion $\rho>0$ (the NB "number of failures"; $\mathrm{Var}=\mu+\mu^2/\rho$). The offset
$\log A_{rt}$ enters the linear predictor but is **not** part of the working response — an earlier
sign/placement error here (offset added to $\eta$ but never subtracted from the Pólya-Gamma working
response) drove a degeneracy ($\rho\to 0$, slopes nuked) that is now fixed and truth-recovery validated.

### 3.2 Pólya-Gamma augmentation

The NB likelihood is rendered conditionally Gaussian by Pólya-Gamma (PG) data augmentation
(Polson, Scott & Windle, 2013). With latent $\omega_{rt}\sim\mathrm{PG}(Y_{rt}+\rho,\,0)$ the count
likelihood contributes a Gaussian pseudo-observation

$$
z_{rt} \;=\; \frac{Y_{rt}-\rho}{2\,\omega_{rt}} \;-\; \log A_{rt},
\qquad
z_{rt}\mid\cdot \;\sim\; \mathcal N\!\left(x_{rt}^\top\beta + u_{c(r)},\; \omega_{rt}^{-1}\right),
$$

so that $\beta$ and the random effects $u$ are drawn from conjugate Gaussian full conditionals. The
dispersion $\rho$ is updated by a Chinese-Restaurant-Table (CRT) Gibbs step. A spike-and-slab prior on a
learnable subset of coefficients performs variable selection; a country random intercept
$u_c\sim\mathcal N(0,\sigma_u^2)$ with a support-preserving prior captures unmodelled national level.

### 3.3 What Stage 1 delivers downstream

For each species and country the sampler returns the posterior mean slope vector
$\beta^{(s)}_{c}=\mathbb E[\beta^{(s)}\mid u_{c}]$ (the `postb_total` array, `[driver × country]`). This is
the **level** term of the plug-in rule (§6). Because the subtype and system splits are formed by
*renormalization within region* (§6.2), the country-level intercept and dispersion cancel and only the
slopes propagate.

---

## 4. Stage 2 — composition models (symmetric multinomial logit)

Both composition stages use the same engine, `mnlogit_rcpp_sym` (`codes/mnlogit_rcpp_sym.R`), a symmetric
zero-sum multinomial-logit Pólya-Gamma Gibbs sampler with horseshoe shrinkage and country random effects.

### 4.1 The symmetric (sum-to-zero) parameterization

For $K$ categories with coefficient matrix $B=[\,b_1,\dots,b_K\,]$, the choice probability is the softmax

$$
P(k\mid x) \;=\; \frac{\exp(x^\top b_k)}{\sum_{j=1}^{K}\exp(x^\top b_j)},
\qquad\text{subject to}\qquad \sum_{k=1}^{K} b_k = 0 .
$$

The zero-sum constraint removes the reference-category indeterminacy of the standard MNL: every category's
coefficient is a **contrast against the category average**, so effects are symmetric and directly
comparable (no baseline is privileged). The composition "delta" reported for category $k$ is exactly
$b_k$ under this constraint. Estimation standardizes $x$ internally and back-transforms draws to the raw
scale, so the returned contrasts are in the same units as the Stage-1 slopes — a prerequisite for the
additive plug-in (§6).

### 4.2 Stage 2a — D/O/F subtype composition (`composition/fit_composition.R`)

For each species, the NUTS2 subtype head counts $(n^{\mathrm D}_r,n^{\mathrm O}_r,n^{\mathrm F}_r)$ follow a
3-category symmetric MNL in the curated drivers:

$$
(n^{\mathrm D}_r,n^{\mathrm O}_r,n^{\mathrm F}_r)\;\sim\;\mathrm{Multinomial}\!\big(N_r,\,\boldsymbol\pi_r\big),
\qquad
\pi_{rg}=\frac{\exp(x_r^\top b_g)}{\sum_{g'}\exp(x_r^\top b_{g'})},\quad \sum_g b_g = 0.
$$

The sum-to-zero contrasts are recovered per posterior draw as
$\delta_{\mathrm D}=b_{\mathrm D}-\bar b,\ \delta_{\mathrm O}=b_{\mathrm O}-\bar b,\ \delta_{\mathrm F}=-\bar b$
with $\bar b=\tfrac13\sum_g b_g$. **Horseshoe** priors on the slopes (`use_horseshoe=TRUE`,
`symmetric_hs=TRUE`) yield either a visibly non-zero effect (posterior CI excludes 0) or a *solid zero*
(tight CI at 0), avoiding the collinear ± compensation that a flat prior produces. A **country random
intercept** absorbs national level. Prior to fitting, drivers are **decorrelated** by greedy pruning until
$\max|\mathrm{corr}|<0.7$, protecting all land-use classes and GDP/Pop.

### 4.3 Stage 2b — organic-system composition (`composition/fit_organic.R`)

The organic split is a 2-category symmetric MNL on the organic-grassland proxy (§2.4), estimated at
**NUTS3**:

$$
(n^{\mathrm{org}}_r,n^{\mathrm{conv}}_r)\;\sim\;\mathrm{Multinomial}\!\big(N_r,\,(\pi_r,1-\pi_r)\big),
\qquad
\pi_r=\frac{\exp(x_r^\top b_{\mathrm{org}})}{\exp(x_r^\top b_{\mathrm{org}})+\exp(x_r^\top b_{\mathrm{conv}})},
$$

with $b_{\mathrm{org}}+b_{\mathrm{conv}}=0$, so the fitted contrast is
$\delta_{\mathrm{org}}=b_{\mathrm{org}}-\tfrac12(b_{\mathrm{org}}+b_{\mathrm{conv}})=b_{\mathrm{org}}$ and the
organic log-odds are $\operatorname{logit}\pi_r = 2\,x_r^\top\delta_{\mathrm{org}}$. Same machinery as §4.2
(horseshoe + country RE + decorrelation).

**Circularity exclusion.** The organic-LU columns (`*_HIO`, `*_LIO`, `*_other_O`) share the map lineage of
the response and are therefore **removed from the driver set** for this stage. The organic pattern is thus
learned from *independent* biophysical/economic drivers (climate, terrain, GDP, conventional land use),
which is what lets it generalize to scenario projections where the organic map is unavailable.

---

## 5. Estimation

- **Sampler.** Pólya-Gamma Gibbs throughout (counts in Stage 1; multinomial via the symmetric MNL in
  Stage 2). Conditionally Gaussian coefficient draws; CRT-Gibbs for NB dispersion.
- **Priors.** Horseshoe (global-local) on composition slopes; spike-and-slab selection + support-preserving
  country-RE prior in the totals model; weakly-informative Gaussian intercepts.
- **Standardization.** Drivers are centered/scaled inside the sampler and coefficients back-transformed to
  raw units, so all three stages' coefficients live in a common space.
- **MCMC settings (reference run).** Totals: `niter = 12000`, `nburn = 4000` (the dispersion and RE
  variances are slow-mixing and need the long chain). Composition: `niter = 4000`, `nburn = 2000`.
- **Convergence.** Split-$\hat R$ and bulk-ESS are computed per parameter; the totals model additionally
  writes a per-target convergence table (`<sampler>_admin_convergence_check_<label>_<target>_<date>.rds`).

---

## 6. Combining the stages — the plug-in parameter table

Implemented in `composition/build_subclass_parameters.R`; output
`output/composition/subclass_country_parameters.csv`.

### 6.1 Additive log-linear combination

The final coefficient for cell $(g,q)$ of species $s$ in country $c$ is the sum of the level and the two
contrasts:

$$
\boxed{\;\gamma^{(s)}_{c,g,q} \;=\; \beta^{(s)}_{c} \;+\; \delta^{(s)}_{g} \;+\; \operatorname{sgn}(q)\,\delta^{(s)}_{\mathrm{org}} \;+\; \tfrac12\,\operatorname{sgn}(q)\,c^{(s)}_{c}\,e_{0}\;}
$$

where $\operatorname{sgn}(\text{org})=+1$, $\operatorname{sgn}(\text{conv})=-1$, $e_0$ is the intercept basis
vector, and $c^{(s)}_{c}$ is the per-country organic anchor of §6.3. Composition/organic contrasts are
mapped onto the totals' driver names ($\texttt{lu\_area\_}X\!\to\!\texttt{log1p\_lu\_area\_}X$, etc.) before
addition, and any driver absent from a stage contributes $0$.

### 6.2 Identification under within-region renormalization

At prediction the cell weight is $w^{(s)}_{r,g,q}\propto\exp\!\big(\gamma^{(s)\top}_{c,g,q}x_r\big)$,
renormalized over the $12$ cells within region $r$. Two useful cancellations follow from the zero-sum
construction:

1. **Level/species cancels the totals baseline.** Adding $\beta^{(s)}_{c}$ to every cell of a species and
   renormalizing leaves subtype/system *shares* unchanged — the totals model sets the herd size, the
   composition contrasts set its split. (The herd size itself is applied separately as the Stage-1
   prediction.)
2. **The organic split depends only on $\delta_{\mathrm{org}}$ (plus the anchor).** For fixed $g$, the
   organic-vs-conventional log-odds within a region are
   $$
   \log\frac{w_{r,g,\text{org}}}{w_{r,g,\text{conv}}}
   = \big(\gamma_{c,g,\text{org}}-\gamma_{c,g,\text{conv}}\big)^\top x_r
   = 2\,\delta_{\mathrm{org}}^\top x_r + c^{(s)}_{c},
   $$
   because $\beta_c$ and $\delta_g$ are common to the two systems of the same subtype and cancel. Hence the
   organic share is $\operatorname{plogis}(2\,\delta_{\mathrm{org}}^\top x_r + c^{(s)}_{c})$ — identical
   across D/O/F — which is why a *single* anchor per (species, country) suffices.

### 6.3 Eurostat level anchoring (moment condition)

The fitted organic pattern reproduces *relative* within-country variation but not the correct national
level (the proxy is grassland area, not certified livestock). We therefore choose the per-country intercept
shift $c^{(s)}_{c}$ so that the area-weighted mean modelled organic share equals the Eurostat national
share:

$$
\frac{\sum_{r\in c} \big(n^{\mathrm{org}}_r+n^{\mathrm{conv}}_r\big)\,\operatorname{plogis}\!\big(2\,\delta_{\mathrm{org}}^\top x_r + c^{(s)}_{c}\big)}{\sum_{r\in c}\big(n^{\mathrm{org}}_r+n^{\mathrm{conv}}_r\big)} \;=\; \pi^{\mathrm{eur}}_{c,s}.
$$

This is one moment condition in one unknown; $c^{(s)}_c$ is obtained by a monotone root-find (`uniroot`,
weights = grassland area). The map thus supplies *where* organic farming concentrates within a country and
Eurostat supplies *how much* it is in aggregate — the level is exact by construction. Validation on the
reference data recovers modelled-share $=$ Eurostat share to three decimals for every country/species.

### 6.4 Output schema

`subclass_country_parameters.csv` has one row per (species, country, subclass, system, driver):

```
species, country, subclass, dof, system, driver, beta_total, delta_dof, delta_org, gamma
```

`subclass` $\in$ {BOVD, BOVO, BOVF, SGTD, SGTO, SGTF}; `system` $\in$ {conventional, organic}. The organic
dimension is an **explicit column**, eliminating the earlier suffix-collision that mislabeled the meat
subclass (`…O`) as organic.

---

## 7. Prediction / downscaling

Given grid-cell drivers $x$ in region $r$ of country $c$:

1. **Herd size.** Predict the species total from Stage 1 (NB mean with the region's offset), then reconcile
   to the region's known/GLOBIOM total.
2. **Composition weights.** For each of the $12$ cells, $w_{g,q}=\exp(\gamma^{(s)\top}_{c,g,q}x)$.
3. **Allocate.** Cell count $=$ herd size $\times\, w_{g,q}/\sum_{g',q'} w_{g',q'}$.

Because every $\gamma$ is a posterior summary of a coherent Bayesian model, the allocation inherits the
totals' spatial fidelity, the composition's driver responses, and — for organic — an exact match to the
Eurostat national accounts.

---

## 8. Reproducibility

Run end-to-end from the repo root (`gamble-core`). The count driver executes Stages 1–2 and the plug-in
in one pass:

```bash
DRIVER_RUN_MODE=production DRIVER_FIT_COMPOSITION=TRUE \
  Rscript run_prior_module_count_model.R
```

which internally runs Section 9:

```
prep/prepare_composition_training.R      # D/O/F (NUTS2) + organic (NUTS3) training tables
composition/fit_composition.R  (D_SPECIES=bov)  # δ_DOF, BOV
composition/fit_composition.R  (D_SPECIES=sgt)  # δ_DOF, SGT
composition/fit_organic.R                       # δ_org (NUTS3, off the organic map)
composition/consolidate_composition.R           # D/O/F report
composition/build_subclass_parameters.R         # β_total + δ_DOF + δ_org (Eurostat-anchored) -> 12 cells
```

**Naming.** Outputs follow the pixel-model convention:
`MODEL_LABEL = <units>_<sampler>_<scheme>_<mode>` (e.g. `NUTS3_count_rcpp_GLOBIOM_production`), saved-model
dirs `saved_model_outputs/<MODEL_LABEL>_<target>_RE_<re_col>`, design `dat_admin_FULL_<MODEL_LABEL>_<ts>.rds`
(cf. `dat_pixel_FULL_…`), flat outputs `<sampler>_admin_beta_median_…` / `<sampler>_admin_convergence_check_…`
(cf. `<sampler>_pixel_…`).

**Key environment switches.** `DRIVER_CLASS_COLS` (classification), `DRIVER_NITER`/`DRIVER_NBURN`,
`DRIVER_YEARS`, `ORG_ANCHOR_YEAR` (Eurostat organic year, default 2010), `D_EXCL_ACCESS`,
`D_COVARIATES`/`D_PROTECT` (composition driver control).

---

## 9. Identifying assumptions and limitations

1. **Organic proxy.** The organic split's *spatial pattern* is identified from organic grassland area, not
   from certified livestock. If organic grazing intensity differs systematically from organic grassland
   share within a country in a way correlated with the drivers, the within-country pattern is biased; the
   *national level* is nonetheless exact (§6.3). Sheep/goats and cattle share the same grassland proxy but
   receive species-specific national levels.
2. **Cut, not full Bayes.** Feedback from composition to totals is severed by construction. This is a
   deliberate robustness choice; a fully joint model would be more efficient but would let a mis-specified
   split perturb the totals and would break Pólya-Gamma conjugacy.
3. **Cross-source resolution.** D/O/F is NUTS2 (Eurostat), organic is NUTS3 (map). The plug-in treats the
   coarser NUTS2 subtype contrasts as constant within their NUTS2 parent — acceptable because subtype
   composition varies smoothly, but it caps the subtype split's spatial resolution.
4. **Static Eurostat organic level.** The 2010 anchor is applied to the reference build; the temporal rule
   (2000 organic $\approx 0$; 2010 backcast; 2018/2020 direct) governs other years.
5. **Confounder exclusions.** Dropping the human-footprint cluster from the composition drivers trades some
   predictive fit for interpretable, non-separated contrasts; re-including it (`D_EXCL_ACCESS=FALSE`) is
   available for sensitivity analysis.

---

## References

- Polson, N., Scott, J., & Windle, J. (2013). *Bayesian inference for logistic models using Pólya–Gamma
  latent variables.* JASA 108(504).
- Plummer, M. (2015). *Cuts in Bayesian graphical models.* Statistics and Computing 25(1).
- Carvalho, C., Polson, N., & Scott, J. (2010). *The horseshoe estimator for sparse signals.* Biometrika
  97(2).
- Zhou, M., & Carin, L. (2015). *Negative binomial process count and mixture modeling.* IEEE PAMI 37(2).

*Companion docs: [`nested_cut_model.md`](nested_cut_model.md), [`livestock_subtype_allocation_spec.md`](livestock_subtype_allocation_spec.md), [`DEVELOPMENT.md`](DEVELOPMENT.md).*
