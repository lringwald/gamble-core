# Livestock subtype allocation — post-processing spec

Carve the spatially-allocated **total** livestock (output of the count model, already
renormalized to GLOBIOM coarse totals) into **12 subtype layers** using observed Eurostat
shares. The count model supplies *all* spatial signal; Eurostat supplies *all* composition.
No land-use proxies, no fitting — just observed regional shares applied as a cascade.

```
species s ∈ {BOV, SGT}
  × production type ∈ {D (dairy), O (meat/other), F (followers)}   ← Eurostat apro_mt categories
    × management ∈ {organic, conventional}                          ← Eurostat org_lstspec
= 12 cells per timestep
```

## Inputs

| input | source | resolution |
|---|---|---|
| `total_{s,i,t}` allocated LSU | count model + downscaling renorm | NUTS3 unit `i`, year `t` |
| BOV categories | `apro_mt_lscatl` | NUTS0 (country); dairy cows also NUTS2 via `agr_r_animal` |
| SGT categories | `apro_mt_lssheep` + `apro_mt_lsgoat` | NUTS0; sheep/goats NUTS2 via `agr_r_animal` |
| organic animals | `org_lstspec` | NUTS0 |
| region crosswalk | NUTS3 → NUTS2 → NUTS0 | — |

Model timesteps: **t ∈ {2000, 2010, 2018/2020}** (outcome years 2000/2010/2020; LU covariates
use the 2018 LUM map for the 2020 tier). Share-year `τ(t)` = the outcome/population year.

## Stage 0 — category → {D, O, F} mapping

**BOV** (`apro_mt_lscatl`):
- **BOVD** = dairy cows
- **BOVO** = non-dairy (suckler) cows
- **BOVF** = followers = heifers + bulls + calves + young (all non-cow categories)

**SGT** (`apro_mt_lssheep` + `apro_mt_lsgoat`, pooled):
- **SGTD** = milk/dairy ewes + dairy she-goats
- **SGTO** = non-dairy (meat) ewes + meat she-goats
- **SGTF** = lambs + kids + other young

> SGT caveat: the dairy/meat split is explicit for cattle ("dairy cows") but patchier for small
> ruminants — use the milk-ewe / dairy-goat categories where reported, else infer SGTD from
> sheep/goat **milk production** (`apro_mk_*`) or hold the country dairy-ratio. It is the
> least-clean number in the scheme; dairy small-ruminants are strongly Mediterranean so the
> country/NUTS2 share carries the signal.

## Stage 1 — production split (best available resolution)

For species `s`, year `t`, region `r`:
1. Country D/O/F shares from the Stage-0 category sums: `pD, pO, pF` (sum to 1).
2. If NUTS2 dairy available (`agr_r_animal`): take `pD` at **NUTS2**, split the remainder
   `(1 − pD)` into O:F by the **country** O:F ratio. (dairy enters at finer resolution.)
3. Apply: `X_{s,type,i,t} = total_{s,i,t} · p_{type, r(i), t}`.

Because shares sum to 1, `BOVD+BOVO+BOVF = total` per unit, and the regional subtype total
= `share × GLOBIOM total` matches Eurostat by construction. **No IPF needed.**

## Stage 2 — organic split (country)

1. `o_{s,r,t}` = organic share = `org_lstspec` / total, per species, country (SGT pools
   organic sheep + goats).
2. Apply within each production type: `X_{s,type,organic,i,t} = X_{s,type,i,t} · o_{s,country(i),t}`;
   conventional = complement.

> Assumption: organic ⊥ production type within a country (same organic fraction on D, O, F)
> unless `org_lstspec` cross-tabs organic by category — it skews toward suckler/extensive, so
> document it. If a cross-tab exists, use type-specific organic shares.

## Temporal handling (timesteps 2000 / 2010 / 2018-2020)

Per (attribute, region, t) pull Eurostat at `τ(t)`; resolve gaps with a documented ladder:
1. exact `τ(t)`,
2. nearest year within ±2,
3. country → EU-aggregate share,
4. earliest/latest available, flagged as extrapolated.

| attribute | 2000 | 2010 | 2018/2020 |
|---|---|---|---|
| production D/O/F (`apro_mt`) | available | available | available |
| organic (`org_lstspec`) | **set to 0** (reporting begins ~2012; assume fully conventional) | **2012 share × backcast factor** | available (direct) |

**Organic pre-2012 is the main temporal gap. Decision (locked):**
- **2000 → organic share = 0** for every species/region. Organic livestock was marginal before
  the CAP organic-support era; treat the 2000 herd as fully conventional.
- **2010 → backcast the 2012 country share**, slightly down-weighted to reflect ~2 years of
  organic growth (organic *grew* into 2012, so 2010 < 2012):
  `o_{c,2010} = o_{c,2012} · ( I_EU,2010 / I_EU,2012 )`
  where `I_EU` is the EU organic index (organic-livestock numbers if available, else organic-area;
  Eurostat `org_lstspec` / `org_cropar` aggregated to EU). This keeps each country's 2012
  cross-sectional structure and scales the level back by the common growth factor.
  *Fallback if no index handy:* a fixed `o_{c,2010} = 0.85 · o_{c,2012}` (≈ 2 yr at ~8%/yr growth).
- **2018/2020 → direct** from `org_lstspec`.

Production shares are safe at all three timesteps. Use `τ(2020)=2020` for population/organic
(both reported in 2020) even though LU covariates use 2018 — the shares describe the livestock
population, not land use.

## Coverage & fallbacks (must run BEFORE allocation)

Model domain = **34 countries / 1380 NUTS3 units** (NUTS3 codes → NUTS2 → NUTS0 by substring).
Eurostat does **not** cover all of them, so a coverage gate + fallback ladder is mandatory.

Known gaps (verify against the actual extracts):
- **Western Balkans — AL, BA, MK, RS/ME:** likely NO `org_lstspec` and NO `agr_r_animal` (NUTS2);
  national `apro_mt` category breakdown may be missing, esp. 2000.
- **UK:** OK for 2000/2010; `org_lstspec`/category data may be **missing for 2020** (post-Brexit).
- **NO, CH (EFTA):** national `apro_mt` usually present; organic/regional partial.
- **EU-27:** national production shares OK all years; NUTS2 dairy EU-mostly; organic via pre-2012 rule.

**Fallback ladder for production D/O/F** — progressive-widening AVERAGE of *shares* (each member
sums to 1; average the full D/O/F vector, renormalize). First rung meeting a minimum support `K`
wins; the rung + pool size are LOGGED per cell (`src`, `n_pool`). Country-specific real data is
exhausted before spatial averaging:
1. **Eurostat exact** (country × year).
2. **Eurostat same country, nearest years** (±1, ±2 … → country time-mean) — temporal widening
   (fills UK-2020, accession-2000-if-reported-later, odd-year SGT).
3. **FAOSTAT country** — terminal *country-specific* source: real cattle/sheep/goat stocks +
   dairy-animal counts, covers AL/BA/MK/RS/NO/CH and pre-2000. Gives D directly; (1−D) split into
   O:F by the borrowed demographic ratio (rung 4 pool). This is what fills the structural
   non-reporters with their ACTUAL composition rather than a neighbour guess.
4. **Eurostat neighbour-cluster mean** (livestock-system groups: Nordic, Baltic, SE-Europe/Balkans,
   Mediterranean, Alpine, Western, Central-East), then broader region — spatial average, only where
   even FAO is blank.
5. **Eurostat published EU aggregate** (`EU27_2020`/`EU28`) — absolute last resort (rarely reached).
6. *(projection years only)* **GLOBIOM regional D/O/F** — native source forward in time.

**NUTS2 dairy missing:** drop to the country D share (already the default).

**Organic missing** (country/year, incl. all Balkans, 2000 everywhere, UK-2020): **organic share = 0**
(assume conventional) — same convention as the 2000 rule; do NOT borrow another country's organic.
Organic has no FAO equivalent, so its ladder is Eurostat-only: nearest-year → 0.

**Validation gate (assert, don't silently fill):** before allocating, check every one of the
1380 units × 3 years resolves to a production share and an organic share (post-fallback). Emit a
coverage table `[country × year × attribute → source/fallback used]`; fail loud if any unit hits
the terminal fallback unexpectedly. Shares are ratios, so Eurostat↔GLOBIOM *level* mismatches are
harmless (we renormalize to GLOBIOM totals) — the only thing that breaks the cascade is a **missing
category**, which the ladder must cover.

## Share-uncertainty propagation (consistent subtype posteriors)

The shares are *estimated from counts*, not exact, so the 12 cells should carry uncertainty —
wider where a share rests on few animals or a fallback rung, tight where it's a large direct read.
Propagate it draw-by-draw alongside the total's posterior.

**Distributions (per region × year):**
- **Production D/O/F** ~ `Dirichlet(α_D, α_O, α_F)`, with `α = effective head counts + ε` (ε ≈ 0.5
  Jeffreys floor). Direct Eurostat: `α` = the actual category head counts (`A2300F`, `A2300G`,
  residual; or milk-ewes etc.) — large counts ⇒ tight. The mean is the share we already compute.
- **Organic** ~ `Beta(α_org, α_conv)`, `α_org = organic heads + ε`, `α_conv = (total − organic) + ε`.
  The 2000 = 0 and missing → 0 cells are **degenerate (point 0)**, not Beta — no organic mass.

**Concentration by fallback rung (this is where `n_pool` earns its keep):** rescale `α` to an
*effective* sample size so borrowed shares are appropriately vague:
| rung | effective `α` |
|---|---|
| Eurostat exact / FAOSTAT | full observed head counts (tight) |
| Eurostat nearest-year / time-mean | observed counts, lightly down-weighted |
| Eurostat cluster mean | `α` rescaled to a small effective N reflecting between-country spread in the cluster (flat ⇒ wide) |
| Eurostat EU aggregate | minimal effective N (very flat) |
Operationally: keep the share **mean** from the resolved rung, set the Dirichlet/Beta **total
concentration** = the rung's effective sample size (logged `n_pool`). Low rung ⇒ low concentration
⇒ honest, wide subtype intervals exactly where the data is thin (Balkans, Norway, sparse SGT).

**Per-draw algorithm** (operates on the saved total posterior + the share tables):
```
for each posterior draw s of the total allocation total[species, i]^(s):
    for each region r, year t:
        π[D,O,F]_r^(s)  ~ Dirichlet(α_DOF, r, t)        # production split
        o_r^(s)         ~ Beta(α_org, α_conv, r, t)      # organic split (0 where degenerate)
        # NUTS2 dairy refinement: draw dairy at NUTS2 from its own Beta, split (1−D) by the
        #   country O:F Dirichlet — each with its own α — then renormalize the 3-vector.
    cell[s,type,mgmt,i]^(s) = total[s,i]^(s) · π[type]_r(i)^(s) · {o or 1−o}_r(i)^(s)
```
This yields a full posterior per cell that, **draw-by-draw**, sums to the total and matches the
regional shares in expectation, with intervals that widen automatically over the fallback regions.

**Data requirement:** the wrangler must emit the underlying **head counts** (not just the
normalized shares) plus the resolved `rung`/`n_pool` per region × year × category, so the
allocation can reconstruct `α`. (Add count columns to `production_shares_nuts0.csv` /
`organic_shares_nuts0.csv`, or a parallel `*_counts.csv`.)

## Driver-conditioned composition (optional) — the merged `γ_j = β + δ_j` parameterization

The flat regional shares above make subtypes differ only at NUTS2/country granularity (no
within-region differentiation — within a region every subtype is the total's pattern × a constant).
To let subtypes differ *on the grid* by drivers — without fine subtype data and without overriding
the Eurostat regional anchors — condition the composition on the count model's drivers, written as
a **shared level + sum-to-zero contrasts**. This is the multinomial–NB factorization, and it
*merges* the count `β` and the share coefficients into one field.

**Parameterization.** Per subtype cell `j`:
```
η_ij = offset_i + X_i · γ_j ,   γ_j = β + δ_j ,   Σ_j δ_j = 0   (per split: D/O/F; and organic 2-way)
```
- **β** = the count model's fitted slopes — the LEVEL (*where the animals are*). **Shared, taken
  as-is from the count model, NOT refit.**
- **δ_j** = sum-to-zero subtype CONTRASTS — the COMPOSITION (*which subtype*). The only new params.
- Total `Σ_j exp(η_ij)` recovers the NB count model (to first order); composition = `softmax(X_i δ_j)`.
- (Same zero-sum/Helmert contrast structure as `mnlogit_rcpp_sym` — reuse that machinery.)

**Fit δ (region level, regularized).** Response = observed regional share vectors (D/O/F at NUTS2
where available, else country; organic at country). Predictors = the same `X`, livestock-weighted-
aggregated to the region. Model = Dirichlet / multinomial-logit (additive-logratio) with a
**horseshoe on δ**; ≈280 NUTS2 obs for dairy, ≈30 country obs for organic, so keep it sparse.
Put `β·X` in as an **offset** so δ captures only the *differential* driver effect beyond the level.

**Predict + pin.** Per grid cell `π_cell = softmax(X_cell · δ)` (fine drivers), then **renormalize
within each region** so `Σ_cell total_cell·π_cell = observed Eurostat regional subtype total`. The
regional level stays exactly the observed number; δ only adds within-region texture.

**Resolution split.** β is identified at fine scale (NUTS3 totals); δ only at region scale (NUTS2
dairy / country organic & F). Organic δ is thinnest → may stay flat (`δ=0` ⇒ flat country share) if
drivers don't help.

**Uncertainty.** β posterior (count model) × δ posterior (contrast fit) ⇒ per-cell composition
posterior; where drivers explain composition it tightens, else it falls back toward the regional
Dirichlet (the flat-share uncertainty). Feeds the same per-draw cell propagation.

**Caveats.** Ecological inference (region→cell relationship; renormalization limits the damage);
thin training ⇒ regularize hard. Note this *estimates* the dairy↔driver association (e.g. pasture
intensity) rather than *assuming* it — reconciling the earlier rejection of the intensity proxy.

**Heavier alternative — full joint.** Estimate `γ_j = β + δ_j` in ONE multi-resolution hierarchical
multinomial-NB (NUTS3 total counts pin `β`; regional shares pin `δ`). Gives a single coherent joint
posterior over level + composition, but heavy and **no extra spatial identifiability** (composition
is still region-limited). Use only if a fully coherent joint posterior is required; the modular
**β-fixed + δ-fit** version gives the same predictions far more cheaply.

**Default ladder:** flat regional shares (simplest) → driver-conditioned δ (when within-region
differentiation matters) → full joint (only if a single joint posterior is needed). Choose per
attribute — e.g. driver-conditioned for dairy (rich NUTS2 data), flat for organic (thin).

## Outputs & coherence checks

12 NUTS3 layers per timestep: `{BOV,SGT} × {D,O,F} × {organic,conventional}`.
- `Σ_{type,mgmt} X = total_{s,i,t}` (± rounding).
- `Σ_{i∈r} X_{type} = ` Eurostat regional type total (by construction).
- every cell `≥ 0` and `≤` its parent (`organic-dairy ≤ dairy ≤ total`).

## Projection (≤ 70 yr)

Same cascade, swapping the data source for the *level* only:
- production D/O/F totals ← **GLOBIOM** (it emits the dairy/other/follower herds directly);
- organic share ← **scenario trajectory** (held, or trended toward a policy target);
- spatial shape ← unchanged (the count-model total + the cascade).
Eurostat calibrates the historical shares; GLOBIOM carries them forward.
