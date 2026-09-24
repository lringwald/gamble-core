# GLOBIOM_subclass — prior module

Flat multinomial logit on the **GLOBIOM_subclass** taxonomy: 27 modelled classes under the curated
6-nest tree (Cropland / Forests / Natural / Pasture / Urban / Waterbodies), 10 km EU grid, country
random effects, focal neighbourhood composition lagged to 2010.

Sibling of `projects/BMLEH_Los1_CAPRI/` — same segmented, resumable mechanics — but carrying the
**standard sampler configuration** adopted 2026-09-16
(`docs/sampler_model_specification.md` §8.0) rather than BMLEH's own settings.

Everything runs **from the gamble-core root**.

---

## Run it

```bash
cd "/Users/leopoldringwald/Library/CloudStorage/GoogleDrive-leopold.ringwald@gmail.com/My Drive/R/gamble-core"

# 0. sanity-check the sampler first (37 checks, ~1 min) -- do this after ANY change to codes/
Rscript codes/test_suite_lu_pixel.R
#    or the whole repo, with a report + figure in output/report/:
Rscript tests/run_all.R

# 1. build the design  (~30 min)  -> output/pixel_model_inputs_GLOBIOM_subclass.rds
projects/GLOBIOM_subclass/gamble_model/run.sh design

# 2. segment 1: 4000 burn-in + 1000 kept, 4 chains  (the expensive one -- carries the burn-in)
projects/GLOBIOM_subclass/gamble_model/run.sh fit

# 3. extend by 500 sweeps, as many times as the diagnostics ask for
projects/GLOBIOM_subclass/gamble_model/run.sh more
projects/GLOBIOM_subclass/gamble_model/run.sh more
```

A **10-minute end-to-end check** before committing to the full run:

```bash
projects/GLOBIOM_subclass/gamble_model/run.sh smoke     # design + 5,000-pixel, 1-chain fit
```

Long runs survive a closed terminal:

```bash
nohup projects/GLOBIOM_subclass/gamble_model/run.sh fit \
  > output/globiom_fit_$(date +%Y%m%d_%H%M).log 2>&1 &
echo "pid $!"; tail -f output/globiom_fit_*.log
```

Kill by **PID**. `pkill -f estimate_prior.R` takes down every project's fit at once, including
BMLEH's.

---

## Batches, and restarting

Every segment writes, under `results/gamble_model/GLOBIOM_subclass/prior_<date>_<time>/`:

| file | what it is |
|---|---|
| `posterior/posterior_batch_*_chain_*.qs` | the streamed draws (the posterior itself) |
| `posterior/final_state_chain_<id>.qs` | each chain's raw internal state — what `more` resumes from |
| `lineage.txt` | segment number, design, parent segment, sweep budget |
| `run_config.rds` | **every sampler switch actually used** — read this before quoting a number |
| `prior_fit.rds`, `per_class_fit.csv` | held-out scores, pooled across all chains and segments |

`more` picks the newest **resumable** segment — one that finished and left chain state — and skips
empty husks from aborted runs. Resuming is near-exact: on BMLEH a resumed chain restarted 91 nats
from where it stopped against an 83-nat draw-to-draw noise floor, versus 2,970 for a cold start.

Scoring pools **every chain and every ancestor segment**, so the posterior is the union of the
segments rather than just the newest one.

```bash
# re-score what is already on disk, without sampling (works on a run killed mid-flight)
projects/GLOBIOM_subclass/gamble_model/run.sh score

# continue a SPECIFIC segment rather than the newest
projects/GLOBIOM_subclass/gamble_model/run.sh more results/gamble_model/GLOBIOM_subclass/prior_2026-09-16_1200
```

**Keep the design fixed across segments.** Chain state is indexed by covariate and class, so
resuming onto a different dump is meaningless. `estimate_prior.R` refuses a dump whose class set
isn't GLOBIOM_subclass (`GB_ALLOW_STALE_DESIGN=TRUE` overrides, diagnostics only).

---

## Measurement arms

The baseline is the standard config. Each arm changes exactly one thing, and each needs its **own
output directory**, which `GB_OUT` provides:

```bash
# BART treatment -- the headline methods claim
GB_BART=TRUE GB_OUT=results/gamble_model/GLOBIOM_subclass/arm_bart \
  projects/GLOBIOM_subclass/gamble_model/run.sh fit

# RE-scale ASIS off -- BMLEH says ASIS is what makes the posterior quotable, an earlier GLOBIOM
# node arm says it HURTS RE-variance ESS (73 -> 27). This settles it on the real design.
GB_RE_ASIS=FALSE GB_OUT=results/gamble_model/GLOBIOM_subclass/arm_noasis \
  projects/GLOBIOM_subclass/gamble_model/run.sh fit

# re_support_strength 0 -- the one knob §8.1b leaves deliberately unreconciled. BMLEH's evidence
# for the standard was produced at 0; flat and nested now run 1.
GB_RE_SUPPORT=0 GB_OUT=results/gamble_model/GLOBIOM_subclass/arm_resupport0 \
  projects/GLOBIOM_subclass/gamble_model/run.sh fit

# RE block size -- the country intercept carries nearly all the RE skill
GB_RE_VARS=intercept GB_OUT=results/gamble_model/GLOBIOM_subclass/arm_icpt \
  projects/GLOBIOM_subclass/gamble_model/run.sh fit

# 5 km instead of 10 km -- MAUP / scale dependence. Needs its own DESIGN (~4x the pixels).
DRIVER_PIXEL_RES=5 projects/GLOBIOM_subclass/gamble_model/run.sh design
```

| env | default | what it does |
|---|---|---|
| `GB_NITER` / `GB_NBURN` / `GB_THIN` | 5000 / 4000 / 1 | sweep budget = 4000 burn + 1000 kept (`more` uses 500 / 0 / 1) |
| `GB_CHAINS` / `GB_CORES` | 4 / cores−2 | chains run in parallel, one process each |
| `GB_NPIX` / `GB_TESTFRAC` | 0 (all) / 0.2 | subsample cap; held-out fraction |
| `GB_SPLIT` | `hybrid` | `hybrid` \| `random` \| `country` \| `block` |
| `GB_BLOCK_MULT` | `5` | block side in units of the 3-cell focal window (5 → 150 km at 10 km) |
| `GB_BUFFER` | `1` | cells around each test block dropped from training |
| `GB_BLOCK_SHARE` | `0.5` | share of the hold-out taken as blocks rather than random |
| `GB_RE_VARS` | `log1p_GDP,log1p_Pop,GHM_HI,CISI` | RE covariates (the intercept is always included) |
| `GB_RE_ASIS` | `TRUE` | RE-scale interweaving |
| `GB_RE_SUPPORT` | `1` | RE sparse-group gate |
| `GB_SLAB_C2` / `GB_SLAB_C2_VAL` | `TRUE` / 4 | sample the RE slab cap, or fix it |
| `GB_BART` | `FALSE` | BART treatment arm |
| `GB_SCORE_POOL` | `lineage` | `self` scores this segment only |

---

## Caveats worth knowing before quoting numbers

- **`thin = 1` is deliberate.** Thinning never improves ESS, it only discards draws, and storage is
  not the constraint (~60 MB per 500 draws × 4 chains).
- **`estimate_slab_c2` defaults to sampled here**, unlike BMLEH, which fixes it. BMLEH's reason is
  design-specific: there the slab never binds (cap σ ≤ 3.8 against σ_re 0.108). Whether that holds
  on 27 GLOBIOM classes is untested — `GB_SLAB_C2=FALSE` is the arm.
- **The hold-out is hybrid, and gives you two numbers.** Half is country-stratified random
  (interpolation), half is whole spatial blocks (transfer), both scored from the same fit. The
  printed `OPTIMISM` line is the gap between them — how much a random split flatters this model, on
  this design. Quote the block number for anything about prediction in unseen areas.
- **Block geometry is derived, not chosen.** The focal covariate is a queen-8 neighbourhood mean
  (hardcoded in `compute_focal_coord`), so a held-out pixel's own predictors are built from its 8
  neighbours' responses. Blocks are `GB_BLOCK_MULT × 3` cells (default 5 → 150 km at 10 km, 75 km
  at 5 km), and a 1-cell ring around every block is **dropped from training entirely** — that
  buffer, not the block size, is what removes the leakage.
- **`GB_BLOCK_MULT = 5` is a judgement, not a measurement.** The right block exceeds the residual
  autocorrelation range, which nobody here has measured — and with the focal in the model the
  residual range may be far shorter than the raw land-use range. Run 3 vs 5 vs 10 as a sensitivity
  arm before leaning on the block number.
- **BART is not in the baseline** on purpose. The GLOBIOM gate found it worth +249 nats held-out at
  ~1.7× fit time, concentrated in abundant terrain-driven classes — so it is a treatment, and the
  linear baseline goes first.
