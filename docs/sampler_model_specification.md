# MNL sampler — full model specification

What `mnlogit_rcpp_sym` fits when every switch is enabled, and which switches you should
actually enable. Written 2026-08-26 against 106 formals across eight subsystems.

**§8 reconciled against the code on 2026-09-16.** Defaults live in the two production callers, not
here; §1-7 describe the model, §8 tells you where the current settings actually are.

Companion to `nested_cut_model.md` (the nesting layer above this sampler) and
`livestock_subtype_methodology.md` (the count/composition side).

---

## 1. Observation model

For pixel `i`, class `k`, group `g(i)`, with `p_all = K` alternatives and baseline class `b`:

    Y_i.  ~  Multinomial( N_i , softmax_k(eta_ik) )

Counts are rescaled by `y_weight` (default `1 / min(rowSums(Y))`).

**The sampler carries a ONE-VS-REST decomposition, not a joint softmax.** Each non-baseline
equation `ip` is Polya-Gamma augmented with predictor `psi_j - c_j`, where

    c_j = log sum_{k != j} exp(psi_k)

`c_j` is not cosmetic. Any term that enters `psi_k` for *every* alternative also moves `c_j`, so a
working response built from a stale `c_j` is inconsistent with the parameter being estimated. That is
exactly what produced the long-standing "alt_spec delta is biased by a factor 2" result: the additive
utility channel was gated off, `c_j` was computed without delta, and the conditional was wrong. See
§6 and `codes/mnlogit_rcpp_sym.R` (guard message on `alt_spec_Z`).

## 2. Linear predictor

    eta_ij = x_i'( mu_j + b_{g(i),j} * gamma_{g(i),j} )
             + alpha * f_j^BART(x_i)
             + sum_b  delta_b  o  z^b_ij

| term | meaning |
|---|---|
| `x_i` | k case-specific covariates: intercept, terrain, climate, socio, plus const-sum blocks (`OC_TOP`, `ROO`, `AWC_TOP`, `VS`, `focal_*`) |
| `mu_j` | pooled fixed effects |
| `b_{g,j}` | country random effects; `gamma` = spike-and-slab inclusion indicators |
| `f^BART` | optional tree ensemble, tempered by `alpha` |
| `delta_b o z^b` | alternative-specific blocks (conditional-logit terms) |

A **case-specific** covariate has one value per pixel and a coefficient per class. An
**alternative-specific** covariate has a value per (pixel, class) and one coefficient (or one per
class). Spatial Y-lags and t-1 land-use state are naturally the latter.

## 3. Fixed effects

    mu_{v,j} ~ N( m_{v,j} , A0 )          A0 = 2 (slopes), 100 (intercept)

Under `use_horseshoe`, a **regularized (Finnish) horseshoe** is added as a harmonic mixture of
precisions:

    prior precision  = 1/(tau_f(v)^2 lambda_v^2) + 1/c^2
    lambda_v^2 | nu  ~ IG( (m+1)/2 , 1/nu_v + ||gamma_v||^2 / (2 tau^2) )
    tau_f^2   | xi   ~ IG( (D_f m + 1)/2 , 1/xi_f + sum ||gamma||^2/lambda^2 / 2 )   per family f
    c^2              ~ IG( (nu_slab + D m)/2 , ... )

Makalic-Schmidt auxiliaries throughout (`nu`, `xi`, `zeta`), so every scale draw is conjugate.

Three modifiers:

- **`symmetric_hs`** — replaces the diagonal penalty with `c_v * M_sym` (`M_sym = I - 11'/p_all`),
  which is exactly baseline-invariant; const-sum blocks get `kron(M_sym, M_b)`. Verified to machine
  precision by suite test C7a: the symmetric quadratic form is constant across all baselines
  (spread 1.3e-15) while the diagonal one moves 3.89x.
- **`hs_groups`** — one `(tau, xi)` per covariate FAMILY and per const-sum block, instead of one
  global scale over every shrunk column. Fitted scales span 620x across families and 80x across
  blocks, so the single-scale assumption was measurably wrong.
- **`fe_support_strength`** — multiplies the horseshoe precision by `(n/PR)^s`. See §8.

## 4. Random effects

    b_{g,j,v}    ~ N( 0 , (sigma_v * r_{v,g} * tau_g)^2 )
    sigma_v      ~ Half-Cauchy(0, A), capped by a Finnish slab (c^2 estimated by default)
    gamma_{v,j,g}~ Bernoulli(pi)                     spike-and-slab inclusion
    b_{.,g}      ~ CAR(W, rho) across groups         if use_car and country_adjacency given

- `r_{v,g}` — **participation-ratio support factor**, computed once from X on within-group-centred
  columns so it measures SLOPE information. Deterministic, therefore identified; the C++ draw and
  the variance update use the same factor. This is the well-conditioned way to stop data-sparse
  groups contributing free random effects.
- `tau_g` — per-group shrinkage factor (count sampler only, experimental; see §8).
- Under `symmetric_hs`, `sigma` is POOLED across categories. The variance must then be estimated on
  UNCENTRED deviations (`re_prec_center = FALSE`, now the automatic default): the draw is uncentred
  and nothing constrains it to the zero-sum subspace, so a centred sufficient statistic zeroes the
  random effects outright.

## 5. Alternative-specific blocks

    delta      ~ N(0, alt_spec_prior_sd^2 I)
    delta | .  ~ N( P^-1 b , P^-1 ),   P = Lambda0 + sum_ip D_ip' W_ip D_ip

Exact conjugate draw: under PG each equation is a Gaussian working regression with known weights
`omega`, so delta is an ordinary weighted-regression coefficient. All blocks are drawn JOINTLY,
because a spatial and a temporal lag are correlated and a one-at-a-time scan would mix slowly along
that ridge.

**The baseline handling differs by coefficient type, and getting it wrong silently attenuates:**

| `coef` | regressor | why |
|---|---|---|
| `"shared"` | `z_ij - z_ib` | `delta*z_ij - delta*z_ib` factorises |
| `"per_class"` | raw `z_ij` | `delta_j z_ij - delta_b z_ib` does NOT factorise; `delta_b` pinned at 0 |

Input forms: a matrix (one shared delta), a named list of matrices (one shared delta each), or a
named list of `list(Z=, coef=)`.

## 6. Sweep order

    PG weights -> coefficients (C++ symmetric or diagonal kernel) -> ASIS interweave
    -> RE variances -> spike-slab -> CAR precision -> horseshoe -> delta (conjugate)
    -> BART trees -> store

## 7. The nesting layer (above this sampler)

`nested_cut` fits this sampler once per tree node, with inclusive-value (`IV_*`) columns carrying
the log-sum-exp from child nodes and CUT/multiple-imputation propagating their uncertainty.
`NCUT_TREE=flat` collapses the tree to a single softmax over all classes.

---

## 8. What to actually switch on

**This section is not the authority on defaults — the callers are.** Two production callers ship
DIFFERENT configurations, both newer than the table below was:

- **flat / pixel** — `run_lu_pixel_model.R:2001-2034` (`sampler_extra$mnlogit_rcpp_sym`)
- **nested** — `codes/nested_cut.R:325-425` (`.ncut_fit_block`)

Read those two blocks before quoting anything here. Where this document and the code disagree, the
code is right and this document is stale.

### 8.0 THE STANDARD CONFIGURATION — BMLEH_Los1 (adopted 2026-09-16)

`projects/BMLEH_Los1_CAPRI/gamble_model/estimate_prior.R` is the reference configuration, because it
is the only fit in this repo with **convergence, calibration and held-out skill evidenced in one
place** — and on the hardest design here (66,861 pixels, 43 classes, 84 covariates, 25 CAPRI groups):

| block | evidence |
|---|---|
| `log_lik` | Rhat 1.036, ESS 379 |
| `mu` (3,612 FE) | Rhat med 1.001 / max 1.525, 2% above tolerance, ESS med 1,623 |
| `sigma_re` | Rhat 1.005, ESS 422 |
| held-out | McFadden 0.2454 on 13,372 pixels (20% random split) |
| calibration | mean absolute gap 0.0019 over 12 predicted-share bins |
| country RE | +8,438 nats (pooled McFadden -0.0009 -> 0.2454) |

    symmetric = TRUE, symmetric_hs = TRUE
    use_horseshoe = TRUE, horseshoe_idx = setdiff(seq_len(ncol(X)), icpt)   # intercept NOT shrunk
    fe_support_strength = 0, re_support_strength = 1
    re_regularize = TRUE, re_asis = TRUE
    estimate_slab_c2 = FALSE, collapse_slab_c2 = 14.69                      # design-dependent, see 8.1b
    use_bart = FALSE                                                        # baseline is linear
    re_idx = intercept + log1p_GDP + log1p_Pop + GHM_HI + CISI
    thin = 1, 4 chains, 2000 burn + 500 kept, extended by resumable segments

**Provenance caveat.** The 2026-09-12 segment predates `run_config.rds` (2026-09-14), so its report
states the switches were not recorded. They are the `run.sh` defaults and almost certainly what ran,
but one further `run.sh more` segment would pin this properly. Do that before quoting these numbers
in print.

### 8.1 What the other two callers now do

Aligned to the standard on 2026-09-16 (`symmetric_hs`, `fe_/re_support_strength`, intercept excluded
from the horseshoe, `re_asis`):

| | flat (`run_lu_pixel_model.R`) | nested (`codes/nested_cut.R`) |
|---|---|---|
| `symmetric_hs` | TRUE | TRUE (`NCUT_SYM_HS`, default flipped) |
| `fe_support_strength` | **0**, explicit | 0, forced |
| `re_support_strength` | **1**, explicit | 1 (`NCUT_RE_SUPPORT`) |
| `horseshoe_idx` | excludes the intercept **by name** | excludes `IV_*` and the intercept (`NCUT_HS_NO_INTERCEPT`) |
| `re_asis` | **TRUE** | **TRUE** |
| `const_sum_blocks` | `"auto"` | `"auto"` |
| BART | full validated block, opt-in | **absent — the nested path has no BART at all** |

**b) Deliberate, documented divergences — set explicitly, not by accident:**

| knob | value | why it differs |
|---|---|---|
| `re_asis` | TRUE everywhere | BMLEH: sigma_re Rhat 1.618 -> 1.007, ESS 7 -> 470, chain spread 0.0220 -> 0.0031. **Counter-measurement kept on the record:** an earlier GLOBIOM-node arm found ASIS hurting RE-variance ESS (73 -> 27). Standard is ON; re-measure before assuming it transfers to a small-class design |
| `estimate_slab_c2` | BMLEH FALSE (c2 = 14.69); flat & nested TRUE | BMLEH: the slab never binds (cap sigma <= 3.8 vs sigma_re 0.108) and sampling it adds the worst-mixing scalar for free. Flat/nested: full Bayes beats the best fixed value held-out, and c2 is identified (starts 4 and 100 both converge to 4.12). Genuinely design-dependent |
| `use_horseshoe` | BMLEH & flat TRUE; nested FALSE | Nested measured it POST-kernel-fix on the real root design at `fe_support_strength = 0`: ridge -2306.1 vs ridge+horseshoe -2369.3, i.e. -63.2 nats. That measurement stands under the new standard, so nested keeps the ridge by choice |
| `re_idx` | BMLEH 5, flat 8, nested `NCUT_RE_COLS` | not yet reconciled; intercept-only measured at parity on skill with better convergence, but only on one node |

**BART (flat only), the validated block:** `bart_symmetric = TRUE`, `bart_base = 0.90`,
`bart_power = 3.0`, `bart_k = 2.0`, `bart_k_prevalence = TRUE` (cap 10),
`store_bart_trees = TRUE`, `do_slim_trees = TRUE`. `use_bart` itself is opt-in
(`DRIVER_USE_BART`, default FALSE, 25 trees) — but when it IS on, that is the configuration.
Symmetric and prevalence-scaled `k` are complementary and super-additive: held-out total
+38.4 -> +75.4, rare classes -51.6 -> -9.0. It stays OUT of the standard deliberately: the GLOBIOM
gate found BART worth +249 nats held-out at ~1.7x fit time, concentrated in abundant terrain-driven
classes — not the right first cut for a target with many rare classes. BART is a treatment arm, not
part of the baseline.

**Store-hash warning.** Flipping a nested default changes behaviour while an unset env var still
hashes the same. The store's `prior` signature therefore records **effective** values, not raw env
strings (2026-09-16). Any future default change must keep that property, or resumed runs will
silently reload draws from the old configuration.

### 8.2 Revised switch status

| switch | status |
|---|---|
| `symmetric_hs` | **the "costs ~23/120/129 nats" gate is VOID.** Those arms ran the pre-2026-08-19 bug (`Mb` applied within-equation instead of `kron(Msym, Mb)`; complement redraw conditioned on post-ASIS `mu_R`), one of them at Rhat 1.91. Post-fix: McFadden -0.018 -> +0.216 on the real root node, and the only head-to-head is IN-SAMPLE at one node (+0.216 sym vs +0.260 diag). Flat ships TRUE, nested ships FALSE, and **no held-out comparison has been run on either** |
| `hs_kernel_live` | **the frozen-kernel bug is FIXED (2026-08-30).** `use_horseshoe = TRUE` used to deliver a plain A0=2 ridge silently; it now reaches the draw. `hs_kernel_live = "frozen"` is retained only to reproduce a pre-fix fit. Measured on the real root design: ridge -2306.1 vs ridge+horseshoe -2369.3 (-63.2 nats) at `fe_support_strength = 0`, and -1273.9 at 2. Nested therefore turns the FE horseshoe OFF by choice |
| `fe_support_strength` | **OPEN DISCREPANCY.** Nested forces 0 and documents `(n/PR)^2` on the symmetric kernel as a median 106x / max 2.8e5x precision inflation costing ~125 nats. The flat driver passes `support_prior_strength = 2` and no override, so it inherits 2 — and with `symmetric_hs = TRUE` that factor IS live, via `c_v` (`mnlogit_rcpp_sym.R:3638`, "feeds the symmetric kernel"). The -125 nats was measured on the nested root, not on the flat design; it has never been measured on the flat path. Do not assume either caller is right |
| `use_country_shrinkage` | gates correctly (tau 0.68 rich vs 0.185 sparse) but `sigma_v x tau_g` is NON-IDENTIFIED and absorbs the fixed effects |
| `alt_spec_Z` | delta recovers (1/1.5/2/3 -> 1.006/1.525/2.021/2.976). Nested exposes it as `NCUT_ALT_BLOCKS`; works flat and at leaves, **fails at a root nest** (`.ncut_node_iv` rebuilds child utilities from the full X and without delta). No equivalent env knob on the flat path — there, a spatial/temporal lag is an ordinary design column (`DRIVER_FOCAL_YEARS`, `DRIVER_PREV_STATE`) |
| `collapse_re_var` | rejected for the pixel model (ESS 20x worse) |
| `re_idx = 1:k` | over-parameterised on every design tested (189 rows / 26 countries gives 6-8 RE parameters per observation) |
| `estimate_c2` (FE slab) | TRUE in both callers, at explicit request. Costs ~1.3 held-out LL vs FIXED on Pasture; kept for consistency with the RE slab, **not** because it improves fit. The inline comment at `run_lu_pixel_model.R:2015` still says it "stays FALSE" — that comment is stale, the code below it sets TRUE |

### 8.3 Where the skill actually is

`use_re = TRUE` with an intercept-carrying RE block is the single largest term: the country RE
carries essentially all of the held-out log-score skill, and the intercept carries essentially all
of the RE. Random slopes tie on skill for 8x the per-group parameters and converge worse
(RE Rhat 1.5 vs 1.1). Start there before touching any shrinkage knob.

**Scoring caveat.** `codes/score_nested_cut.R` / `run/score.R` is IN-SAMPLE — it measures
reproduction, not generalisation, and mechanically favours the richer RE block. The only genuine
held-out harness in the repo is `experiments/mixing/re_idx_tradeoff.R` (country-stratified, with
predictions averaged over posterior draws rather than plugged in at the mean).

## 9. Known traps

- **`support_prior_strength` used to drive two unrelated mechanisms** (the RE sparse-group gate and
  the FE horseshoe multiplier). Split into `re_support_strength` / `fe_support_strength`; the old
  name still works as a fallback for both.
- **Inert arguments.** MNL: `a_pi`, `b_pi`. Count: `tau0_dev`, `p0_re`, `bart_batch_size`,
  `symmetric` — all signature-only. `tau0_dev` is passed in production and ignored; the count
  sampler now warns.
- **`tau0` / `p0` are not tuning knobs.** The xi update uses `rate = 1 + 1/tau^2`, i.e.
  Half-Cauchy(0,1); `tau0_pooled` is computed, printed, and never enters the sampler.
- **The pixel suite has no `alt_spec` coverage.** That is why `tests/test_altspec_recovery.R` and
  `tests/test_altspec_multiblock.R` exist as separate files.
- **A disconnected parameter still "covers truth"**, because its posterior IS its wide prior. Gate
  recovery tests on posterior sd vs PRIOR sd, not on coverage.
