# MNL sampler — full model specification

What `mnlogit_rcpp_sym` fits when every switch is enabled, and which switches you should
actually enable. Written 2026-08-26 against 106 formals across eight subsystems.

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

**With all switches on this is not a model you should run.** Several are individually measured as
costly or not yet identified:

| switch | status |
|---|---|
| `symmetric_hs` | works; costs ~23 nats vs diagonal (RE-engaged, one node). Use if baseline invariance is the goal, not for skill |
| `hs_kernel_live` | correct, but -41 nats; **catastrophic (-1274) unless `fe_support_strength = 0`** |
| `use_country_shrinkage` | gates correctly (tau 0.68 rich vs 0.185 sparse) but `sigma_v x tau_g` is NON-IDENTIFIED and absorbs the fixed effects |
| `alt_spec_Z` | delta now recovers (1/1.5/2/3 -> 1.006/1.525/2.021/2.976); still behind `allow_unvalidated`, never run on a real design |
| `collapse_re_var` | rejected for the pixel model (ESS 20x worse) |
| `re_idx = 1:k` | over-parameterised on every design tested (189 rows / 26 countries gives 6-8 RE parameters per observation) |

**A coherent "everything sensible on" configuration:**

    use_horseshoe    = TRUE
    hs_groups        = <per-family / per-block tau>
    symmetric_hs     = only if baseline invariance is the goal
    fe_support_strength = 0        # never pair a live FE horseshoe with this
    re_support_strength = 1        # the sparse-group gate; keep it
    use_re           = TRUE, re_idx = intercept (+ a few measured slopes)
    re_regularize    = TRUE, estimate_slab_c2 = TRUE
    use_spike_slab   = TRUE
    alt_spec_Z       = <spatial and/or temporal lag blocks>
    use_bart         = FALSE
    use_car          = FALSE

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
