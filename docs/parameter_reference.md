# Parameter reference — priors, full conditionals, draw algorithms

Every sampled quantity in `mnlogit_rcpp_sym` / `mncount_rcpp`: what it is, its prior, its full
conditional, how it is drawn, and where. Companion to `sampler_model_specification.md` (the model as
a whole) and `nested_cut_model.md` (the tree above it).

Notation: `i` pixel, `j`/`k` alternative (class), `ip` equation (baseline-removed, `j = pp[ip]`),
`v` covariate, `g` group (country), `b` block. `K = p_all` alternatives, `p = K-1` equations,
`k` covariates, `G` groups.

---

## 0. Likelihood and augmentation

    Y_i.  ~  Multinomial( N_i , softmax_k(eta_ik) )

**One-vs-rest, not a joint softmax.** Each equation is reduced to a binary logit against the rest:

    psi_ij = eta_ij ,   c_ij = log sum_{k != j} exp(psi_ik)
    logit  = psi_ij - c_ij

**Polya-Gamma augmentation.** With `omega ~ PG(N, psi - c)` the binary logit becomes Gaussian:

    omega_i,ip     ~  PG( N_i , psi_ij - c_ij )
    y~_i,ip        =  kappa_i,j / omega_i,ip + c_ij          (working response)
    y~_i,ip | .    ~  N( eta_ij , 1 / omega_i,ip )

This is why nearly every conditional below is conjugate: conditional on `omega`, each equation is a
weighted Gaussian regression. `kappa = Y - N/2` is the PG-centred count.

> **`c_ij` depends on every parameter that enters any `psi_ik`.** A conditional derived while holding
> a STALE `c` is wrong. This is not pedantry: the alt-spec `delta` read exactly 2x low for months
> because its channel was gated off, so `c` was computed without `delta` while `delta` was estimated
> against it.

---

## 1. `beta` / `mu` — coefficients            [C++ `gibbs_step_*`]

    eta_ij = x_i' beta_.j                 (pooled)   or   x_i' (mu_.j + b_.j,g(i))   (with REs)

**Prior** `beta_v,j ~ N(m_vj, A0)`, `A0 = 2` slopes / `100` intercept, plus the horseshoe precision
(§3) when active.

**Full conditional** — Gaussian, per equation:

    P  = prior_P + hs_prec + X' diag(omega_.ip) X
    Pb = prior_Pb + X' ( kappa_.j + omega_.ip * c_.ip - omega_.ip * f_.ip )
    beta_.ip | .  ~  N( P^-1 Pb , P^-1 )

`f` is the additive channel (BART + alt-spec). **Draw**: Cholesky of the precision
(`chol_sample_precision_cpp`), never an explicit inverse. Under `symmetric_hs` the K equations are
drawn JOINTLY per covariate with precision `c_v * M_sym`, `M_sym = I - 11'/K`.

---

## 2. `b_g` — random effects              [C++ `gibbs_step_re_ncp`]

    b_v,j,g  ~  N( 0 , (sigma_v * r_v,g * tau_g)^2 )

Non-centred: `b = sigma * r * tau * z`, `z ~ N(0,I)`, drawn as a Gaussian conditional exactly as §1
with `mu` as the offset. ASIS interweaves a centred redraw of `mu` each sweep.

- `r_v,g` — participation-ratio support factor, **fixed** (computed once from X), hence identified
- `tau_g` — optional learned per-group factor (count model only; see §8)

---

## 3. `lambda_v`, `tau_f`, `c2` — the regularised (Finnish) horseshoe

**Prior**, as a harmonic mixture of precisions — this identity is why it stays conjugate:

    1/sigma~_v^2  =  1/(tau_f^2 lambda_v^2)  +  1/c^2
    lambda_v ~ C+(0,1),   tau_f ~ C+(0,tau0),   c^2 ~ InvGamma(nu/2, nu s^2/2)

**Makalic-Schmidt auxiliaries** turn every half-Cauchy into two inverse-gamma draws:

    nu_v  | lambda_v^2  ~ IG( 1 , 1 + 1/lambda_v^2 )
    lambda_v^2 | .      ~ IG( (m+1)/2 , 1/nu_v + ||gamma_v||^2 / (2 tau_f^2) )
    xi_f  | tau_f^2     ~ IG( 1 , 1 + 1/tau_f^2 )
    tau_f^2 | .         ~ IG( (D_f m + 1)/2 , 1/xi_f + sum_v ||gamma_v||^2/lambda_v^2 / 2 )
    zeta  | c^2         ~ IG( 1 , 1 + 1/c^2 )
    c^2   | .           ~ IG( (nu + D m)/2 , nu s^2/2 + sum ||gamma||^2 / 2 )

`m = p` (or the block rank `Kb-1`), `D_f` = covariates in family `f`. All exact IG draws — no MH, no
slice. `hs_groups` gives each covariate family and each const-sum block its own `(tau_f, xi_f)`.

> **`tau0` is NOT a tuning knob.** The `xi` update uses `rate = 1 + 1/tau^2`, i.e. Half-Cauchy(0,1);
> `tau0_pooled` is computed, printed, and never enters the sampler.

> **Under `symmetric_hs` the statistic is `||gamma_v||^2 = b' M_sym b`** (the zero-sum magnitude),
> computed on CLR-centred coefficients, so it is rotation- and baseline-invariant.

---

## 4. `sigma_v` — RE scale        [C++ `update_re_precision_hc(_sym)`]

    sigma_v ~ C+(0, A),  optionally capped:  tau_eff = tau_raw + 1/c2_re

**Unregularised** (`re_regularize = FALSE`) — exact Gamma:

    a_v   | .  ~ IG( 1 , tau_v + 1/A^2 )
    tau_v | .  ~ Gamma( 1/2 + n_v/2 , a_v + ss_v/2 )        ss_v = sum_{g,ip} (b - mu)^2 / r^2

**Regularised** — `tau_eff` is not Gamma, so a **1-D slice on `log tau_raw`**:

    lp(l) = (n/2) log(tau_raw + 1/c2) - (tau_raw + 1/c2) ss/2 - l/2 - a tau_raw + l

**`ss` must be in the SAME coordinates as the draw.** It is whitened by `1/r^2` because the draw
scales by `r`. The symmetric variant additionally CENTRED `ss` across categories while the draw
stayed uncentred — that zeroed the random effects outright (RE sd 0.0000, -51.6 nats). Hence
`re_prec_center` now resolves to FALSE automatically whenever the symmetric updater is used.

---

## 5. `delta_b` — alternative-specific coefficients

    eta_ij  +=  sum_b delta_b o z^b_ij

**Prior** `delta ~ N(0, s^2 I)`. **Full conditional — exact multivariate Gaussian**, because under PG
each equation is a weighted regression and `delta` is an ordinary coefficient on `z`:

    P = Lambda0 + sum_ip D_ip' W_ip D_ip ,  W_ip = diag(omega_.ip)
    b = sum_ip D_ip' W_ip r_ip ,            r_ip = y~_ip - X beta_.ip - f^BART_ip
    delta | .  ~  N( P^-1 b , P^-1 )

All blocks drawn JOINTLY (a spatial and a temporal lag are correlated; a coordinate scan crawls).
The per-equation design `D_ip` depends on the coefficient type:

| `coef` | constraint | `D_ip` column k | reads as |
|---|---|---|---|
| `shared` | one delta | `z_ij - z_ib` | one stickiness for all classes |
| `per_class` | `delta_b = 0` | `z_ij` on `k = ip` only | class j vs a reference class |
| `symmetric` | `sum_j delta_j = 0` | `z_ib + 1{k=ip} z_ij` | class j vs the AVERAGE (baseline-invariant) |

Derivation for `symmetric`: with `delta_b = -sum_{k!=b} delta_k`,
`eta_ij - eta_ib = sum_k delta_k [ 1{k=j} z_ij + z_ib ]`, hence the full p-column design.

**`scale = "sd"` (default).** `delta` multiplies a SHARE whose spread varies 10-30x across classes,
so raw deltas are not comparable. Four classes simulated with an IDENTICAL standardised effect of
0.800 return raw deltas 6.8 / 6.7 / 31.7 / 31.0 — a pure units artifact. Scaling is applied in
`nested_cut` (not the sampler) so the FIT-TIME factor can be REPLAYED when the block is rebuilt for
inclusive values or prediction; recomputing it on another sample would silently change what delta
means.

**Interactions.** A block may be `w_m * z_ij`, giving `delta_i = delta_0 + sum_m gamma_m w_im`:
`gamma_m < 0` = the covariate makes a transition MORE likely; `> 0` = it locks the pixel in.

---

## 6. `gamma` (inclusion), `pi` — spike-and-slab on RE cells

    gamma_v,ip,g ~ Bernoulli(pi_v,ip) ,  pi ~ Beta(a_pi, b_pi)

Conditional inclusion odds = prior odds x the marginal likelihood ratio of the cell with and without
its RE. Conjugate Bernoulli/Beta draws.

## 7. `tau_spatial`, CAR

    b_.,g ~ CAR(W, rho):   precision  tau_sp * (D - rho W)      over GROUPS
    tau_sp | . ~ Gamma( a + G/2 , b + b' (D - rho W) b / 2 )

`rho` is fixed (`car_rho`), clamped below the stability limit.

## 8. `tau_g` — per-group RE shrinkage (count model, EXPERIMENTAL)

    tau_g ~ regularised-C+(0, tau0_country, slab_c2_country)

Multiplied into `r_v,g`, so one scale gates ALL of a group's REs. **Gates correctly** (0.68 for
data-rich groups vs 0.185 for sparse ones) but **`sigma_v * tau_g` is NOT identified** — only the
product is, so `tau` shrinks while `sigma` inflates until the RE space absorbs the fixed effects
(b[x1] 0.731 -> 0.256 against truth 0.700). Anchoring must be imposed INSIDE the draw; a post-hoc
rescale breaks the prior's own cap. Default OFF.

## 3b. FE horseshoe: the lambda/tau scaling ridge (real, but largely cosmetic)

`post_kappa_pooled` stores the SHRINKAGE FACTOR `kappa = 1/(1 + tau^2 lambda^2)`, not lambda.
(The convergence harness labels this block "lambda (HS)"; it is kappa.)

**The scaling ridge is real.** `(tau, lambda) -> (c*tau, lambda/c)` leaves the model invariant, a
continuous 1-D orbit. Chains do sit at different points on it — 4 chains, 5k px:

    final tau2      0.0673  0.0981  0.0927  0.0350    ratio 2.80
    median lambda2  0.4726  0.3386  0.3929  1.1180    ratio 3.30
    product         0.0318  0.0332  0.0364  0.0391    ratio 1.23   <- agrees

tau2 and lambda2 move in OPPOSITE directions; the product is far more stable. Textbook.

**But kappa's Rhat ~1.75 / ESS 6 / B/W 15 overstates the problem.** kappa is invariant to the
scaling, and chains agree on it to three decimals (mean 0.826 / 0.824 / 0.824 / 0.824; median
across-chain range 0.023 on a [0,1] scale). B/W is large because kappa barely moves WITHIN a
chain — a small denominator, not a large numerator. Rhat is scale-free and so flags a
near-constant quantity on tiny absolute differences.

Not entirely cosmetic: a 0.023 range near kappa = 0.82 moves `1-kappa` from 0.174 to 0.151,
about 13% relative variation in retained signal. Small, real, probably not worth chasing.

**Do NOT impose an ordering constraint (lambda < tau) to fix this.** Two reasons: (1) an
inequality is a positive-measure truncation of a continuous orbit — it fences the ridge without
pinning it, and the chain still slides inside the fence; continuous degeneracies need a
MEASURE-ZERO normalisation (fix geometric-mean lambda = 1, or pin tau) or a move ALONG the orbit
(a joint MH proposal `(tau,lambda) -> (c*tau, lambda/c)`, likelihood-invariant so the ratio is
prior + Jacobian only — same pattern as the `re_mean_shift` interweave). (2) Large `lambda_j` IS
the horseshoe's mechanism for letting signals escape shrinkage; capping it below the global scale
converts the horseshoe into an expensive ridge.

**Consequence for reparameterisation work:** NCP / Makalic-Schmidt would help traverse this orbit
but cannot identify lambda separately from tau, and the identified combination already agrees
across chains. Expect cleaner tau/lambda diagnostics, not better inference.

**Does the horseshoe earn its keep? UNRESOLVED.** Held-out LL, chains pooled, kernel live,
`fe_support_strength = 0`:

| split | HS | ridge (A0=2) | diff | note |
|---|---|---|---|---|
| 1 | -3724.1 | -3696.8 | -27.3 | ridge better |
| 2 | -3718.8 | **-4004.1** | +285.3 | ridge BLEW UP |
| 3 | -3712.5 | -3690.7 | -21.8 | both arms ESS 9 — unusable |

mean +78.7, sd 178.9 -> not distinguishable. Typical case mildly favours the RIDGE (~20-27 nats);
the mean is carried entirely by one ridge failure. HS held-out spans 12 nats across splits, the
ridge spans 313 — consistent with the horseshoe buying VARIANCE REDUCTION rather than mean
accuracy, but the split-2 blow-up has NOT been traced and one diverged chain among the four
pooled would explain it with no robustness story.

## 7b. RE slab `c2` — weakly identified, but KEEP IT ON (measured)

`re_regularize` caps the RE variance via `tau_eff = tau_raw + 1/c2`; `estimate_slab_c2` samples
`c2` by a 1-D slice (InvGamma, `slab_df_re = 10`). Its diagnostics look pathological and the
pathology is real, but removing it is WORSE.

**c2 is weakly identified on the pixel design.** The slab only bites when RE variances press
against the cap, and they do not: of 1560 RE cells only 130 are active (92% pinned/masked), and
their SDs are tiny (median 0.038), so `tau_eff` is ~680 while `1/c2` is ~0.63 — the slab supplies
**0.1% of the precision for the median active cell**, >10% for only 20% of them. So c2 is informed
by a handful of cells and otherwise follows its prior: 4 chains land at 1.03 / 2.90 / 0.97 / 1.14,
Rhat 1.34, ESS 10, between/within variance ratio **107** (stuck, not merely autocorrelated).

**It is a ridge, not slow exploration.** `slab_c2_n_slice = 10` (ten extra slice sweeps per Gibbs
iteration, same conditional) changes nothing: ESS 10 -> 9, Rhat 1.339 -> 1.427. Same aliasing
family as `tau_g` and `kappa_v` — c2 enters ONLY through `tau_eff`.

**But it earns its keep.** Held-out log-lik, chains pooled, 3 splits (8k px, 4 chains x 1200):

| config | held-out LL | log_lik ESS |
|---|---|---|
| regularize + **estimate c2** | **best** (+7.3 nats vs no-slab, 3/3 splits) | 12 |
| regularize + fixed c2 | worst | 20 |
| no slab | middle | 35 |

Removing the slab improves joint mixing 2.9x and costs ~7 nats (0.2%) consistently. A weakly
identified c2 still averages over regularization strengths, which beats any single fixed value —
it cannot be learned, but it is doing work. **Do not "fix" this by turning it off.** Treat c2's
Rhat as a nuisance-hyperparameter diagnostic, not a convergence failure of the estimands, and buy
back the mixing with iterations.

## 8a. `re_mean_shift` — sum-to-zero RE identification (default OFF, costs sigma_re ESS)

MH interweave moving `mean_g(b_{v,j,g})` into `mu_{v,j}`. The likelihood is exactly invariant, so
only the prior ratio gates the move; the RE-variance draw then charges `G-1` free dimensions
(`dof_mean_pinned`, which exists ONLY in the `_sym` updater — see the note below). `mu` becomes
the population-averaged effect by construction, and it recovers that identity on synthetic data
(leak 0.0108 -> 0.0001).

**It reliably costs `sigma_re` ESS on the real design, and buys nothing at scale.** ESS_bulk
median, GLOBIOM pixel, baseline -> +shift:

| run | `log_lik` | `sigma_re` | `slab_c2` |
|---|---|---|---|
| 6k px, 4ch x 1500, diagonal | 46 -> 38 | 142 -> **14** | 79 -> **14** |
| 6k px, 4ch x 1500, symmetric | 292 -> 491 | 85 -> 92 | 8 -> 55 |
| **15k px, 5ch x 3000, symmetric** | 749 -> 659 | 275 -> **55** | 11 -> 16 |

The middle row once looked like an endorsement (`log_lik` +68%, `slab_c2` +618%). The larger run
does not reproduce it: `sigma_re` drops 80% and `log_lik` is slightly worse. Those 6k numbers were
small-sample noise in an under-converged regime — `slab_c2` ESS of 8 vs 55 is not a comparison.
The signal that holds across ALL THREE runs is the `sigma_re` penalty. Keep it OFF.

**`slab_c2` is a standing defect of the symmetric path, not something the shift causes or cures.**
At 15k it is Rhat 1.395 / ESS 11 in the BASELINE arm. The full-Bayes slab barely moves between
chains under the symmetric updater; that is an open problem.

Note `symmetric_hs` is a SEPARATE argument from `symmetric`: the latter sets the zero-sum coding
of the response, the former selects the symmetric shrinkage prior AND the `_sym` RE-variance
updater (the only one accepting `dof_mean_pinned`). Passing only `symmetric = TRUE` silently
measures the diagonal horseshoe with the dof half inert. `re_mean_shift_dof` keeps the two
separable so they can never again be confounded.

## 8b. `kappa_v` — joint FE/RE gate (MNL, EXPERIMENTAL, default OFF)

    mu_{v,j}  ~ N(0, kappa_v^2 * sigma~_{v,j}^2)      (FE, horseshoe kernel)
    b_{v,j,g} ~ N(0, kappa_v^2 * sigma_{v,j}^2)       (RE, per-covariate scale)
    kappa_v^2 ~ C+(0, 1)   [Makalic-Schmidt: kappa^2 | nu ~ IG(1/2, 1/nu + ss/2)]

Flag `joint_fe_re_shrink`. One scale per covariate on BOTH channels, so a covariate cannot
survive by hiding in the block the shrinkage does not see. `kappa` reaches the RE draw and the
RE variance update through the SAME matrix (`re_supp_kap`), which is what stops `sigma` from
absorbing `1/kappa` (the `tau_g` failure in section 8).

**Back-compat when OFF is exact** — verified bit-identical on all 10 numeric outputs against a
sampler copy with the blocks physically removed. Every path is inside `if (joint_shrink_on)`,
so no extra RNG draws are consumed.

**But `kappa` is weakly identified, for the same reason as `tau_g`.** It is aliased with
`sigma_v` in the RE block, which supplies `p*G` of the `p*G + p` effective observations in its
conditional; those deviations are whitened by `sigma`, so `ss/n_eff ~ 1` and `kappa` sits near 1
(posterior mean 1.0, sd 0.16, min over draws 0.61 — never near the ~1e-2 a real gate would reach).
Measured against a 3-seed noise floor it moves a pure-null covariate's RE sd (0.048 vs OFF
0.121/0.160/0.195) but leaves a signal-carrying one inside the seed spread. Treat it as a mild
extra prior, not a variable-selection gate, and never report it without the noise floor.

Identified alternative, if wanted: let `kappa` REPLACE `sigma_v` instead of multiplying it — the
`mu`/`b` ratio stays expressible through `tau*lambda`, at the cost of collapsing the per-class RE
variance to per-covariate.

## 9. `r` — negative-binomial dispersion (count model)

    y ~ NB(mu, r) ,  r ~ Gamma(a, b)

Drawn by **CRT-Gibbs**: augment `L_i ~ CRT(y_i, r)` (a sum of Bernoullis), then
`r | . ~ Gamma(a + sum L_i, b - sum log(1 - p_i))`. Exact, and replaces a slow log-RW Metropolis.

## 10. BART `f`

Standard BART back-fitting on the PG working residual `y~ - X beta`, with weights `omega`. Tempered
by `bart_alpha` during burn-in. Enters the same additive channel as `delta`, and each nets the other
out of its own residual.

---

## Sweep order

    omega -> beta/mu (C++) -> ASIS -> sigma_v -> spike-slab -> CAR -> horseshoe -> delta -> BART -> store

## Which draws are exact

| exact conjugate | Gaussian: `beta`, `mu`, `b_g`, `delta`;  IG/Gamma: all horseshoe scales, `sigma_v` (unregularised), `pi`, `tau_sp`;  CRT-Gibbs: `r` |
|---|---|
| **1-D slice** | `sigma_v` regularised, `c2_re`, `tau_g` |
| **Metropolis** | none in the MNL path (the `delta` MH was replaced by its exact conditional) |
