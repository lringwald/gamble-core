# Cut / IV-imputation nested multinomial logit — model & computation

**Code:** `codes/nested_cut.R` · **Test:** `test_nested_cut.R` · **Reuses:** `codes/mnlogit_rcpp_sym.R` (unchanged)

*Living document — kept in step with the implementation. Equation-by-equation model spec (§1–8) then a technical appendix on the algorithm (A–C).*

---

## 1. Data and notation

- Pixels $i = 1,\dots,n$. Covariates $x_i \in \mathbb{R}^{p}$ (row of design $X$; **must include an `intercept` column**).
- Fine land-use classes $j = 1,\dots,J$. Observed **compositional response** $y_i \in \Delta^{J-1}$ (shares, $\sum_j y_{ij}=1$); the sampler treats $y_i$ as expected multinomial counts.
- Random-effect groups (regions/countries) $g(i)\in\{1,\dots,G\}$.
- A **nesting tree** $\mathcal{T}$ partitions the fine classes hierarchically (§2).

## 2. The nesting tree

$\mathcal{T}$ is a rooted tree whose **leaves are fine classes** and whose **internal nodes are nests**. A node is either

- a **terminal nest** — a set $c$ of fine classes $\{j\in c\}$ (a singleton if $|c|=1$), or
- an **internal nest** — with ordered children $\mathrm{ch}(c)=\{c_1,\dots,c_{K_c}\}$, each itself a node.

Let $\mathcal{F}(c)$ = the fine classes under node $c$. The path from the root to a leaf $j$ is $c^{(0)}\!=\!\text{root}\supset c^{(1)}\supset\cdots\supset c^{(L_j)}=\{j\}$.

Example (AgMIP): `root = {Cropland={Arable={Wheat,…}, Permanent={Grapes,…}}, Forest={…}, Grassland, Built_up, …}`.

## 3. Choice probabilities and inclusive values

The joint fine-class probability factorises along the path (sequential/recursive nested logit):

$$
P(j \mid x_i) \;=\; \prod_{\ell=1}^{L_j} P\!\left(c^{(\ell)} \,\middle|\, c^{(\ell-1)}, x_i\right). \tag{3.1}
$$

At an internal node $c$ with children $\mathrm{ch}(c)$, the child choice is a softmax over child utilities:

$$
P(c_k \mid c, x_i) \;=\; \frac{\exp\!\big(v_{i,c_k}\big)}{\sum_{k'=1}^{K_c}\exp\!\big(v_{i,c_{k'}}\big)},
\qquad
v_{i,c_k} \;=\; x_i^{\top}\gamma_{c_k} \;+\; \lambda_{c_k}\, \mathrm{IV}_{i,c_k}. \tag{3.2}
$$

The **inclusive value** (McFadden logsum) of a child $c_k$ is the expected maximum utility of *its own* choice set, i.e. the log-partition of the subtree rooted at $c_k$:

$$
\mathrm{IV}_{i,c_k} \;=\;
\begin{cases}
\displaystyle \log\!\sum_{j\in c_k} \exp\!\big(x_i^{\top}\beta_{j}\big), & c_k \text{ a terminal (multi-class) nest},\\[2mm]
\displaystyle \log\!\sum_{c'\in\mathrm{ch}(c_k)} \exp\!\big(v_{i,c'}\big), & c_k \text{ an internal nest},\\[1mm]
0, & c_k \text{ a singleton / degenerate (no signal)}.
\end{cases} \tag{3.3}
$$

Here $\beta_j$ are the terminal sub-utilities and $v_{i,c'}$ are the child utilities (3.2) one level down — so (3.3) is **recursive**. The **dissimilarity parameter** $\lambda_{c_k}\in(0,1]$ measures within-nest correlation: $\lambda=1$ collapses the nest to flat MNL; $\lambda\to 0$ makes the nest act as a single alternative. Dropping the IV term ($\lambda_{c_k}\equiv 0$) gives the **factorized** variant (§8).

**Behavioural reading.** A pixel whose best sub-alternative (e.g. wheat) scores high gets a high $\mathrm{IV}$, which raises $P(\text{parent nest})$ (e.g. Cropland) — the "attractiveness of the best crop feeds the farm-or-not decision" coupling that flat factorization severs.

## 4. Symmetric zero-sum parameterisation and the λ rescaling

Each node's choice model (3.2) is fit by the **symmetric (zero-sum) multinomial sampler** `mnlogit_rcpp_sym`, which parameterises coefficients on the $K$-simplex with the identifiability constraint

$$
\sum_{k=1}^{K}\gamma_{c_k}=0,\qquad \sum_{k=1}^{K}(\text{coef on any covariate})=0. \tag{4.1}
$$

Consequence for **λ**: a *true* dissimilarity $\lambda_{c_k}$ enters only child $c_k$'s utility, but under the zero-sum constraint the estimated coefficient of the shared column $\mathrm{IV}_{c_k}$ is spread across categories and appears on $c_k$ as $\lambda_{c_k}\,(K-1)/K$. We therefore recover

$$
\boxed{\;\hat\lambda_{c_k} \;=\; \big(\text{zero-sum coef of }\mathrm{IV}_{c_k}\text{ on child }c_k\big)\times \frac{K}{K-1}\;}\qquad (K=K_c). \tag{4.2}
$$

*(Omitting the $K/(K-1)$ factor biases $\hat\lambda$ low — e.g. reads 0.29 instead of 0.58 for $K=2$.)*

## 5. Random effects and priors (per node)

Every node fit reuses the production MNL machinery: random effects on a chosen covariate subset $\mathcal{R}$ (country intercept + socio-economic + terrain by default), a regularised half-Cauchy/horseshoe RE-variance, and a support-aware prior:

$$
\gamma_{c_k} = \bar\gamma_{c_k} + u^{(g)}_{c_k},\quad
u^{(g)}_{c_k}\sim \mathcal N\!\big(0,\ \Sigma^{\mathrm{RE}}_{\text{reg}}\big),\quad
\Sigma^{\mathrm{RE}}_{\text{reg}} = \text{Finnish-HS cap } \big(\sigma^2 \le c_2\big). \tag{5.1}
$$

RE act **within each level independently** — the country-baseline RE that carries most skill in the flat model is present at every node. (Full cross-level RE coupling would require the joint sampler, out of scope here.)

## 6. Estimation — the cut / multiple-imputation posterior

### 6.1 The cut (modular Bayes)

Let $\theta_c=(\gamma_{\cdot},\lambda_{\cdot})$ be a node's parameters and $\beta_{\downarrow c}$ the parameters of its subtree. Full-information Bayes targets the joint $p(\theta_c,\beta_{\downarrow c}\mid \text{all data})$, in which $\mathrm{IV}$ makes $\beta_{\downarrow c}$ enter $\theta_c$'s likelihood **non-linearly** (breaks Pólya-Gamma conjugacy). We instead target the **cut distribution** (Plummer 2015; Jacob et al. 2017):

$$
\pi_{\mathrm{cut}}(\beta_{\downarrow c},\theta_c)
\;=\;
\underbrace{p\!\big(\beta_{\downarrow c}\mid \mathcal D_{\downarrow c}\big)}_{\text{subtree fit, leaf data only}}\;\cdot\;
\underbrace{p\!\big(\theta_c \mid \mathcal D_{c},\, \mathrm{IV}(\beta_{\downarrow c})\big)}_{\text{parent fit given imputed IV}}. \tag{6.1}
$$

Information flows **up only** — the parent never re-informs the subtree ("feedback is cut"). This is the modular, limited-information posterior; it is exact for the sub-models and consistent for $\theta_c$, trading full efficiency for tractability and robustness to sub-model misspecification.

### 6.2 The imputation

The intractable piece is $\mathrm{IV}(\beta_{\downarrow c})$, a function of the subtree posterior. We Monte-Carlo it with $M$ **imputations** drawn from the subtree posterior:

$$
\beta_{\downarrow c}^{(m)} \sim p(\beta_{\downarrow c}\mid \mathcal D_{\downarrow c}),\qquad
\mathrm{IV}^{(m)}_{i,c_k} = \text{(3.3) evaluated at } \beta^{(m)},\qquad m=1,\dots,M. \tag{6.2}
$$

For each $m$ we append the imputed columns $\{\mathrm{IV}^{(m)}_{c_k}\}$ to the node design, $X^{(m)}_c=[\,X \,\|\, \mathrm{IV}^{(m)}\,]$, and draw a parent posterior conditional on that imputation:

$$
\big\{\theta_c^{(m,r)}\big\}_{r=1}^{R} \sim p\!\big(\theta_c \mid \mathcal D_c,\ X^{(m)}_c\big),\qquad R=\texttt{draws\_per\_impute}. \tag{6.3}
$$

### 6.3 Pooling (Rubin's rule)

The cut posterior for $\theta_c$ is the mixture over imputations — approximated by **pooling** all draws:

$$
\pi_{\mathrm{cut}}(\theta_c) \;\approx\; \frac{1}{M}\sum_{m=1}^{M} p\!\big(\theta_c\mid \mathcal D_c, X^{(m)}_c\big)
\;\;\Longleftrightarrow\;\;
\big\{\theta_c^{(m,r)}\big\}_{m\le M,\,r\le R}\ \text{pooled}. \tag{6.4}
$$

The pooled spread automatically contains **between-imputation variance** (leaf uncertainty, Rubin's $B$) + **within variance** (parent sampling, $\bar W$): $\widehat{\mathrm{Var}} = \bar W + (1+1/M)\,B$. $\lambda_{c_k}$'s credible interval (via 4.2 on the pooled draws) is therefore **honestly wider** than the mean-plug-in point.

### 6.4 Why $M$ separate fits, not one chain that randomises the IV

A tempting shortcut is a *single* parent chain that, each iteration, resamples one stored draw $IV^{(m)}$ and updates $\theta$. This is the **naive cut sampler** and it is **biased**: resampling $m$ each sweep gives the mixture *kernel* $K=\frac1M\sum_m P_m$ (with $P_m$ the parent Gibbs kernel of stationary $\pi_m=p(\theta\mid IV^{(m)})$), whose stationary distribution is **not** $\frac1M\sum_m\pi_m$ unless the $\pi_m$ coincide (they don't). The $\theta$-chain *lags* the freshly-swapped $IV$ — it never equilibrates to any one $\pi_m$ before the design changes — and the resulting error typically **under-propagates** the leaf uncertainty (CIs too narrow), defeating the purpose. It would be exact only if $\theta\mid IV$ could be drawn exactly in one shot, but that requires the PG latent $\omega$ and hence more than one sweep (Plummer 2015; Jacob et al. 2017).

### 6.5 The moment carrier — `iv_mode = "moments"` (light path)

Because $IV$ enters the parent utility **linearly** ($v=x'\gamma+\lambda\,IV$), the parent needs only the first two moments of the inclusive value, not all $M$ draws. The $M$ coherent IV fields are also **low-rank** — each is one smooth function of a leaf draw, so their variation lives in a few principal directions. We exploit both with an **unscented (sigma-point)** representation: take the field mean $\mu$ and, for the top $Q$ principal 1-SD perturbation fields $d_q$, the sigma points $\mu,\ \mu\pm\sqrt3\,d_q$. Fit the parent at these **$1+2Q$** fields (hot-started) and pool with weights $\{1-Q/3,\ \tfrac16,\dots\}$ that reproduce the mean and leading-direction variances exactly (3-pt Gauss–Hermite per direction). This propagates the leaf/IV uncertainty **correct to second order** at **$1+2Q$ fits instead of $M$** (e.g. 3 vs 30). It is an approximation (Gaussian-in-the-leading-directions), tight precisely because $IV$ is linear; validated against the exact `"draws"` cut ($\lambda$ and CI agree, ~4–8× faster). $Q$ = `moment_rank`.

**The Gaussian assumption is *checked*, not assumed — `iv_mode = "auto"`.** Quasi-separation in sparse leaves (e.g. rare crop × small-country cells in the 14-class Arable nest) breeds **skewed** IV posteriors, which the sigma-point (mean/variance) representation cannot capture. Two structural facts soften this: the log-sum-exp **down-weights exactly those ill-behaved cells** ($e^{x'\beta}\!\approx\!0$ for a rare class), so $IV$ is dominated by the well-identified common classes; and the low-rank carrier keeps the coherent (Gaussian) modes and drops the local skewed ones. But it is not safe to *assume* — so `auto` **diagnoses each node for free** (the $M$ IV fields are already computed): it tests the skewness / excess-kurtosis of the leading principal-direction scores ($\sqrt{\lambda_q}\,v_q$) against tolerances (`moment_skew_tol`, `moment_kurt_tol`; defaults 1.0 / 2.0 ≈ 2 sampling-SD at $M\!\approx\!30$). Pass ⇒ `moments`; fail ⇒ **automatic fallback to the exact `draws` cut for that node**, logged. The sparse Arable parent thus routes to `draws` while the well-behaved macro nests use `moments`. Skew (the separation signature) is caught reliably; symmetric heavy tails are harder to detect at small $M$ but also less damaging to a mean/variance match. `nested_cut_modes(fit)` reports the mode used per node.

We otherwise draw each imputation's parent posterior with its **own** conditionally-correct chain (nested MCMC), which is exact for the cut. To recover the efficiency the single-chain idea was after, we **hot-start** imputation $m$ from imputation $m{-}1$'s `final_state`: since $IV^{(m)}\!\approx\!IV^{(m-1)}$ (same $X$, neighbouring leaf draws), the chain is already near $\pi_m$, so $m>1$ pay only a short warm burn-in `nburn_warm` (default $\approx\text{nburn}/4$) while $m{=}1$ pays the full cold burn-in. Each fit still conditions on its **own fixed** $IV^{(m)}$ and re-equilibrates → the pooled draws stay the correct cut posterior; the induced serial dependence across $m$ costs only a little effective sample size. **Regime note:** the win is real when the fit is *sweep-dominated* (large-$n$ multi-class leaves); small macro fits are *setup-dominated* (per-fit overhead $\gg$ burn-in sweeps), where the complementary lever is to run the embarrassingly-parallel imputations across cores instead.

## 7. Prediction — coherent cut draws

A predictive draw $d$ **threads one imputation index $m(d)$ through the entire tree** and picks one parent draw per node, so leaf and macro pieces are mutually consistent:

$$
\widehat P^{(d)}(j\mid x) \;=\; \prod_{\ell=1}^{L_j} \mathrm{softmax}_k\!\Big(x^{\top}\gamma^{(d)}_{c^{(\ell)}} + \lambda^{(d)}_{c^{(\ell)}}\,\mathrm{IV}^{(m(d))}_{c^{(\ell)}}(x)\Big), \tag{7.1}
$$

with $\mathrm{IV}^{(m(d))}$ recomputed on the prediction grid from the **same** subtree draw $m(d)$ (recursion 3.3). The set $\{\widehat P^{(d)}\}_{d=1}^{D}$ is the cut predictive: its mean is the point prediction, its quantiles give per-pixel / per-class-total credible intervals.

## 8. Variants (one framework, `(tree, use_iv)`)

| Variant | Call | Structure |
|---|---|---|
| **Flat MNL** | `tree = list(all = <all J classes>)`, `use_iv=FALSE` | single softmax over all fine classes (= current pixel model) |
| **Factorized nested** | hierarchical `tree`, `use_iv=FALSE` | levels independent; $\lambda\equiv 0$, no IV columns |
| **IV-nested (cut)** | hierarchical `tree`, `use_iv=TRUE` | inclusive-value coupling + $\lambda$ with CI (this document) |

All three run through `nested_cut_fit()` and are directly comparable on held-out joint log-likelihood.

## 9. Parameterisation and the fine-class effective symmetric table

**Native parameters are symmetric per nest.** Every node is fit by `mnlogit_rcpp_sym` in the zero-sum parameterisation (4.1), so each node returns a coefficient matrix whose columns sum to zero over *that node's* choice set — including the inclusive-value columns (whence $\lambda$ via 4.2). A fine class $j$ therefore owns a **path** of symmetric coefficients (e.g. Wheat = root's `Cropland` column ⊕ Cropland's `Arable` column ⊕ Arable's `Wheat` column); singletons own the root column directly. `summary_nested_cut()` returns one zero-sum matrix per node + the $\lambda$ table.

**Effective fine-class table (derived).** To compare with the *flat* symmetric heatplot, collapse the per-nest matrices into a single $[J\times P]$ table of **average marginal effects**

$$
\tilde\beta_{jp} \;=\; \Big\langle \frac{\partial \log P(j\mid x_i)}{\partial x_p}\Big\rangle_i \;-\; \frac1J\sum_{j'}\Big\langle\frac{\partial \log P(j'\mid x_i)}{\partial x_p}\Big\rangle_i, \tag{9.1}
$$

centred over fine classes ⇒ **zero-sum**, hence directly comparable to the flat coefficients. The derivative is the chain rule over the path (3.1): at each internal node,
$\partial_p\log P(c_k\mid c)=\partial_p v_{c_k}-\sum_{k'}P(c_{k'})\,\partial_p v_{c_{k'}}$ with
$\partial_p v_{c_k}=\gamma_{c_k,p}+\sum_{c'}b_{IV_{c'},c_k}\,\partial_p IV_{c'}$, and the inclusive-value derivative is itself recursive: $\partial_p IV_{\text{leaf}}=\sum_m P(m)\beta_{m,p}$, $\partial_p IV_{\text{nest}}=\sum_{c'}P(c')\partial_p v_{c'}$. `effective_symmetric_table(fit, X)` computes this per imputation → **median + 95% CI** (leaf/IV uncertainty propagated). *It is a derived summary — the model's native, exact parameters remain the per-nest symmetric matrices.* Note the marginal effect is $x$-dependent (nested-logit reality); (9.1) reports its sample average.

---

# Appendix A — Computation algorithm

**A.1 Recursive cut fit** `.ncut_node(node, X, Y, in_pixels, …)`

```
fit_node(node, in_pixels):
  if node is terminal:
     if singleton -> return {singleton}
     keep <- in_pixels & (rowSums(Y[,fine])>0);  if too few -> {even-split}
     fit  <- mnlogit_rcpp_sym( X[keep], normalize(Y[keep,fine]) )      # PG-Gibbs, symmetric
     return {leaf, beta_draws = M subsampled posterior draws of postb_pooled}   # (6.2)
  else (internal, children c_1..c_K):
     for each child: children[c] <- fit_node(child, in_pixels & child-present)   # bottom-up
     Ynode <- column-bind child memberships (aggregated shares)                  # (3.1) counts
     warm <- NULL                                                                # hot-start state
     for m in 1..M:                                                              # (6.2)-(6.3)
        IV^{(m)} <- [ node_iv(children[c], X, m) for c in iv_children ]          # recursion (3.3)
        X^{(m)}  <- cbind(X, IV^{(m)})            (skipped if use_iv=FALSE)
        (nb,ni)  <- (m==1 ? (nburn, niter) : (nburn_warm, nburn_warm+niter-nburn))   # §6.4
        fit_m    <- mnlogit_rcpp_sym( X^{(m)}[keep], normalize(Ynode[keep]), init_state=warm, nburn=nb, niter=ni )
        warm     <- fit_m$final_state                                            # thread into m+1
        parent_draws[[m]] <- R draws of fit_m$postb_pooled
        lambda_draws[c]   <- append( coef(IV_c on c) * K/(K-1) )                 # (4.2)
     return {nest, children, parent_draws, lambda_draws}
```

**A.2 Inclusive-value recursion** `.ncut_node_iv(nd, X_new, m)` — evaluates (3.3) for imputation $m$: leaf ⇒ `lse(X_new %*% beta_draws[[m]])`; internal ⇒ rebuild `[X_new ‖ child-IV^{(m)}]`, take one parent draw, `lse(·)`; singleton/degenerate ⇒ 0.

**A.3 Coherent predictive** `.ncut_predict_draw(nd, X_new, share, m, P)` — descends the tree distributing probability mass `share`; at each internal node computes the child softmax from `[X_new ‖ IV^{(m)}]` and one parent draw, recurses on `share × child-prob`; at a leaf distributes `share × sub-softmax`. `predict_nested_cut` loops $D$ draws (each an imputation index $m$) → array `[n × J × D]`.

**A.4 Summary** `summary_nested_cut` — per node: posterior-mean zero-sum coefficient matrix (rows = covariates incl. IV, cols = choices) + $\lambda$ table (median, 2.5%, 97.5%) from the pooled draws.

# Appendix B — Complexity, cost, and the PG reuse

- **Fits performed:** one sub-model fit per terminal multi-class nest (**once**), plus $M$ parent fits per internal node. Depth does **not** multiply $M$ — the imputation index is shared down the tree, so cost is $O\!\big(\sum_{\text{internal }c} M\big)$ parent fits, not $M^{\text{depth}}$.
- **Per fit:** a standard symmetric Pólya-Gamma MNL Gibbs update — $O(\text{niter}\cdot n_c\cdot K_c\cdot p)$ with conjugate Gaussian $\beta$ draws (no MH, no HMC). Parent fits are cheap (macro level, few classes); the multi-class leaves (e.g. 14 arable crops) dominate but run once.
- **Memory:** each node stores $M$ leaf draws (leaf) or $M\times R$ parent draws (internal) of small $[k\times K]$ matrices.
- **Why PG survives:** conditional on an imputed $\mathrm{IV}^{(m)}$, every node model is an *ordinary* multinomial logit in its own parameters → full PG conjugacy is retained. The non-conjugacy of the true joint is confined to, and resolved by, the imputation/cut.

# Appendix C — Identification, consistency, caveats

- **Consistency, not full efficiency.** The cut estimator is consistent for $\theta_c$ and $\lambda$; it is limited-information, so credible intervals are wider than (but honest, unlike) the plug-in, and narrower than the — unavailable here — full joint.
- **λ range.** $\lambda\in(0,1]$ is *not* imposed; a fitted $\hat\lambda>1$ or $<0$ flags a nest inconsistent with random-utility nesting (diagnostic, not an error).
- **Generated-regressor attenuation** is *reduced* (uncertainty propagated) but not eliminated — the point $\hat\lambda$ can still attenuate slightly toward 0 when the shared covariates already carry the signal; the IV payoff is empirical.
- **Degenerate nests** (singletons, or nests too sparse in a fold) contribute $\mathrm{IV}=0$ and split mass evenly — they carry no coupling.
- **No feedback by design.** If the macro data are highly informative about the sub-utilities, the cut discards that information; this is the deliberate robustness/tractability trade (§6.1).

*References: McFadden (1978) nested logit & inclusive values; Plummer (2015) "Cuts in Bayesian graphical models"; Jacob, Murray, Holmes, Robert (2017) "Better together? Statistical learning in models made of modules"; Rubin (1987) multiple imputation.*
