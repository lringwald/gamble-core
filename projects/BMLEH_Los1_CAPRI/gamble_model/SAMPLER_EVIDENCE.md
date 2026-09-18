# Why the BMLEH sampler is configured the way it is

The runs these numbers came from were deleted on 2026-09-11: they were fitted on a design that was
wrong in three ways (AgMIP class set, GLOB_country RE keys, lagged focal with no crop types). The
SAMPLER conclusions survive the design being wrong, so they are kept here rather than re-derived.

## RE-scale ASIS is what makes the posterior quotable — `BM_RE_ASIS=TRUE`

Controlled comparison: identical config, identical hot-start states, only `re_asis` differed.

| block      | no ASIS            | ASIS               |
|------------|--------------------|--------------------|
| log_lik    | Rhat 1.878, ESS 6  | **1.037, ESS 89**  |
| sigma_re   | Rhat 1.618, ESS 7  | **1.007, ESS 470** |
| mu (FE)    | 1.038 / 1.809 max  | 1.018 / 1.351      |
| b_g (RE)   | 1.030 / 1.820 max  | 1.015 / 1.347      |
| kappa (HS) | 1.129 / 1.895 max  | 1.046 / 1.465      |

Per-chain `sigma_re` without ASIS: 0.0576 / 0.0779 / 0.0559 / 0.0773 (spread 0.0220).
With ASIS: 0.1091 / 0.1080 / 0.1094 / 0.1063 (spread 0.0031).

Cost: on equal draws ASIS scored **0.2391** held-out against the un-mixed control's **0.2455**.
Better mixing finds a LARGER RE variance, so the stuck chains were over-pooled and scored marginally
better while being unconverged. Converged is worth 0.006 McFadden for a prior consumed with intervals.
The only other cost is ~14% median ESS on mu/b_g, the expected price of an extra interweaving move.

## The slab was never the problem — `BM_SLAB_C2=FALSE`

`c2` sat at ~9-26 with Rhat 2.59, which looked like the cause. It is not: the slab caps variance at
sqrt(c2) ~ 3-5 while sigma_re is 0.108, so it never binds. Fixing it changed nothing
(log_lik Rhat 1.878 fixed vs 1.706 sampled). Sampling it only adds the worst-mixing scalar for free.

Two further hypotheses tested and rejected:
- FE<->RE aliasing: sigma_re vs FE size correlates **+0.99**, not negative; the larger-sigma chains
  have LARGER fixed effects and fit WORSE (-0.996 with log_lik). Not a trade-off.
- horseshoe: kappa_pooled differs 4% between chain groups, c2 identical.

What it actually was: sigma_re still drifting upward after 5000 sweeps (p 1e-26 to 1e-40) with
log_lik flat -- four chains at different points on one slow climb, not separate modes.

## More sweeps alone does not fix it

Segment 1 -> segment 2 added 3000 sweeps at thin 4 (2.5x the draws):
log_lik Rhat 1.903 -> 1.706, ESS 6 -> 6. Within-chain sampling improved 4x (mu ESS 400 -> 1612)
while the joint chain did not converge. That is why ASIS, not budget, was the answer.

## Resume is exact — `run.sh more`

A resumed chain restarts **91 nats** from where it stopped, against an **83-nat** draw-to-draw noise
floor; a cold start from the same point is 2970 away. Per-chain identity preserved 4/4 on the
parameter state (diagonal 0.7789 vs off-diagonal 0.9334, against a 0.7536 within-chain floor).
`final_state` carries `a_re`, so ASIS continues rather than restarting.
