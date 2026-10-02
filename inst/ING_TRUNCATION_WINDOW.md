# The ING tau^2 truncation window

How `Prior_Setup_GLMM()` chooses the default `disp_lower` / `disp_upper`
bounds for `dIndependent_Normal_Gamma` (ING) Block-2 dispersion priors, and
why they are quantiles of the *limiting posterior* rather than of the prior.

## Setup

For each random-effect component `k`, the Block-2 hyper-regression treats the
`J` group-level coefficients as observations (`n_w = J`) with `p_k` Block-2
predictors and classical RE-variance estimate `tau2_k`.  The ING prior on the
precision `1/tau2` follows the glmbayesCore default calibration
(`compute_gaussian_prior()` with `k = 1`; see glmbayesCore Chapter A12):

```
shape = (n0 + 1 + p_k) / 2          # shape_ING convention
rate  = tau2_k * (n0 + p_k - 1) / 2 # = tau2_k * (shape - 1)
```

where `n0 = n_prior_dispersion` (from `pwt_dispersion`, by default derived
from the coefficient `pwt`; weak values are fine precisely because the
truncation window below does not depend on `n0`).  Because
`rate = tau2_k * (shape - 1)`, the implied inverse-Gamma
prior on `tau2` has mean exactly `tau2_k` for every `n0` and `p_k`.

The sampler (`glmbayesCore::two_block_rNormal_reg_v2`) requires *both*
truncation bounds: each `tau2_k` draw is hard-truncated to
`[disp_lower, disp_upper]` by renormalized inverse-CDF sampling.  A fixed,
two-sided window serves three purposes:

1. `disp_lower > 0` is the conservative plug-in for the eigenvalue / total
   variation calibration of `m_convergence` (smaller `tau2` means stronger
   block coupling, so the bound computed at `disp_lower` dominates the whole
   truncated support).
2. A finite window guarantees geometric ergodicity of the two-block chain.
3. Fixing the window across all inner Gibbs sweeps keeps the truncation
   state-independent (the alternative -- per-sweep posterior-derived bounds,
   as glmbayes uses for one-shot fits -- would make the invariant
   distribution sweep-dependent).

## Why not prior quantiles?

The first implementation used the central 98% *prior*-mass interval
(0.01/0.99 quantiles of the inverse-Gamma prior).  That window degenerates as
the prior weakens: the prior shape is roughly `n0/2`, so `n0 -> 0` makes the
prior heavy-tailed and the window stretches without bound -- while the
conditional posterior, whose Gamma shape gains `J/2` from the data
*regardless of* `n0`, stays concentrated near `tau2_k`.  Consequences (schools
example, `J = 47`, intercept component):

| `pwt_dispersion` | prior window  | width ratio | candidates per accepted draw |
|------------------|---------------|-------------|------------------------------|
| 0.4              | [150, 469]    | 3.1         | ~4.6                         |
| 0.2              | [110, 656]    | 5.9         | ~28                          |
| 0.1              | [80, 934]     | 11.7        | ~159                         |

The envelope accept-reject sampler must cover the whole window, so its
acceptance rate collapses as the prior weakens -- paying an exploding
computational price for tail regions the posterior never visits (the
window's posterior coverage tends to 100% as `n0 -> 0`).

## The limiting-posterior window

glmbayesCore Chapter A12 (Theorem 2) gives the weak-prior limit of the
posterior Gamma for the precision as `n_prior -> 0`.  Applied to the Block-2
hyper-regression with the mean-matched plug-in `tau2_k`:

```
a_inf = (J + 1) / 2
b_inf = tau2_k * (J - 1) / 2     # so b_inf / (a_inf - 1) = tau2_k
```

The default window is the central 98% mass of this limiting law, inverted to
the `tau2` scale:

```
disp_lower = 1 / qgamma(0.99, shape = a_inf, rate = b_inf)
disp_upper = 1 / qgamma(0.01, shape = a_inf, rate = b_inf)
```

Properties:

- **Independent of `n0` and `p_k`.**  One window per component, computable at
  prior-setup time from `J` and `tau2_k` alone -- so it remains a fixed
  truncation across all Gibbs sweeps.
- **Mean-matched.**  The limiting law's implied dispersion mean is exactly
  `tau2_k`, like the prior, so the window always brackets the classical
  estimate.
- **Asymptotically exact coverage.**  As `n0 -> 0` the exact posterior
  converges to the limiting law, so the bounds become genuine 0.01/0.99
  posterior percentiles.  For any `n0 > 0` the exact posterior has shape
  `(n0 + 1 + J)/2 > a_inf` (strictly more concentrated), so coverage exceeds
  98% -- the window is conservative for every prior strength.  Schools
  example: coverage 99.7% / 99.1% / 98.6% / 98.1% at `pwt_dispersion`
  0.4 / 0.2 / 0.1 / 0.02.
- **Stable sampling cost.**  Width ratio ~2.6 for `J = 47` (vs. 12+ for weak
  prior-based windows), keeping the envelope's candidates-per-draw roughly
  constant as priors weaken.

## Caveats

- The truncation is part of the model: ~1% of posterior mass is clipped in
  each tail (asymptotically).  This slightly shortens the reported extreme
  `tau2` quantiles.  Diagnose with the `Cand/draw` column of
  `summary()` and the printed window in `print(Prior_Setup_GLMM(...))`.
- The per-sweep conditional posterior's rate fluctuates with the current
  group-effect draws; the coverage statements above use the mean-matched
  plug-in rate.  The shape -- which controls the width -- is exact.
- Users can always override the window by constructing
  `dIndependent_Normal_Gamma()` pfamilies manually with their own
  `disp_lower` / `disp_upper` (both required for sampling, and the
  `n_prior <= J` guard still applies).
