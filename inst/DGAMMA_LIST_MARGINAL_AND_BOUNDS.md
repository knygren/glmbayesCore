# Per-group `dGamma()` measurement dispersion: the marginal Gamma and its bounds

How `dGamma_list()` (`R/dGamma_list_Prior_Setup_GLMM.R`) builds a
per-group Independent Normal--Gamma (ING) prior on Block~1 measurement
dispersion `sigma^2_j`, and how it derives that prior's truncation window.
Companion to `inst/TAU2_ING_FORMULAS.md` (the same conjugate machinery for
`tau^2_k`) and to `lmebayes`'s `inst/GROUP_DISPERSION_HYPERPRIOR.md` /
`inst/ING_TRUNCATION_WINDOW.md`. Reproducible with
`data-raw/group_dgamma_bounds_derivation_check.R`.

**Status (current default):** the Gamma density fed to `dGamma()` is
`g$rate` -- the Chapter A12 **?3.3.4 marginal ING rate** (`beta` integrated
out) -- not `g$rate_gamma` (**?3.3.5**'s fixed-`beta` rate). `rate_gamma` is
still computed and stored for diagnostic comparison only (see Part I and the
dev-only print in `Prior_Setup_GLMM()`). This is a change from an
earlier version of this document/code, which fed `rate_gamma` to the
sampler; Part I below explains why `rate` is the theoretically correct
choice for this particular sampler, and Part II shows why it also improves
alignment with the (unchanged) truncation window.

**2026 update -- Part VI is now the permanent default, and the bounds moved
upstream.** Two things changed after Parts II/III/V below were written; see
**Part VII** for the current implementation:

1. Part VI's model-derived `Omega_j` fold-in ("also integrating out the
   prior mean `mu_j`", i.e. the fixed effects `gamma`, not just `b_j`) is no
   longer an unimplemented extension -- it is now applied unconditionally
   inside `Prior_Setup_GLMM()`/`.lmebayes_calibrate_ing_prior_measurement_group()`.
   There is no opt-out.
2. `disp_lower`/`disp_upper` are no longer built by `dGamma_list()` from a
   separate, `n_combined,j`-based `shape_w`/`rate_w` (Part II's construction,
   with the `disp_center = "dispersion2"` alternative from Part III). They
   are now computed once inside `Prior_Setup_GLMM()`, as literal
   `(1-max_disp_perc_measurement)`/`max_disp_perc_measurement` quantiles
   (scalar, or a per-group vector -- see `Prior_Setup_GLMM()`'s
   `max_disp_perc_measurement` argument) of the **same**
   `Gamma(shape_ING,j, rate,j)` that is the actual sampling prior (the
   Part-VI-extended one) -- i.e. a truncated version of the real prior,
   not a separately-mean-matched window. `dGamma_list()` now just reads
   `disp_lower`/`disp_upper` off `ing_prior_measurement_group` verbatim.

Parts II/III's `shape_w`/`rate_w`/BLUP-inflation/`dispersion2`-envelope-
centering machinery (and the `disp_center`, `disp_upper_anchor`, and
`n_rss_iter` arguments that controlled it) has been removed from
`dGamma_list()`; those sections are kept below **for historical record and
worked-example intuition only** -- they no longer describe shipped
behavior.

---

## 0. Setup and notation

For each group (school) `j = 1..J`:

| Symbol | Code | Meaning |
|---|---|---|
| `n_j` | `g$n_j` | Observations in group `j` |
| `X_j` | `model.matrix(block_formula, dat_j)` | `n_j x p_re` design (RE-relevant predictors only; see `.lmebayes_block_formula_from_re()`) |
| `Y_j` | ��� | Response, group `j` |
| `bhat_j` | `coef(glm.fit(X_j, Y_j))` | Group `j`'s own OLS/WLS estimate |
| `mu_j` | `.lmebayes_block_formula_prior_mu()` | Prior mean: null-model intercept, zero slopes (defaults) |
| `V0_j` | `vcov(glm_full)` | Sampling covariance of `bhat_j` |
| `sd_tau` | `sqrt(diag(Sigma_ranef))` | Population RE sd, shared across groups |
| `w_j` | `pwt_measurement` (scalar or length-`J`) | Group `j`'s Block~1 prior weight |
| `n_prior,j` | `w_j/(1-w_j) * n_j` | Group `j`'s effective prior sample size |
| `n_combined,j` | `n_prior,j + n_j` | Total effective sample size |
| `p_re` | `length(design$re_coef_names)` | Number of RE components (2 in the worked example: intercept + `distracted_ppvt`) |

`.lmebayes_calibrate_ing_prior_measurement_group()` (`R/mixed_rmerb_helpers.R`)
builds, **independently for each group**, a coefficient-scale prior covariance

```
pwt_j    <- diag(V0_j) / (diag(V0_j) + sd_tau^2)          # per-coefficient shrinkage weight
Sigma_j  <- V0_j * outer(sqrt((1-pwt_j)/pwt_j), sqrt((1-pwt_j)/pwt_j))   # coefficient-scale prior cov
Sigma_0j <- Sigma_j / dispersion_classical_j               # dispersion-free (Zellner form)
```

and calls `compute_gaussian_prior(X_j, Y_j, weights=1, offset=0, dispersion=NULL,
n_effective=n_j, bhat=bhat_j, mu=mu_j, Sigma_0=Sigma_0j, Sigma=Sigma_j,
n_prior=n_prior,j, k=1)`.

---

## Part I -- The marginal Gamma distribution (per-group prior)

This section specializes `vignettes/Chapter-A12.Rmd` Sections 3.1, 3.3.1, and
3.3.4-3.3.5 (the general Normal--Gamma / ING / `dGamma()` calibration) to one
group's Block~1 design; see that vignette for the full derivation and proof
that this is the correct marginal (i.e. `beta` integrated out) Gamma law on
the residual precision.

**Marginal quadratic term** (Chapter A12 ?3.1, specialized to group `j`,
weights `=1`):

```
S_marg,j = RSS_ols,j + (bhat_j - mu_j)' * Mi_j * (bhat_j - mu_j)
Mi_j     = (Sigma_0j + (X_j'X_j)^{-1})^{-1}
RSS_ols,j = sum((Y_j - X_j %*% bhat_j)^2)
```

The second term is a **penalty for how far the group's own OLS fit sits from
its prior mean `mu_j`** (relative to how confident the prior is, via
`Sigma_0j`), on top of the group's raw residual sum of squares.

**Calibrated point estimate** (`g$sigma2_hat`, `cal$dispersion`):

```
sigma2_hat,j = S_marg,j / (n_j - p_re)
```

Because `S_marg,j >= RSS_ols,j`, **`sigma2_hat,j >= sigma2_ols,j` whenever
`bhat_j != mu_j`** -- the quad-penalty term always inflates the calibrated
dispersion above the raw OLS residual variance. This gap is largest exactly
for groups whose own fit deviates most from the shared null-model prior mean
(see Part IV: this is also where `blup_infl,j` tends to be largest).

### Two candidate Gamma rates from the same `compute_gaussian_prior()` call

`.lmebayes_compute_ing_prior_cal_from_sigma()` calls
`compute_gaussian_prior()` once per group and gets back **two** rates that
answer different questions (`R/compute_gaussian_prior.R`):

```
shape_ING,j = (n_prior,j + 1)/2 + p_re/2                                    # k = 1; shared by both rates below

## ?3.3.4 -- marginal rate (beta integrated out against N(mu_j, Sigma_j)):
rate,j       = 0.5 * S_marg,j * (n_prior,j + p_re - 1) / (n_j - p_re)

## ?3.3.5 -- fixed-beta rate (RSS at the point-blend beta_star,j):
beta_star,j  = (1 - w_j) * bhat_j + w_j * mu_j                              # Zellner blend at group j's own w_j
RSS_star,j   = sum((Y_j - X_j %*% beta_star,j)^2)
rate_gamma,j = 0.5 * RSS_star,j * (n_prior,j + p_re - 1) / (n_j - p_re)
```

Both are calibrated **per group**, using group `j`'s own `bhat_j`, `mu_j`,
and `w_j` -- neither is one shared prior across all `J` groups (contrast with
the "consistent design" `lmebayes`'s `inst/GROUP_DISPERSION_HYPERPRIOR.md`
?8 argues for; see Part V).

**Key algebraic fact:** `rate,j / (shape_ING,j - 1) = S_marg,j / (n_j - p_re)
= sigma2_hat,j` **exactly**, for every `n_prior,j`, `k`, and `p_re` (the same
cancellation `compute_gaussian_prior()`'s own Step E comment documents for
the un-inflated `shape`/`rate` pair). `rate_gamma,j` has no such identity:
`RSS_star,j` need not equal `S_marg,j` (they differ whenever `w_j != 0`, and
empirically `S_marg,j > RSS_star,j` for essentially every group on the
39-school fixture below), so `rate_gamma,j / (shape_ING,j - 1) != sigma2_hat,j`
in general.

### The Gamma actually fed to the sampler: `rate` (?3.3.4), not `rate_gamma` (?3.3.5)

`dGamma_list.Prior_Setup_GLMM()` passes `g$rate` -- **not**
`g$rate_gamma` -- to each group's `dGamma(shape = g$shape_ING, rate = g$rate,
...)`. `rate_gamma` remains on `object$ing_prior_measurement_group` (and is
printed alongside `rate` by the dev-only comparison table
`Prior_Setup_GLMM()` emits whenever `dispformula` requests per-group
dispersion) purely for diagnostic comparison; it is not consumed downstream.

**Why `rate`, not `rate_gamma` -- the choice is forced by what the sampler
does with `beta`.** The two rates are the Gamma parameters of two different
conditional statements about `sigma^2_j`:

- `rate` (?3.3.4): the rate of the **marginal** law for `1/sigma^2_j`, with
  `beta_j` integrated out against its own prior `N(mu_j, Sigma_j)`. Correct
  when `beta_j` is itself a random quantity yet to be drawn.
- `rate_gamma` (?3.3.5): the rate **conditional on `beta_j` fixed** at the
  specific point `beta_star,j`. Correct when the consumer of the prior treats
  `beta` as a known plug-in and never updates it (`?dGamma`'s general
  description: a Gibbs step "where the beta and dispersion parameters are
  updated separately", with `beta` "held fixed" -- the intended use for a
  one-shot `Prior_Setup()` + `dGamma()` illustration).

Which one is right depends entirely on what the actual consumer of
`(shape_ING, rate)` does with `beta_j`. Trace it for the Block~1 ING
sampler `dGamma_list()` output feeds:

1. `dGamma_list(ps)` becomes `lmerb(..., dispersion_ranef = dGamma_list(ps))`'s
   per-group `dispersion_ranef`.
2. `.two_block_block1_ing_group_draw_one_chain()` (`R/rLMM_reg.R`) calls
   `.rLMM_ing_one_group_draw()` once per group, per sweep.
3. On the (typical) full-rank path, that calls `rindepNormalGamma_reg()`
   (`R/simfunction.R`), which performs a **joint** Independent Normal-Gamma
   draw: `sigma2_j` is drawn first from the *marginal* Gamma (`beta_j`
   integrated out -- exactly the ?3.3.4 law), and **then** `b_j | sigma2_j`
   is drawn from its Normal conditional.

`beta_j` (the group's random effect, `b_j`) is **not** fixed at any point
value during this step -- it is simulated fresh, every sweep, for every
group. Because the sampler's own derivation integrates `beta_j` out before
conditioning on `sigma2_j`, the prior handed to it must already be that same
marginal law -- i.e. `(shape_ING, rate)` from ?3.3.4. Feeding it
`rate_gamma` (?3.3.5) instead silently substitutes a Gamma calibrated
*conditional on a specific point `beta_star,j`* for the *marginal* Gamma the
joint draw's own derivation assumes -- understating `E[sigma^2_j]` by
however much the quadratic-penalty gap (`S_marg,j - RSS_ols,j`) exceeds the
Zellner-blend RSS gap (`RSS_star,j - RSS_ols,j`). This is exactly the
`pct_rate` column the dev-only comparison print reports per group (see
worked numbers below); it is not a rounding-level discrepancy -- on the
39-school fixture it reaches 40-55% for the most BLUP-inflated schools.

**"Conditional on Block 2" qualifier.** Both `mu_j` (the prior mean for
`b_j`) and `Sigma_j` (built from `sd_tau = sqrt(diag(Sigma_ranef))`) are
supplied by the *current* Block~2 hyperparameters -- Block~2's own
random-effect mean/covariance is treated as **given, fixed input** to this
Block~1 calibration step, not re-derived from Block~1's own data. So the
?3.3.4 marginal used here is the correct precision prior for the `(b_j,
sigma2_j)` joint draw **given Block~2's current RE covariance and mean** --
it is not a full joint marginal that additionally integrates over
uncertainty in `sd_tau`/`Sigma_ranef` themselves. (Part III's `dispersion2`
construction is the closest existing approximation to that additional
layer, and it is applied only to the truncation *bounds*, not to the prior
rate itself; see Part III.) Within a two-block Gibbs sampler this is the
standard and correct convention -- each block's own conditional prior is
calibrated treating the *other* block's current state as fixed -- and it is
exactly the same convention `inst/TAU2_ING_FORMULAS.md` / `ING_TRUNCATION_WINDOW.md`
use for Block~2's own `tau^2_k` priors (conditional on Block~1's `bhat`).

### Worked numbers: `rate` vs. `rate_gamma` (39-school `big_word_club` fixture)

From the dev-only comparison table `Prior_Setup_GLMM(..., dispformula =
~school_id)` prints (`.lmebayes_print_ing_prior_measurement_group_compare()`,
`R/mixed_rmerb_helpers.R`); `pct_rate = 100*(rate - rate_gamma)/rate_gamma`:

```
group  shape_ING  rate_gamma    rate    pct_rate
  33     2.1667      246.90    383.90    +55.5%
  41     2.1111      316.90    476.61    +50.4%
   6     2.1111      306.64    435.33    +42.0%
   3     2.2222     1217.71   1223.49     +0.5%
  20     2.3333      420.81    420.84     +0.0%
```

Groups 6, 33, and 41 are exactly the three groups flagged
`asymmetric_window = TRUE` in Part II below -- their own fit deviates most
from the null-model prior mean, so both the quad-penalty gap (driving
`rate`'s advantage over `rate_gamma`) and `blup_infl,j` (driving the
window's upper-tail widening) are largest for the same underlying reason.
Groups where `beta_star,j approx bhat_j` (small `w_j`, well-behaved fit)
show `pct_rate approx 0`, as expected -- both rates agree when the
fixed-point RSS and the marginal RSS nearly coincide.

---

## Part II -- The truncation window (HISTORICAL -- superseded by Part VII)

**No longer shipped behavior.** `dGamma_list()` no longer builds this
`shape_w`/`rate_w`/`n_combined,j` window, and the `disp_upper_anchor`
argument described below has been removed. Kept for the worked-example
intuition and as a record of the design this superseded; see Part VII for
what `disp_lower`/`disp_upper` are computed from now.

### Why a window at all

Same three reasons as `lmebayes`'s `inst/ING_TRUNCATION_WINDOW.md` for
`tau^2_k`: a `disp_lower > 0` plug-in for convergence calibration, guaranteed
geometric ergodicity of the envelope accept/reject sampler, and a
sweep-independent (fixed) truncation so the invariant distribution does not
drift across Gibbs sweeps.

### Mean-matched construction

`dGamma_list.Prior_Setup_GLMM()` (`R/dGamma_list_Prior_Setup_GLMM.R`)
does **not** reuse either Part I rate for the window. Instead it builds a
fresh Gamma whose **mean is pinned to `sigma2_hat,j` by construction**, with
`n_combined,j` controlling only the spread:

```
shape_w,j <- (n_combined,j + 1)/2 + p_re/2
rate_w,j  <- sigma2_hat,j * (n_combined,j + p_re - 1)/2
```

Algebraically, `rate_w,j / (shape_w,j - 1) = sigma2_hat,j` **for every
`n_combined,j`** -- mean and spread are deliberately decoupled, so widening or
narrowing the window (by changing what feeds `n_combined,j`) never moves its
center.

### Alignment with Part I's prior mean

Because `rate,j / (shape_ING,j - 1) = sigma2_hat,j` exactly (Part I), the
window built here is mean-matched to **the same point** the ?3.3.4 prior
itself asserts for `E[sigma^2_j]` -- prior and window now share a center by
construction, for every group, regardless of `pwt_measurement`,
`max_disp_perc`, or `disp_upper_anchor`.

This was **not** true while `rate_gamma` fed the sampler: whenever the
quad-penalty gap exceeded the Zellner-blend RSS gap (Part I), the prior's
own implied mean (`rate_gamma,j/(shape_ING,j-1)`) could sit measurably below
`sigma2_hat,j` -- on the 39-school fixture, groups 33 and 41 had their prior
mean fall *below* `disp_lower,j` entirely (e.g. group 33: prior mean 211.6 vs.
`disp_lower = 226.2`), which is visible in `lmebayes::summary_sigma2()` as
`Pr(Prior_tail) approx 0` and elevated `Cand/draw` for exactly those groups.
Switching to `rate` moves each group's prior mean to `sigma2_hat,j`, which by
construction always lies strictly inside `[disp_lower,j, disp_upper,j]`
(the window is centered there). This does not change `disp_lower`,
`disp_upper`, `blup_infl`, `R_lo`, or `R_hi` at all -- only where the prior
density sits *relative to* that unchanged window.

### Asymmetric upper tail (`disp_upper_anchor = "blup"`, the default)

`lmer`'s BLUP for group `j`'s random effect can sit well away from group
`j`'s own OLS fit. `.lmebayes_group_blup_rss_inflation()` computes

```
blup_infl,j = RSS_blup,j / RSS_ols,j >= 1
```

(`RSS_blup,j` = residual sum of squares at the reference `lmer` fit's BLUP
coefficients for group `j`; `RSS_ols,j` at group `j`'s own OLS fit). The upper
tail is then anchored at an inflated rate:

```
rate_u,j <- rate_w,j * blup_infl,j        (rate_u,j == rate_w,j when disp_upper_anchor = "symmetric")

disp_lower,j <- 1 / qgamma(max_disp_perc,     shape = shape_w,j, rate = rate_w,j)
disp_upper,j <- 1 / qgamma(1 - max_disp_perc, shape = shape_w,j, rate = rate_u,j)
```

`.lmebayes_dgamma_window_cross_percentiles()` additionally reports cross
percentiles `R_lo`, `R_hi` diagnosing how asymmetric the resulting window is
(flagged via `asymmetric_R_lo` / `asymmetric_R_hi`, printed with
`warn_asymmetric = TRUE`).

### Worked numbers (39-school `big_word_club` fixture, `p_re = 2`, `max_disp_perc = 0.8`)

```
school 41: n_j=11  sigma2_hat=428.6  blup_infl=1.52  ->  disp_lower= 289.6  disp_upper= 823.3
school  9: n_j=11  sigma2_hat=137.1  blup_infl=1.26  ->  disp_lower=  92.6  disp_upper= 217.5
```

---

## Part III -- Integrating over random-effect uncertainty (HISTORICAL -- superseded by Part VII)

**No longer shipped behavior.** The `disp_center`/`n_rss_iter` arguments and
the `EnvelopeCentering()`-style trace correction described below have been
removed from `dGamma_list()`. Kept for the worked-example intuition; see
Part VII for the current bounds construction.

Was implemented as an opt-in argument: `dGamma_list(ps, disp_center = "dispersion2")`
(default remains `disp_center = "sigma2_hat"`, i.e. Part II unchanged). See
Part V for status.

### Motivation

`lmer`'s BLUP shrinkage assumes **one pooled residual variance across all
groups**. When a group's *true* `sigma^2_j` differs substantially from the
pooled value, `lmer` mis-states that group's sampling uncertainty and
shrinks its BLUP too much (if `sigma^2_j` is smaller than pooled) or too
little (if larger) -- inflating or deflating `RSS_blup,j` relative to
`RSS_ols,j` for reasons that have nothing to do with what the *bounds* should
actually cover. `blup_infl,j` is a reasonable, cheap proxy for "how much
wider does the upper tail need to be", but it borrows its magnitude from a
model (`lmer`, homogeneous variance) that disagrees with the one being
fit (heterogeneous per-group variance) -- see Part IV for direct evidence
this heterogeneity is real.

### `EnvelopeCentering`'s trace correction (`src/EnvelopeCentering.cpp`)

At a fixed working dispersion `dispersion2`, the posterior mean of the random
effect `b_j` under prior `b_j ~ N(mu_j, Sigma_ranef)` and Gaussian likelihood
is the ridge/GLS estimator

```
b2_j = (X_j'X_j / dispersion2 + P)^{-1} (X_j'Y_j / dispersion2 + P %*% mu_j),   P = Sigma_ranef^{-1}
Cov(b2_j) = (X_j'X_j / dispersion2 + P)^{-1}
```

The **expected** residual sum of squares under this posterior -- not just the
RSS *at* the posterior mean -- follows from the standard bias--variance
(law-of-total-variance) decomposition of `E[(Y_j - X_j b_j)'(Y_j - X_j b_j)]`
over `b_j`'s posterior:

```
RSS_precomputed,j = ||Y_j - X_j %*% b2_j||^2  +  tr(X_j'X_j %*% Cov(b2_j))
                     ^^^^^^^^^^^^^^^^^^^^^^^     ^^^^^^^^^^^^^^^^^^^^^^^^^^
                     RSS at the posterior mean    correction for b_j's own
                                                   remaining posterior uncertainty
```

`dispersion2` and `RSS_precomputed,j` are then updated to a fixed point (10
iterations in both `EnvelopeCentering.cpp` and the reproduction script):

```
dispersion2 <- RSS_precomputed,j / (n_j - p_re)     # data-only path -- see pitfall below
```

### Mean-matching the refined estimate

Exactly the Part II construction, with `sigma2_hat,j` replaced by the
converged `dispersion2,j`, and **no separate BLUP widening** -- the trace term
already prices in `b_j`'s uncertainty directly, symmetrically:

```
shape_w_prop,j <- (n_combined,j + 1)/2 + p_re/2      # unchanged
rate_w_prop,j  <- dispersion2,j * (n_combined,j + p_re - 1)/2   # mean-matched at dispersion2,j, not sigma2_hat,j

disp_lower_prop,j <- 1 / qgamma(max_disp_perc,     shape = shape_w_prop,j, rate = rate_w_prop,j)
disp_upper_prop,j <- 1 / qgamma(1 - max_disp_perc, shape = shape_w_prop,j, rate = rate_w_prop,j)
```

### The double-counting pitfall

**Do not** seed the `dispersion2` fixed point (or the final mean-match) with
group `j`'s own `shape_ING,j` / `rate,j` (or `rate_gamma,j`) as if they were
a fresh prior to be updated with `n_j/2`, `RSS_precomputed,j/2` on top.
`shape_ING,j` / `rate,j` (Part I) are **already** calibrated from group
`j`'s full `n_j` observations (via `S_marg,j`, itself derived from `bhat_j`)
-- just scaled
down to an equivalent `n_prior,j`-sized weight for a different purpose (the
prior object handed to the sampler). Running the same `n_j` observations
through the update a second time inflates the effective weight of the
group's own data well beyond what `n_prior,j` is meant to represent, and
produces artificially tight windows. Concretely, on the 39-school fixture this
mistake pushed school 9's BLUP point estimate (`sigma2_blup = 160.4`) outside
its own proposed upper bound (`disp_upper_prop = 157.2`) -- a symptom that a
group's own reference-model point estimate should essentially never land
outside an 80%-mass window built for exactly that group. Mean-matching
directly at `dispersion2,j` with the *same* `n_combined,j` convention Part II
already uses (rather than re-deriving shape/rate through `shape_ING,j` +
`n_j/2`) avoids this: `n_j` enters through `RSS_precomputed,j`'s point
estimate exactly once.

### Worked numbers (same fixture and setting)

| group | n_j | sigma2_ols | sigma2_blup | blup_infl | dispersion2 (proposed center) | disp_lower (current) | disp_upper (current) | disp_lower (proposed) | disp_upper (proposed) | %delta lower | %delta upper |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 41 | 11 | 283.6 | 432.3 | 1.52 | 472.6 | 289.6 | 823.3 | 319.3 | 595.4 | +10.3% | -27.7% |
| 33 | 12 | 209.4 | 292.5 | 1.40 | 339.7 | 225.3 | 573.7 | 233.5 | 425.8 | +3.7% | -25.8% |
| 6  | 11 | 272.9 | 373.3 | 1.37 | 404.4 | 263.5 | 672.2 | 273.3 | 509.6 | +3.7% | -24.2% |
| 9  | 11 | 127.4 | 160.4 | 1.26 | 156.6 | 92.6 | 217.5 | 105.8 | 197.3 | +14.2% | -9.3% |

With the pitfall fixed, the BLUP point estimate falls inside the proposed
window for **all 39 groups** (verified by
`data-raw/group_dgamma_bounds_derivation_check.R`). Across all 39 groups,
median `%delta lower = +11.2%` (proposed lower bound is looser than current),
median `%delta upper = +6.0%`, but the three highest-`blup_infl` groups (6,
33, 41) get **narrower** upper bounds than the current BLUP-inflated ones
(-24% to -28%) while most other groups get modestly wider ones (up to +23%).

---

## Part IV -- Empirical validation against an independent classical estimator

Both `sigma2_hat,j` (Part I's quad-penalized point estimate) and the
`blup_infl,j`-driven asymmetry (Part II) rest on the premise that residual
variance genuinely differs across groups -- `lmer` assumes otherwise. This is
directly testable: fit the *same* random-slope structure with `nlme::lme()`,
adding a per-group residual-variance term via `varIdent()`:

```r
fit_hetero <- nlme::lme(
  score_ppvt ~ 1 + distracted_ppvt,
  random  = list(school_id = nlme::pdDiag(~ 1 + distracted_ppvt)),
  weights = nlme::varIdent(form = ~ 1 | school_id),
  data    = dat
)
```

**Heteroscedasticity is statistically supported.** LRT of `lme4` (pooled
`sigma^2`) vs. `nlme` (per-group `sigma^2`), same RE structure: `57.36` on
`38` df, `p = 0.023`.

**`nlme`'s REML per-group variances track both our point estimates closely**
(`cor(sigma2_ols, sigma2_nlme) = 0.985`, `cor(sigma2_hat, sigma2_nlme) =
0.995`), confirming the per-group heterogeneity direction and rough magnitude
independently of this package's own calibration.

**The quad penalty (`sigma2_hat,j > sigma2_ols,j`) is validated exactly where
it is largest.** `nlme`'s estimate falls strictly between `sigma2_ols,j` and
`sigma2_hat,j` for 16 of 39 groups; restricting to the 11 groups with the
largest `blup_infl,j` (equivalently, the largest quad-penalty gap
`sigma2_hat,j - sigma2_ols,j`), it holds for **all 11** (correlation between
`blup_infl,j` and "nlme falls between" = `0.69`; between the quad-penalty gap
and "nlme falls between" = `0.73`). For groups with `blup_infl,j ~ 1` (quad
penalty negligible), `nlme` instead lands slightly below both -- its own
REML shrinkage-toward-pooled effect on noisy small-`n_j` variance estimates,
unrelated to the quad penalty.

### `glmmTMB` as the calibration reference itself, not just an external check

The `nlme::lme()` fit above is an *external* validation: it never feeds back
into `dGamma_list()`'s numbers. `Prior_Setup_GLMM(..., dispformula =
~<group_name>)` instead makes an equivalent `glmmTMB::glmmTMB()` fit --
`glmmTMB(formula, data, dispformula = dispformula, REML = TRUE)` -- the
*source* of every calibration quantity that would otherwise come from the
pooled `lmer` fit: `fixef`, the RE variances `tau^2_k` (hence `sd_tau`), and
the BLUP coefficients used for `blup_infl,j` (Part II). `dispformula = ~1`
is unaffected (still pooled `lmer`/`glmer`); see `fit_ref`, `mer_fit`, and
`calibration_source` on the returned object.

Concretely, this changes *which* model's random-slope shrinkage produces
`sd_tau` and the BLUP residuals in Part I/II above -- previously always
`lmer`'s pooled-`sigma^2` fit, now `glmmTMB`'s per-group-`sigma^2` fit, which
is the model whose assumptions actually match a heteroscedastic
`dGamma_list()`. On the 39-school fixture used throughout this document the
two references are close enough that Part I/II/IV's worked numbers are
essentially unchanged (each group has a reasonably large `n_j`, so the
random-slope variance is well identified by either optimizer). `fit_ref`
also exposes `sigma2_group` -- `glmmTMB`'s own per-group observation-level
dispersion (`predict(fit_ref, type = "disp")`, aggregated by group) -- as a
further diagnostic cross-check against `sigma2_hat,j`, independent of the
`nlme` comparison in Part IV.

**Caution on small or weakly-identified groups.** Unlike `nlme::lme()`,
which is only ever an external check here, `glmmTMB`'s fit directly drives
`dGamma_list()`'s bounds when `dispformula` requests per-group dispersion.
`glmmTMB` and `lmer` solve the same REML problem with different optimizers
and can converge to different local optima for a nearly-singular RE
covariance; on small fixtures (few groups or few observations per group) we
have observed `glmmTMB` collapse a random-slope variance to a boundary value
(`Std.Dev. approx 0`) where `lmer` finds a well-separated-from-zero estimate,
with both fits reporting a converged, positive-definite Hessian (so
`.lmebayes_glmmtmb_convergence_issues()` cannot flag it as a fitting
failure -- it is a legitimate, if inconvenient, local optimum). A
near-zero `tau^2_k` inflates that component's prior precision and can stop
the `dIndependent_Normal_Gamma` sampler from converging. There is currently
no automatic fallback for this; if `Prior_Setup_GLMM(..., dispformula =
~<group_name>)` produces an implausibly small `sd_tau` component, compare
against the pooled (`dispformula = ~1`) fit's `Sigma_ranef` and consider
whether the per-group model is over-parameterized for the available data.

---

## Part V -- Status and open questions (as of the original Parts I-IV; see Part VII for what shipped after)

- **The sampler now uses the ?3.3.4 marginal `rate`, not the ?3.3.5
  `rate_gamma`** (Part I). This is the current default and is not
  configurable via a public argument; `rate_gamma` is retained on
  `ing_prior_measurement_group` and shown alongside `rate` by
  `Prior_Setup_GLMM()`'s dev-only comparison print for anyone who wants
  to inspect the gap, but nothing downstream consumes it. This does not
  change `sigma2_hat` or any other calibration quantity described elsewhere
  in this document -- see Part II's "Alignment with Part I's prior mean"
  (historical) for the one place the two used to interact via the window.
- **Part III was implemented as an opt-in `dGamma_list(ps, disp_center =
  "dispersion2")` argument** (an R-level reproduction of
  `EnvelopeCentering()`'s trace correction, run per group with the
  double-counting pitfall above avoided). `disp_center`,
  `disp_upper_anchor`, and the underlying
  `.lmebayes_group_dispersion2_envelope_centering()`/
  `.lmebayes_group_blup_rss_inflation()` helpers have since been **removed**
  (Part VII); `data-raw/group_dgamma_bounds_derivation_check.R`'s
  independent reproduction of both is kept as a historical companion script
  and is no longer exercised by the shipped `dGamma_list()`.
- **Tension with `lmebayes`'s `inst/GROUP_DISPERSION_HYPERPRIOR.md` ?8.**
  That note argues the "consistent design" keeps **one shared prior** across
  all groups, varying only the truncation window per group, to avoid
  double-counting each group's data in its own prior. `dGamma_list()` instead
  calibrates `shape_ING,j` / `rate,j` **per group** (Part I). The
  default `pwt_measurement = 0.01` keeps `n_prior,j << n_j`, so the resulting
  double-count in the *prior itself* is proportionally small -- but it is not
  zero.

---

## Part VI -- Extension: also integrating out the prior mean `mu_j` (IMPLEMENTED, permanent default)

Part I's marginal Gamma treats `mu_j` as a known constant. This section
derives the change needed to also account for **fixed, known uncertainty
about `mu_j` itself** -- i.e. to integrate out not just `b_j` (already
done, Part I) but a further, fixed budget of uncertainty about where `b_j`'s
prior mean actually is. As of the 2026 migration (Part VII) this is
implemented and is the unconditional default -- see "Status" below for
exactly what was built and how it maps onto the derivation.

### Setup: the `u_j = b_j - mu_j` substitution

Write `b_j = u_j + \bar\mu`, with `u_j \sim N(0, \Sigma_j)` as in Part 0 and
`\bar\mu` a fixed anchor (in place of the current per-group `mu_j`).
Substituting into the likelihood,

```
Y_j = X_j b_j + e_j = X_j u_j + X_j bar_mu + e_j
```

To additionally budget for **fixed, known** uncertainty about `bar_mu`
itself -- call it `Omega_j`, a fixed coefficient-scale covariance matrix,
*not* a further hyper-prior to be integrated (there is no third level of
the hierarchy here; `Omega_j` is given, exactly the way `Sigma_j` already
is) -- the natural definition is a **widened fixed prior covariance** for
`b_j` around `bar_mu`:

```
Sigma_j' = Sigma_j + Omega_j        # both fixed, absolute (coefficient-scale) units
```

### Why `Omega_j` does *not* get rescaled by `sigma2_j`

It is tempting to ask whether `Omega_j` needs to be "dispersion-scaled"
before combining it with `Sigma_0j`. It does not, for the same reason `Sigma_j`
itself is never rescaled by the live `sigma2_j` the ING sampler draws each
sweep: **the ING prior decouples `b_j`'s covariance from the precision
parameter entirely** --

```939:941:vignettes/Chapter-A12.Rmd
The independent Normal���Gamma (ING) prior replaces the conjugate covariance structure
tau^{-1}*Sigma_0 with a fixed coefficient-scale covariance Sigma, ...
```

`Sigma_0j` (Part 0) is "dispersion-free" only in the bookkeeping sense that
it is `Sigma_j` divided by a **fixed, already-known classical plug-in**,
`dispersion_classical,j` (`.lmebayes_ing_prior_measurement_group_glm_inputs()`,
`R/mixed_rmerb_helpers.R`: `dispersion_classical <- rss / (n_j - nvar)`,
computed once from group `j`'s own OLS/WLS residuals, before any
calibration or sampling). It is *not* divided by the unknown `sigma2_j` the
sampler will later draw -- that division exists purely so
`compute_gaussian_prior()` can reuse the classical conjugate
marginal-likelihood algebra as a **calibration device** for `(shape_ING,
rate)`; it makes no claim that `b_j`'s actual runtime prior covariance
scales with dispersion. `Omega_j`, being fixed for the same reason `Sigma_j`
is fixed, belongs on exactly the same footing:

```
Omega_0,j = Omega_j / dispersion_classical,j     # same fixed reference Sigma_0j already divides by
Sigma_0j' = Sigma_0j + Omega_0,j = (Sigma_j + Omega_j) / dispersion_classical,j
```

This is a **single, non-iterative substitution** -- no fixed point, no
proportionality-to-`sigma2_j` question to resolve, because
`dispersion_classical,j` is fixed and known before calibration starts (unlike
Part III's `dispersion2`, which genuinely does need a fixed point, because
it is a *posterior* quantity that depends on the working dispersion).

### Resulting formulas

```
Mi_j'         = (Sigma_0j + Omega_0,j + (X_j'X_j)^{-1})^{-1}         # was Mi_j = (Sigma_0j + (X_j'X_j)^{-1})^{-1}
S_marg,j'     = RSS_ols,j + (bhat_j - bar_mu)' Mi_j' (bhat_j - bar_mu)
sigma2_hat,j' = S_marg,j' / (n_j - p_re)                              # same denominator as Part I
rate,j'       = 0.5 * S_marg,j' * (n_prior,j + p_re - 1) / (n_j - p_re)
shape_ING,j   unchanged                                              # depends only on n_prior_j, p_re
```

### Direction of the effect

`Omega_0,j \succeq 0` implies `Sigma_0j + Omega_0,j + \mathrm{Ginv} \ge
Sigma_0j + \mathrm{Ginv}` in the Loewner order, hence `Mi_j' \le Mi_j`. For
the *same* deviation vector, the quadratic penalty shrinks -- widening the
fixed budget of uncertainty about the anchor attenuates how much weight the
calibration places on any single group's apparent deviation from it. If
`bar_mu` also differs from the current per-group `mu_j` (e.g. sourced
externally rather than from group `j`'s own null-model fit, see below), the
net change in `sigma2_hat,j` is not guaranteed monotone, since the
deviation vector itself shifts too.

### What stays unaffected

- `shape_ING,j`, `RSS_ols,j`, and the `n_j - p_re` denominator (Part I).
- Part II's window construction is purely mechanical in `sigma2_hat,j` -- it
  automatically re-centers at `sigma2_hat,j'`, preserving the "prior mean =
  window center" identity Part II's "Alignment with Part I's prior mean"
  already establishes.
- `rate_gamma,j` (?3.3.5) has no marginalization to extend -- it is a
  fixed-point RSS evaluated at `beta_star,j`, not an integral. The only
  sensible analog is `mu_j \to \bar\mu` inside the Zellner blend
  `beta_star,j = (1 - w_j) bhat_j + w_j \bar\mu`. This sharpens, rather than
  changes, Part I's existing argument for why `rate` (not `rate_gamma`) is
  the theoretically forced choice for the joint ING draw.

### Consistency with the runtime `Sigma`

`compute_gaussian_prior()` returns its `Sigma` argument unchanged when one is
supplied (`.lmebayes_compute_ing_prior_cal_from_sigma()` passes `Sigma =
Sigma_j` as such an override) -- that is the *same* matrix the ING sampler
conditions `b_j | sigma2_j` on at runtime. If the goal is only to calibrate a
better `(shape_ING, rate)` while leaving the runtime `b_j` draw unchanged,
only `Sigma_0` needs widening in the calibration call. If the intent is for
the runtime prior on `b_j` to *also* reflect the wider uncertainty, both
arguments must be widened consistently -- `Sigma = Sigma_j + Omega_j` and
`Sigma_0 = Sigma_0j + Omega_0,j` -- otherwise the calibrated Gamma and the
coefficient covariance the sampler actually draws from would silently
disagree about how much is known about `b_j`.

### Where `Omega_j` and `bar_mu` would have to come from

Two caveats, both already present elsewhere in this document in a different
guise:

1. **Independence from group `j`'s own data.** Under the current default
   (`intercept_source = "null_model"`, Part 0), `mu_j`'s intercept already
   comes from `lm(y ~ 1, data = dat_j)` -- **group `j`'s own data**, the same
   data that produces `bhat_j`. The additive-covariance step above assumes
   `bar_mu`'s uncertainty is independent of `bhat_j`'s sampling error; that
   only holds if `bar_mu`/`Omega_j` are sourced from information that does
   *not* reuse group `j`'s own `y_j` -- e.g. genuinely external, cross-group
   (Block~2) information, the same convention `sd_tau`/`Sigma_ranef` already
   follow in Part 0.
2. **The `O(1/J)` empirical-Bayes caveat.** If `bar_mu`/`Omega_j` are
   themselves estimated by pooling across all `J` groups, group `j`'s own
   contribution to that pooled estimate is a partial, `O(1/J)` double-count
   -- the same caveat Part V already notes for `shape_ING,j`/`rate,j`, and
   discussed at length in `lmebayes`'s `inst/GROUP_DISPERSION_HYPERPRIOR.md`
   ?8. Standard, generally accepted empirical-Bayes practice (Efron, 2010),
   but not a zero double-count.

### Status

**Implemented (permanent default, no opt-out) in
`.lmebayes_calibrate_ing_prior_measurement_group()`
(`R/mixed_rmerb_helpers.R`), wired from `Prior_Setup_GLMM()`.** With one
simplification relative to the general derivation above: `mu_j` itself is
**not** changed (no `mu_j \to \bar\mu` substitution) -- only `Sigma_j \to
Sigma_j + \Omega_j` is applied, i.e. the *same* per-group prior mean from
Part 0 (`.lmebayes_block_formula_prior_mu()`, still group `j`'s own
null-model intercept / zero slopes) with a wider covariance around it. This
is exactly the "Resulting formulas" substitution above with `bar_mu = mu_j`,
so `S_marg,j' = RSS_ols,j + (bhat_j - mu_j)' Mi_j' (bhat_j - mu_j)` and the
rest of the "Resulting formulas" block applies verbatim.

`Omega_j` is **model-derived** (not invented): it is diagonal across RE
components `k` (each `gamma_k` calibrated independently, Part~2 of
`vignettes/Chapter-A12.Rmd`'s Block~2 structure) with

```
Omega_j[k, k] = W_j[[k]] %*% Sigma_fixef_k %*% t(W_j[[k]])
```

where `W_j[[k]]` is group `j`'s own row of `design$W[[k]]` (the Block~2
hyper-design for RE component `k`) and `Sigma_fixef_k =
prior_list[[k]]$Sigma_fixef` is that component's own calibrated
fixed-effect/`gamma_k` covariance (already computed earlier in
`Prior_Setup_GLMM()`, the same object `pfamily_list()` uses to build
`dNormal()`/Block~2 priors) -- i.e. `mu_j[k]`'s uncertainty stand-in is
exactly the model's own Block~2 hyper-regression uncertainty about
`W_j[[k]] %*% gamma_k`, propagated through group `j`'s own design row.

Of the two source caveats above: caveat 1 (independence from group `j`'s
own data) does not apply here, because `mu_j` itself is unchanged -- only
`Sigma_j`'s width changes, and `Omega_j` never enters `mu_j`. Caveat 2 (the
`O(1/J)` empirical-Bayes double-count) does apply: `Sigma_fixef_k` comes
from `vcov(fit_ref)` fit on the pooled data (all `J` groups, including
group `j`'s own observations), so group `j` makes a partial, `O(1/J)`
contribution to its own `Omega_j`. This is the same standard,
generally-accepted empirical-Bayes caveat already noted for
`shape_ING,j`/`rate,j` themselves (Part V), not a new or larger one.
`disp_lower,j`/`disp_upper,j` are no longer built from Part II's
`shape_w`/`rate_w`; see Part VII.

---

## Part VII -- Current implementation (2026 migration): bounds moved upstream

Summarizes what actually ships today, after folding Part VI's `Omega_j` in
as the permanent default and moving `disp_lower`/`disp_upper` out of
`dGamma_list()` and into `Prior_Setup_GLMM()`.

### What changed, in one place

`.lmebayes_calibrate_ing_prior_measurement_group()` (`R/mixed_rmerb_helpers.R`)
is now the single source of truth for per-group `sigma^2_j` calibration.
For each group `j`, in order:

```
Sigma_j          <- Part 0's shrinkage formula (unchanged)          # coefficient-scale, dispersion-free-eligible
Omega_j          <- Part VI's model-derived fold-in (see above)     # diagonal across RE components k
cal              <- compute_gaussian_prior(..., Sigma_0 = (Sigma_j + Omega_j)/dispersion_classical,j,
                                                 Sigma   = Sigma_j + Omega_j, ...)
                     # -> cal$shape_ING, cal$rate, cal$dispersion (== sigma2_hat,j)
disp_lower,j     <- 1 / qgamma(max_disp_perc,     shape = cal$shape_ING, rate = cal$rate)
disp_upper,j     <- 1 / qgamma(1 - max_disp_perc, shape = cal$shape_ING, rate = cal$rate)
```

via `.lmebayes_ing_prior_quantile_window()` -- the same construction already
used for the pooled `ing_prior_measurement` case (Block~1 dispersion when
`dispformula = ~1`) and for `ing_prior`'s `tau^2_k` window (Chapter A12,
Theorem 2). `shape_ING,j`, `rate,j`, `sigma2_hat,j`, `disp_lower,j`,
`disp_upper,j`, `max_disp_perc`, and `omega_j` (the matrix itself, for
diagnostics) are all stored directly on
`ps$ing_prior_measurement_group[[j]]`.

### Why literal quantiles of the actual prior, not a separately-matched window

Parts II/III built a *different* Gamma (`shape_w,j`/`rate_w,j`, mean-matched
at `sigma2_hat,j` via `n_combined,j`) purely to control the window's spread
independently of the actual `Gamma(shape_ING,j, rate,j)` fed to the sampler.
That decoupling is no longer needed once the marginal itself already
integrates out both `b_j` and `gamma` (Part VI) -- the actual sampling prior
is now the right object to truncate directly. This guarantees `sigma2_hat,j`
(the prior's own implied mean-ish point estimate) always lies inside the
window by construction (`shape_ING,j > 1` is required of any valid
calibration for the mean to exist, and the window brackets the whole central
mass around it), with no separate spread parameter (`n_combined,j`) to keep
in sync.

**Amended in Part VIII below**: the window's *shape/rate* are not literally
`(shape_ING,j, rate,j)` -- they are `(shape_ING,j, rate,j)` inflated by the
group's own `n_j` observations, i.e. the *posterior* (not prior) spread.
`shape_ING,j`/`rate,j` themselves, the values actually fed to the sampler as
the Gamma prior, are unaffected either way.

### `dGamma_list()` after the migration

`dGamma_list.Prior_Setup_GLMM()` (`R/dGamma_list_Prior_Setup_GLMM.R`)
is now a thin consumer: for each group it reads
`g$shape_ING`/`g$rate`/`g$disp_lower`/`g$disp_upper` off
`object$ing_prior_measurement_group` and calls `dGamma()` directly. Its
override argument is `max_disp_perc_measurement` (`NULL` by default), a
scalar or a named/positional length-`J` vector (one value per group level,
via the same `.lmebayes_expand_scalar_or_vector()` resolver used inside
`Prior_Setup_GLMM()`). `NULL` (the default) reuses every group's own
stored `g$max_disp_perc` bounds verbatim; supplying a scalar or vector
recomputes fresh quantiles of the same `Gamma(shape_ING,j, rate,j)` for
every group whose resolved override differs from its own stored
`g$max_disp_perc` (a **per-group** comparison, not a single top-level
`object$max_disp_perc` field -- `Prior_Setup_GLMM()` itself now accepts
`max_disp_perc_measurement` as a scalar or per-group vector too, so
different groups can already carry different stored windows before
`dGamma_list()` is even called). The `disp_center`, `disp_upper_anchor`, and
`n_rss_iter` arguments (Parts II/III) have been **removed** from the
function's formal signature -- callers who still pass them get a silent
no-op (absorbed by `...`), not an error, since `dGamma_list()` keeps a `...`
parameter. The `warn_asymmetric`/`print_asymmetric`/`asymmetric_R_lo`/
`asymmetric_R_hi` diagnostic (which flagged windows where the BLUP-anchored
upper tail and the OLS-anchored lower tail crossed very different
percentiles of the underlying Gamma) is also gone -- there is no longer a
separate upper/lower anchor to be asymmetric between. The
`"window_diagnostics"` attribute is still attached (now just `group`, `n_j`,
`sigma2_hat`, `shape_ING`, `rate`, `disp_lower`, `disp_upper` -- no
`asymmetric_window`/`blup_infl`/`R_lo`/`R_hi` columns).

### Verified equivalence

`data-raw/_scratch_compare_partvi_handrolled.R` confirms
`Prior_Setup_GLMM()`'s default `shape_ING,j`/`rate,j`/`sigma2_hat,j` now
match, to floating-point precision, the demo `Ex_13b`'s pre-migration
hand-rolled reproduction of this same Part VI derivation (which called
`.lmebayes_ing_prior_measurement_group_glm_inputs()` and
`.lmebayes_compute_ing_prior_cal_from_sigma()` directly). `data-raw/
_scratch_validate_disp_bounds_migration.R` additionally checks: every
group's `sigma2_hat,j` lands inside its own `[disp_lower,j, disp_upper,j]`;
`dGamma_list(ps)`'s returned `pfamily` objects carry the stored
shape/rate/bounds verbatim; a different `max_disp_perc_measurement` triggers
a recompute (narrower window at a looser `max_disp_perc_measurement`, as
expected); and the now-removed `disp_center`/`disp_upper_anchor`/`n_rss_iter`/
`warn_asymmetric` arguments are silently accepted without error.

---

## Part VIII -- Correction: the window uses the posterior shape (`+n_j/2`), not the prior shape

Part VII's initial migration truncated the *prior* `Gamma(shape_ING,j,
rate,j)` directly:

```
disp_lower,j  <- 1 / qgamma(max_disp_perc,     shape = shape_ING,j, rate = rate,j)
disp_upper,j  <- 1 / qgamma(1 - max_disp_perc, shape = shape_ING,j, rate = rate,j)
```

That is too wide relative to what the sampler itself actually needs the
window for. `src/EnvelopeDispersionBuild.cpp` builds its own truncation
window internally, as a *fallback*, whenever `disp_lower`/`disp_upper` are
not supplied -- and it does so from the **posterior**-adjusted shape/rate,
not the prior's:

```
shape2  <- Shape + n_w / 2   # Shape = the input (prior) shape_ING,j; n_w = n_j
rate3   <- rate  + RSS_post / 2
```

(see `inst/BLOCK_ING_RINDEPNORMALGAMMA_REG.md`, "Lemma 2 -- centering"). This
`+n_j/2` term is exactly the historical Part II `shape_w,j` construction's
own `n_combined,j = n_prior,j + n_j` (vs. Part VII's `n_prior,j`-only
`shape_ING,j`) -- Part II was right about *what* to inflate by, Part VII's
first migration pass dropped that inflation when it unified the window with
the sampling prior.

### The fix, scoped to the window only

`.lmebayes_calibrate_ing_prior_measurement_group()`
(`R/mixed_rmerb_helpers.R`) now computes the window from a **separate,
posterior-shape** Gamma, while leaving `shape_ING,j`/`rate,j` -- the values
actually stored on `ps$ing_prior_measurement_group[[j]]$shape_ING`/`$rate`
and fed to the sampler as the Gamma prior -- untouched:

```
shape_post,j  <- shape_ING,j + n_j / 2
rate_post,j   <- sigma2_hat,j * (shape_post,j - 1)   # mean-matched at the SAME sigma2_hat,j
disp_lower,j  <- 1 / qgamma(max_disp_perc,     shape = shape_post,j, rate = rate_post,j)
disp_upper,j  <- 1 / qgamma(1 - max_disp_perc, shape = shape_post,j, rate = rate_post,j)
```

`rate_post,j` is mean-matched (`sigma2_hat,j = rate_post,j /
(shape_post,j - 1)`) rather than reusing `rate,j` unchanged, for the same
reason Part II mean-matched `rate_w,j`: inflating only the shape while
keeping the rate fixed would shift the window's *center* away from
`sigma2_hat,j`, not just tighten its spread. Algebraically,
`(shape_post,j, rate_post,j)` reproduces Part II's `(shape_w,j, rate_w,j)`
exactly, since both equal `((n_prior,j + n_j + 1)/2 + p_re/2,
sigma2_hat,j * (n_prior,j + n_j + p_re - 1)/2)` when `n_combined,j =
n_prior,j + n_j` -- confirmed to floating-point precision by
`data-raw/_scratch_check_window_njfix.R`.

`dGamma_list.Prior_Setup_GLMM()`'s own recompute path (triggered when a
caller overrides `max_disp_perc_measurement` for a group away from its
stored value) applies the identical `shape_ING + n_j/2` / mean-matched-rate
construction, so overriding the window's `max_disp_perc_measurement` at
`dGamma_list()` time stays consistent with what `Prior_Setup_GLMM()`
would have computed had it been called with that same value in the first
place.

### What does *not* change

- `shape_ING,j`/`rate,j` (the actual Gamma prior handed to the sampler) --
  unchanged; this correction only widens/narrows the *window*.
- The pooled Block~1 **window** still uses prior quantiles of the pooled
  shape/rate (no `+n/2` inflation); that gap is unchanged here.
- **Centers (updated separately):** pooled `dispformula = ~1` now aggregates
  per-group \(S_{\mathrm{marg},j}\) (Part I / Part VI) into
  \(\hat\sigma^2_{\mathrm{pool}}=(\sum_j S_{\mathrm{marg},j})/(n-J p_{\mathrm{re}})\),
  and Block~2 `pop.ing_prior` centers \(\tau^2_k\) on the same A12 §3.3.4
  hyper-regression marginal (see `Prior_Setup_GLMM` /
  `inst/TAU2_ING_FORMULAS.md`). Only the per-group
  `ing_prior_measurement_group` **window** construction is the subject of
  this Part VII note.

---

## References

- `vignettes/Chapter-A12.Rmd` -- the general Normal--Gamma / ING / `dGamma()`
  marginal-Gamma derivation this document specializes (Part I, Part VI).
- `inst/TAU2_ING_FORMULAS.md` -- the same conjugate Gamma--Normal update
  mechanics applied to `tau^2_k` instead of measurement dispersion.
- `lmebayes`'s `inst/ING_TRUNCATION_WINDOW.md` -- the truncation-window
  rationale this document's Part II specializes to per-group bounds.
- `lmebayes`'s `inst/GROUP_DISPERSION_HYPERPRIOR.md` -- precursor exploration
  of group-level dispersion heterogeneity and the double-counting hazard
  (?8), on the same `big_word_club` fixture; also the source of the
  `O(1/J)` empirical-Bayes caveat cited in Part VI.
- Pinheiro, J.C. and Bates, D.M. (2000). *Mixed-Effects Models in S and
  S-PLUS.* Springer. -- `varIdent()` heteroscedastic variance structures used
  in Part IV.
- Efron, B. (2010). *Large-Scale Inference: Empirical Bayes Methods for
  Estimation, Testing, and Prediction.* Cambridge University Press. --
  justification for the `O(1/J)` empirical-Bayes double-count in Part VI /
  `GROUP_DISPERSION_HYPERPRIOR.md` ?8.

## Runnable check

```r
## Historical (Parts I-IV's pre-migration window construction; not exercised
## by the shipped dGamma_list() -- see Part VII):
Rscript data-raw/group_dgamma_bounds_derivation_check.R

## Current implementation (Part VII):
Rscript data-raw/_scratch_compare_partvi_handrolled.R
Rscript data-raw/_scratch_validate_disp_bounds_migration.R
```

Reproduces the Part II/III worked-numbers table (all 39 groups) and the Part
IV `nlme` validation statistics from the live `big_word_club` fixture.
