# Guide: Adding a new pfamily / link combination to glmbayes

This document is the canonical checklist for extending `glmbayes` with a new
conjugate or semi-conjugate prior family (`pfamily`) and its associated
simulation function.  Follow every step in order; each step depends on the
previous one being correct.

---

## Background: the pfamily contract

Every `pfamily` object is a list with exactly these named slots:

| Slot | Type | Role |
|---|---|---|
| `pfamily` | character | Name tag (e.g. `"dBeta"`) |
| `prior_list` | list | All prior hyperparameters, plus `mu`, `Sigma` surrogates |
| `okfamilies` | character vector | GLM family names this prior supports |
| `plinks` | function(family) | Returns allowed link names for a given family; NULL if unsupported |
| `simfun` | function | Low-level sampler; called by `rglmb()` dispatch |
| `call` | call | Matched call from the constructor |

---

## Step 1 — Derive the conjugate posterior update

Before writing any code:

1. Write down the likelihood kernel `p(y | θ)` for the new family/link.
2. Write down the prior kernel `p(θ)`.
3. Derive the posterior `p(θ | y) ∝ p(y | θ) p(θ)`.
4. Identify which parameters are **fixed** (known constants that enable
   closed-form conjugacy) vs **estimated** (drawn each iteration).
5. Write down the closed-form posterior hyperparameter update rules.
6. Identify any constraints (e.g. intercept-only, non-negative response,
   response in [0,1]).

### Implemented conjugate pairs for reference

| Likelihood | Link | Prior on | Fixed | Update |
|---|---|---|---|---|
| Poisson | identity | rate λ | — | Gamma(α + Σy, β + n) |
| Gamma | identity | rate β | shape k (`lik_shape`) | Gamma(α + nk, β + Σy) |
| Binomial | identity | prob θ | — | Beta(α + Σ successes, β + Σ failures) |
| Gaussian | identity | precision τ | mean μ | Gamma(α + n/2, β + RSS/2) |

---

## Step 2 — `R/pfamily.R`: write the constructor

Add a new exported function at the bottom of `pfamily.R`.

**Required elements:**
- Validate all hyperparameter inputs (type, length, positivity/range).
- Build the Normal-style surrogate `mu` and `Sigma` matrices from the prior
  moments (used by `glmb()` pre-simulation and `summary()` reporting).
- Set `okfamilies` (character vector of supported family names).
- Define `plinks(family)` closure that returns allowed links per family.
- Set `simfun` to the new simulation function (written in Step 3).
- Build `prior_list` storing all hyperparameters plus `mu` and `Sigma`.
- Set `attr(prior_list, "Prior Type")` and `attr(outlist, "Prior Type")`.
- Set `class(outlist) <- "pfamily"` and `outlist$call <- match.call()`.

**Roxygen tags needed:**
```r
#' @export
#' @rdname pfamily
#' @order N        # next available integer
```

All `@param` entries must be documented — `R CMD check` will flag any
undocumented argument.

---

## Step 3 — `R/simfunction.R`: write `rXxx_reg()`

Model the new function on `rGamma_Conjugate_reg()`.  Key structural requirements:

### 3a. Function signature
```r
rXxx_reg <- function(n, y, x, prior_list, offset = NULL,
                     weights = 1, family = gaussian(),
                     Gridtype = 2, n_envopt = NULL,
                     use_parallel = TRUE, use_opencl = FALSE,
                     verbose = FALSE, progbar = FALSE)
```

### 3b. Standard preamble (copy from `rGamma_Conjugate_reg`)
- `call <- match.call()`
- Rename `weights → wt`, `offset → alpha`
- Validate and normalize `y`, `x`, `wt`, `alpha`
- Extract prior hyperparameters from `prior_list`

### 3c. Family / link guard
```r
okfamilies <- c(...)
# Check family$family %in% okfamilies and family$link %in% oklinks
```

### 3d. Scalar design guard (if conjugate, intercept-only)
Call `.check_gamma_conjugate_scalar_design()` **or** add an inline guard if
the error messages specific to your prior name matter.  The guard enforces:
- `ncol(x) == 1L` (intercept-only)
- `prior_list$beta` is 1×1
- All offsets are zero
- Weights are constant (unless the conjugate update handles heterogeneous
  weights explicitly — verify this from the math)

### 3e. Response validation
Enforce any response constraints (e.g. y ∈ [0,1] for binomial, y > 0 for
Gamma/exponential).

### 3f. Conjugate draw
```r
# Compute posterior hyperparameters
# ...
# Draw n_draw IID samples
coef_out <- matrix(rXxx(n_draw, ...), nrow = n_draw, ncol = 1L)
disp_out <- rep(dispersion_value, n_draw)
draws_out <- matrix(1L, nrow = n_draw, ncol = 1L)  # each draw accepted in 1 step
```

### 3g. Posterior mode
Provide an analytic posterior mode where it exists; fall back to the posterior
mean otherwise.

### 3h. Assemble output list (must match `rGamma_Conjugate_reg` structure exactly)
```r
outlist <- list(
  coefficients  = coef_out,
  coef.mode     = coef_mode_out,
  dispersion    = disp_out,
  Prior         = list(mean = ..., Precision = ..., ...),
  prior.weights = wt,
  y             = y,
  x             = x,
  famfunc       = glmbfamfunc(family),  # requires Step 4
  iters         = draws_out,
  Envelope      = NULL
)
# then append: family, offset2, formula, model, data, call
class(outlist) <- c(outlist$class, "rglmb", "glmb", "glm", "lm")
```

### 3i. Roxygen tag
```r
#' @family simfuncs
```

---

## Step 4 — `R/simulationpipeline.R`: add a `glmbfamfunc` branch

`glmbfamfunc` must have a branch for every family/link that a simulation
function uses — without it, DIC, `logLik()`, and `directional_tail()` will
fail with *"object 'f1' not found"*.

Add the new branch just before `out = list(f1=f1, ...)` at the bottom of the
function.  Each branch must define **all five** closures:

| Closure | Role | Signature |
|---|---|---|
| `f1` | Negative log-likelihood | `f1(b, y, x, alpha=0, wt=1)` |
| `f2` | Negative log-posterior | `f2(b, y, x, mu, P, alpha=0, wt=1)` |
| `f3` | Gradient of f2 w.r.t. b | `f3(b, y, x, mu, P, alpha=0, wt=1)` |
| `f4` | Deviance (2 × NLL contrast vs saturated) | `f4(b, y, x, alpha=0, wt=1, dispersion=1)` |
| `f7` | Fisher information matrix | `f7(b, y, x, mu, P, alpha=0, wt=1)` |

**Common derivations:**
- f3 = `−∂L/∂b + P(b−μ)` where `∂L/∂b = X' diag(wt) ∂ℓ/∂η`
- f7 = `X' diag(wt × w_i) X` where `w_i = (∂μ/∂η)² / V(μ_i)` (GLM weight)
- f4 = `2 × f1(b,y,x,α,wt/φ) + 2 × Σ log p(y_i | ŷ_i=y_i)`  (saturated NLL)

If `glmbfamfunc` needs an extra parameter (e.g. `lik_shape`), add it as an
optional argument with a sensible default.

Update the **documentation table** in the roxygen header listing all
implemented family/link combinations.

---

## Step 5 — `R/prior.R`: add `conj_Xxx` to `Prior_Setup()`

If the new prior supports calibration via `Prior_Setup()`:

1. Add a `conj_Xxx <- NULL` initialization near the Poisson block.
2. Determine the conditions (family, link, `ncol(x) == 1L`, scalar `pwt`).
3. Compute the effective prior sample size:
   `n_prior = (pwt / (1 - pwt)) * n_eff`
4. Derive prior hyperparameters from `n_prior` and the data summary statistic.
5. Build `conj_Xxx <- list(...)` with at minimum `beta`, plus the
   hyperparameters needed by the `pfamily` constructor.
6. Include `conj_Xxx` in the returned `prior_list`.
7. Add a `conj_Xxx` branch to `print.PriorSetup()` (search for the
   `conj_poisson` print block to find the right location).
8. Update the `@return` roxygen docs for `Prior_Setup` to document
   `conj_Xxx`.

---

## Step 6 — Documentation

- Run `devtools::document()` after writing all roxygen blocks.
- Run `spelling::spell_check_package()` and add any new legitimate words to
  `inst/WORDLIST`.
- Verify `R CMD check` produces no warnings about undocumented arguments.

---

## Step 7 — Tests (`tests/testthat/`)

Add at least one test file covering:

1. Constructor validation (bad inputs produce informative errors).
2. Conjugate draw mean ≈ analytic posterior mean (within Monte Carlo noise).
3. Conjugate draw SD ≈ analytic posterior SD.
4. `summary()` runs without error.
5. `logLik()` and `DIC_Info()` run without error (validates `glmbfamfunc`).
6. `Prior_Setup()` produces non-NULL `conj_Xxx` under the right conditions.

---

## Step 8 — Example / vignette

- Add a minimal example to `inst/examples/Ex_dXxx.R` (for `@example` in
  roxygen).
- Add a subsection to `vignettes/Chapter-02.Rmd` following the pattern of
  the Beta–Binomial (§2) or Gamma–Gamma (§4) sections.

---

## Quick cross-reference: where each prior name appears

When you name a new prior `dXxx`, search for every occurrence of
`dGamma_Conjugate` (or the nearest analogue) in the files below and create
the corresponding `dXxx` entry:

| File | What to add |
|---|---|
| `R/pfamily.R` | Constructor `dXxx()` |
| `R/simfunction.R` | Simulation function `rXxx_reg()` |
| `R/simulationpipeline.R` | `glmbfamfunc` branch + doc table row |
| `R/prior.R` | `conj_Xxx` block in `Prior_Setup()` + print method |
| `inst/WORDLIST` | Any new technical terms flagged by spell check |
| `tests/testthat/` | `test-dXxx.R` |
| `inst/examples/` | `Ex_dXxx.R` |
| `vignettes/Chapter-02.Rmd` | New section |
