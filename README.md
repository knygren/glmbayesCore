# glmbayesCore

![GitHub release (latest by date)](https://img.shields.io/github/v/release/knygren/glmbayesCore?label=version)
![License: GPL-2](https://img.shields.io/badge/license-GPL--2-blue.svg)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/knygren/glmbayesCore/R-CMD-check.yaml?label=R%20CMD%20Check)

**glmbayesCore** is the compiled sampling engine that powers the glmbayes
ecosystem. It holds the C++/OpenCL envelope samplers, the family-function
infrastructure, and the R-level prior and simulation interfaces that
downstream packages depend on. End users should install
[glmbayes](https://github.com/knygren/glmbayes) rather than this package
directly.

The relationship to the broader ecosystem parallels how `StanHeaders` /
`rstan` serve as the compiled backbone for `rstanarm`: **glmbayesCore** is
the infrastructure layer; **glmbayes** and the in-development **lmebayes**
are the user-facing packages built on top of it.

**Current staging note.** This tree currently ships the iid GLM/LM envelope
engine used by **glmbayes**. Mixed-model (LMM/GLMM / two-block) engines are
still part of the long-term **glmbayesCore** API and are under active
development in the temporary
[lmebayesCore](https://github.com/knygren/lmebayesCore) fork (consumed by
[lmebayes](https://github.com/knygren/lmebayes)). Think of **lmebayesCore**
as a development holding package: features return here gradually once the
iid backend is stable (CRAN / **glmbayes** re-import path).

---

## Package Ecosystem

**Target architecture** (after mixed-model reintegration):

```
                ┌─────────────────────────────────────────┐
                │           End-user packages             │
                │   glmbayes  ·  lmebayes  ·  (others)    │
                └──────────────────┬──────────────────────┘
                                   │ Imports / LinkingTo
                ┌──────────────────▼──────────────────────┐
                │              glmbayesCore               │
                │  iid GLM/LM · LMM/GLMM · OpenCL         │
                │  pfamily · simfunctions · rglmb/rlmb    │
                └──────────────────┬──────────────────────┘
                                   │ Imports
                ┌──────────────────▼──────────────────────┐
                │   opencltools  ·  nmathopencl            │
                │   Rcpp · RcppArmadillo · RcppParallel   │
                └─────────────────────────────────────────┘
```

**Temporary staging** (today): **glmbayes** → **glmbayesCore** (iid);
**lmebayes** → **lmebayesCore** (full fork including mixed-model stack).
The fork collapses back into **glmbayesCore** as features are merged.

**glmbayes** adds the formula interface (`glmb()`, `lmb()`), MCMC diagnostics,
and the full suite of S3 methods that mirror base-R's `lm()` / `glm()`.

**lmebayes** (in development) extends the engine to linear / generalized
linear mixed-effects models (`lmerb()`, `glmerb()`).

---

## What Is Inside glmbayesCore

### C++ sampling engine (`src/`)

The core is organized under the `glmbayes::` namespace:

| Sub-namespace | Key files | Role |
|---|---|---|
| `glmbayes::fam` | `famfuncs.h`, `famfuncs_*.cpp` | Negative log-posterior (`f2`) and gradient (`f3`) for gaussian, poisson, binomial, Gamma |
| `glmbayes::env` | `EnvelopeBuild*.cpp`, `EnvelopeEval.cpp`, `EnvelopeSort.cpp`, `EnvelopeSize.cpp`, `Set_Grid.cpp`, `Set_LogP.cpp` | Piecewise-exponential envelope construction (Nygren & Nygren, 2006) |
| `glmbayes::sim` | `rNormalGLM.cpp`, `rIndepNormalGammaReg.cpp`, `rNormalGammaReg.cpp`, `rNormalReg.cpp`, `rGammaGamma.cpp`, `rGammaGaussian.cpp` | Posterior samplers |
| `glmbayes::rng` | `rng_utils.cpp` | Thread-safe RNG wrappers for parallel sampling |
| `glmbayes::progress` | `progress_utils.cpp` | Optional progress bar support |

Export wrappers in `export_wrappers.cpp` and `kernel_wrappers.cpp` expose
selected entry points to R via Rcpp.

### OpenCL kernels (`inst/cl/`)

For systems with an OpenCL-capable device, envelope construction can be
offloaded to the GPU. The `inst/cl/` tree contains family/link `f2`/`f3`
kernels, an OpenCL port of R Mathlib probability functions, and shim headers.
Kernel loading for exploration uses **opencltools**; runtime GPU assembly uses
`kernel_loader.cpp` and `kernel_runners.cpp`.

### R-level infrastructure (`R/`)

| File | Role |
|---|---|
| `pfamily.R` | Prior-family constructors (`dNormal`, `dNormal_Gamma`, `dIndependent_Normal_Gamma`, `dGamma`, `dBeta`) and the `pfamily()` generic |
| `prior.R` | `Prior_Setup()`, `Prior_Check()`, and helper utilities for default hyperparameters |
| `simfunction.R` | Low-level simulation functions (`rNormal_reg`, `rNormalGamma_reg`, `rindepNormalGamma_reg`, `rGamma_reg`, …) and the `simfunction()` introspection generic |
| `simulationpipeline.R` | `glmbfamfunc()`, envelope R exports, standardized samplers |
| `rglmb.R` / `rlmb.R` | Matrix-input samplers — the primary R-level interface for **glmbayes** |
| `envelopeorchestrator.R` | R orchestration of multi-step envelope building and optional GPU dispatch |
| `compute_gaussian_prior.R` | Gaussian-specific prior calibration utilities |

---

## Architecture: How pfamilies Route to Simulation Functions

A `pfamily` object is a self-contained prior specification. Every constructor
bundles the hyperparameters into a `prior_list` **and** embeds a `simfun`
function pointer. When `rglmb()` draws samples, it calls
`pfamily$simfun(y, x, prior_list, family, ...)` — there is no internal
`switch` on prior type.

```
rglmb(y, x, pfamily = dNormal(...), family = poisson())
          │
          └─► pfamily$simfun  ──►  rNormal_reg()
                                       │
                              family == gaussian?
                              ├── Yes ──► conjugate multivariate normal draw
                              └── No  ──► envelope sampling (Nygren & Nygren, 2006)
                                              │
                                              └──► rNormalGLM (C++)
```

| pfamily constructor | Embedded `simfun` | Posterior path |
|---|---|---|
| `dNormal()` | `rNormal_reg()` | Conjugate MVN draw (Gaussian); subgradient envelope sampling (other families) |
| `dNormal_Gamma()` | `rNormalGamma_reg()` | Conjugate Normal-Gamma draw (Gaussian only) |
| `dIndependent_Normal_Gamma()` | `rindepNormalGamma_reg()` | Joint coefficient + dispersion envelope (Gaussian; non-conjugate) |
| `dGamma(Inv_Dispersion = TRUE)` | `rGamma_reg()` | Gamma prior on inverse dispersion |
| `dGamma(Inv_Dispersion = FALSE)` | `rGamma_Conjugate_reg()` | Conjugate Gamma–Poisson or Gamma–Gamma (intercept-only, identity link) |
| `dBeta()` | `rBeta_reg()` | Conjugate Beta–Binomial (intercept-only, identity link) |

`Prior_Setup()` fits an auxiliary GLM and returns calibrated hyperparameters
on the same scale as the design matrix.

---

## Architecture: How Simulation Functions Route to C++ Samplers

### `rNormal_reg()`

```
rNormal_reg(y, x, prior_list, family, ...)
       │
  family$family == "gaussian"?
  ├── Yes ──► direct MVN draw via backsolve / Cholesky
  └── No  ──► EnvelopeOrchestrator (R)
                   ├── EnvelopeBuild (C++)
                   └── rNormalGLM (C++)   [accept-reject; optional OpenCL envelope]
```

### `rindepNormalGamma_reg()`

```
rindepNormalGamma_reg(y, x, prior_list, ...)
       │
       └──► rIndepNormalGammaReg (C++)
                   ├── EnvelopeBuild_Ind_Normal_Gamma per dispersion grid point
                   └── joint accept-reject over (beta, dispersion)
```

### `rGamma_reg()`

```
rGamma_reg(y, x, prior_list, family, ...)
       │
  family$family == "gaussian"?
  ├── Yes ──► rGammaGaussian (C++)
  └── No  ──► rGammaGamma (C++)
```

---

## Architecture: How `rglmb()` Orchestrates a Draw

`rglmb()` validates the `family × pfamily` combination and delegates sampling
to the `simfun` embedded in the `pfamily` object. In **glmbayes**, `glmb()` and
`lmb()` wrap `rglmb()` / `rlmb()` with formula parsing.

```
rglmb(y, x, family = poisson(), pfamily = dNormal(mu, Sigma), n = 1000)
  │
  ├─ 1. Resolve family
  ├─ 2. Unpack pfamily (okfamilies, plinks, prior_list, simfun)
  ├─ 3. Validate combination
  ├─ 4. outlist ← simfun(...)
  └─ 5. Post-process → class c("rglmb", "glmb", "glm", "lm")
```

Adding a new prior family requires a new pfamily constructor and simulation
function — not changes to `rglmb()` itself.

---

## Function overview

Symbols below are exported from **glmbayesCore** today (iid path). End users
typically load **glmbayes** (or **lmebayes** for mixed models). Mixed-model
exports temporarily ship from **lmebayesCore** and will return here.

**Maintainers:** current slim inventories live in
[inst/R_FUNCTION_INVENTORY.md](inst/R_FUNCTION_INVENTORY.md)
([exports](inst/R_EXPORTED_AND_DOCUMENTED.md),
[Core-only by type](inst/R_CORE_ONLY_EXPORTS.md),
[reachability](inst/R_EXPORT_REACHABILITY.md),
[internal helpers](inst/R_INTERNAL_HELPERS.md)).

### Shared with **glmbayes** (iid GLM / LM)

#### Retain as **glmbayes** re-exports

| Function | Role |
|----------|------|
| `Prior_Setup()`, `Prior_Check()` | Default prior calibration and prior predictive checks |
| `pfamily()`, `dNormal()`, `dNormal_Gamma()`, `dIndependent_Normal_Gamma()`, `dGamma()`, `dBeta()` | Prior-family constructors |
| `multi_prior_setup()`, `multi_rlmb()` | Multi-response Gaussian prior setup / LM sampler |
| `rglmb()`, `rlmb()` | Matrix-level Bayesian GLM / LM samplers |
| `diagnose_glmbayes()` | OpenCL / GPU diagnostic report |

#### Phase out of **glmbayes** (stay in **glmbayesCore**)

| Function | Role |
|----------|------|
| `compute_gaussian_prior()` | Internal Gaussian calibration used inside `Prior_Setup()` |
| `simfunction()`, `glmbfamfunc()` | Simulation registry and GLM family pipeline helpers |
| `rNormal_reg()`, `rNormalGamma_reg()`, `rindepNormalGamma_reg()`, `rGamma_reg()`, `rBeta_reg()`, … | Low-level `simfunction` samplers |
| `rNormalGLM_std()`, `rIndepNormalGammaReg_std()`, `glmb.wfit()`, `glmb_Standardize_Model()` | Standardized envelope path and fitter hooks |
| `EnvelopeBuild()`, `EnvelopeOrchestrator()`, `EnvelopeSize()`, … | Accept–reject envelope machinery |
| `pnorm_ct()`, `rnorm_ct()`, `pinvgamma_ct()`, `rgamma_ct()`, … | Truncated-distribution C++ callbacks |

### Planned mixed-model API (temporary **lmebayesCore**; returns here)

These are part of the long-term **glmbayesCore** surface for **lmebayes**.
They are not exported from this tree today.

| Area | Examples |
|------|----------|
| Setup | `model_setup()`, `Prior_Setup_lmebayes()`, `pfamily_list()` |
| Matrix drivers | `rlmerb()`, `rglmerb()` |
| Two-block / sweep | `rGLMM_reg*`, `rLMM_reg*`, `rGLMM_sweep()`, `two_block_*`, `plot_sweep_history_diag()` |
| Block helpers | `build_mu_all()`, ICM helpers, `block_rNormalReg()` / `block_rNormalGLM()` |

Typical **lmebayes** workflow (via **lmebayesCore** for now):
`model_setup()` → `Prior_Setup_lmebayes()` → `pfamily_list(ps)` →
`lmerb()` / `glmerb()`.

---

## Developer Interface Levels

### Level 1 — C++ (via `LinkingTo`)

```cpp
#include "glmbayesCore/famfuncs.h"
#include "glmbayesCore/Envelopefuncs.h"
#include "glmbayesCore/simfuncs.h"
#include "glmbayesCore/R_interface.h"
```

### Level 2 — R simulation functions

```r
library(glmbayesCore)
fit <- rindepNormalGamma_reg(
  y = y, x = X, n = 2000,
  prior_list = dIndependent_Normal_Gamma(mu, Sigma, shape, rate)$prior_list,
  family = gaussian()
)
```

### Level 3 — `rglmb()` / `rlmb()` with pfamily objects

```r
ps  <- Prior_Setup(y, X, family = poisson())
fit <- rglmb(y = y, x = X, n = 1000,
             pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma),
             family  = poisson())
```

---

## Installation

**GitHub / R-Universe** (recommended for developers):

```r
install.packages("glmbayesCore",
                 repos = c("https://cloud.r-project.org",
                           "https://knygren.r-universe.dev"))
```

**From source** (required for OpenCL GPU support):

```r
install.packages("glmbayesCore", type = "source",
                 repos = "https://knygren.r-universe.dev")
```

See [Chapter 16 — Large models: GPU acceleration using OpenCL](https://knygren.r-universe.dev/articles/glmbayes/Chapter-16.html)
for system-level setup instructions.

**Dependencies that must be installed first:**

```r
install.packages(c("Rcpp", "RcppArmadillo", "RcppParallel", "MASS", "Rdpack"))
install.packages(c("opencltools", "nmathopencl"),
                 repos = "https://knygren.r-universe.dev")
```

---

## Extending glmbayesCore

### Adding a new pfamily

See `inst/ADDING_PFAMILY.md`. In summary:

1. Write a constructor in `pfamily.R` that builds `prior_list` and sets `simfun`.
2. Implement or reuse a simulation function in `simfunction.R`.
3. If a new C++ sampler is needed, add it under `src/`, register via
   `Rcpp::compileAttributes()`, and expose it through `export_wrappers.cpp`.
4. For GPU support, add the corresponding `f2`/`f3` OpenCL kernel under
   `inst/cl/src/` and register it in `kernel_loader.cpp`.

### Block Gibbs / mixed-model engines

The orchestrator pattern is intentionally generic: validate a model
specification, unpack a routing object, call the embedded `simfun`,
post-process. Mixed-effects drivers in **lmebayes** (`lmerb()`, `glmerb()`)
build on two-block Gibbs engines that will return to this package. While
those engines are staged in **lmebayesCore**, architecture notes
(ergodicity, `rGLMM_sweep` / Block~1–Block~2 call chains, C++ migration
plans) live there and will move back with the code.

### `R/` symbol inventory

[inst/R_FUNCTION_INVENTORY.md](inst/R_FUNCTION_INVENTORY.md)

---

## Key References

- Nygren, K.N. and Nygren, L.M. (2006). Likelihood Subgradient Densities. *Journal of the American Statistical Association*, 101(475), 1144–1156. — The accept-reject envelope method at the heart of the non-Gaussian samplers.
- Lindley, D.V. and Smith, A.F.M. (1972). Bayes estimates for the linear model. *Journal of the Royal Statistical Society B*, 34, 1–41. — Conjugate Normal-Gamma foundations.
- Gelman, A. et al. (2013). *Bayesian Data Analysis*, 3rd ed. — Reference for prior specifications and dispersion modeling.

A complete bibliography is in `inst/REFERENCES.bib`.

---

## Future plans

- **Reintegrate mixed-model stack from temporary lmebayesCore:** Gradually
  merge LMM/GLMM setup (`model_setup`, `Prior_Setup_lmebayes`,
  `pfamily_list`), matrix drivers (`rlmerb` / `rglmerb`), and two-block /
  block-ING engines back into **glmbayesCore**, then point **lmebayes** at
  this package again and retire **lmebayesCore**.
- **Sweep-outer drivers and `sweep_history` on all two-block paths:**
  Mixed-model sampling should use a **sweep-outer** loop (all chains
  complete inner sweep `m`, then `m+1`, …) on every route, for consistency
  with `rGLMM_sweep()` / `rGLMM_reg_*`. Each stored draw should attach
  **`sweep_history`** (class `two_block_sweep_history`) so `print()` and
  `plot_sweep_history_diag()` can diagnose inner-Gibbs convergence.
  Sweep-outer R drivers and ING pilot/main paths already capture history in
  the **lmebayesCore** fork; gaps remain on some Gaussian fixed-τ² / fixed-σ²
  routes that still use a chain-outer C++ driver without per-sweep
  cross-chain stats.
- **C++ inner-chain loops and within-block parallel sampling:** Per-sweep
  **inner chain** loops (Block~1 + Block~2 updates across replicate chains)
  should migrate from R orchestration into **`src/*.cpp`** drivers.
  Parallelism should be **within-block, across chains** at fixed inner
  sweep `m` (not parallel inner sweeps): **Block~1** random-effect updates
  over replicate chains first (**higher priority**); **Block~2** fixed-effect
  / hyperparameter updates over chains where safe (**ideal follow-on**).
  Use native threading (e.g. `RcppParallel`) so large `n` does not pay full
  R-loop overhead.
- **OpenCL loader alignment with glmbayes:** `LinkingTo: opencltools`, thin
  manifest-based loader.
- **CRAN path for the iid engine:** Submit the slim backend, point
  **glmbayes** at **glmbayesCore**, and strip duplicated backend code from
  **glmbayes** — then grow mixed-model features back into Core.

---

## License

GPL-2. See the `LICENSE` file and `inst/COPYRIGHTS` for attribution of
incorporated R Mathlib sources.
