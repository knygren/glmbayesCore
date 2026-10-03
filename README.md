# glmbayesCore

![GitHub release (latest by date)](https://img.shields.io/github/v/release/knygren/glmbayesCore?label=version)
![License: GPL-2](https://img.shields.io/badge/license-GPL--2-blue.svg)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/knygren/glmbayesCore/R-CMD-check.yaml?label=R%20CMD%20Check)

**glmbayesCore** is the compiled sampling engine that powers the glmbayes
ecosystem. It holds the C++/OpenCL envelope samplers, the family-function
infrastructure, and the R-level prior and simulation interfaces that
downstream packages depend on. End users should install
[glmbayes](https://github.com/knygren/glmbayes) rather than this package
directly. Backend documentation is in the package vignettes (see
[Vignettes](#vignettes)).

The relationship to the broader ecosystem parallels how `StanHeaders` /
`rstan` serve as the compiled backbone for `rstanarm`: **glmbayesCore** is
the infrastructure layer; **glmbayes** and the in-development **lmebayes**
are the user-facing packages built on top of it.

**Current staging note.** This tree ships the iid GLM/LM envelope engine used
by **glmbayes**, plus the mixed-model **setup** layer used by **lmerb()**:
design and identifiability (`model_setup()`, `check_identifiability()`),
prior calibration (`Prior_Setup_GLMM()`, `Prior_SetupGroup()`), `pfamily`
builders (`pfamily_list()`, `dGamma_list()`) and row-group samplers
(`rNormal_reg_group()`, `rNormalGLM_reg_group()`). The two-block Gibbs
**engines** (`rlmerb()`, `rglmerb()`, `rLMM_reg*`, `rGLMM_reg*`,
`two_block_*`) are still developed in the temporary
[lmebayesCore](https://github.com/knygren/lmebayesCore) fork, which
[lmebayes](https://github.com/knygren/lmebayes) uses. **lmebayesCore** is a
holding package: features move back here once they are stable.

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
| `glmbayes::sim::group` | `simfuncs_groups.h`, `group_utils.cpp`, `rNormalRegGroups.cpp`, `rNormalGLMGroups.cpp` | Row-group partition, prior layout, groupwise draws |
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
| `model_setup.R` / `check_identifiability.R` | Mixed-model design extraction, reference `lmer`/`glmer` fit, two-level identifiability |
| `Prior_Setup_GLMM.R` | Block 1 / Block 2 prior calibration for mixed models |
| `pfamily_list*.R` / `dGamma_list*.R` | Calibrated setups → named lists of `pfamily` objects |
| `Prior_SetupGroup.R` / `normalize_group.R` / `simfunction_group.R` | Row-group priors, partitions and groupwise samplers |

C++ → R callback inventory (both packages): see **glmbayes**
[`data-raw/CPP_R_CALLBACK_INVENTORY.md`](https://github.com/knygren/glmbayes/blob/main/data-raw/CPP_R_CALLBACK_INVENTORY.md).
Re-scan with `Rscript data-raw/cpp_r_callback_inventory.R` in either package tree.

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

Symbols below are exported from **glmbayesCore** today. End users typically
load **glmbayes** (or **lmebayes** for mixed models). The mixed-model setup
functions are exported here; the two-block Gibbs engines still ship from
**lmebayesCore** and will return here.

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

### Mixed-model setup (exported; used by **lmebayes**)

These functions prepare a single-grouping-factor mixed model for the
two-block Gibbs samplers. They are documented in vignette Part 4
([Core 20](#part-4-mixed-model-design-and-priors-enabling-lmerb)).

| Function | Role |
|----------|------|
| `model_setup()` | `lmer`-style formula → sampler design (`y`, `D`, `group`, `W`), reference `lmer`/`glmer` (and optional per-group-dispersion `glmmTMB`) fit |
| `check_identifiability()` | Within-group (Block 1) and across-group (Block 2) identifiability checks |
| `Prior_Setup_GLMM()` | Calibrated Block 2 (`pop.*`) and Block 1 (`group.*`) priors from the reference fit |
| `pfamily_list()` | Block 2 `dNormal` / `dIndependent_Normal_Gamma` priors, one per random-effect coefficient (also a `Prior_SetupGroup` method) |
| `dGamma_list()` | Per-group measurement-dispersion `dGamma` priors (Gaussian, `dispformula = ~group`) |
| `Prior_SetupGroup()`, `normalize_group()` | Independent per-group `Prior_Setup()` and row-group partitions |
| `rNormal_reg_group()`, `rNormalGLM_reg_group()` | Conditionally independent groupwise draws (Block 1 building block) |

### Mixed-model engines (temporary **lmebayesCore**; return here)

| Area | Examples |
|------|----------|
| Matrix drivers | `rlmerb()`, `rglmerb()` |
| Two-block / sweep | `rGLMM_reg*`, `rLMM_reg*`, `rGLMM_sweep()`, `two_block_*`, `plot_sweep_history_diag()` |

Typical **lmebayes** workflow:
`model_setup()` → `Prior_Setup_GLMM()` → `pfamily_list(ps)` (and
`dGamma_list(ps)` for per-group dispersion) → `lmerb()` / `glmerb()`.

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

OpenCL requires a source install of **glmbayesCore** with GPU-capable drivers; see
[Core 11 — GPU acceleration using OpenCL](https://knygren.r-universe.dev/articles/glmbayesCore/Core-11.html)
and [Core 08](https://knygren.r-universe.dev/articles/glmbayesCore/Core-08.html)
in [Vignettes](#vignettes) below.

**Dependencies that must be installed first:**

```r
install.packages(c("Rcpp", "RcppArmadillo", "RcppParallel", "MASS", "Rdpack"))
install.packages(c("opencltools", "nmathopencl"),
                 repos = "https://knygren.r-universe.dev")
```

---

## Vignettes

The **glmbayesCore** vignettes are grouped into **Parts**, following the
same scheme as **glmbayes**. Vignette numbers (`Core-NN`) are permanent
identifiers and are not renumbered, so a Part need not use consecutive
numbers. Subsections use the `Core-NN-S0M` pattern, as in
**glmbayes** `Chapter-02-S0M`. Appendix-style **`Chapter-A*`** names belong
to **glmbayes** only.

After installing from [R-Universe](https://knygren.r-universe.dev/glmbayesCore), use the
links below or in R: `vignette("Core-01", package = "glmbayesCore")`,
`browseVignettes("glmbayesCore")`.

### Part 1: Orientation

- **Core 01 — Overview of the glmbayesCore package**  
  https://knygren.r-universe.dev/articles/glmbayesCore/Core-01.html

### Part 2: Estimation and accept–reject sampling

How the iid samplers draw from each posterior, and how `Prior_Setup()`
calibrates the priors they use.

- **Core 02 — Overview of Estimation Procedures**  
  https://knygren.r-universe.dev/articles/glmbayesCore/Core-02.html

- **Core 03 — Simulation Methods - Likelihood Subgradient Densities**  
  https://knygren.r-universe.dev/articles/glmbayesCore/Core-03.html

- **Core 04 — Accept–Reject Sampling for Dispersion in Gamma Regression**  
  https://knygren.r-universe.dev/articles/glmbayesCore/Core-04.html

- **Core 05 — Accept–Reject Sampling for gaussian Regression models with independent normal-gamma priors**  
  https://knygren.r-universe.dev/articles/glmbayesCore/Core-05.html

- **Core 09 — Implementation Companion for Independent Normal-Gamma**  
  https://knygren.r-universe.dev/articles/glmbayesCore/Core-09.html

- **Core 10 — Technical Derivations for Priors Returned by `Prior_Setup()`**  
  https://knygren.r-universe.dev/articles/glmbayesCore/Core-10.html

### Part 3: Envelopes, parallel CPU, and OpenCL

How envelopes are built and how that work is spread across CPU threads or a
GPU.

- **Core 06 — Overview of Envelope Related Functions**  
  https://knygren.r-universe.dev/articles/glmbayesCore/Core-06.html

- **Core 07 — Parallel Sampling Implementation using RcppParallel**  
  https://knygren.r-universe.dev/articles/glmbayesCore/Core-07.html

- **Core 08 — Accelerated EnvelopeBuild Implementation using OpenCL**  
  https://knygren.r-universe.dev/articles/glmbayesCore/Core-08.html

- **Core 11 — Large models: GPU acceleration using OpenCL (backend installation)**  
  https://knygren.r-universe.dev/articles/glmbayesCore/Core-11.html

### Part 4: Mixed-model design and priors (enabling lmerb)

The setup layer behind `lmerb()` / `glmerb()` in **lmebayes**: design,
identifiability, two-block prior calibration, and conversion to `pfamily`
objects. The user-level tutorials are **glmbayes** Chapters 17 and 18.

- **Core 20 — Mixed-model design and priors in glmbayesCore: overview**  
  https://knygren.r-universe.dev/articles/glmbayesCore/Core-20.html

  - **Core 20-S01 — Model setup and identifiability** (`model_setup()`, `check_identifiability()`)  
    https://knygren.r-universe.dev/articles/glmbayesCore/Core-20-S01.html

  - **Core 20-S02 — Calibrating mixed-model priors** (`Prior_Setup_GLMM()`)  
    https://knygren.r-universe.dev/articles/glmbayesCore/Core-20-S02.html

  - **Core 20-S03 — From calibrated priors to pfamily objects** (`pfamily_list()`, `dGamma_list()`)  
    https://knygren.r-universe.dev/articles/glmbayesCore/Core-20-S03.html

  - **Core 20-S04 — Row-group priors and samplers** (`Prior_SetupGroup()`, `normalize_group()`, `rNormal_reg_group()`, `rNormalGLM_reg_group()`)  
    https://knygren.r-universe.dev/articles/glmbayesCore/Core-20-S04.html

### Part 5: Two-block Gibbs engines (reserved)

Numbers **Core 25–29** are reserved for the two-block Gibbs engines
(`rlmerb()`, `rglmerb()`, `rLMM_reg*`, `rGLMM_reg*`, `two_block_*`). Until
those engines move back from
[lmebayesCore](https://github.com/knygren/lmebayesCore), their
documentation lives in that package. **Core 12** is reserved for a
developer guide on adding a new `pfamily` (today in `inst/ADDING_PFAMILY.md`).

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

---

## Key References

- Nygren, K.N. and Nygren, L.M. (2006). Likelihood Subgradient Densities. *Journal of the American Statistical Association*, 101(475), 1144–1156. — The accept-reject envelope method at the heart of the non-Gaussian samplers.
- Lindley, D.V. and Smith, A.F.M. (1972). Bayes estimates for the linear model. *Journal of the Royal Statistical Society B*, 34, 1–41. — Conjugate Normal-Gamma foundations.
- Gelman, A. et al. (2013). *Bayesian Data Analysis*, 3rd ed. — Reference for prior specifications and dispersion modeling.

A complete bibliography is in `inst/REFERENCES.bib`.

---

## Future plans

- **Reintegrate mixed-model engines from temporary lmebayesCore:** The
  LMM/GLMM setup layer (`model_setup`, `Prior_Setup_GLMM`, `pfamily_list`,
  `dGamma_list`) is now in **glmbayesCore**. Next, merge the matrix drivers
  (`rlmerb` / `rglmerb`) and the two-block / block-ING engines (vignette
  Part 5), then point **lmebayes** at this package again and retire
  **lmebayesCore**.
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
