# glmbayes vs glmbayesCore — comparison assessment

**Generated:** 2026-07-16  
**glmbayesCore:** `C:\Rpackages\glmbayesCore` @ `90ddf14`  
**glmbayes:** `C:\Rpackages\glmbayes` @ `94d0262`  

Regenerate: `source("data-raw/make_compare_glmbayes_glmbayesCore.R")`

---

## 1. Purpose and method

This document inventories and compares the **iid GLM/LM + OpenCL** surfaces of
**glmbayes** and **glmbayesCore** before Batch 7b (OpenCL loader alignment).
Mixed-model code temporarily staged in **lmebayesCore** is out of scope.

### Identity rules

- **Byte-identical:** SHA256 of raw file bytes.
- **Normalized-identical:** SHA256 after UTF-8 decode, CRLF→LF, strip trailing
  whitespace per line (catches Windows line-ending-only drift).
- **Differing:** neither byte nor normalized hash matches.
- **Package-rename noise:** when comparing function bodies, a secondary check
  replaces `glmbayes`/`glmbayesCore` and `_glmbayes_`/`_glmbayesCore_` with a
  placeholder; matches after that are labeled **rename-only**.

### Scope

| Tree | Included |
|------|----------|
| `R/` | all `*.R` |
| `src/` | top-level `*.cpp` / `*.h` / `*.c` (excl. `RcppExports.cpp` from engine tables) |
| `src/nmath/` | Core vendored Mathlib (absent in glmbayes) |
| `inst/cl/` | `*.cl`, dependency manifests (`.tsv`/`.rds`), `README.md` |

`RcppExports.R` / `RcppExports.cpp` are listed where present but treated as
generated noise for engine-parity narrative.

---

## 2. File inventory by type

### 2.1 `R/`

| Category | Count |
|----------|------:|
| glmbayesCore `R/*.R` | 34 |
| glmbayes `R/*.R` | 46 |
| Shared basenames | 26 |
| Byte-identical | 4 |
| Normalized-identical only | 2 |
| Differing | 20 |
| glmbayesCore-only | 8 |
| glmbayes-only | 20 |

#### Identical (byte)

- `compute_gaussian_prior.R`
- `fitter_functions.R`
- `internal_rcppparallel.R`
- `summary.rgamma_reg.R`

#### Identical after normalize only

- `data-BikeSharing.R`
- `data-carinsca.R`

#### Differing (shared basename)

- `RcppExports.R`
- `data-AMI.R`
- `data-Boston_centered.R`
- `data-Cleveland.R`
- `envelopeorchestrator.R`
- `formula.summary.rglmb.R`
- `gamma_ct.R`
- `globals.R`
- `gpu_diagnostics.R`
- `invgamma_ct.R`
- `normal_ct.R`
- `pfamily.R`
- `prior.R`
- `rcpp_wrappers.R`
- `rglmb.R`
- `rlmb.R`
- `simfunction.R`
- `simulationpipeline.R`
- `summary.rglmb.R`
- `zzz.R`

#### glmbayesCore-only

- `dic_info.R`
- `glmbayesCore-package.R`
- `ing_prior_guard.R`
- `multi_prior_setup.R`
- `multi_rNormal_reg.R`
- `multi_rlmb.R`
- `residuals.rglmb.R`
- `summary.mrglmb.R`

#### glmbayes-only

- `anova.glmb.R`
- `case.names.glmb.R`
- `confint.glmb.R`
- `deviance.rglmb.R`
- `directional_tail.R`
- `dummy.coef.glmb.R`
- `extractDIC.R`
- `get_opencl_core_count.R`
- `glmb.R`
- `glmbayes-package.R`
- `influence.glmb.R`
- `lmb.R`
- `logLik.glmb.R`
- `plot.glmb.R`
- `predict.glmb.R`
- `residuals.glmb.R`
- `simulate.glmb.R`
- `summary.glmb.R`
- `summary.mlmb.R`
- `vcov.glmb.R`

### 2.2 `src/` (top-level, excl. `RcppExports.cpp`)

| Category | Count |
|----------|------:|
| glmbayesCore top-level | 43 |
| glmbayes top-level | 43 |
| Shared basenames | 42 |
| Byte-identical | 14 |
| Normalized-identical only | 3 |
| Differing | 25 |
| glmbayesCore-only | 1 |
| glmbayes-only | 1 |

#### Identical (byte)

- `EnvelopeEval.cpp`
- `EnvelopeSort.cpp`
- `Set_Grid.cpp`
- `Set_LogP.cpp`
- `configure_OpenCL.cpp`
- `cuda_probe.cpp`
- `famfuncs.h`
- `famfuncs_poisson.cpp`
- `invgamma_ct.cpp`
- `kernel_runners.cpp`
- `opencl_detect.cpp`
- `rNormalReg.cpp`
- `rng_utils.h`
- `rnorm_ct.cpp`

#### Identical after normalize only

- `EnvelopeCentering.cpp`
- `glmbayes_getRegisteredNamespace.cpp`
- `glmbayes_getRegisteredNamespace.h`

#### Differing (shared basename)

- `EnvelopeBuild.cpp`
- `EnvelopeBuild_Ind_Normal_Gamma.cpp`
- `EnvelopeDispersionBuild.cpp`
- `EnvelopeOrchestrator.cpp`
- `EnvelopeSize.cpp`
- `Envelopefuncs.h`
- `OpenCL_helper.cpp`
- `R_interface.h`
- `export_wrappers.cpp`
- `famfuncs_Gamma.cpp`
- `famfuncs_binomial.cpp`
- `famfuncs_gaussian.cpp`
- `kernel_loader.cpp`
- `kernel_wrappers.cpp`
- `opencl.h`
- `openclPort.h`
- `progress_utils.cpp`
- `progress_utils.h`
- `rGammaGamma.cpp`
- `rGammaGaussian.cpp`
- `rIndepNormalGammaReg.cpp`
- `rNormalGLM.cpp`
- `rNormalGammaReg.cpp`
- `rng_utils.cpp`
- `simfuncs.h`

#### glmbayesCore-only

- `package_ns.h`

#### glmbayes-only

- `rNormalGLMBlocks.cpp`

### 2.3 `src/nmath/` (Core-only vendored Mathlib)

| Category | Count |
|----------|------:|
| glmbayesCore `src/nmath` files | 127 |
| glmbayes `src/nmath` files | 0 |

**glmbayes** does not vendor `src/nmath/`; it relies on **nmathopencl**
(Imports) for OpenCL Mathlib pieces and does not ship the R Mathlib `.c`
tree under `src/`.

<details><summary>glmbayesCore src/nmath file list</summary>

- `nmath/bd0.c`
- `nmath/bessel.h`
- `nmath/bessel_i.c`
- `nmath/bessel_j.c`
- `nmath/bessel_k.c`
- `nmath/bessel_y.c`
- `nmath/beta.c`
- `nmath/chebyshev.c`
- `nmath/choose.c`
- `nmath/cospi.c`
- `nmath/d1mach.c`
- `nmath/dbeta.c`
- `nmath/dbinom.c`
- `nmath/dcauchy.c`
- `nmath/dchisq.c`
- `nmath/dexp.c`
- `nmath/df.c`
- `nmath/dgamma.c`
- `nmath/dgeom.c`
- `nmath/dhyper.c`
- `nmath/dlnorm.c`
- `nmath/dlogis.c`
- `nmath/dnbeta.c`
- `nmath/dnbinom.c`
- `nmath/dnchisq.c`
- `nmath/dnf.c`
- `nmath/dnorm.c`
- `nmath/dnt.c`
- `nmath/dpois.c`
- `nmath/dpq.h`
- `nmath/dt.c`
- `nmath/dunif.c`
- `nmath/dweibull.c`
- `nmath/fmax2.c`
- `nmath/fmin2.c`
- `nmath/fprec.c`
- `nmath/fround.c`
- `nmath/fsign.c`
- `nmath/ftrunc.c`
- `nmath/gamma.c`
- `nmath/gamma_cody.c`
- `nmath/gammalims.c`
- `nmath/i1mach.c`
- `nmath/imax2.c`
- `nmath/imin2.c`
- `nmath/lbeta.c`
- `nmath/lgamma.c`
- `nmath/lgammacor.c`
- `nmath/log1p.c`
- `nmath/mlutils.c`
- `nmath/nmath.h`
- `nmath/nmath2.h`
- `nmath/pbeta.c`
- `nmath/pbinom.c`
- `nmath/pcauchy.c`
- `nmath/pchisq.c`
- `nmath/pexp.c`
- `nmath/pf.c`
- `nmath/pgamma.c`
- `nmath/pgeom.c`
- `nmath/phyper.c`
- `nmath/plnorm.c`
- `nmath/plogis.c`
- `nmath/pnbeta.c`
- `nmath/pnbinom.c`
- `nmath/pnchisq.c`
- `nmath/pnf.c`
- `nmath/pnorm.c`
- `nmath/pnt.c`
- `nmath/polygamma.c`
- `nmath/ppois.c`
- `nmath/pt.c`
- `nmath/ptukey.c`
- `nmath/punif.c`
- `nmath/pweibull.c`
- `nmath/qDiscrete_search.h`
- `nmath/qbeta.c`
- `nmath/qbinom.c`
- `nmath/qcauchy.c`
- `nmath/qchisq.c`
- `nmath/qexp.c`
- `nmath/qf.c`
- `nmath/qgamma.c`
- `nmath/qgeom.c`
- `nmath/qhyper.c`
- `nmath/qlnorm.c`
- `nmath/qlogis.c`
- `nmath/qnbeta.c`
- `nmath/qnbinom.c`
- `nmath/qnbinom_mu.c`
- `nmath/qnchisq.c`
- `nmath/qnf.c`
- `nmath/qnorm.c`
- `nmath/qnt.c`
- `nmath/qpois.c`
- `nmath/qt.c`
- `nmath/qtukey.c`
- `nmath/qunif.c`
- `nmath/qweibull.c`
- `nmath/rbeta.c`
- `nmath/rbinom.c`
- `nmath/rcauchy.c`
- `nmath/rchisq.c`
- `nmath/rexp.c`
- `nmath/rf.c`
- `nmath/rgamma.c`
- `nmath/rgeom.c`
- `nmath/rhyper.c`
- `nmath/rlnorm.c`
- `nmath/rlogis.c`
- `nmath/rmultinom.c`
- `nmath/rnbinom.c`
- `nmath/rnchisq.c`
- `nmath/rnorm.c`
- `nmath/rpois.c`
- `nmath/rt.c`
- `nmath/runif.c`
- `nmath/rweibull.c`
- `nmath/sexp.c`
- `nmath/sign.c`
- `nmath/signrank.c`
- `nmath/snorm.c`
- `nmath/standalone/sunif.c`
- `nmath/standalone/test.c`
- `nmath/stirlerr.c`
- `nmath/toms708.c`
- `nmath/wilcox.c`

</details>

### 2.4 `inst/cl/`

| Category | Count |
|----------|------:|
| glmbayesCore relevant paths | 54 |
| glmbayes relevant paths | 54 |
| Shared relative paths | 54 |
| Byte-identical | 53 |
| Normalized-identical only | 0 |
| Differing | 1 |
| glmbayesCore-only | 0 |
| glmbayes-only | 0 |

#### Identical (byte)

All **51** shared `*.cl` files are byte-identical. Total byte-identical shared paths (incl. manifests): **53**.

<details><summary>Byte-identical path list</summary>

- `OPENCL.cl`
- `R_ext_internals/Parse.cl`
- `R_ext_internals/R_ext_internals.cl`
- `R_ext_internals/stats_package.cl`
- `R_ext_runtime/Arith.cl`
- `R_ext_runtime/Error.cl`
- `R_ext_runtime/MathThreads.cl`
- `R_ext_runtime/Memory.cl`
- `R_ext_runtime/Print.cl`
- `R_ext_runtime/RS.cl`
- `R_ext_runtime/Random.cl`
- `R_ext_runtime/Riconv.cl`
- `R_ext_runtime/Utils.cl`
- `R_ext_types/Boolean.cl`
- `R_ext_types/Complex.cl`
- `R_ext_types/Constants.cl`
- `R_ext_types/Visibility.cl`
- `R_ext_types/libextern.cl`
- `R_shims/Rconfig.cl`
- `R_shims/Rdefines.cl`
- `R_shims/Rinternals.cl`
- `System/stdint.cl`
- `libR_shims/libR.cl`
- `nmath/Rmath.cl`
- `nmath/bd0.cl`
- `nmath/chebyshev.cl`
- `nmath/cospi.cl`
- `nmath/dbinom.cl`
- `nmath/dgamma.cl`
- `nmath/dnorm.cl`
- `nmath/dpois.cl`
- `nmath/dpq.cl`
- `nmath/fmax2.cl`
- `nmath/gamma.cl`
- `nmath/gammalims.cl`
- `nmath/kernel_dependency_index.rds`
- `nmath/kernel_dependency_index.tsv`
- `nmath/lgamma.cl`
- `nmath/lgammacor.cl`
- `nmath/log1p.cl`
- `nmath/nmath.cl`
- `nmath/pgamma_utils.cl`
- `nmath/pnorm.cl`
- `nmath/refactored.cl`
- `nmath/stirlerr.cl`
- `nmath/stirlerr_cycle_dependent.cl`
- `nmath/stirlerr_cycle_free.cl`
- `src/f2_f3_binomial_cloglog.cl`
- `src/f2_f3_binomial_logit.cl`
- `src/f2_f3_binomial_probit.cl`
- `src/f2_f3_gamma.cl`
- `src/f2_f3_gaussian.cl`
- `src/f2_f3_poisson.cl`

</details>

#### Identical after normalize only

*(none)*

#### Differing

- `README.md`

#### glmbayesCore-only

*(none)*

#### glmbayes-only

*(none)*

#### Shared f2/f3 kernel entry points (reference)

- `src/f2_f3_binomial_cloglog.cl`: `f2_f3_binomial_cloglog`
- `src/f2_f3_binomial_logit.cl`: `f2_f3_binomial_logit`
- `src/f2_f3_binomial_probit.cl`: `f2_f3_binomial_probit`
- `src/f2_f3_gamma.cl`: `f2_f3_gamma`
- `src/f2_f3_gaussian.cl`: `f2_f3_gaussian`
- `src/f2_f3_poisson.cl`: `f2_f3_poisson`

### 2.5 OpenCL-related DESCRIPTION fields

| Field | glmbayesCore | glmbayes |
|-------|:------------:|:--------:|
| Imports `opencltools` | True | True |
| Imports `nmathopencl` | True | True |
| LinkingTo `opencltools` | False | True |

---

## 3. Function-level report (differing files)

For each shared file that is not byte- or normalize-identical, top-level
functions / kernels are parsed and compared.

### 3.1 Differing `R/` files

#### `R/RcppExports.R`

Generated by Rcpp::compileAttributes(); skip detailed body compare.

#### `R/data-AMI.R`

Parsed symbols — Core: 0, glmbayes: 0.

| Category | Count |
|----------|------:|
| Identical body | 0 |
| Rename-only (package string) | 0 |
| Different body | 0 |
| Core-only | 0 |
| glmbayes-only | 0 |

#### `R/data-Boston_centered.R`

Parsed symbols — Core: 0, glmbayes: 0.

| Category | Count |
|----------|------:|
| Identical body | 0 |
| Rename-only (package string) | 0 |
| Different body | 0 |
| Core-only | 0 |
| glmbayes-only | 0 |

#### `R/data-Cleveland.R`

Parsed symbols — Core: 0, glmbayes: 0.

| Category | Count |
|----------|------:|
| Identical body | 0 |
| Rename-only (package string) | 0 |
| Different body | 0 |
| Core-only | 0 |
| glmbayes-only | 0 |

#### `R/envelopeorchestrator.R`

Parsed symbols — Core: 2, glmbayes: 2.

| Category | Count |
|----------|------:|
| Identical body | 2 |
| Rename-only (package string) | 0 |
| Different body | 0 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Identical:** `EnvelopeCentering`, `EnvelopeOrchestrator`

#### `R/formula.summary.rglmb.R`

Parsed symbols — Core: 1, glmbayes: 1.

| Category | Count |
|----------|------:|
| Identical body | 0 |
| Rename-only (package string) | 0 |
| Different body | 1 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Different:** `formula.summary.rglmb`

#### `R/gamma_ct.R`

Parsed symbols — Core: 1, glmbayes: 1.

| Category | Count |
|----------|------:|
| Identical body | 1 |
| Rename-only (package string) | 0 |
| Different body | 0 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Identical:** `rgamma_ct`

#### `R/globals.R`

Parsed symbols — Core: 0, glmbayes: 0.

| Category | Count |
|----------|------:|
| Identical body | 0 |
| Rename-only (package string) | 0 |
| Different body | 0 |
| Core-only | 0 |
| glmbayes-only | 0 |

#### `R/gpu_diagnostics.R`

Parsed symbols — Core: 2, glmbayes: 2.

| Category | Count |
|----------|------:|
| Identical body | 0 |
| Rename-only (package string) | 0 |
| Different body | 1 |
| Core-only | 1 |
| glmbayes-only | 1 |

**Different:** `diagnose_glmbayes`

**Core-only:** `glmbayesCore_has_opencl`

**glmbayes-only:** `has_opencl`

#### `R/invgamma_ct.R`

Parsed symbols — Core: 3, glmbayes: 3.

| Category | Count |
|----------|------:|
| Identical body | 3 |
| Rename-only (package string) | 0 |
| Different body | 0 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Identical:** `pinvgamma_ct`, `qinvgamma_ct`, `rinvgamma_ct`

#### `R/normal_ct.R`

Parsed symbols — Core: 2, glmbayes: 2.

| Category | Count |
|----------|------:|
| Identical body | 2 |
| Rename-only (package string) | 0 |
| Different body | 0 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Identical:** `pnorm_ct`, `rnorm_ct`

#### `R/pfamily.R`

Parsed symbols — Core: 7, glmbayes: 7.

| Category | Count |
|----------|------:|
| Identical body | 6 |
| Rename-only (package string) | 0 |
| Different body | 1 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Identical:** `dBeta`, `dGamma`, `dNormal`, `dNormal_Gamma`, `pfamily`, `print.pfamily`

**Different:** `dIndependent_Normal_Gamma`

#### `R/prior.R`

Parsed symbols — Core: 3, glmbayes: 3.

| Category | Count |
|----------|------:|
| Identical body | 3 |
| Rename-only (package string) | 0 |
| Different body | 0 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Identical:** `Prior_Check`, `Prior_Setup`, `print.PriorSetup`

#### `R/rcpp_wrappers.R`

Parsed symbols — Core: 22, glmbayes: 23.

| Category | Count |
|----------|------:|
| Identical body | 0 |
| Rename-only (package string) | 17 |
| Different body | 3 |
| Core-only | 2 |
| glmbayes-only | 3 |

**Rename-only:** `.EnvelopeBuild_Ind_Normal_Gamma_cpp`, `.EnvelopeBuild_cpp`, `.EnvelopeCentering_cpp`, `.EnvelopeEval_cpp`, `.EnvelopeOrchestrator_cpp`, `.EnvelopeSet_Grid_cpp`, `.EnvelopeSet_LogP_cpp`, `.EnvelopeSize_cpp`, `.glmb_Standardize_Model_cpp`, `.rGammaGamma_cpp`, `.rGammaGaussian_cpp`, `.rIndepNormalGammaReg_cpp`, `.rIndepNormalGammaReg_std_cpp`, `.rIndepNormalGammaReg_std_parallel_cpp`, `.rNormalGLM_cpp`, `.rNormalGLM_std_cpp`, `.rNormalGammaReg_cpp`

**Different:** `.EnvelopeDispersionBuild_cpp`, `.gpu_names_cpp`, `.rNormalReg_cpp`

**Core-only:** `.glmbayesCore_has_opencl_cpp`, `.rIndepNormalGammaReg_with_envelope_cpp`

**glmbayes-only:** `.get_opencl_core_count_cpp`, `.has_opencl_cpp`, `.rNormalGLMBlocks_cpp`

#### `R/rglmb.R`

Parsed symbols — Core: 2, glmbayes: 2.

| Category | Count |
|----------|------:|
| Identical body | 2 |
| Rename-only (package string) | 0 |
| Different body | 0 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Identical:** `print.rglmb`, `rglmb`

#### `R/rlmb.R`

Parsed symbols — Core: 2, glmbayes: 2.

| Category | Count |
|----------|------:|
| Identical body | 2 |
| Rename-only (package string) | 0 |
| Different body | 0 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Identical:** `print.rlmb`, `rlmb`

#### `R/simfunction.R`

Parsed symbols — Core: 15, glmbayes: 12.

| Category | Count |
|----------|------:|
| Identical body | 10 |
| Rename-only (package string) | 0 |
| Different body | 2 |
| Core-only | 3 |
| glmbayes-only | 0 |

**Identical:** `.check_gamma_conjugate_scalar_design`, `logdiffexp`, `print.rGamma_reg`, `print.simfunction`, `rBeta_reg`, `rGamma_Conjugate_reg`, `rGamma_reg`, `rNormalGamma_reg`, `simfunction`, `simfunction.default`

**Different:** `rNormal_reg`, `rindepNormalGamma_reg`

**Core-only:** `.lmebayes_check_disp_bounds_or_stop`, `.rindepNormalGamma_reg_impl`, `rindepNormalGamma_reg_with_envelope`

#### `R/simulationpipeline.R`

Parsed symbols — Core: 14, glmbayes: 14.

| Category | Count |
|----------|------:|
| Identical body | 14 |
| Rename-only (package string) | 0 |
| Different body | 0 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Identical:** `EnvelopeBuild`, `EnvelopeDispersionBuild`, `EnvelopeEval`, `EnvelopeOpt`, `EnvelopeSetGrid`, `EnvelopeSetLogP`, `EnvelopeSize`, `EnvelopeSort`, `dpois2`, `glmb_Standardize_Model`, `glmbfamfunc`, `print.glmbfamfunc`, `rIndepNormalGammaReg_std`, `rNormalGLM_std`

#### `R/summary.rglmb.R`

Parsed symbols — Core: 5, glmbayes: 2.

| Category | Count |
|----------|------:|
| Identical body | 0 |
| Rename-only (package string) | 0 |
| Different body | 2 |
| Core-only | 3 |
| glmbayes-only | 0 |

**Different:** `print.summary.rglmb`, `summary.rglmb`

**Core-only:** `.rglmb_get_offset`, `.rglmb_prior_precision`, `summary.rlmb`

#### `R/zzz.R`

Parsed symbols — Core: 4, glmbayes: 4.

| Category | Count |
|----------|------:|
| Identical body | 3 |
| Rename-only (package string) | 0 |
| Different body | 1 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Identical:** `.onAttach`, `.opencl_runtime_sniff`, `.opencl_startup_quiet`

**Different:** `.opencl_startup_message`

### 3.2 Differing `src/` files

#### `src/EnvelopeBuild.cpp`

Parsed symbols — Core: 1, glmbayes: 1.

| Category | Count |
|----------|------:|
| Identical body | 1 |
| Rename-only (package string) | 0 |
| Different body | 0 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Identical:** `EnvelopeBuild`

#### `src/EnvelopeBuild_Ind_Normal_Gamma.cpp`

Parsed symbols — Core: 1, glmbayes: 1.

| Category | Count |
|----------|------:|
| Identical body | 0 |
| Rename-only (package string) | 0 |
| Different body | 1 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Different:** `EnvelopeBuild_Ind_Normal_Gamma`

#### `src/EnvelopeDispersionBuild.cpp`

Parsed symbols — Core: 20, glmbayes: 20.

| Category | Count |
|----------|------:|
| Identical body | 18 |
| Rename-only (package string) | 0 |
| Different body | 2 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Identical:** `EnvBuildLinBound_cpp`, `RSS`, `UB2`, `bisection_root`, `bound_rss_over_dispersion`, `compute_envelope_geometry_cpp`, `compute_mixture_and_outputs_cpp`, `g_of_t`, `hprime_of_t`, `max_vec`, `minimize_ub2_over_dispersion`, `rss_face_at_disp`, `rss_face_bound_from_cache_cpp`, `rss_face_quadratic_sum_internal`, `run_ub2_pilot_block`, `thetabar_const_cpp`, `ub2_min_exact_1d`, `ub2_reduced`

**Different:** `EnvelopeDispersionBuild`, `bound_ub2_over_dispersion`

#### `src/EnvelopeOrchestrator.cpp`

Parsed symbols — Core: 1, glmbayes: 1.

| Category | Count |
|----------|------:|
| Identical body | 0 |
| Rename-only (package string) | 0 |
| Different body | 1 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Different:** `EnvelopeOrchestrator`

#### `src/EnvelopeSize.cpp`

Parsed symbols — Core: 1, glmbayes: 1.

| Category | Count |
|----------|------:|
| Identical body | 0 |
| Rename-only (package string) | 0 |
| Different body | 1 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Different:** `EnvelopeSize`

#### `src/Envelopefuncs.h`

Parsed symbols — Core: 21, glmbayes: 16.

| Category | Count |
|----------|------:|
| Identical body | 8 |
| Rename-only (package string) | 0 |
| Different body | 8 |
| Core-only | 5 |
| glmbayes-only | 0 |

**Identical:** `decl:EnvelopeBuild_Ind_Normal_Gamma`, `decl:EnvelopeCentering`, `decl:EnvelopeOrchestrator`, `decl:EnvelopeSet_Grid_C2`, `decl:EnvelopeSet_Grid_C2_pointwise`, `decl:EnvelopeSort_cpp`, `decl:UB2`, `decl:setlogP_C2`

**Different:** `decl:EnvelopeBuild`, `decl:EnvelopeDispersionBuild`, `decl:EnvelopeEval`, `decl:EnvelopeSet_Grid`, `decl:EnvelopeSet_LogP`, `decl:EnvelopeSize`, `decl:RSS`, `decl:rss_face_at_disp`

**Core-only:** `check_disp_bounds_or_stop`, `check_disp_bounds_or_stop(Rcpp::Nullable<Rcpp::NumericVector> disp_lower, Rcpp::Nullab)`, `check_disp_bounds_or_stop(Rcpp::Nullable<double> disp_lower, Rcpp::Nullable<double> di)`, `decl:bound_ub2_over_dispersion`, `decl:check_disp_bounds_or_stop`

#### `src/OpenCL_helper.cpp`

Parsed symbols — Core: 4, glmbayes: 3.

| Category | Count |
|----------|------:|
| Identical body | 2 |
| Rename-only (package string) | 0 |
| Different body | 0 |
| Core-only | 2 |
| glmbayes-only | 1 |

**Identical:** `copyVector`, `flattenMatrix`

**Core-only:** `glmbayesCore_has_opencl`, `opencl_core_count_for_scaling`

**glmbayes-only:** `has_opencl`

#### `src/R_interface.h`

Parsed symbols — Core: 39, glmbayes: 40.

| Category | Count |
|----------|------:|
| Identical body | 32 |
| Rename-only (package string) | 0 |
| Different body | 4 |
| Core-only | 3 |
| glmbayes-only | 4 |

**Identical:** `decl:r_as_matrix`, `decl:r_as_numeric`, `decl:r_as_vector`, `decl:r_expand_grid`, `decl:r_format`, `decl:r_gaussian`, `decl:r_interactive`, `decl:r_lm_fit`, `decl:r_lm_wfit`, `decl:r_optim`, `decl:r_qgamma`, `decl:r_readline`, `decl:r_runif`, `decl:r_sys_time`, `decl:r_system_file`, `decl:r_try`, `r_as_matrix`, `r_as_numeric`, `r_as_vector`, `r_expand_grid`, `r_format`, `r_gaussian`, `r_interactive`, `r_lm_fit`, `r_lm_wfit`, `r_optim`, `r_qgamma`, `r_readline`, `r_runif`, `r_sys_time`, `r_system_file`, `r_try`

**Different:** `r_envelope_opt`, `r_envelope_sort`, `r_rNormal_reg_wfit`, `r_rgamma_ct`

**Core-only:** `decl:pkg_env`, `pkg_env`, `r_glmbfamfunc`

**glmbayes-only:** `decl:r_envelope_opt`, `decl:r_envelope_sort`, `decl:r_rNormal_reg_wfit`, `decl:r_rgamma_ct`

#### `src/export_wrappers.cpp`

Parsed symbols — Core: 48, glmbayes: 50.

| Category | Count |
|----------|------:|
| Identical body | 42 |
| Rename-only (package string) | 0 |
| Different body | 2 |
| Core-only | 4 |
| glmbayes-only | 6 |

**Identical:** `EnvelopeBuild_Ind_Normal_Gamma_cpp_export`, `EnvelopeBuild_Ind_Normal_Gamma_cpp_export(const Rcpp::NumericVector& bStar, const Rcpp::NumericMatrix&)`, `EnvelopeBuild_cpp_export`, `EnvelopeBuild_cpp_export(Rcpp::NumericVector bStar, Rcpp::NumericMatrix A, Rcpp::Nume)`, `EnvelopeCentering_cpp_export`, `EnvelopeCentering_cpp_export(const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& x, )`, `EnvelopeDispersionBuild_cpp_export`, `EnvelopeDispersionBuild_cpp_export(const Rcpp::List& Env, double Shape, double Rate, const Rcpp)`, `EnvelopeEval_cpp_export`, `EnvelopeEval_cpp_export(const Rcpp::NumericMatrix& G4, const Rcpp::NumericVector& y,)`, `EnvelopeOrchestrator_cpp_export`, `EnvelopeOrchestrator_cpp_export(const Rcpp::NumericVector& bstar2, const Rcpp::NumericMatrix)`, `EnvelopeSet_Grid_cpp_export`, `EnvelopeSet_Grid_cpp_export(const Rcpp::NumericMatrix& GIndex, const Rcpp::NumericMatrix)`, `EnvelopeSet_LogP_cpp_export`, `EnvelopeSet_LogP_cpp_export(const Rcpp::NumericMatrix& logP, const Rcpp::NumericVector& )`, `EnvelopeSize_cpp_export`, `EnvelopeSize_cpp_export(const arma::vec& a, const Rcpp::NumericMatrix& G1, int Gridt)`, `UB2_cpp_export`, `UB2_cpp_export(double dispersion, const Rcpp::List& cache, const Rcpp::Nume)`, `glmb_Standardize_Model_cpp_export`, `glmb_Standardize_Model_cpp_export(const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& x, )`, `gpu_names_cpp_export`, `gpu_names_cpp_export()`, `rGammaGamma_cpp_export`, `rGammaGamma_cpp_export(int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatr)`, `rGammaGaussian_cpp_export`, `rGammaGaussian_cpp_export(int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatr)`, `rIndepNormalGammaReg_std_cpp_export`, `rIndepNormalGammaReg_std_cpp_export(int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatr)`, `rIndepNormalGammaReg_std_parallel_cpp_export`, `rIndepNormalGammaReg_std_parallel_cpp_export(int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatr)`, `rNormalGLM_cpp_export`, `rNormalGLM_cpp_export(int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatr)`, `rNormalGLM_std_cpp_export`, `rNormalGLM_std_cpp_export(int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatr)`, `rNormalGammaReg_cpp_export`, `rNormalGammaReg_cpp_export(int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatr)`, `rNormalReg_cpp_export`, `rNormalReg_cpp_export(int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatr)`, `rss_face_at_disp_cpp_export`, `rss_face_at_disp_cpp_export(double dispersion, const Rcpp::List& cache, const Rcpp::Nume)`

**Different:** `rIndepNormalGammaReg_cpp_export`, `rIndepNormalGammaReg_cpp_export(int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatr)`

**Core-only:** `glmbayesCore_has_opencl_cpp_export`, `glmbayesCore_has_opencl_cpp_export()`, `rIndepNormalGammaReg_with_envelope_cpp_export`, `rIndepNormalGammaReg_with_envelope_cpp_export(int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatr)`

**glmbayes-only:** `get_opencl_core_count_cpp_export`, `get_opencl_core_count_cpp_export()`, `has_opencl_cpp_export`, `has_opencl_cpp_export()`, `rNormalGLMBlocks_cpp_export`, `rNormalGLMBlocks_cpp_export(int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatr)`

#### `src/famfuncs_Gamma.cpp`

Parsed symbols — Core: 6, glmbayes: 6.

| Category | Count |
|----------|------:|
| Identical body | 5 |
| Rename-only (package string) | 0 |
| Different body | 1 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Identical:** `f1_gamma`, `f2_f3_gamma`, `f2_gamma`, `f2_gamma_rmat`, `f3_gamma`

**Different:** `dgamma_glmb`

#### `src/famfuncs_binomial.cpp`

Parsed symbols — Core: 16, glmbayes: 16.

| Category | Count |
|----------|------:|
| Identical body | 13 |
| Rename-only (package string) | 0 |
| Different body | 3 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Identical:** `f1_binomial_cloglog`, `f1_binomial_logit`, `f1_binomial_probit`, `f2_binomial_cloglog`, `f2_binomial_cloglog_rmat`, `f2_binomial_logit`, `f2_binomial_logit_rmat`, `f2_binomial_probit`, `f2_f3_binomial_cloglog`, `f2_f3_binomial_logit`, `f3_binomial_cloglog`, `f3_binomial_logit`, `f3_binomial_probit`

**Different:** `dbinom_glmb`, `f2_binomial_probit_rmat`, `f2_f3_binomial_probit`

#### `src/famfuncs_gaussian.cpp`

Parsed symbols — Core: 10, glmbayes: 10.

| Category | Count |
|----------|------:|
| Identical body | 10 |
| Rename-only (package string) | 0 |
| Different body | 0 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Identical:** `Inv_f3_precompute_disp`, `Inv_f3_with_disp`, `Inv_f3_with_disp_rmat`, `dnorm_glmb`, `f1_gaussian`, `f2_f3_gaussian`, `f2_gaussian`, `f2_gaussian_rmat`, `f2_gaussian_rmat_mat`, `f3_gaussian`

#### `src/kernel_loader.cpp`

Parsed symbols — Core: 9, glmbayes: 9.

| Category | Count |
|----------|------:|
| Identical body | 1 |
| Rename-only (package string) | 0 |
| Different body | 5 |
| Core-only | 3 |
| glmbayes-only | 3 |

**Identical:** `resolve_kernel_path`

**Different:** `load_kernel_library`, `load_kernel_source`, `load_library_for_kernel`, `load_library_for_kernel_cross_package`, `load_likelihood_subgradient_program`

**Core-only:** `load_likelihood_subgradient_program_v2`, `parse_cl_tag`, `read_tsv_index`

**glmbayes-only:** `get_opencl_core_count`, `load_program_preload`, `opencltools_take_cstr`

#### `src/kernel_wrappers.cpp`

Parsed symbols — Core: 1, glmbayes: 1.

| Category | Count |
|----------|------:|
| Identical body | 0 |
| Rename-only (package string) | 0 |
| Different body | 1 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Different:** `f2_f3_opencl`

#### `src/opencl.h`

Parsed symbols — Core: 4, glmbayes: 3.

| Category | Count |
|----------|------:|
| Identical body | 1 |
| Rename-only (package string) | 0 |
| Different body | 2 |
| Core-only | 1 |
| glmbayes-only | 0 |

**Identical:** `decl:f2_f3_opencl`

**Different:** `decl:load_likelihood_subgradient_program`, `decl:openclPort`

**Core-only:** `decl:load_likelihood_subgradient_program_v2`

#### `src/openclPort.h`

Parsed symbols — Core: 9, glmbayes: 12.

| Category | Count |
|----------|------:|
| Identical body | 5 |
| Rename-only (package string) | 0 |
| Different body | 2 |
| Core-only | 2 |
| glmbayes-only | 5 |

**Identical:** `decl:configureOpenCL`, `decl:copyVector`, `decl:detect_num_gpus_internal`, `decl:flattenMatrix`, `decl:gpu_names`

**Different:** `decl:load_kernel_library`, `decl:load_kernel_source`

**Core-only:** `decl:glmbayesCore_has_opencl`, `decl:opencl_core_count_for_scaling`

**glmbayes-only:** `decl:get_opencl_core_count`, `decl:has_opencl`, `decl:load_library_for_kernel`, `decl:load_library_for_kernel_cross_package`, `decl:load_program_preload`

#### `src/progress_utils.cpp`

Parsed symbols — Core: 2, glmbayes: 1.

| Category | Count |
|----------|------:|
| Identical body | 0 |
| Rename-only (package string) | 0 |
| Different body | 1 |
| Core-only | 1 |
| glmbayes-only | 0 |

**Different:** `progress_bar`

**Core-only:** `progress_bar_finish`

#### `src/progress_utils.h`

Parsed symbols — Core: 14, glmbayes: 13.

| Category | Count |
|----------|------:|
| Identical body | 12 |
| Rename-only (package string) | 0 |
| Different body | 1 |
| Core-only | 1 |
| glmbayes-only | 0 |

**Identical:** `begin`, `decl:begin`, `decl:localtime_r`, `decl:localtime_s`, `decl:now_hms`, `decl:print_completed`, `format_hms`, `format_hms(int h, int m, int s)`, `format_int_with_commas`, `now_hms`, `print_completed`, `timestamp_cpp`

**Different:** `decl:progress_bar`

**Core-only:** `decl:progress_bar_finish`

#### `src/rGammaGamma.cpp`

Parsed symbols — Core: 1, glmbayes: 1.

| Category | Count |
|----------|------:|
| Identical body | 0 |
| Rename-only (package string) | 0 |
| Different body | 1 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Different:** `rGammaGamma`

#### `src/rGammaGaussian.cpp`

Parsed symbols — Core: 1, glmbayes: 1.

| Category | Count |
|----------|------:|
| Identical body | 0 |
| Rename-only (package string) | 0 |
| Different body | 1 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Different:** `rGammaGaussian`

#### `src/rIndepNormalGammaReg.cpp`

Parsed symbols — Core: 6, glmbayes: 6.

| Category | Count |
|----------|------:|
| Identical body | 4 |
| Rename-only (package string) | 0 |
| Different body | 2 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Identical:** `g1_face_at_disp`, `g2_face_at_disp`, `rIndepNormalGammaReg_std_parallel`, `rIndepNormalGammaReg_worker`

**Different:** `rIndepNormalGammaReg`, `rIndepNormalGammaReg_std`

#### `src/rNormalGLM.cpp`

Parsed symbols — Core: 7, glmbayes: 7.

| Category | Count |
|----------|------:|
| Identical body | 6 |
| Rename-only (package string) | 0 |
| Different body | 1 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Identical:** `make_rcppparallel_pilot_result`, `rNormalGLM`, `rNormalGLM_std`, `rNormalGLM_std_parallel`, `rNormalGLM_worker`, `run_rcppparallel_pilot`

**Different:** `glmb_Standardize_Model`

#### `src/rNormalGammaReg.cpp`

Parsed symbols — Core: 1, glmbayes: 1.

| Category | Count |
|----------|------:|
| Identical body | 0 |
| Rename-only (package string) | 0 |
| Different body | 1 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Different:** `rNormalGammaReg`

#### `src/rng_utils.cpp`

Parsed symbols — Core: 6, glmbayes: 6.

| Category | Count |
|----------|------:|
| Identical body | 2 |
| Rename-only (package string) | 0 |
| Different body | 4 |
| Core-only | 0 |
| glmbayes-only | 0 |

**Identical:** `log_p_inv_gamma_ct_safe`, `runif_safe`

**Different:** `log_p_inv_gamma_safe`, `p_inv_gamma_safe`, `q_inv_gamma_safe`, `rinvgamma_ct_safe`

#### `src/simfuncs.h`

Parsed symbols — Core: 11, glmbayes: 11.

| Category | Count |
|----------|------:|
| Identical body | 6 |
| Rename-only (package string) | 0 |
| Different body | 4 |
| Core-only | 1 |
| glmbayes-only | 1 |

**Identical:** `decl:glmb_Standardize_Model`, `decl:rGammaGaussian`, `decl:rIndepNormalGammaReg_std_parallel`, `decl:rNormalGLM`, `decl:rNormalGLM_std`, `decl:rNormalReg`

**Different:** `decl:rGammaGamma`, `decl:rIndepNormalGammaReg`, `decl:rIndepNormalGammaReg_std`, `decl:rNormalGammaReg`

**Core-only:** `decl:rNormalGLM_optim_poisson_log`

**glmbayes-only:** `decl:rNormalGLMBlocks`

### 3.3 Differing `inst/cl/` files

#### `inst/cl/README.md`

Non-.cl differing file; hash-only.

---

## 4. Summary overview

### Snapshot

| Layer | Shared | Byte-identical | Norm-only | Differing |
|-------|-------:|---------------:|----------:|----------:|
| `R/` | 26 | 4 | 2 | 20 |
| `src/` (top-level) | 42 | 14 | 3 | 25 |
| `inst/cl/` | 54 | 53 | 0 | 1 |

### What is already the same

- **OpenCL program sources:** all **51** shared `*.cl` files
  are **byte-identical** (entry kernels under `src/f2_f3_*.cl` plus
  nmath/shim prelude). Manifests (`.tsv`/`.rds`) also match.
  The only shared `inst/cl` text that differs is **`README.md`**
  (documentation wording, not kernels).
- **C++ already in lockstep (byte):** 14 files,
  notably **`kernel_runners.cpp`**, `configure_OpenCL.cpp`,
  `opencl_detect.cpp`, `rNormalReg.cpp`, `rnorm_ct.cpp`, `invgamma_ct.cpp`,
  `EnvelopeEval.cpp`, `EnvelopeSort.cpp`, `Set_Grid.cpp`, `Set_LogP.cpp`,
  `famfuncs.h`, `famfuncs_poisson.cpp`, `cuda_probe.cpp`, `rng_utils.h`.
- **C++ normalize-only (line endings):** `EnvelopeCentering.cpp`, `glmbayes_getRegisteredNamespace.cpp`, `glmbayes_getRegisteredNamespace.h`.
- **R already in lockstep (byte):** `compute_gaussian_prior.R`, `fitter_functions.R`, `internal_rcppparallel.R`, `summary.rgamma_reg.R`.
- **R normalize-only:** `data-BikeSharing.R`, `data-carinsca.R`.

### Structural / packaging differences (not logic forks)

| Topic | glmbayesCore | glmbayes |
|-------|--------------|----------|
| Package identity | `package_ns.h`, Core symbol prefixes | No `package_ns.h`; glmbayes prefixes |
| Formula / S3 UX | Matrix API focus (`rglmb`/`rlmb`) | Full `glmb()`/`lmb()` + many S3 methods |
| Multi-response helpers | `multi_*` R files present | Not in this tree |
| Vendored CPU nmath under `src/` | Yes (127 files) | No |
| OpenCL Mathlib at runtime | Imports `nmathopencl`; still vendors `inst/cl/nmath` | Same Imports; thin loader via opencltools C API |
| `LinkingTo: opencltools` | **False** | **True** |
| `kernel_loader.cpp` | Fat in-tree `system.file` / dependency assembly | Thin opencltools C-API wrapper |
| Blocks leftover | Stripped | Still has `rNormalGLMBlocks.cpp` |

### Function-level deltas on differing files (§3 aggregates)

| Layer | Identical bodies | Rename-only | Different bodies |
|-------|-----------------:|------------:|-----------------:|
| Differing `R/` | 48 | 17 | 11 |
| Differing `src/` | 168 | 0 | 49 |

Rename-only means the function text matches after substituting package
name / dynlib prefixes. **Different** is the residual that needs a human
merge decision when Core becomes the shared backend.

### Substantive engine diffs (high level)

- **OpenCL loader (Batch 7b target):** `kernel_loader.cpp` is a real fork. Parsed symbols — identical 1, different 5, Core-only 3, glmbayes-only 3.
  Meanwhile **`kernel_runners.cpp` is already byte-identical**, so the
  enqueue/runtime path does not need a port — only the program-assembly
  loader and `LinkingTo: opencltools`.
- **Envelope / famfuncs / samplers:** several large `src/` files still
  differ beyond rename (see §3.2 `different` lists for
  `EnvelopeBuild*.cpp`, `rNormalGLM.cpp`, `rIndepNormalGammaReg.cpp`,
  `export_wrappers.cpp`, …). Treat these as a separate sync track from
  OpenCL loader alignment.
- **Core-only R capabilities:** `rindepNormalGamma_reg_with_envelope`,
  multi-response samplers (`multi_*`), `ing_prior_guard` helpers.
- **glmbayes-only R layer:** formula modelling + diagnostics S3
  (`glmb.R`, `lmb.R`, `summary.glmb.R`, `predict.glmb.R`, …) — stays in
  glmbayes when Core is the engine.

### Recommended Batch 7b actions (code not done here)

1. Port glmbayes’ thin `kernel_loader.cpp` into Core; archive the fat
   loader under `src/backup/` (same pattern as glmbayes).
2. Add `LinkingTo: opencltools` to Core `DESCRIPTION`.
3. Keep shared `*.cl` trees as-is (already identical); optionally sync
   `inst/cl/README.md` wording. Later: decide whether to stop shipping
   vendored nmath/shims and rely solely on **nmathopencl**.
4. Do **not** port `rNormalGLMBlocks.cpp` back into Core as part of OpenCL
   alignment (mixed-model path belongs with gradual lmebayesCore merge).

---

## 5. Appendix

### Regenerate

```r
# From glmbayesCore package root
source("data-raw/make_compare_glmbayes_glmbayesCore.R")
# Optional override:
# Sys.setenv(GLMBAYES_COMPARE_ROOT = "C:/path/to/glmbayes")
```

### Machine-readable sidecar

`data-raw/compare_glmbayes_glmbayesCore.json` — file hashes and function
bucket lists for the differing files.

