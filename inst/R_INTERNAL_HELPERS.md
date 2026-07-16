# `R/` — internal helpers (`@noRd`)

Undocumented helpers in **glmbayesCore** today. Mixed-model / two-block /
lme4-design helpers are temporarily in **lmebayesCore** (Batches 2–5
staging) and will return with those engines.

Index: [R_FUNCTION_INVENTORY.md](R_FUNCTION_INVENTORY.md).

---

## Simulation internals

| Symbol | File | Role |
|--------|------|------|
| `.rindepNormalGamma_reg_impl` | `simfunction.R` | Shared ING implementation behind public wrappers |
| `.lmebayes_check_disp_bounds_or_stop` | `simfunction.R` | Dispersion-bound validation for ING |
| `simfunction.default` | `simfunction.R` | Default method for `simfunction()` |
| `dpois2` | `simulationpipeline.R` | Poisson density helper for pipeline |

## Prior guards / DIC

| Symbol | File | Role |
|--------|------|------|
| `.ing_n_prior_from_shape` | `ing_prior_guard.R` | Map ING shape to effective prior sample size |
| `.ing_stop_if_prior_exceeds_data` | `ing_prior_guard.R` | Guard against prior dominating data |
| `DIC_Info` | `dic_info.R` | DIC / information criteria helper |

## Parallelism

| Symbol | File | Role |
|--------|------|------|
| `use_RcppParallel` | `internal_rcppparallel.R` | Toggle / probe RcppParallel usage |

## Rcpp positional wrappers (`rcpp_wrappers.R`)

Thin `.Call` bridges (not public API):

`.rNormalGLM_cpp`, `.rNormalReg_cpp`, `.rIndepNormalGammaReg_cpp`,
`.rIndepNormalGammaReg_with_envelope_cpp`, `.rNormalGammaReg_cpp`,
`.rGammaGaussian_cpp`, `.rGammaGamma_cpp`, `.rNormalGLM_std_cpp`,
`.rIndepNormalGammaReg_std_cpp`, `.rIndepNormalGammaReg_std_parallel_cpp`,
`.EnvelopeCentering_cpp`, `.EnvelopeSize_cpp`, `.EnvelopeBuild_cpp`,
`.EnvelopeBuild_Ind_Normal_Gamma_cpp`, `.EnvelopeEval_cpp`,
`.EnvelopeDispersionBuild_cpp`, `.EnvelopeOrchestrator_cpp`,
`.EnvelopeSet_Grid_cpp`, `.EnvelopeSet_LogP_cpp`,
`.glmb_Standardize_Model_cpp`, `.glmbayesCore_has_opencl_cpp`,
`.gpu_names_cpp`

Generated `.Call` stubs also live in `R/RcppExports.R` (do not edit by hand).
