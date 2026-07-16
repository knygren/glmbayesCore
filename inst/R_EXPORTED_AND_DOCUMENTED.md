# `R/` — exported and documented functions

Catalog of symbols in **`NAMESPACE`** (or with a `man/*.Rd` page) **in this
tree today**. Mixed-model exports (`rlmerb`, `rGLMM_*`, `two_block_*`,
`model_setup`, …) are temporarily in **lmebayesCore** pending gradual
reintegration into **glmbayesCore**.

Index: [R_FUNCTION_INVENTORY.md](R_FUNCTION_INVENTORY.md).

---

## Also exported from **glmbayes** (shared API)

Keep signatures aligned while both packages export a symbol.

### Retain as **glmbayes** re-exports

| Function | File | Role |
|----------|------|------|
| `Prior_Setup()` | `prior.R` | Default prior calibration |
| `Prior_Check()` | `prior.R` | Prior predictive checks |
| `pfamily()` | `pfamily.R` | Prior-family generic |
| `dNormal()` | `pfamily.R` | Multivariate Normal prior |
| `dGamma()` | `pfamily.R` | Gamma prior |
| `dBeta()` | `pfamily.R` | Beta prior |
| `dNormal_Gamma()` | `pfamily.R` | Normal–Gamma prior |
| `dIndependent_Normal_Gamma()` | `pfamily.R` | Independent Normal–Gamma prior |
| `multi_prior_setup()` | `multi_prior_setup.R` | Multi-response Gaussian prior setup |
| `multi_rlmb()` | `multi_rlmb.R` | Multi-response LM draws |
| `rglmb()` | `rglmb.R` | Matrix-level Bayesian GLM sampler |
| `rlmb()` | `rlmb.R` | Matrix-level Bayesian LM sampler |
| `diagnose_glmbayes()` | `gpu_diagnostics.R` | OpenCL / GPU diagnostic report |

### Phase out of **glmbayes** (stay in **glmbayesCore**)

| Function | File | Role |
|----------|------|------|
| `compute_gaussian_prior()` | `compute_gaussian_prior.R` | Gaussian calibration inside `Prior_Setup()` |
| `simfunction()` | `simfunction.R` | Simulation registry generic |
| `glmbfamfunc()` | `simulationpipeline.R` | GLM family pipeline helpers |
| `rNormal_reg()` | `simfunction.R` | Normal / envelope coefficient sampler |
| `rNormalGamma_reg()` | `simfunction.R` | Conjugate Normal–Gamma sampler |
| `rindepNormalGamma_reg()` | `simfunction.R` | Independent Normal–Gamma sampler |
| `rindepNormalGamma_reg_with_envelope()` | `simfunction.R` | ING with returned envelope |
| `rGamma_reg()` | `simfunction.R` | Gamma dispersion sampler |
| `rGamma_Conjugate_reg()` | `simfunction.R` | Conjugate Gamma rate sampler |
| `rBeta_reg()` | `simfunction.R` | Conjugate Beta–Binomial sampler |
| `multi_rNormal_reg()` | `multi_rNormal_reg.R` | Multi-response Normal sampler |
| `multi_rNormalGamma_reg()` | `multi_rlmb.R` | Multi-response Normal–Gamma sampler |
| `multi_rindepNormalGamma_reg()` | `multi_rlmb.R` | Multi-response ING sampler |
| `rNormalGLM_std()` | `simulationpipeline.R` | Standardized Normal GLM sampler |
| `rIndepNormalGammaReg_std()` | `simulationpipeline.R` | Standardized ING sampler |
| `glmb.wfit()` / `rNormal_reg.wfit()` | `fitter_functions.R` | Weighted fitter hooks |
| `glmb_Standardize_Model()` | `simulationpipeline.R` | Design / prior standardization |
| `EnvelopeBuild()`, `EnvelopeEval()`, `EnvelopeSize()`, `EnvelopeSort()`, `EnvelopeOpt()`, `EnvelopeSetGrid()`, `EnvelopeSetLogP()`, `EnvelopeDispersionBuild()` | `simulationpipeline.R` | Envelope machinery |
| `EnvelopeOrchestrator()`, `EnvelopeCentering()` | `envelopeorchestrator.R` | Envelope orchestration |
| `glmbayesCore_has_opencl()` | `gpu_diagnostics.R` | OpenCL availability probe |
| `pnorm_ct()`, `rnorm_ct()` | `normal_ct.R` | Truncated normal callbacks |
| `pinvgamma_ct()`, `qinvgamma_ct()`, `rinvgamma_ct()` | `invgamma_ct.R` | Inverse-gamma callbacks |
| `rgamma_ct()` | `gamma_ct.R` | Truncated gamma callback |

---

## S3 methods (Core)

Registered in this package's `NAMESPACE`:

`print` / `summary` / `residuals` / `formula` methods for `rglmb`, `rlmb`,
`rGamma_reg`, `mrglmb`, `pfamily`, `simfunction`, `PriorSetup`, `glmbfamfunc`.
