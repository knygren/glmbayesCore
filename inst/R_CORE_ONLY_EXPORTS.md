# **glmbayesCore** exports by function type

Catalog of exports **in this tree today** (iid GLM/LM). Mixed-model /
two-block / block-ING exports are temporarily in **lmebayesCore** and will
return here.

Companion overlap matrix: [R_EXPORTED_AND_DOCUMENTED.md](R_EXPORTED_AND_DOCUMENTED.md).

---

## Prior / pfamily

`Prior_Setup`, `Prior_Check`, `compute_gaussian_prior`, `pfamily`, `dNormal`,
`dGamma`, `dBeta`, `dNormal_Gamma`, `dIndependent_Normal_Gamma`,
`multi_prior_setup`

## Matrix samplers

`rglmb`, `rlmb`, `multi_rlmb`

## Simulation functions (`simfunction`)

`simfunction`, `rNormal_reg`, `rNormalGamma_reg`, `rindepNormalGamma_reg`,
`rindepNormalGamma_reg_with_envelope`, `rGamma_reg`, `rGamma_Conjugate_reg`,
`rBeta_reg`, `multi_rNormal_reg`, `multi_rNormalGamma_reg`,
`multi_rindepNormalGamma_reg`

## Standardized / fitter hooks

`rNormalGLM_std`, `rIndepNormalGammaReg_std`, `glmb_Standardize_Model`,
`glmb.wfit`, `rNormal_reg.wfit`, `glmbfamfunc`

## Envelope machinery

`EnvelopeBuild`, `EnvelopeEval`, `EnvelopeSize`, `EnvelopeSort`, `EnvelopeOpt`,
`EnvelopeSetGrid`, `EnvelopeSetLogP`, `EnvelopeDispersionBuild`,
`EnvelopeOrchestrator`, `EnvelopeCentering`

## OpenCL / diagnostics

`diagnose_glmbayes`, `glmbayesCore_has_opencl`

## Truncated-distribution callbacks

`pnorm_ct`, `rnorm_ct`, `pinvgamma_ct`, `qinvgamma_ct`, `rinvgamma_ct`,
`rgamma_ct`
