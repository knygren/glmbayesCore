# glmbayesCore (development version)

* **`pfamily_list()` generic:** New S3 generic for building a named list of
  pfamily objects from a prior-specification container.  Downstream
  packages provide methods (e.g. `lmebayes` for `Prior_Setup_lmebayes()`
  objects, mapping each random-effect component to `dNormal()` or
  `dIndependent_Normal_Gamma()`).

* **Convergence rate for the two-block sampler:** New **`two_block_rate()`**
  computes the eigenvalues of
  `A = P11^{-1/2} P12 P22^{-1} P21 P11^{-1/2}` (Nygren 2020, Remark 8) for
  the joint Gaussian posterior targeted by **`two_block_rNormal_reg()`**,
  without ever forming the `J*p_re x J*p_re` Block 1 precision: the cross
  moment is accumulated per group with `p_re x p_re` solves followed by a
  single `q x q` symmetric eigendecomposition. The maximal eigenvalue
  `lambda*` is the geometric TV contraction rate of the sampler;
  `m_for_tol(tol)` returns the implied number of inner Gibbs sweeps. For
  non-Gaussian families explicit IRLS-style `weights` give a local-Gaussian
  heuristic. Validated against a dense brute-force construction of the joint
  precision and against the observed contraction of the ICM mean recursion
  (`lmerb_posterior_mean()`), which contracts at exactly `lambda*`.

* **Likelihood precision at the posterior mode:** New
  **`two_block_mode_weights()`** evaluates per-observation likelihood
  precisions (IRLS/Fisher weights) at a supplied random-effects value -
  typically the joint posterior mode from `glmerb_posterior_mode()` - and
  assembles the per-group likelihood precision blocks `Z_j' W_j Z_j`.
  Weights are computed generically from the family object
  (`w_i = wt_i mu'(eta_i)^2 / (V(mu_i) phi)`): exact observed Hessian for
  canonical links (gaussian, poisson-log, binomial-logit), expected (Fisher)
  information otherwise - including correct probit/cloglog/Gamma-log weights
  where `glmbfamfunc()$f7` carries copy-pasted logistic weights.  The
  `weights` component feeds `two_block_rate(weights = )` directly, providing
  the local-Gaussian heuristic input for extending the TV-rate analysis to
  non-Gaussian `glmerb` models.  Validated against `f7` on its correct
  branches and against the exact Gaussian rate path.

* **Explicit TV convergence bounds:** New **`two_block_tv_bound()`**
  evaluates the total-variation bound between the `l`-step kernel and the
  target (Nygren 2020) from the `two_block_rate()` spectrum, two ways:
  `method = "theorem3"` computes the exact per-eigendirection terms
  `d_i^(l)` using the closed form `erf_n(x) = pchisq(2 x^2, n)` with
  `r_i^(l) = (1 - a_{i-1}^{2l})/(1 - a_i^{2l})`; `method = "corollary1"`
  evaluates the looser geometric envelope with explicit constants. With the
  chain started at the exact posterior mean (as `lmerb` does), the mean term
  vanishes identically (`D0 = 0` default) and only the variance-convergence
  sum remains, which decays like `lambda*^{2l}` - twice the exponent of the
  crude `(lambda*)^m` proxy. **`two_block_l_for_tv()`** inverts the bound to
  give the number of inner Gibbs sweeps required for a target tolerance, and
  `print.two_block_rate()` now tabulates proxy vs Theorem 3 vs Corollary 1
  sweeps. On the lmerb big_word_club example (`lambda* = 0.839`): TV <= 1e-3
  needs 16 sweeps (Theorem 3) / 23 (Corollary 1) vs 40 for the proxy.

* **Two-block Gibbs loop in C++:** The main loop of
  **`two_block_rNormal_reg()`** (Block 1 random-effects update, Block 2
  hyperparameter update, `m_convergence` inner steps, replicate sampling) now
  runs entirely in C++ (`two_block_rNormal_reg_cpp_export` in
  `src/twoBlockGibbs.cpp`), eliminating per-iteration R/C++ round trips. This
  is a port-only change: the R wrapper still performs input validation,
  `glmbfamfunc()` resolution, and output assembly, and the C++ driver calls
  the same per-block samplers (`rNormalGLM` envelope sampler, `rNormalReg`)
  in the same order as the previous R loop. Draws are statistically
  equivalent but not bit-reproducible against the old R loop because the C++
  rejection sampler uses its own RNG stream (compare averages over many
  draws, not individual draws).

* **Faster GLM block sampling:** **`block_rNormalGLM()`** now performs block
  partitioning and prior payload assembly in C++
  (`block_rNormalGLM_cpp_export`), removing per-call R overhead in block
  Gibbs loops (e.g. Block 1 of the **lmebayes** two-block sampler). The
  sampling algorithm itself is unchanged: each block still calls the existing
  `rNormalGLM()` envelope sampler serially. Posterior modes are numerically
  identical to the previous R-prep path; individual draws follow the same
  distribution but are not bit-reproducible against the old path (compare
  means over longer runs). Present-but-`NULL` prior elements (e.g.
  `dispersion = NULL`) are treated as absent, matching R `is.null()`
  semantics.

# glmbayes 0.9.6

## Highlights

* **Multi-response `lmb()`:** **`lmb()`** now handles both univariate and
  multivariate responses with a single unified interface, mirroring the behaviour
  of R's **`lm()`**. When the response has a single column the result is an
  **`lmb`** object (unchanged from prior releases). When the formula specifies
  multiple response columns (e.g. `cbind(y1, y2) ~ x`), **`lmb()`** fits a
  separate Bayesian linear model per response column and returns a named list
  with class **`mlmb`**. For the multi-response case, **`pfamily`** must be a
  list of **`pfamily`** objects with exactly one entry per response column;
  passing a single **`pfamily`** object is an error. Summary, print, and
  coefficient methods for **`mlmb`** objects are included.

* **Conjugate GLM priors (Poisson, binomial, Gamma):** New closed-form IID
  sampling paths for intercept-only models with identity links. **`dBeta()`**
  with **`rBeta_reg()`** supports Beta–Binomial(identity) conjugate updates;
  **`dGamma(Inv_Dispersion = FALSE)`** with **`rGamma_Conjugate_reg()`**
  supports Gamma–Poisson(identity) and Gamma–Gamma(identity) rate priors.
  **`Prior_Setup()`** can calibrate conjugate hyperparameters for these
  families (weighted Poisson rate and binomial probability defaults). See
  **`?dBeta`**, **`?dGamma`**, and the Chapter 02 / Chapter 07–11 vignettes.

* **Vignette structure:** Reworked **Chapter 00** as a roadmap across five
  main parts plus technical appendices. **Chapter 02** is now a conceptual
  introduction to single-parameter conjugacy; worked examples move to
  **Chapter 02-S01** through **Chapter 02-S05** (Beta–Binomial, Normal–Normal,
  Gamma–Poisson, exposure-weighted Poisson, and related topics). A **Companion
  textbooks** section in Chapter 00 indexes optional Bayes Rules! and `LearnBayes`
  appendices tied to the main GLM chapters.

* **`opencltools` import:** Core host/runtime OpenCL discovery and diagnostics
  (`detect_*`, PATH helpers, environment checks) now live in the **`opencltools`**
  package (`Imports`, >= 0.8.0). **glmbayes** keeps package-specific entry
  points (`glmbayesCore_has_opencl()`, `diagnose_glmbayes()`) that report compile-time
  OpenCL status for this build while delegating shared GPU/runtime checks—reducing
  duplicated maintenance in **glmbayes**.

* **Bayes Rules! companion examples:** Optional vignette appendices reproduce
  book datasets and published posterior summaries using **`lmb()`**, **`glmb()`**,
  **`Prior_Setup()`**, and **`dNormal()`** (suggested package **`bayesrules`** for
  data only). Coverage includes **`bikes`** (Ch. 03), **`weather_perth`** (Ch. 08–09),
  **`equality_index`** (Ch. 10), Gamma–Poisson conjugacy (Ch. 02-S04), and a
  scope note for Gamma regression (Ch. 11). Comparison tables use **printed book
  values**, not live **`rstanarm`** fits. See **Chapter 00** § Companion textbooks.

* **`LearnBayes` examples:** **Chapter 02-S04**, Appendix A, maps the
  **`hearttransplants`** example from Albert (2009) / `LearnBayes` (exposure-weighted
  Gamma–Poisson conjugacy) to **`glmb()`** with analytic Albert posteriors for
  verification (suggested package **`LearnBayes`**).

## Other changes

* Expanded **testthat** coverage for **`dBeta()`** / binomial(identity) conjugate
  paths and related **`glmb()`** integration.

# glmbayes 0.9.5

* **Tests / CRAN:** All **OpenCL**-specific **testthat** blocks now call
  **`skip_on_cran()`** (in addition to **`skip_if_no_opencl()`**), consistent
  with existing Boston/Cleveland OpenCL tests. OpenCL coverage remains for local
  checks and source builds with OpenCL; CRAN checks avoid parallel/GPU-heavy
  tests that could trigger **CPU time vs elapsed time** NOTES.

# glmbayes 0.9.4

* **Vignettes:** A vignette that previously used the `notangle` engine now
  uses the standard R Markdown vignette machinery (`knitr` /
  `rmarkdown::html_vignette`), so builds align with CRAN expectations and
  vignette index ordering should be consistent with the rest of the package.

* **OpenCL sources (`inst/cl`):** Removed unused or superseded material,
  consolidated kernels and library fragments, and aligned `.cl` layout and
  dependency tagging with the conventions used in 'openclport' and
  'nmathopencl' (prelude, shims, `nmath/` stems, family kernels under
  `src/`). See `inst/cl/README.md` for how the assembled program is stitched.

* **OpenCL program assembly:** Reworked loading so the full OpenCL program is
  built from explicit fragments (global header, `nmath` closure, family/link
  kernels) rather than ad hoc concatenation—clearer ownership of what enters
  GPU compilation and easier parity with CPU paths.

* **Tests:** Added and expanded **testthat** coverage aimed at OpenCL code
  paths (including binomial examples that exercise GPU envelope evaluation),
  complementing existing Cleveland-style checks.

* **Bug fix — binomial OpenCL:** Binomial `f2_f3` OpenCL kernels now evaluate
  the data log-likelihood with the same **proportion × trial-count**
  semantics as **`dbinom_glmb`** on the CPU (`round` successes and trials,
  clamped probability). This fixes envelope / PLSD failures for aggregated
  binomial data (e.g. `cbind(successes, failures)` / `MASS::menarche`) where
  the previous kernels treated **`y`** like a raw success count.

# glmbayes 0.9.3

* Published on CRAN.
* Version bump in response to CRAN resubmission feedback.

# glmbayes 0.9.2

* Version bump in preparation for resubmission incorporating CRAN review feedback.

# glmbayes 0.9.1

* Wrapped OpenCL-dependent examples in `\donttest{}` for CRAN compliance.
* Reduced iteration counts in rlmb Gibbs sampler example to stay within
  CRAN example time limits on slower check machines.

# glmbayes 0.9.0

First CRAN submission. This release is a stable pre-release with a
near-complete feature set relative to earlier development builds.

## Highlights

### Bayesian Generalized Linear (glmb) and Linear (lmb) modeling functions:

  `glmb()` is a Bayesian analog for the classical `glm()` function while
  `lmb()` covers Gaussian models. Calls largely mirror those for the 
  classical functions but leverage pfamilies for prior specifications.
  Method functions largely mirror those for the classical functions. 
  Samples generated by the functions are largely iid samples 
  (no MCMC convergence dignostics are needed).

### Implemented Likelihood families/ link functions:
   
  Most of the families implemented in the `glm()` function are also implemented 
  in the `glmb()` function (the `lmb()` function covers only gaussian() families). 
  Link functions that lead to log-concave likelihood functions are generally 
  implemented.  Specifically, we have the following:
  
  **Supported likelihoods:** gaussian (identity), Poisson / quasi-Poisson
  (log), binomial / quasi-binomial (logit, probit, cloglog), Gamma (log).

### Prior Family functions:

 `pfamily` constructors are used to specify priors and play the same
  kind of role for the prior specifications as `family` constructors 
  and `link` functions play for the likelihoods. Specifically, we
  have the following:

  **Supported Priors:** Normal (all families/links), Normal–Gamma and 
  independent Normal–Gamma (gaussian families), and Gamma-on-precision 
  (gaussian and Gamma families).
  
### Prior_Setup function:
 
  The package comes with a convenient `Prior_Setup()` function that provides 
  default prior input parameters for each of the implemented models. Basic calls
  (without tailoring) mirror traditional calls to the `glmb()` and `lmb()`
  functions respectively and only require the user to provide the model formula
  and (if not the gaussian family) the family/link function. 
  
  The function can also be used to easily adjust prior specifications 
  (see documentation for details).
  
### Extensive Method functions:
  
  The package comes with extensive method functions that mirror those 
  for the classical functions.  These include dedicated `print()`,
  `summary()`, `predict()` and `simulate()` functions.

### Lower Level Modeling functions:

  The package comes with lower level modeling/simulation functions
  that advanced users can use to implement block Gibbs samplers. These
  generally come with less overhead than the `glmb()` and `lmb()` functions 
  and are called internally by the the higher level modeling functions.

### RcppParallel and OpenCL GPU Acceleration Implementations
  
  Some of the simulation functions comes with use_parallel and use_opencl options
  that speed up simulation for higher dimensional models.
  
### Extensive help files, vignettes, examples and demos

  The package also comes with extensive help files for the varios functions 
  that are complemented with a rich set of vignettes. A large number of 
  examples and demos are also availabel (see the READM.md file for a sample).

---

## Earlier development history (0.1.x series)

The notes below summarize major work during the initial development series
before the 0.9.0 pre-release.

### OpenCL and GPU acceleration

- Completed the OpenCL-based grid construction framework for large models.
- Added GPU-aware envelope sizing and improved OpenCL failure handling.
- Introduced diagnostic utilities to assess OpenCL availability and
  performance.
- Improved configure scripts to detect OpenCL and provide informative
  messages.
- Expanded OpenCL documentation and added a dedicated vignette chapter.

### Parallel CPU sampling (RcppParallel)

- Enabled parallel envelope construction and parallel iid sampling.
- Added pilot functions for large-dimension grid estimation.
- Implemented thread-safe parallel sampling for independent normal-gamma
  models.

### Core statistical improvements

- Migrated to an improved independent normal-gamma simulation algorithm.
- Added theoretical derivations for independent normal-gamma regression.
- Improved UB2 and RSS minimization routines, including scaling corrections.
- Enhanced `Prior_Setup()` to support family-specific prior construction.
- Added dedicated envelope evaluation and sizing functions.

### Package infrastructure

- Significant cleanup to remove NOTES and improve CRAN readiness.
- Improved configure and Makevars files for portability.
- Added testthat tests, including OpenCL-specific tests.
- Consolidated envelope-building functions into a cleaner structure.

### Documentation

- Major updates to README and package-level documentation.
- Added multiple new vignettes and expanded existing ones.
- Improved examples for `lmb()`, `rlmb()`, and OpenCL models.

### Bug fixes (0.1.x era)

- Corrected scaling in UB2 minimization.
- Improved error handling for missing OpenCL functionality.
- Fixed various small issues uncovered during parallelization work.
