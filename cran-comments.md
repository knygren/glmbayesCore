# CRAN submission comments — glmbayesCore 0.5.3

## Summary

This is a resubmission of the new package glmbayesCore (0.5.3).

### Changes in response to reviewer comments (0.5.2)

* Added `\value` documentation for `diagnose_glmbayes()` /
  `glmbayesCore_has_opencl()` in `gpu_diagnostics.Rd` (via `@return` in
  roxygen).
* Replaced ungated `print()` / `cat()` console output with `message()` /
  `warning()` in `Prior_Check()`, `rglmb()`, and `rlmb()`.
* Refactored `diagnose_glmbayes()` to return a `"diagnose_glmbayes"` object
  and print the human-readable report via `print.diagnose_glmbayes()`.

### Additional changes in 0.5.3

* Fixed Windows linking against **RcppParallel** / TBB in `configure.win`
  (`-DRCPP_PARALLEL_USE_TBB=1` and `RcppParallel::RcppParallelLibs()`),
  addressing undefined reference errors to `tbb::detail::r1::*` on current
  Windows toolchains.
* Limited RcppParallel to 2 threads in the `pfamily`, `Prior_Setup`, and
  `simfuncs` examples only (`\dontshow{setThreadOptions(numThreads = 2)}`),
  addressing the NOTE “Examples with CPU time > 2.5 times elapsed time”.
  Package default parallelism for users is unchanged.

## Test environments

* local Windows, `R CMD check --as-cran`: 0 errors | 0 warnings | 0 notes
  (aside from the expected “New submission” NOTE)

* win-builder (CRAN): 0 errors | 0 warnings | 0 notes on **r-release**,
  **r-devel**, and **r-oldrel** 

---
_This file is listed in `.Rbuildignore` and is not included in the built source
tarball. When submitting, paste the content above into the “Optional comments”
field on the CRAN submission form at_
https://cran.r-project.org/submit.html
