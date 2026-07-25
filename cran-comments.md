# CRAN submission comments — glmbayesCore 0.5.2

## Summary

This is a resubmission of the new package glmbayesCore (0.5.2).

Changes in response to reviewer comments:

* Added `\value` documentation for `diagnose_glmbayes()` /
  `glmbayesCore_has_opencl()` in `gpu_diagnostics.Rd` (via `@return` in
  roxygen).
* Replaced ungated `print()` / `cat()` console output with `message()` /
  `warning()` in `Prior_Check()`, `rglmb()`, and `rlmb()`.
* Refactored `diagnose_glmbayes()` to return a `"diagnose_glmbayes"` object
  and print the human-readable report via `print.diagnose_glmbayes()`.

## Test environments

* local Windows, `R CMD check --as-cran`: 0 errors | 0 warnings | 0 notes
  (aside from the expected “New submission” NOTE)

* win-builder (CRAN): 0 errors | 0 warnings | 0 notes on **r-release**,
  **r-devel**, and **r-oldrel** (aside from the expected “New submission”
  NOTE on each).

---
_This file is listed in `.Rbuildignore` and is not included in the built source
tarball. When submitting, paste the content above into the “Optional comments”
field on the CRAN submission form at_
https://cran.r-project.org/submit.html
