# CRAN submission comments — glmbayesCore 0.5.3

## Summary

This is a resubmission of the new package glmbayesCore (0.5.3).

Changes since 0.5.2:

* Fixed Windows linking against **RcppParallel** / TBB in `configure.win`
  (`-DRCPP_PARALLEL_USE_TBB=1` and `RcppParallel::RcppParallelLibs()`),
  addressing undefined reference errors to `tbb::detail::r1::*` on current
  Windows toolchains.

Prior reviewer feedback (0.5.2) remains addressed:

* `\value` for `diagnose_glmbayes()` / `glmbayesCore_has_opencl()`
* Ungated `print()` / `cat()` replaced with `message()` / `warning()` /
  S3 `print.diagnose_glmbayes()`

## Test environments

* local Windows, `R CMD check --as-cran`: (update after local check)
* win-builder / R-universe Windows: (update after rebuild confirms TBB link)

---
_This file is listed in `.Rbuildignore` and is not included in the built source
tarball. When submitting, paste the content above into the “Optional comments”
field on the CRAN submission form at_
https://cran.r-project.org/submit.html
