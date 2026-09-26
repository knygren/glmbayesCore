# CRAN submission comments — glmbayesCore 0.5.4

## Summary

This is an update of glmbayesCore from 0.5.3 (currently on CRAN) to 0.5.4.

### Changes in 0.5.4

* **Configure (Linux/macOS):** `-DUSE_OPENCL` is set only when a **non-PoCL**
  OpenCL platform exposes at least one **GPU** device (same policy as
  **glmbayes**), avoiding PoCL cache NOTEs on CRAN debian-gcc.
* **Configure policy:** Removed `tools/rcpp_include.R` /
  `tools/patch_rcpp_function_h.R` and related Function.h / registered-namespace
  probing from `configure` and `configure.win`. Builds rely on standard
  **`LinkingTo: Rcpp`** and **Rcpp (>= 1.1.1)** instead of recommending a
  GitHub install of Rcpp (same CRAN policy fix as **glmbayes**).
* **`residuals.rglmb()` / `residuals.rlmb()` / `residuals.summary.rglmb()`:**
  Fixed `ysim` to substitute for the observed response (fitted values held
  fixed), aligning with **glmbayes** `residuals.glmb()` semantics. Default
  (`ysim = NULL`) behavior is unchanged.

## Test environments

* local Windows, `R CMD check --as-cran`: (update after local check of 0.5.4)
* win-builder (CRAN): (update after 0.5.4 rebuild)

---
_This file is listed in `.Rbuildignore` and is not included in the built source
tarball. When submitting, paste the content above into the “Optional comments”
field on the CRAN submission form at_
https://cran.r-project.org/submit.html
