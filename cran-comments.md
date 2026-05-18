# CRAN submission comments — glmbayes 0.9.5

## Summary

Resubmission after CRAN feedback. **OpenCL-specific tests** now use **`skip_on_cran()`**
(they already used **`skip_if_no_opencl()`**). CRAN’s default **CPU-only** builds are
unchanged for end users; this release only avoids running parallel/OpenCL **testthat**
code on CRAN (which could produce **CPU time vs elapsed time** NOTES on some
configurations). OpenCL tests still run off CRAN when the package is built with
OpenCL support.

Other content relative to **0.9.4** remains as before (OpenCL layout, vignette
machinery, binomial OpenCL fix, expanded coverage)—see **NEWS.md**.

## Test environments

- Local (**Windows**): `R CMD check --as-cran` — OK  
- **Win-builder** (devel, release, and oldrel): OK *(re-run after changes)*  
- **macOS builder**: OK *(re-run after changes)*  
- **R-hub**: OK *(optional re-run)*  
- **R-Universe**: OK  
- OpenCL build/run smoke-tested on Linux + NVIDIA (Vast.ai): OK  

---
_This file is listed in `.Rbuildignore` and is not included in the built source
tarball. When submitting, paste the content above into the “Optional comments”
field on the CRAN submission form at_
https://cran.r-project.org/submit.html
