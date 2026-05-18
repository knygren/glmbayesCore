# CRAN submission comments — glmbayes 0.9.4

## Summary

Maintenance release. Most changes are internal (OpenCL source layout and program
assembly, extra **testthat** coverage where OpenCL is enabled at compile time, one
vignette switched to standard **R Markdown** machinery). **CRAN’s default
builds remain CPU-only** (no OpenCL), so typical platform checks exercise the same
code paths as before for end users installing binaries.

Optional GPU/OpenCL paths are unchanged from a packaging perspective: still
conditional on building from source with OpenCL available.

## Test environments

- Local (**Windows**): `R CMD check` — OK  
- **Win-builder** (devel, release, and oldrel): OK  
- **macOS builder**: OK  
- **R-hub**: OK  
- **R-Universe**: OK
- OpenCL build/run smoke-tested on Linux + NVIDIA (Vast.ai): OK

---
_This file is listed in `.Rbuildignore` and is not included in the built source
tarball. When submitting, paste the content above into the “Optional comments”
field on the CRAN submission form at_
https://cran.r-project.org/submit.html
