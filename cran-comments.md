# CRAN submission comments — glmbayes 0.9.0

## Package summary

glmbayes provides iid sampling for Bayesian Generalized Linear Models
(Gaussian, Poisson, Binomial, Gamma) via accept-reject methods based on
likelihood subgradients (Nygren & Nygren, 2006). It mirrors the interface
of base R's glm() and lm(), and optionally accelerates envelope
construction via OpenCL for high-dimensional models. OpenCL is an optional
capability; the package detects its absence at build time and disables that
code path gracefully — all checks pass on platforms without OpenCL.

## Test environments

### Local (developer machine)
- Windows 11, ASUS TUF F16, GeForce RTX GPU, OpenCL installed
- R 4.5.1, glmbayes built with OpenCL enabled
- Command: `devtools::check(vignettes = TRUE, args = "--as-cran", remote = TRUE, manual = TRUE)`
- 0 errors, 0 warnings, 3 notes
  1. New submission (see Notes)
  2. Rcpp workaround (see Notes)
  3. Long-running examples on OpenCL-enabled machine (see Notes)
   
### Win-builder
- R-devel (4.6.0 RC, 2026-04-20 r89921 ucrt):    0 errors, 0 warnings, 2 notes
- R-release (4.5.3, 2026-03-11 ucrt):             0 errors, 0 warnings, 2 notes
- R-oldrelease (4.4.3 patched, 2026-02-12 ucrt):  0 errors, 0 warnings, 3 notes

The additional note on R-oldrelease relates to the Rcpp version workaround
(see Notes).

### Mac-builder
- macOS release (mac.R-project.org): 0 errors, 0 warnings, N notes
- macOS devel  (mac.R-project.org): 0 errors, 0 warnings, N notes

### R-universe
- All platforms pass except wasm (WebAssembly), which is expected:
  the package includes compiled C/C++ code that is not compatible
  with the wasm toolchain.

### rhub (via rhub::rhub_check())
- linux, macos-arm64, windows, m1-san, atlas, c23,
  clang16–clang22, gcc13–gcc16, intel, lto, mkl,
  nold, noremap, ubuntu-clang, ubuntu-gcc12,
  ubuntu-release, donttest:
  0 errors, 0 warnings, N notes
  [Note: Rcpp was pre-installed manually on some rhub platforms —
  see Rcpp note below]
- valgrind, clang-asan, clang-ubsan, gcc-asan:
  0 errors, 0 warnings, N notes
- rchk: [describe outcome and explain here]

### GPU / OpenCL on Linux (Vast.ai virtual machine)
- Ubuntu [version], OpenCL enabled, R [version]
- Confirms OpenCL code path builds and runs correctly outside Windows
- Result: 0 errors, 0 warnings, N notes


## Notes

All checks produced 0 errors and 0 warnings. The following 3 notes were
observed on the local Windows machine (R 4.5.3, OpenCL enabled):

1. **New submission**

       Maintainer: 'Kjell Nygren <kjell.a.nygren@gmail.com>'
       New submission

   Expected for an initial CRAN submission. No action required.

2. **Rcpp listed in more than one field**

       Package listed in more than one of Depends, Imports, Suggests, Enhances:
         'Rcpp'
       A package should be listed in only one of these fields.

   Rcpp 1.1.1 (current on R-release and R-oldrelease) and Rcpp 1.1.1-1
   (current on R-devel) differ in a header that glmbayes links against.
   `Imports: Rcpp (>= 1.1.1)` ensures compatibility on release platforms;
   `Suggests: Rcpp (>= 1.1.1-1)` captures the R-devel requirement without
   breaking release builds. This is a temporary workaround for an upstream
   Rcpp incompatibility across the current R release/devel boundary.

3. **Examples with long CPU or elapsed time**

       Examples with CPU (user + system) or elapsed time > 5s
                        user  system elapsed
       Boston_centered 150.89  16.16  105.20
       Cleveland        42.25   3.00   29.34
       rlmb              8.57   1.79    8.30

   Boston_centered and Cleveland are GPU/OpenCL examples guarded by
   `has_opencl()` and do not execute on machines without OpenCL installed.
   They will not appear on CRAN check servers. The `rlmb` example reflects
   genuine sampling time for a Bayesian linear model and cannot be
   meaningfully reduced without undermining the demonstration.
   


### Note 1: New submission
This is the first submission of glmbayes to CRAN.

### Note 2: [Rcpp versioning note]
[Your explanation here — e.g. whether this is a known upstream issue,
whether it affects functionality, and any relevant Rcpp version details]

### Note 3: [New note you mentioned]
[Explanation here]

### Note on rchk
[rchk checks for PROTECT issues in C code. Describe what rchk flagged,
whether it is a false positive, and what you did to investigate or
mitigate it. If the flag is in Rcpp-generated code rather than your
own C, say so explicitly.]

---
_This file is listed in `.Rbuildignore` and is not included in the built
source tarball. When submitting, paste the content above into the
"Optional comments" field on the CRAN submission form at
https://cran.r-project.org/submit.html._