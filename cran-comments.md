# CRAN submission comments — glmbayes 0.9.0

## Package summary

glmbayes provides iid sampling for Bayesian Generalized Linear Models
(Gaussian, Poisson, Binomial, Gamma) via accept-reject methods based on
likelihood subgradients (Nygren & Nygren, 2006). It mirrors the interface
of base R's glm() and lm(), and optionally accelerates envelope
construction via OpenCL for high-dimensional models. OpenCL is an optional
capability; the package detects its absence at build time and disables that
code path gracefully — all checks pass on platforms without OpenCL. See 
README.md and NEWS.md for additional package details.


## Test environments

### Win-builder

- R release 
    -R 4.6.0 (2026-04-24 ucrt)
    -Platform: x86_64-w64-mingw32 (Windows Server 2022 x64)
    -Toolchain: GCC 14.3.0 (C/C++) and GNU Fortran 14.3.0
    -Rcpp 1.1.1-1.1
    -0 errors, 0 warnings, 1 note (New submission)

    
- R-oldrelease 
    -R 4.5.3 (2026-03-11 ucrt)
    -Platform: x86_64-w64-mingw32 (Windows Server 2022 x64)
    -Toolchain: GCC 14.3.0 (C/C++) and GNU Fortran 14.3.0
    -Rcpp 1.1.1
    -0 errors, 0 warnings, 2 notes (New submission; long-running non-OpenCL example timing)

  1. New submission (see Notes)
  2. Long-running non-OpenCL (see Notes - oldrelease only)

A version of package prior to a change to the Rcpp Include/Suggest setup ran as below

- R-devel   
    -R Under development (unstable) (2026-04-25 r89962 ucrt)
    -Platform: x86_64-w64-mingw32 (Windows Server 2022 x64)
    -Toolchain: GCC 14.3.0 (C/C++) and GNU Fortran 14.3.0
    -Rcpp 1.1.1-1.1
    -0 errors, 0 warnings, 3 notes (New submission, )

where the additional notes were

a) Package listed in more than one of Depends, Imports, Suggests, Enhances

b) Skipping checking math rendering: package 'V8' unavailable

A rerun after removing these currently runs into an ERROR: failed to lock directory 'd:/RCompile/CRANguest/R-devel/lib' for modifying
(likely due to re-run being too close to the previous one).

### Mac-builder

- macOS release (mac.R-project.org)
    - Build system: `r-release-macosx-arm64|4.6.0|macosx|macOS 26.2 (25C56)|Mac mini|Apple M1||en_US.UTF-8|macOS 14.4|clang-1700.6.3.2|GNU Fortran (GCC) 14.2.0`
    - R Under development (unstable) (2026-03-22 r89674)
    - Platform: aarch64-apple-darwin23
    - Toolchain: Apple clang 17.0.0 (clang-1700.6.3.2), GNU Fortran (GCC) 14.2.0
    - Rcpp 1.1.1-1.1 
    - Status: OK (macbuilder reports final status instead of standard E/W/N summary)

- macOS devel (mac.R-project.org)
    - Build system: `r-release-macosx-arm64|4.6.0|macosx|macOS 26.2 (25C56)|Mac mini|Apple M1||en_US.UTF-8|macOS 14.4|clang-1700.6.3.2|GNU Fortran (GCC) 14.2.0`
    - R Under development (unstable) (2026-03-22 r89674)
    - Platform: aarch64-apple-darwin23
    - Toolchain: Apple clang 17.0.0 (clang-1700.6.3.2), GNU Fortran (GCC) 14.2.0
    - Rcpp 1.1.1-1.1
    - Status: OK (macbuilder reports final status instead of standard E/W/N summary)

Both macOS release and devel submissions currently report the same effective build profile
on arm64 and now finish with `Status: OK`.

### R-universe
- R-universe: all non-wasm platforms pass with 0 errors, 0 warnings, and 0 notes.
- wasm (WebAssembly) remains expected to fail because the package includes compiled
  C/C++ code that is not compatible with the wasm toolchain.

### rhub (via rhub::rhub_check())

**Platforms with regular checks:**

| Platform              | R version (svn)   | Rcpp version | E/W/N        |
|-----------------------|-------------------|--------------|--------------|
| atlas                 | R 4.7.0 (r89961)  | 1.1.1-1.1    | OK           |
| c23*                  | R 4.6.0 (r89623)  | 1.1.1        | OK           |
| clang16*              | R 4.6.0 (r89629)  | 1.1.1        | OK           |
| clang17*              | R 4.6.0 (r89629)  | 1.1.1        | 1 NOTE       |
| clang18*              | R 4.6.0 (r89623)  | 1.1.1        | 1 NOTE       |
| clang19*              | R 4.6.0 (r89629)  | 1.1.1        | 1 NOTE       |
| clang20*              | R 4.6.0 (r89623)  | 1.1.1        | 1 NOTE       |
| clang21               | R 4.7.0 (r89961)  | 1.1.1-1.1    | 1 NOTE       |
| clang22               | R 4.7.0 (r89961)  | 1.1.1-1.1    | OK           |
| donttest              | R 4.7.0 (r89961)  | 1.1.1-1.1    | OK           |
| gcc13*                | R 4.6.0 (r89629)  | 1.1.1        | OK           |
| gcc14*                | R 4.6.0 (r89629)  | 1.1.1        | OK           |
| gcc15*                | R 4.6.0 (r89629)  | 1.1.1        | OK           |
| gcc16                 | R 4.7.0 (r89961)  | 1.1.1-1.1    | OK           |
| intel*                | R 4.6.0 (r89439)  | 1.1.1        | OK           |
| linux (R-devel)       | R 4.7.0 (r89961)  | 1.1.1-1.1    | OK           |
| lto                   | R 4.6.0 (r89956)  | 1.1.1.1      | OK           |
| m1-san (R-devel)      | R 4.6.0 (r89961)  | 1.1.1-1.1    | OK           |
| macos-arm64 (R-devel) | R 4.6.0 (r89961)  | 1.1.1-1.1    | OK           |
| mkl                   | R 4.7.0 (r89955)  | 1.1.1-1.1    | OK           |
| nold                  | R 4.7.0 (r89961)  | 1.1.1-1.1    | OK           |
| noremap*              | R 4.6.0 (r89623)  | 1.1.1        | OK           |
| ubuntu-clang          | R 4.7.0 (r89961)  | 1.1.1-1.1    | OK           |
| ubuntu-gcc12          | R 4.7.0 (r89874)  | 1.1.1-1.1    | OK           |
| ubuntu-next           | R 4.6.0 (r89961)  | 1.1.1.1      | OK           |
| ubuntu-release        | R 4.6.0 (r89956)  | 1.1.1.1      | OK           |
| windows (R-devel)     | R 4.7.0 (r89962)  | 1.1.1-1.1    | OK           |

`*` Platforms where R/Rcpp version inconsistencies prevent installation of 
Rcpp 1.1.1-1 or later. Rcpp 1.1.1 installs correctly on these platforms if custom installed. 
The boundary appears to be r89746 (i.e., R 4.6.0 below r89746 requires Rcpp 1.1.1 
instead of Rcpp 1.1.1-1 or later). If the R versions on these systems get updated to the 
release version or later, these should migrate to the latest CRAN Rcpp version.

**Platforms with special checks:**

| Platform    | R version (svn)   | Rcpp version | E/W/N        |
|-------------|-------------------|--------------|--------------|
| clang-asan  | R 4.7.0 (r89961)  | 1.1.1-1.1    | 1 NOTE       |
| clang-ubsan | R 4.7.0 (r89961)  | 1.1.1-1.1    | 1 NOTE       |
| gcc-asan    | R 4.7.0 (r89961)  | 1.1.1-1.1    | OK           |
| valgrind    | R 4.7.0 (r89961)  | 1.1.1-1.1    | OK           |

- Remaining NOTE on some clang-based rhub platforms is environment/toolchain-provided:
  non-portable compile flag `-Wp,-D_FORTIFY_SOURCE=3` reported by `R CMD check`.
  This flag is injected by the build image, not by glmbayes Makevars.

- Sanitizer/valgrind diagnostic summary for special-check platforms:
  - `clang-asan` and `gcc-asan` focus on memory safety (out-of-bounds access, use-after-free,
    double/invalid free); no sanitizer findings attributable to package code were reported.
  - `clang-ubsan` focuses on undefined behavior (e.g., invalid shifts/casts, misalignment,
    overflow-related UB); no UBSAN findings attributable to package code were reported.
  - `valgrind` focuses on runtime memory errors/leaks; no invalid read/write or leak diagnostics
    attributable to package code were reported.
  - These diagnostics apply to code paths exercised during the rhub check workload
    (examples/tests/vignettes run on those platforms).

- rchk: [describe outcome and explain here]


## Comments Related to Notes appearing on various systems

All checks produced 0 errors and 0 warnings. The following 3 notes were
observed on select systems.

### Note: **New submission** 

       Maintainer: 'Kjell Nygren <kjell.a.nygren@gmail.com>'
       New submission

   Expected for an initial CRAN submission. No action required.

### Note: **Non-OpenCL Examples with long CPU or elapsed time**

       Examples with CPU (user + system) or elapsed time > 5s
                user  system elapsed
       rlmb    12.60    0.45   10.61

This appears only on select platforms/machines. On many, this note is never
triggered as elapsed time falls below the 5-second threshold. Example compares 
a two-block gibbs sampler (MCMC) to the main lmb implementation (MC). Sufficently long 
runs and needed to make comparisons valid statistically.


### Note on rchk
[rchk checks for PROTECT issues in C code. Describe what rchk flagged,
whether it is a false positive, and what you did to investigate or
mitigate it. If the flag is in Rcpp-generated code rather than your
own C, say so explicitly.]

## GPU/OpenCL test environments

### Local (developer machine)

#### Local release/near-release (`R 4.6.0 RC`)
- Environment: Windows 11 x64 (build 26200), ASUS TUF F16, NVIDIA GeForce RTX GPU, OpenCL available
- R version: 4.6.0 RC (2026-04-22 r89945 ucrt), platform `x86_64-w64-mingw32`
- Toolchain:
  - R compiled by `gcc.exe (GCC) 14.3.0`, `GNU Fortran (GCC) 14.3.0`
  - Package install/check used C compiler `gcc.exe (GCC) 14.2.0`, C++ compiler `G++ (GCC) 14.2.0`
- Command: `devtools::check(vignettes = TRUE, args = "--as-cran", remote = TRUE, manual = TRUE)`
- Result: `0 errors | 0 warnings | 2 notes`
  1. New submission (see Notes above)
  2. Long-running OpenCL examples (see Notes below)
- OpenCL/GPU context: this run uses the local OpenCL-capable NVIDIA GPU (`has_opencl()` true), so the GPU code path is exercised.
- Example timing note (from `R CMD check`):
  - `Boston_centered`: user 51.22s, system 4.64s, elapsed 31.44s
  - `Cleveland`: user 11.09s, system 3.25s, elapsed 7.92s

#### Local oldrel (`R 4.5.3`)
- Environment: Windows 11 x64 (build 26200), ASUS TUF F16, NVIDIA GeForce RTX GPU, OpenCL available
- R version: 4.5.3 (2026-03-11 ucrt), platform `x86_64-w64-mingw32`
- Toolchain:
  - R compiled by `gcc.exe (GCC) 14.3.0`, `GNU Fortran (GCC) 14.3.0`
  - Package install/check used C compiler `gcc.exe (GCC) 14.2.0`, C++ compiler `G++ (GCC) 14.2.0`
- Command: `devtools::check(vignettes = TRUE, args = "--as-cran", remote = TRUE, manual = TRUE)`
- Result: `0 errors | 0 warnings | 2 notes`
  1. New submission
  2. Long-running OpenCL examples
- OpenCL/GPU context: this run uses the local OpenCL-capable NVIDIA GPU (`has_opencl()` true), so the GPU code path is exercised.
- Example timing note (from `R CMD check`):
  - `Boston_centered`: user 78.87s, system 7.53s, elapsed 59.11s
  - `Cleveland`: user 16.86s, system 2.12s, elapsed 10.47s

### GPU / OpenCL on Linux (Vast.ai virtual machine)
- Environment: Ubuntu 22.04 (gcc/g++ 11.4.0), NVIDIA GeForce RTX 5060 Ti, OpenCL runtime detected

#### Vast.ai oldrel (`R 4.5.3`)
- R version: 4.5.3 (svn r89597), invoked explicitly via `/opt/R/4.5.3/bin/Rscript`
- Install sources: glmbayes from r-universe (`https://knygren.r-universe.dev/src/contrib/...`), dependencies from CRAN source tarballs
- Dependencies built from source in this run: `rbibutils`, `coda`, `Rcpp`, `RcppParallel`, `Rdpack`, `RcppArmadillo`
- Rcpp: `1.1.1-1.1` (configure reports normalized `1.1.1.1.1`), simulated `Function.h` branch = 2 (`R_VERSION < 4.6.0 || R_SVN_REVISION < 89746`)
- OpenCL detection in configure:
  - headers found in `/usr/include`
  - OpenCL library found in `/usr/lib/x86_64-linux-gnu`
  - runtime probe succeeded
- Build/install result: source install completed successfully (`* DONE (glmbayes)`), Step 1 elapsed `301.18s`
- Repos used by install script in this run: `https://knygren.r-universe.dev`, `https://cloud.r-project.org`
- OpenCL examples (`run_opencl_examples()`) result: `OK`, Step 2 elapsed `80.06s`
  - examples exceeding 5s:
    - `example:Boston_centered`: user 138s, system 19.2s, elapsed 67s
    - `example:Cleveland`: user 33s, system 3.9s, elapsed 13s
- Combined remote run timing (install + examples): `381.24s`

#### Vast.ai release (`R 4.6.0`)
- R version: 4.6.0 (svn r89956), invoked explicitly via `/opt/R/4.6.0/bin/Rscript`
- Install sources: glmbayes from r-universe (`https://knygren.r-universe.dev/src/contrib/...`), dependencies from CRAN source tarballs
- Dependencies built from source in this run: `rbibutils`, `coda`, `Rcpp`, `RcppParallel`, `Rdpack`, `RcppArmadillo`
- Rcpp: `1.1.1-1.1` (configure reports normalized `1.1.1.1.1`), simulated `Function.h` branch = 3 (`R_VERSION < 4.6.0 || R_SVN_REVISION < 89746` is false)
- OpenCL detection in configure:
  - headers found in `/usr/include`
  - OpenCL library found in `/usr/lib/x86_64-linux-gnu`
  - runtime probe succeeded
- Build/install result: source install completed successfully (`* DONE (glmbayes)`), Step 1 elapsed `341.55s`
- Repos used by install script in this run: `https://knygren.r-universe.dev`, `https://cloud.r-project.org`
- OpenCL examples (`run_opencl_examples()`) result: `OK`, Step 2 elapsed `79.40s`
  - examples exceeding 5s:
    - `example:Boston_centered`: user 140s, system 19.5s, elapsed 67s
    - `example:Cleveland`: user 24s, system 3.9s, elapsed 12s
- Combined remote run timing (install + examples): `420.97s`

- These VM runs validate source build + install and OpenCL runtime execution on Linux; no `R CMD check` run on these VMs.

### Note: **OpenCL Examples with long CPU or elapsed time**

       Examples with CPU (user + system) or elapsed time > 5s
                        user  system elapsed
       Boston_centered 150.89  16.16  105.20
       Cleveland        42.25   3.00   29.34

`Boston_centered` and `Cleveland` are GPU/OpenCL examples: the heavy path runs only when
`has_opencl()` is true. On typical CRAN check machines there is no usable OpenCL device, so
that path is skipped and these timings are **not** what CRAN sees. The numbers above are
from an OpenCL-capable maintainer machine and are included to document behavior and
performance validation.

**Timing columns:** `elapsed` is wall-clock time (what CRAN’s example-timing checks care about).
`user` and `system` are CPU time for the R process; on parallel workloads (host threads,
OpenMP/BLAS, OpenCL/device overlap), CPU time can be summed across cores so **user + system
can exceed elapsed** even though wall-clock time is lower.

**Why this matters for users:** The main motivation for the OpenCL path is speed on **large**
problems (many covariates and/or observations). **Wall-clock gains from parallel GPU/OpenCL
execution typically grow with problem size**—often dramatically—because work maps well to the
device and amortizes setup. The examples above use **modest** data sizes so they stay
runnable in documentation; they **understate** the speedup you should expect on realistically
large models, where the advantage over a serial CPU-only path is likely to be **much**
larger.

---
_This file is listed in `.Rbuildignore` and is not included in the built
source tarball. When submitting, paste the content above into the
"Optional comments" field on the CRAN submission form at
https://cran.r-project.org/submit.html