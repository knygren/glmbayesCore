# glmbayesCore 0.5.4

## Bug fixes

* **`residuals.rglmb()` / `residuals.rlmb()` / `residuals.summary.rglmb()`:**
  Fixed the `ysim` argument to substitute for the **observed response** `y`
  (fitted value held fixed at each draw), not the reverse. The previous
  behavior substituted `ysim` for the fitted value while holding `y` fixed,
  which does not correspond to a standard posterior-predictive residual check
  and contradicted this function's own `ysim` ("simulated responses")
  documentation. The corrected convention matches **glmbayes**'s
  `residuals.glmb()` (simulate new response data, recompute residuals against
  the same fitted values, and compare to the actual residuals for outlier
  diagnostics) — the two packages' `ysim` semantics were inconsistent since
  this function was added post-split; they are now aligned. Behavior for the
  default case (`ysim = NULL`) is unchanged.

* **Configure (Linux/macOS):** `-DUSE_OPENCL` is set only when a **non-PoCL**
  OpenCL platform exposes at least one **GPU** device (same policy as
  **glmbayes** 0.9.75). PoCL-only ICD stacks (typical on CRAN debian-gcc
  incoming) no longer enable OpenCL at compile time, so kernel builds do not
  touch `~/.cache/pocl` (avoids incoming NOTE on `tempfile_*` there). Real GPU
  installs (NVIDIA/AMD/Intel) are unchanged. `configure.win` is unaffected
  (Windows never ran a runtime OpenCL probe).

* **Configure policy:** Removed `tools/rcpp_include.R` / `tools/patch_rcpp_function_h.R`
  and the `glmbayes_getRegisteredNamespace` compatibility shim (`src/glmbayes_getRegisteredNamespace.{cpp,h}`)
  from `configure` and `configure.win`. Builds now rely on standard
  **`LinkingTo: Rcpp`** and the CRAN **Rcpp (>= 1.1.1)** requirement instead of
  probing/patching `Rcpp/Function.h` or recommending a GitHub install of Rcpp
  (same policy fix as **glmbayes** 0.9.73/0.9.75).

* **CI:** Replaced `.github/workflows/rhub.yaml` with the current
  **glmbayes** version, removing dead matrix/report logic that referenced the
  now-removed `tools/patch_rcpp_function_h.R` script.

# glmbayesCore 0.5.3

## Bug fixes

* Windows builds now link against the TBB libraries provided by
  **RcppParallel** (`RcppParallel::RcppParallelLibs()` and
  `-DRCPP_PARALLEL_USE_TBB=1` in `configure.win`). This fixes undefined
  reference errors to `tbb::detail::r1::wait_on_address` /
  `notify_by_address_one` seen on current Windows toolchains
  (win-builder / R-universe).

* Cap RcppParallel to 2 threads in the `pfamily`, `Prior_Setup`, and
  `simfuncs` examples (`\dontshow{setThreadOptions(numThreads = 2)}`)
  so CRAN checks do not report “CPU time > 2.5 times elapsed time”.
  Default parallel sampling for users is unchanged.

# glmbayesCore 0.5.2

## Initial CRAN submission

* First CRAN release of **glmbayesCore**. This package provides the compiled
  sampling engine and related “core” R interfaces (prior families, matrix-level
  samplers such as `rglmb` / `rlmb`, envelope construction, and optional OpenCL
  acceleration) that previously lived inside **glmbayes** and are now packaged
  separately for reuse.

* End users of Bayesian GLM/LM modelling should continue to use **glmbayes**
  for the formula interface and S3 methods. **glmbayesCore** is the developer
  / backend layer those packages build on.

* Future work includes pointing **glmbayes** (and other downstream packages)
  at this Core package as their shared sampling engine, so that the iid
  envelope stack is maintained in one place.

* Version **0.5.2** (rather than 0.1.0) reflects that the engine API is already
  substantially complete from its prior life inside **glmbayes**, while leaving
  room for integration-driven releases before a 1.0.0 freeze.
