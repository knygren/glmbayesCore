# glmbayesCore 0.5.4

## Bug fixes

* **Configure (Linux/macOS):** `-DUSE_OPENCL` is set only when a **non-PoCL**
  OpenCL platform exposes at least one **GPU** device (same policy as
  **glmbayes** 0.9.75). PoCL-only ICD stacks (typical on CRAN debian-gcc
  incoming) no longer enable OpenCL at compile time, so kernel builds do not
  touch `~/.cache/pocl` (avoids incoming NOTE on `tempfile_*` there). Real GPU
  installs (NVIDIA/AMD/Intel) are unchanged. `configure.win` is unaffected
  (Windows never ran a runtime OpenCL probe).

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
