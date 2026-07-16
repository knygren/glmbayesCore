# Archived OpenCL program loaders (not compiled)

Reference copies of superseded GPU program assembly code. These files are
**not** linked into the glmbayesCore shared library (R only compiles `*.cpp`
directly under `src/`, not under subdirectories).

## `kernel_loader_fat_pre_opencltools.cpp`

Pre-opencltools / in-tree loader: prelude, shims, and selective `nmath/` were
all read from **glmbayesCore** `inst/cl/` via `system.file` and local TSV
parsing. Replaced by the production loader in `src/kernel_loader.cpp`, which
delegates to the **opencltools** C API and loads prelude/nmath from
**nmathopencl** (same assembly recipe as **glmbayes**).

Archived: 2026-07-16.
