# Headers for `LinkingTo: glmbayesCore`

Installed as `include/glmbayesCore/` when **glmbayesCore** is installed.
Downstream packages that call the glmbayesCore samplers from C/C++ use the
**C API**:

```cpp
#include <glmbayesCore/glmbayesCore_capi.h>

SEXP out = glmbayesCore_rNormalReg(n, y, x, mu, P, offset, wt, dispersion,
                                   f2, f3, start, family, link, Gridtype);
```

## Requirements

- **`Imports: glmbayesCore`** -- the DLL must be loaded; its `.onLoad`
  registers the entry points.
- **`LinkingTo: glmbayesCore`** -- supplies `-I` for this header only.
- **No `PKG_LIBS`** (symbols resolve via `R_GetCCallable`).
- **`glmbayesCore_api_version()`** -- call for ABI compatibility checks
  (currently `1`).

## Entry points

| C API | Wraps |
|---|---|
| `glmbayesCore_rNormalGLMGroups` | `glmbayes::sim::group::rNormalGLMGroups` |
| `glmbayesCore_rNormalReg` | `glmbayes::sim::rNormalReg` |

Each takes the wrapped function's arguments as `SEXP`, in the same order, and
returns the same list. Errors are raised as R errors.

## Adding an entry point

Add the `extern "C"` wrapper and its `R_RegisterCCallable` line in
`src/glmbayesCore_ccallables.cpp`, then the typedef, accessor, and wrapper in
`glmbayesCore_capi.h`.
