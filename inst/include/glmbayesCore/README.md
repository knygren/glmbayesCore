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
| `glmbayesCore_group_rNormalGLM` | `glmbayes::sim::group::group_rNormalGLM_cpp_export` |
| `glmbayesCore_group_rNormalReg` | `glmbayes::sim::group::group_rNormalReg_cpp_export` |
| `glmbayesCore_two_block_build_mu_all` | `glmbayes::sim::twoblock::two_block_build_mu_all` |
| `glmbayesCore_two_block_groups_prior_with_tau2` | `glmbayes::sim::twoblock::two_block_groups_prior_with_tau2` |
| `glmbayesCore_two_block_groups_one_chain_from_mu_P` | `glmbayes::sim::twoblock::two_block_groups_one_chain_from_mu_P_impl` |
| `glmbayesCore_fixef_list_from_batch_chain` | `glmbayes::sim::twoblock::fixef_list_from_batch_chain` |
| `glmbayesCore_batch_tau2_chain_row` | `glmbayes::sim::twoblock::batch_tau2_chain_row` |
| `glmbayesCore_batch_b_assign_slice` | `glmbayes::sim::twoblock::batch_b_assign_slice` (in place) |
| `glmbayesCore_batch_b_n_chains` | `glmbayes::sim::twoblock::batch_b_n_chains` |
| `glmbayesCore_ensure_batch_b_dimnames` | `glmbayes::sim::twoblock::ensure_batch_b_dimnames` (in place) |
| `glmbayesCore_batch_iters_ranef_add` | `glmbayes::sim::twoblock::batch_iters_ranef_add` (in place) |
| `glmbayesCore_x_hyper_as_matrix` | `glmbayes::sim::twoblock::x_hyper_as_matrix` |
| `glmbayesCore_rIndepNormalGammaReg` | `glmbayes::sim::rIndepNormalGammaReg` |
| `glmbayesCore_progress_bar` | `glmbayes::progress::progress_bar` |
| `glmbayesCore_progress_bar_finish` | `glmbayes::progress::progress_bar_finish` |

Each takes the wrapped function's arguments as `SEXP`, in the same order, and
returns the same value (`R_NilValue` for functions that return nothing; scalar
results come back as length-1 vectors). The "in place" helpers modify their
first argument, as the wrapped functions do. Errors are raised as R errors.

## Adding an entry point

Add the `extern "C"` wrapper and its `R_RegisterCCallable` line in
`src/glmbayesCore_ccallables.cpp`, then the typedef, accessor, and wrapper in
`glmbayesCore_capi.h`.
