# `R/` function inventory (index)

Maintainer-only index (under `data-raw/inventories/`; not shipped in the
package tarball) for symbols under **`R/`** in **glmbayesCore** today
(iid GLM/LM envelope engine). Mixed-model APIs are temporarily developed in
**lmebayesCore** and will be reintegrated here gradually.

| Document | Contents |
|----------|----------|
| **[R_EXPORTED_AND_DOCUMENTED.md](R_EXPORTED_AND_DOCUMENTED.md)** | `NAMESPACE` exports grouped by **glmbayes** overlap and Core-only |
| **[R_CORE_ONLY_EXPORTS.md](R_CORE_ONLY_EXPORTS.md)** | Same exports by function type (simulation, envelopes, …) |
| **[R_EXPORT_REACHABILITY.md](R_EXPORT_REACHABILITY.md)** | Notes on export roles / phase-out from **glmbayes** |
| **[R_INTERNAL_HELPERS.md](R_INTERNAL_HELPERS.md)** | `@noRd` helpers and C++ glue |

Scratch checks live in `data-raw/` (not run by `test_check`).
