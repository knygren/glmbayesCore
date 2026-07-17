# Export reachability

Notes on how **glmbayesCore** exports are reached today (slim iid tree).

Index: [R_FUNCTION_INVENTORY.md](R_FUNCTION_INVENTORY.md).

---

## Primary paths (this tree)

| Entry | Reaches |
|-------|---------|
| **glmbayes** `glmb()` / `lmb()` | `rglmb` / `rlmb` → `pfamily$simfun` → `rNormal_reg` / … → C++ |
| Direct `rglmb()` / `rlmb()` | Same as above |
| Direct `simfunction` calls | Thin R wrappers → `rcpp_wrappers.R` → C++ |
| Envelope exports | Used by `EnvelopeOrchestrator` and advanced callers; also re-exported historically from **glmbayes** (phase-out candidates) |

## Mixed-model paths (temporary **lmebayesCore**)

`rlmerb` / `rglmerb`, `rGLMM_*` / `rLMM_*`, `two_block_*`, `model_setup`,
`Prior_Setup_lmebayes`, `pfamily_list`, block ING, `build_mu_all`, ICM
helpers. **lmebayes** currently imports these from **lmebayesCore**; they
are planned to return to **glmbayesCore**.

## Phase-out from **glmbayes**

Low-level simulation, envelope, and truncated-distribution exports should
eventually be **glmbayesCore**-only. User-facing retain set:
`Prior_Setup`, `Prior_Check`, `pfamily` constructors, `rglmb` / `rlmb`,
`multi_prior_setup`, `multi_rlmb`, `diagnose_glmbayes`.
