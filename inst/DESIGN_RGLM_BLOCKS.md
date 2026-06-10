# Design: conditionally independent block sampling (`rNormalGLM`)

Maintainer-facing implementation plan for sampling GLM posterior blocks when
the full conditional factorizes across observation groups. Canonical doc;
excluded from the source tarball via `.Rbuildignore`.

---

## 1. Purpose and scope

### Goal

Speed up and structure **binomial, Poisson, and Gamma** (and quasi variants)
sampling when the model is **conditionally independent across blocks**:

\[
p(\beta_1,\ldots,\beta_k \mid y)
  = \prod_{b=1}^{k} p(\beta_b \mid y_b).
\]

### In scope

- Same number of coefficients per block: **`l1 = ncol(x)`** fixed across blocks.
- Varying observations per block: **`l2_b`**, with **`l2 = sum_b l2_b`**.
- Pipeline per block: existing **`rNormalGLM`** steps in C++ (`optim` →
  `glmb_Standardize_Model` → `EnvelopeBuild` → `rNormalGLM_std` /
  `rNormalGLM_std_parallel`).
- Typical use: **block Gibbs** with **`n = 1`** per block; also **`n > 1`**
  iid draws per block.

### Out of scope (v1)

- **Gaussian** path (`rNormal_reg` / `rNormalReg` / `lm.fit`) — not this doc.
- Joint sampling when blocks are **not** conditionally independent (coupled
  prior, shared \(\beta\) across rows, cross-block columns on the same rows).
- Replacing R **`optim`** with native C++ BFGS (planned later; see Phase 4).

### Reference implementation today

- Non-blocked GLM draws: **`rNormal_reg()`** → **`.rNormalGLM_cpp()`** →
  **`rNormalGLM()`** in `src/rNormalGLM.cpp`.
- Use case: Eight Schools Gibbs — `data-raw/make_Chapter13_Eight_Schools_gibbs_output.R`.

---

## 2. Naming conventions

Aligned with existing simfuncs (`rNormal_reg`, `rNormalGamma_reg`,
`rindepNormalGamma_reg`) and C++ (`rNormalGLM`, `rNormalReg`,
`rIndepNormalGammaReg`).

| Layer | Name | File (planned) |
|-------|------|----------------|
| **R (public)** | `rNormalGLM_reg_block()` | **`R/simfunction_block.R`** |
| **R (docs)** | `@name simfuncs_block` / `@family simfuncs_block` | same file; `@seealso` **`simfuncs`** |
| **R (internal)** | `.rNormalGLMBlocks_cpp()` | `R/rcpp_wrappers.R` |
| **R helpers** | `normalize_block()`, `normalize_prior_for_blocks()` | **`R/simfunction_block_utils.R`** (or top of `simfunction_block.R` until small) |
| **C++ (export)** | `rNormalGLMBlocks_cpp_export()` | `src/export_wrappers.cpp` |
| **C++ (core)** | `rNormalGLMBlocks()` | `src/rNormalGLMBlocks.cpp` — loops blocks, calls `rNormalGLM()` each time |

**R file pairing:** **`R/simfunction.R`** holds **`pfamily$simfun`** IID samplers
(`rNormal_reg`, …). **`R/simfunction_block.R`** holds block **full conditionals**
only — adjacent in the directory and in documentation.

**Note:** The public R name **`rNormalGLM_reg_block`** mirrors `rNormal_reg` /
`rindepNormalGamma_reg` (`*_reg_*` on R only). The C++ names follow **`rNormalGLM*`** without `_reg_`, matching `rNormalGLM` / `rNormalGLM_std`. Plural **`Blocks`** on the C++ multi-block entry only.

### Future block samplers (same file family)

| R (later) | C++ (later) | Notes |
|-----------|-------------|--------|
| `rNormal_reg_block` | `rNormalRegBlock(s)` | Gaussian WLS; optional |
| `rNormalGamma_reg_block` | `rNormalGammaRegBlock(s)` | Per-block Normal–Gamma |
| `rindepNormalGamma_reg_block` | `rIndepNormalGammaRegBlock(s)` | Heavy; envelope per block |
| `rGamma_reg_block` | TBD | Scalar conjugate units |

Registered symbol (generated): `_glmbayes_rNormalGLMBlocks_cpp_export`.

---

## 3. Statistical assumptions

### Conditional independence

All of the following must hold:

1. **Row partition:** each observation belongs to exactly one block.
2. **Likelihood factorization:** block \(b\) depends only on \(\beta_b\) and
   data \((y_b, X_b)\).
3. **Prior factorization:** \(p(\beta_1,\ldots,\beta_k)=\prod_b p(\beta_b\mid\mu_b,P_b)\)
   (block-diagonal precision on stacked \(\beta\), or explicit per-block priors).

Hyperparameters fixed at a higher Gibbs level are fine: conditional on them,
blocks are independent.

### Observation-level vectors

When supplied, these share the **same blocking as `y` and `x`**:

- **`offset`:** `length(offset) == l2` (or `1`, recycled to `l2`), sliced by block.
- **`weights`:** same as `offset`.

### Design matrix layout

- Stacked: `y` length `l2`, `x` is `l2 × l1`, block-major row order if using
  `l2_blocks` counts.
- Block \(b\) rows only “see” \(\beta_b\) (other columns zero, or separate
  `x_b` with `l1` columns).

---

## 4. R user API

### Primary function

```r
rNormalGLM_reg_block(
  n,
  y,
  x,
  block,
  prior_list = NULL,
  prior_lists = NULL,
  offset = NULL,
  weights = 1,
  family = gaussian(),   # error if family$family == "gaussian"
  Gridtype = 2L,
  n_envopt = NULL,
  use_parallel = TRUE,     # within-block draw parallelism (n > 1)
  use_opencl = FALSE,
  verbose = FALSE,
  progbar = FALSE
)
```

### `block` argument (user-facing)

One of:

| Form | Length | Meaning |
|------|--------|---------|
| `factor` / `integer` | `l2` | Observation membership |
| `integer` | `k` | `l2_blocks`: contiguous counts, `sum == l2` |
| `list` | `k` | Named or unnamed row index vectors into `y` / `x` |

### Prior inputs

**Extraction and `Sigma` → `P` conversion now run in C++** (Phase 2b):
`normalize_block_cpp` / `normalize_prior_for_blocks_cpp` /
`prior_payload_from_blocks` in `src/block_utils.cpp`, with the legacy R
helpers (`normalize_block`, `normalize_prior_for_blocks`) retained for tests
and external callers. The C++ versions mirror R semantics exactly, including
`is.null()` behaviour: a present-but-`NULL` element (e.g.
`list(dispersion = NULL)` from the two-block Gibbs prior builder) is treated
as absent. `Sigma` → `P` uses `arma::inv_sympd` (LAPACK `dpotrf` + `dpotri`),
matching R's `chol2inv(chol(Sigma))` bit-for-bit.

Accepted shapes (normalize to canonical):

1. **Single `prior_list`** — recycled to every block (`mu` length `l1`,
   one `P` or `Sigma`).
2. **`prior_list$mu`** as `l1 × k` matrix; **`prior_list$P`** or **`Sigma`**
   as `list` length `k` or `l1 × l1 × k` array.
3. **`prior_list$blocks`** — named list aligned with `levels(block)`.
4. **`prior_lists`** — length `k` (or `1` to recycle) of `prior_list` objects.

Canonical output of `normalize_prior_for_blocks()`:

```r
list(
  mu = <numeric l1 or matrix l1 x k>,
  P_blocks = <list of length 1 or k, each l1 x l1>,
  prior_by_block = <logical>,
  dispersion = <scalar or length-k>,
  # optional passthrough: ddef, shape, rate if needed later
)
```

Validation (R): PD checks on `P`, dim match `ncol(x)`, symmetric `P`,
`length(weights)`, `length(offset)` in `c(1, l2)`.

### Return object (Phase 1 R implementation)

- **`coefficients`:** matrix **`k × l1`** — **rows = blocks**, **columns = parameters** (one Gibbs draw, `n = 1`).
- **`coef.mode`:** matrix **`k × l1`**, same layout.
- **`dispersion`:** length-`k` vector (one value per block from `rNormal_reg`).
- **`y`, `x`, `offset`:** full stacked inputs (as `rNormal_reg` does).
- **`Prior`:** list describing shared or per-block priors.
- **`blocks`:** optional `list` of per-block summaries.
- **`attr(block_info)`** or element **`block_info`:** `k`, `l2_blocks`, `ids`.
- **`Envelope`:** `list` length `k` (or `NULL` / omitted in v1 if size is an issue).
- Classes: `"rglmb"`, `"glmb"`, etc. (optional; match `rNormal_reg` non-Gaussian branch where helpful for utilities — not required for Gibbs loops).

---

## 4b. Integration policy (package boundaries)

### v1: building blocks only

- **Do not** wire **`rNormalGLM_reg_block`** into **`rglmb()`**, **`rlmb()`**, or
  **`pfamily$simfun`** by default. Those paths keep the **`simfunction.R`**
  contract: `simfun(n, y, x, prior_list, offset, weights, family, …)` with no
  `block` argument.
- **Intended use:** user-defined **block Gibbs** (and similar) outer loops, as in
  Eight Schools scripts — one call draws one block conditional (often **`n = 1`**).
- **Guarantee:** correct draw from \(p(\beta_b \mid y_b,\ldots)\) when conditional
  independence holds (Section 3).
- **Not guaranteed:** convergence of the full chain; geometric ergodicity; bounds on
  total-variation distance to the target. **Convergence diagnostics are the user’s
  responsibility** (e.g. `coda`, custom monitors, vignette patterns).

### Optional later wrappers (not v1)

| Layer | Purpose |
|-------|---------|
| **`rglmb_block()`**, **`rlmb_block()`** | Convenience: formula / `block` / multi-block sweep — only if a clear API is needed; separate from **`rglmb`/`rlmb`**. |
| **`lmerb()`**, **`glmerb()`** | Bayesian **`lmer`/`glmer`-style** entry points — **deferred** until reliable convergence theory and diagnostics (e.g. geometric ergodicity, computable TV bound to target) are available. |

Until then, document in **`?rNormalGLM_reg_block`**: *one block full conditional, not a complete MCMC scheme.*

---

## 5. R internal helpers

### `normalize_block(block, l2)`

Returns:

```r
list(
  k = k,
  ids = <character or NULL>,
  l2_blocks = <integer length k>,
  starts = <integer length k>,   # 1-based row starts if contiguous
  rows = <list of index vectors>  # always populated for C++
)
```

Rules:

- `length(block) == l2` → `split(seq_len(l2), block)`.
- `length(block) == k` → treat as `l2_blocks` (contiguous).
- `is.list(block)` → validate disjoint cover of `1:l2`.

### `normalize_prior_for_blocks(prior_list, prior_lists, block_info, l1, k)`

Precedence:

1. `prior_lists` if not `NULL`.
2. else `prior_list$blocks` matched to `block_info$ids`.
3. else matrix `prior_list$mu` with `ncol == k`.
4. else shared single `prior_list`.

Convert `Sigma` → `P` as in `rNormal_reg()`.

---

## 6. R–C++ boundary

### Wrapper rule

**`.rNormalGLMBlocks_cpp()`** in `R/rcpp_wrappers.R`: positional `.Call` only,
no preprocessing (same contract as `.rNormalGLM_cpp()`).

### High-level export (Phase 2b)

**`block_rNormalGLM()`** (R) now does validation + family/link resolution +
`glmbfamfunc()` only, then calls **`.block_rNormalGLM_cpp()`** →
**`block_rNormalGLM_cpp_export()`** (`src/block_utils.cpp`), which performs
block partition and prior payload assembly in C++ and delegates to the
unchanged `rNormalGLMBlocks()` → `rNormalGLM()` pipeline. `f2`/`f3` are still
R closures used by R `optim` per block (Phase 4 unchanged).

**Equivalence policy:** posterior modes (`coef.mode`) are deterministic and
must match the legacy R-prep payload tightly; individual rejection-sampler
draws are *not* comparable across payload assembly paths — compare draw
means over longer runs (see `data-raw/test_block_rNormalGLM_cpp.R`).

### Proposed `.Call` arguments

```r
.rNormalGLMBlocks_cpp(
  n, y, x, offset, wt,
  dispersion,
  mu,              # numeric matrix l1 x 1 or l1 x k
  P_blocks,        # list of l1 x l1 matrices (length 1 or k)
  prior_by_block,  # logical
  row_blocks,      # list of integer vectors (length k), 1-based
  f2, f3, start,
  family, link,
  Gridtype, n_envopt,
  use_parallel,
  use_opencl, verbose
)
```

Alternative: pass `l2_blocks` + `starts` instead of `row_blocks` when
contiguous-only v1 is desired; document that v1 may require contiguous blocks
for minimal SEXP payload.

### `f2` / `f3`

Still required for **R `optim`** inside each per-block `rNormalGLM()` call until Phase 4.
Passed once from `glmbfamfunc(family)` on the R side.

---

## 7. C++ structure

### Implementation

**`rNormalGLMBlocks()`** only (no separate single-block C++ wrapper):

- loop `b = 0 .. k-1`
- slice `y`, `x`, `offset`, `wt` via `row_blocks[b]`
- `mu_b`, `P_b` from `mu.col(b)` or shared
- `out_b = rNormalGLM(...)` — the existing non-block sampler
- stack `coefficients` / `coef.mode` as `k × l1` matrices

**`rNormalGLM()`** is unchanged; a single-block call is just `k == 1`.

### Parallelism

| Flag | Effect |
|------|--------|
| `use_parallel && n > 1` | `rNormalGLM_std_parallel` inside a block only |
| Block loop | **Serial** over `b = 0 .. k-1` (no across-block parallelism) |

Across-block parallelism is intentionally omitted for Gibbs-friendly, reproducible chains.

### Files to touch

| File | Change |
|------|--------|
| `src/rNormalGLMBlocks.cpp` | `rNormalGLMBlocks()` loop → `rNormalGLM()` |
| `src/simfuncs.h` | Declarations |
| `src/export_wrappers.cpp` | `rNormalGLMBlocks_cpp_export` |
| `R/RcppExports.R` | Generated |
| `R/rcpp_wrappers.R` | `.rNormalGLMBlocks_cpp` |
| `R/simfunction_block_utils.R` | `normalize_*` (`@keywords internal`) |
| `R/simfunction_block.R` | `rNormalGLM_reg_block`, `@family simfuncs_block` |
| `R/simfunction.R` | **No change** (pfamily simfuns only) |
| `tests/testthat/test-rNormalGLM_reg_block.R` | New |

Run **`Rcpp::compileAttributes()`** after adding exports.

---

## 8. Phased implementation

| Phase | Deliverable |
|-------|-------------|
| **0** | This document (done) |
| **1** | **Done (R):** `normalize_block`, `normalize_prior_for_blocks`, `rNormalGLM_reg_block` → **`rNormal_reg(n=1)`** per block; return **`k × l1`** matrices |
| **2** | C++: `rNormalGLMBlocks`, wire `.Call` |
| **3** | Tests: `k == 1` vs `rNormal_reg`; two-block toy; offset slicing; shared vs per-block prior |
| **4** | Perf: optional C++ mode finder (`f2_f3_*`) replacing per-block `optim` |
| **5** | Docs: roxygen, `NEWS.md`, vignette cross-link (Chapter 17 / Eight Schools) |

---

## 9. Testing plan

- **Equivalence:** `k == 1`, same data as `rNormal_reg(..., family = binomial())`
  (up to tolerance / seed).
- **Two blocks:** different `l2_b`, shared `mu`/`P`.
- **Per-block priors:** `prior_lists`, `prior_list$blocks`, `mu` matrix `l1 × k`.
- **`offset`:** `NULL`, scalar `1`, length `l2`, correct slice per block.
- **Gibbs:** `n = 1`, `use_parallel = FALSE`.
- **IID:** `n = 100`, one block vs stacked call (if model equivalent).
- **Errors:** `family = gaussian()`; mismatched `length(block)`; non-PD `P`.

---

## 10. Open questions

- Return per-block **`Envelope`** lists (large) vs strip from default output.
- v1: contiguous **`l2_blocks` only** vs general **`row_blocks`** in C++.
- Column order in **`coefficients`**: block-major vs interleaved (document).

**Resolved:** **`pfamily$simfun`** / **`rglmb`/`rlmb`** — **no** default hook; block
samplers live in **`simfunction_block.R`** only (Section 4b).

---

## 11. Related files

| Topic | Path |
|-------|------|
| GLM sampler | `src/rNormalGLM.cpp` |
| IID simfuns / `pfamily` | `R/simfunction.R` (`rNormal_reg`, …) |
| Block simfuns (planned) | `R/simfunction_block.R` |
| C++ wrappers | `R/rcpp_wrappers.R`, `src/export_wrappers.cpp` |
| Family f2/f3 | `src/famfuncs_*.cpp`, `R/glmbfamfunc` |
| Prior setup | `R/prior.R` (`Prior_Setup`) |
| pfamily checklist style | `inst/ADDING_PFAMILY.md` |
| Gibbs example | `data-raw/make_Chapter13_Eight_Schools_gibbs_output.R` |
| BikeSharing Block 2 benchmark | `data-raw/benchmark_BikeSharing_rNormalGLM_reg_block.R` |

---

## 12. Changelog (design doc)

| Date | Note |
|------|------|
| 2026-05-28 | Initial design; R: `rNormalGLM_reg_block`; C++: `rNormalGLMBlocks` (loop → `rNormalGLM`). |
| 2026-05-28 | R files: `simfunction_block.R`, `simfunction_block_utils.R`; integration policy (no `rglmb`/`rlmb` v1); future `rglmb_block` / `lmerb` noted. |
| 2026-06-09 | Phase 2b: block partition + prior payload ported to C++ (`block_rNormalGLM_cpp_export` in `src/block_utils.cpp`); `block_rNormalGLM()` slimmed to validation + `glmbfamfunc`; `Sigma`→`P` via `inv_sympd`; present-but-NULL prior elements treated as absent (two-block Gibbs compatibility). |
| 2026-06-10 | Phase 3 (port-only): full `two_block_rNormal_reg()` Gibbs loop moved to C++ (`two_block_rNormal_reg_cpp_export` in `src/twoBlockGibbs.cpp`): mu_all assembly, Block 1 via existing block exports, Block 2 via `rNormalReg` core, `m_convergence` inner steps, replicate sampling, progress bar. R wrapper retains validation + `glmbfamfunc` + output assembly. Equivalence: averages over many draws (C++ rejection sampler uses its own RNG stream); regression test `data-raw/test_two_block_cpp.R`. |
