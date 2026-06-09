# Design: conditionally independent Gaussian block sampling (`block_rNormalReg`)

Maintainer-facing implementation plan for the **Gaussian** block Gibbs path.
Companion to `inst/DESIGN_RGLM_BLOCKS.md` (GLM envelope blocks).
Excluded from the source tarball via `.Rbuildignore`.

---

## 1. Purpose and scope

### Goal

Enable **direct C++ calls** for Gaussian conditionally independent block sampling—the path used by Block 1 of `two_block_rNormal_reg` / `lmerb` / `glmerb` when `family = gaussian()`.

### In scope (Phase 1 — Gaussian only)

- `block_rNormalReg()` public API unchanged in signature and return shape
- C++ helpers: block partition, prior payload, block loop
- Native Gaussian `rNormalReg` kernel (no R `lm.fit` callback per block)
- Drop dead-weight `f2`/`f3` R `Function` threading on the Gaussian block path
- Tests against pre-migration behavior and `big_word_club` school-level draws

### Out of scope (later phases)

- `block_rNormalGLM` / GLM envelope path (`DESIGN_RGLM_BLOCKS.md`)
- `two_block_rNormal_reg` inner Gibbs loop in C++
- Block 2 (`multi_rNormal_reg`) C++
- OpenCL / parallel envelope paths

---

## 2. Current architecture (pre-migration)

```
block_rNormalReg (R)
  -> normalize_block, normalize_prior_for_blocks
  -> .prior_payload_for_rNormalGLMBlocks_cpp
  -> glmbfamfunc(gaussian()) -> f2, f3
  -> .rNormalRegBlocks_cpp
  -> rNormalRegBlocks (C++)
       -> rNormalReg (C++) -> R lm.fit  [per block]
```

**Bottleneck:** For `J` groups, each `block_rNormalReg(n=1)` runs `J` times `lm.fit` from C++.

---

## 3. Target architecture (Phase 1)

```
block_rNormalReg (R thin wrapper: dimnames, class, call)
  -> .block_rNormalReg_cpp
  -> block_rNormalReg_cpp_export
       -> normalize_block_cpp
       -> normalize_prior_for_blocks_cpp + prior_payload_blocks_cpp
       -> rNormalRegBlocks
            -> rNormalRegGaussian (pure Armadillo, no lm.fit)
```

---

## 4. Implementation bites

| Bite | Item | File(s) |
|------|------|---------|
| 1 | This design doc | `inst/DESIGN_RNORMALREG_BLOCKS.md` |
| 2 | `rNormalRegGaussian` (no `lm.fit`) | `src/rNormalReg.cpp`, `simfuncs.h` |
| 3 | `normalize_block_cpp` | `src/block_utils.cpp` |
| 4 | `prior_payload_blocks_cpp` | `src/block_utils.cpp` |
| 5 | `block_rNormalReg_cpp_export` | `export_wrappers.cpp`, `rcpp_wrappers.R` |
| 6 | Slim `block_rNormalReg()` | `R/simfunction_block.R` |
| 7 | Tests | `data-raw/test_block_rNormalReg_cpp.R` |

---

## 5. API layers

| Caller | Entry |
|--------|-------|
| End user / lmebayes | `block_rNormalReg(...)` |
| Prebuilt partition | `.rNormalRegBlocks_cpp(...)` |
| Internal | `.block_rNormalReg_cpp(...)` |

---

## 6. Success criteria

- Numerical equivalence: max abs(coef.mode diff) < 1e-10 on fixed-seed examples
- No `lm.fit` in hot path for Gaussian blocks
- `devtools::check()` passes

---

## 7. Future phases

| Phase | Item |
|-------|------|
| 2 | `block_rNormalGLM_cpp_export` + native `glmbayes::fam` dispatch |
| 3 | Block 2 conjugate hyper draw C++ |
| 4 | `mu_all` C++ helper |
| 5 | `two_block_rNormal_reg` inner loop in C++ |
