#ifndef GLMBAYESCORE_CAPI_H
#define GLMBAYESCORE_CAPI_H

#include <R_ext/Rdynload.h>
#include <R_ext/Error.h>
#include <Rinternals.h>

/*
 * glmbayesCore C-callable API
 *
 * Usage:
 *   #include <glmbayesCore/glmbayesCore_capi.h>
 *   SEXP out = glmbayesCore_rNormalReg(n, y, x, mu, P, offset, wt,
 *                                      dispersion, f2, f3, start, family,
 *                                      link, Gridtype);
 *
 * Notes:
 * - Arguments and return values are SEXP for C ABI portability. Each entry
 *   point has the same arguments, in the same order, as the glmbayesCore C++
 *   function it wraps, and returns the same list.
 * - Errors are raised as R errors.
 * - Use glmbayesCore_api_version() for compatibility checks.
 *
 * Layout (same as opencltools_capi.h):
 *   1) function pointer typedefs
 *   2) internal resolver (R_GetCCallable with one clear error message)
 *   3) typed accessors (resolve once, then cache the pointer)
 *   4) convenience wrappers (what downstream code calls)
 */

#ifdef __cplusplus
extern "C" {
#endif

  /* -------- function pointer typedefs -------- */

  typedef int (*glmbayesCore_api_version_t)(void);

  typedef SEXP (*glmbayesCore_rNormalGLMGroups_t)(
      SEXP n, SEXP y, SEXP x, SEXP offset, SEXP wt, SEXP dispersion,
      SEXP mu, SEXP P_blocks, SEXP prior_by_block, SEXP row_blocks,
      SEXP f2, SEXP f3, SEXP family, SEXP link, SEXP Gridtype,
      SEXP n_envopt, SEXP use_parallel, SEXP use_opencl, SEXP verbose);

  typedef SEXP (*glmbayesCore_rNormalReg_t)(
      SEXP n, SEXP y, SEXP x, SEXP mu, SEXP P, SEXP offset, SEXP wt,
      SEXP dispersion, SEXP f2, SEXP f3, SEXP start, SEXP family,
      SEXP link, SEXP Gridtype);

  /* -------- internal resolver -------- */

  static inline DL_FUNC glmbayesCore_capi_resolve_(const char* name) {
    DL_FUNC p = R_GetCCallable("glmbayesCore", name);
    if (!p) {
      Rf_error(
        "glmbayesCore C-callable '%s' not found; check package version/load order.",
        name
      );
    }
    return p;
  }

  /* -------- typed accessors -------- */

  static inline glmbayesCore_api_version_t glmbayesCore_api_version_fn(void) {
    static glmbayesCore_api_version_t fn = NULL;
    if (!fn) fn = (glmbayesCore_api_version_t) glmbayesCore_capi_resolve_("glmbayesCore_api_version");
    return fn;
  }

  static inline glmbayesCore_rNormalGLMGroups_t glmbayesCore_rNormalGLMGroups_fn(void) {
    static glmbayesCore_rNormalGLMGroups_t fn = NULL;
    if (!fn) fn = (glmbayesCore_rNormalGLMGroups_t) glmbayesCore_capi_resolve_("glmbayesCore_rNormalGLMGroups");
    return fn;
  }

  static inline glmbayesCore_rNormalReg_t glmbayesCore_rNormalReg_fn(void) {
    static glmbayesCore_rNormalReg_t fn = NULL;
    if (!fn) fn = (glmbayesCore_rNormalReg_t) glmbayesCore_capi_resolve_("glmbayesCore_rNormalReg");
    return fn;
  }

  /* -------- convenience wrappers -------- */

  static inline int glmbayesCore_api_version(void) {
    return glmbayesCore_api_version_fn()();
  }

  /* Wraps glmbayes::sim::group::rNormalGLMGroups. */
  static inline SEXP glmbayesCore_rNormalGLMGroups(
      SEXP n, SEXP y, SEXP x, SEXP offset, SEXP wt, SEXP dispersion,
      SEXP mu, SEXP P_blocks, SEXP prior_by_block, SEXP row_blocks,
      SEXP f2, SEXP f3, SEXP family, SEXP link, SEXP Gridtype,
      SEXP n_envopt, SEXP use_parallel, SEXP use_opencl, SEXP verbose) {
    return glmbayesCore_rNormalGLMGroups_fn()(
        n, y, x, offset, wt, dispersion, mu, P_blocks, prior_by_block,
        row_blocks, f2, f3, family, link, Gridtype, n_envopt, use_parallel,
        use_opencl, verbose);
  }

  /* Wraps glmbayes::sim::rNormalReg. */
  static inline SEXP glmbayesCore_rNormalReg(
      SEXP n, SEXP y, SEXP x, SEXP mu, SEXP P, SEXP offset, SEXP wt,
      SEXP dispersion, SEXP f2, SEXP f3, SEXP start, SEXP family,
      SEXP link, SEXP Gridtype) {
    return glmbayesCore_rNormalReg_fn()(
        n, y, x, mu, P, offset, wt, dispersion, f2, f3, start, family, link,
        Gridtype);
  }

#ifdef __cplusplus
}
#endif

#endif /* GLMBAYESCORE_CAPI_H */
