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
 *   function it wraps, and returns the same value (R_NilValue for functions
 *   that return nothing).
 * - The batch_* / ensure_* helpers modify their first argument in place, as
 *   the wrapped functions do; it must be a double array.
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

  typedef SEXP (*glmbayesCore_rNormalGLMGroupChains_t)(
      SEXP y, SEXP x, SEXP offset, SEXP wt, SEXP dispersion, SEXP mu,
      SEXP P, SEXP f2, SEXP f3, SEXP family, SEXP link, SEXP Gridtype,
      SEXP n_envopt, SEXP use_opencl, SEXP check, SEXP verbose,
      SEXP mode_method, SEXP setup_method, SEXP start, SEXP warm_method,
      SEXP hess_method);

  typedef SEXP (*glmbayesCore_rNormalReg_t)(
      SEXP n, SEXP y, SEXP x, SEXP mu, SEXP P, SEXP offset, SEXP wt,
      SEXP dispersion, SEXP f2, SEXP f3, SEXP start, SEXP family,
      SEXP link, SEXP Gridtype);

  typedef SEXP (*glmbayesCore_group_rNormalGLM_t)(
      SEXP n, SEXP y, SEXP x, SEXP group, SEXP prior_list, SEXP prior_lists,
      SEXP offset, SEXP wt, SEXP f2, SEXP f3, SEXP family, SEXP link,
      SEXP Gridtype, SEXP n_envopt, SEXP use_parallel, SEXP use_opencl,
      SEXP verbose);

  typedef SEXP (*glmbayesCore_group_rNormalReg_t)(
      SEXP n, SEXP y, SEXP x, SEXP group, SEXP prior_list, SEXP prior_lists,
      SEXP offset, SEXP wt, SEXP f2, SEXP f3, SEXP Gridtype);

  typedef SEXP (*glmbayesCore_two_block_build_mu_all_t)(
      SEXP x_hyper, SEXP fixef, SEXP re_names, SEXP group_levels);

  typedef SEXP (*glmbayesCore_two_block_groups_prior_with_tau2_t)(
      SEXP base_prior, SEXP tau2_vec, SEXP ptypes, SEXP re_names,
      SEXP mu_all);

  typedef SEXP (*glmbayesCore_two_block_groups_one_chain_from_mu_P_t)(
      SEXP mu_all, SEXP P, SEXP dispersion, SEXP ddef, SEXP design,
      SEXP family, SEXP re_names, SEXP group_levels, SEXP f2, SEXP f3,
      SEXP f2_gauss, SEXP f3_gauss);

  typedef SEXP (*glmbayesCore_fixef_list_from_batch_chain_t)(
      SEXP batch_fixef, SEXP chain_i, SEXP re_names);

  typedef SEXP (*glmbayesCore_batch_tau2_chain_row_t)(
      SEXP batch_tau2, SEXP chain_i);

  typedef SEXP (*glmbayesCore_batch_b_assign_slice_t)(
      SEXP b_store, SEXP chain_i, SEXP b_draw);

  typedef SEXP (*glmbayesCore_batch_b_n_chains_t)(SEXP batch_b);

  typedef SEXP (*glmbayesCore_ensure_batch_b_dimnames_t)(
      SEXP batch_b, SEXP group_levels, SEXP re_names, SEXP n_chains);

  typedef SEXP (*glmbayesCore_batch_iters_ranef_add_t)(
      SEXP iters_ranef, SEXP chain_i, SEXP iters_mean);

  typedef SEXP (*glmbayesCore_x_hyper_as_matrix_t)(SEXP x_hyper);

  typedef SEXP (*glmbayesCore_rIndepNormalGammaReg_t)(
      SEXP n, SEXP y, SEXP x, SEXP mu, SEXP P, SEXP offset, SEXP wt,
      SEXP shape, SEXP rate, SEXP max_disp_perc, SEXP disp_lower,
      SEXP disp_upper, SEXP Gridtype, SEXP n_envopt, SEXP use_parallel,
      SEXP use_opencl, SEXP verbose, SEXP progbar);

  typedef SEXP (*glmbayesCore_progress_bar_t)(SEXP x, SEXP N, SEXP prefix);

  typedef SEXP (*glmbayesCore_progress_bar_finish_t)(SEXP newline);

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

#define GLMBAYESCORE_CAPI_ACCESSOR_(name)                                  \
  static inline name##_t name##_fn(void) {                                 \
    static name##_t fn = NULL;                                             \
    if (!fn) fn = (name##_t) glmbayesCore_capi_resolve_(#name);            \
    return fn;                                                             \
  }

  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_api_version)
  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_rNormalGLMGroups)
  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_rNormalGLMGroupChains)
  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_rNormalReg)
  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_group_rNormalGLM)
  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_group_rNormalReg)
  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_two_block_build_mu_all)
  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_two_block_groups_prior_with_tau2)
  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_two_block_groups_one_chain_from_mu_P)
  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_fixef_list_from_batch_chain)
  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_batch_tau2_chain_row)
  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_batch_b_assign_slice)
  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_batch_b_n_chains)
  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_ensure_batch_b_dimnames)
  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_batch_iters_ranef_add)
  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_x_hyper_as_matrix)
  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_rIndepNormalGammaReg)
  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_progress_bar)
  GLMBAYESCORE_CAPI_ACCESSOR_(glmbayesCore_progress_bar_finish)

#undef GLMBAYESCORE_CAPI_ACCESSOR_

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

  /* Wraps glmbayes::sim::group::rNormalGLMGroupChains: one group, one draw
   * per chain (column of mu). mode_method: 0 optim, 1 vmmin with C++ f2/f3,
   * 2 optim plus comparison against vmmin (result field mode_compare),
   * 3 Newton on the C++ f2/f3 and exact Hessian, 4 vmmin plus comparison
   * against 3, 5 Newton on the allocation-free fam f2_f3_f4 kernels, 6 vmmin
   * plus comparison against 5 (result field newton_compare; Newton falls
   * back to vmmin when it fails).
   * setup_method: 0 R grid setup, 1 C++, 2 both with an exact comparison
   * (R kept; result field setup_compare). start: R_NilValue or an
   * ncol(x) x ncol(mu) matrix of starting modes (NA column = cold);
   * warm_method: 0 cold, 1 warm start, 2 both with comparison (cold kept;
   * result field warm_compare). hess_method: 0 finite differences, 1 exact
   * (finite differences if not finite or not PD), 2 both with comparison
   * (finite differences kept; result field hessian_compare). The result
   * also has mode_hessian. */
  static inline SEXP glmbayesCore_rNormalGLMGroupChains(
      SEXP y, SEXP x, SEXP offset, SEXP wt, SEXP dispersion, SEXP mu,
      SEXP P, SEXP f2, SEXP f3, SEXP family, SEXP link, SEXP Gridtype,
      SEXP n_envopt, SEXP use_opencl, SEXP check, SEXP verbose,
      SEXP mode_method, SEXP setup_method, SEXP start, SEXP warm_method,
      SEXP hess_method) {
    return glmbayesCore_rNormalGLMGroupChains_fn()(
        y, x, offset, wt, dispersion, mu, P, f2, f3, family, link, Gridtype,
        n_envopt, use_opencl, check, verbose, mode_method, setup_method,
        start, warm_method, hess_method);
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

  /* Wraps glmbayes::sim::group::group_rNormalGLM_cpp_export. */
  static inline SEXP glmbayesCore_group_rNormalGLM(
      SEXP n, SEXP y, SEXP x, SEXP group, SEXP prior_list, SEXP prior_lists,
      SEXP offset, SEXP wt, SEXP f2, SEXP f3, SEXP family, SEXP link,
      SEXP Gridtype, SEXP n_envopt, SEXP use_parallel, SEXP use_opencl,
      SEXP verbose) {
    return glmbayesCore_group_rNormalGLM_fn()(
        n, y, x, group, prior_list, prior_lists, offset, wt, f2, f3, family,
        link, Gridtype, n_envopt, use_parallel, use_opencl, verbose);
  }

  /* Wraps glmbayes::sim::group::group_rNormalReg_cpp_export. */
  static inline SEXP glmbayesCore_group_rNormalReg(
      SEXP n, SEXP y, SEXP x, SEXP group, SEXP prior_list, SEXP prior_lists,
      SEXP offset, SEXP wt, SEXP f2, SEXP f3, SEXP Gridtype) {
    return glmbayesCore_group_rNormalReg_fn()(
        n, y, x, group, prior_list, prior_lists, offset, wt, f2, f3,
        Gridtype);
  }

  /* Wraps glmbayes::sim::twoblock::two_block_build_mu_all. */
  static inline SEXP glmbayesCore_two_block_build_mu_all(
      SEXP x_hyper, SEXP fixef, SEXP re_names, SEXP group_levels) {
    return glmbayesCore_two_block_build_mu_all_fn()(
        x_hyper, fixef, re_names, group_levels);
  }

  /* Wraps glmbayes::sim::twoblock::two_block_groups_prior_with_tau2. */
  static inline SEXP glmbayesCore_two_block_groups_prior_with_tau2(
      SEXP base_prior, SEXP tau2_vec, SEXP ptypes, SEXP re_names,
      SEXP mu_all) {
    return glmbayesCore_two_block_groups_prior_with_tau2_fn()(
        base_prior, tau2_vec, ptypes, re_names, mu_all);
  }

  /* Wraps glmbayes::sim::twoblock::two_block_groups_one_chain_from_mu_P_impl. */
  static inline SEXP glmbayesCore_two_block_groups_one_chain_from_mu_P(
      SEXP mu_all, SEXP P, SEXP dispersion, SEXP ddef, SEXP design,
      SEXP family, SEXP re_names, SEXP group_levels, SEXP f2, SEXP f3,
      SEXP f2_gauss, SEXP f3_gauss) {
    return glmbayesCore_two_block_groups_one_chain_from_mu_P_fn()(
        mu_all, P, dispersion, ddef, design, family, re_names, group_levels,
        f2, f3, f2_gauss, f3_gauss);
  }

  /* Wraps glmbayes::sim::twoblock::fixef_list_from_batch_chain. */
  static inline SEXP glmbayesCore_fixef_list_from_batch_chain(
      SEXP batch_fixef, SEXP chain_i, SEXP re_names) {
    return glmbayesCore_fixef_list_from_batch_chain_fn()(
        batch_fixef, chain_i, re_names);
  }

  /* Wraps glmbayes::sim::twoblock::batch_tau2_chain_row. */
  static inline SEXP glmbayesCore_batch_tau2_chain_row(
      SEXP batch_tau2, SEXP chain_i) {
    return glmbayesCore_batch_tau2_chain_row_fn()(batch_tau2, chain_i);
  }

  /* Wraps glmbayes::sim::twoblock::batch_b_assign_slice (in place). */
  static inline SEXP glmbayesCore_batch_b_assign_slice(
      SEXP b_store, SEXP chain_i, SEXP b_draw) {
    return glmbayesCore_batch_b_assign_slice_fn()(b_store, chain_i, b_draw);
  }

  /* Wraps glmbayes::sim::twoblock::batch_b_n_chains. */
  static inline SEXP glmbayesCore_batch_b_n_chains(SEXP batch_b) {
    return glmbayesCore_batch_b_n_chains_fn()(batch_b);
  }

  /* Wraps glmbayes::sim::twoblock::ensure_batch_b_dimnames (in place). */
  static inline SEXP glmbayesCore_ensure_batch_b_dimnames(
      SEXP batch_b, SEXP group_levels, SEXP re_names, SEXP n_chains) {
    return glmbayesCore_ensure_batch_b_dimnames_fn()(
        batch_b, group_levels, re_names, n_chains);
  }

  /* Wraps glmbayes::sim::twoblock::batch_iters_ranef_add (in place). */
  static inline SEXP glmbayesCore_batch_iters_ranef_add(
      SEXP iters_ranef, SEXP chain_i, SEXP iters_mean) {
    return glmbayesCore_batch_iters_ranef_add_fn()(
        iters_ranef, chain_i, iters_mean);
  }

  /* Wraps glmbayes::sim::twoblock::x_hyper_as_matrix. */
  static inline SEXP glmbayesCore_x_hyper_as_matrix(SEXP x_hyper) {
    return glmbayesCore_x_hyper_as_matrix_fn()(x_hyper);
  }

  /* Wraps glmbayes::sim::rIndepNormalGammaReg. */
  static inline SEXP glmbayesCore_rIndepNormalGammaReg(
      SEXP n, SEXP y, SEXP x, SEXP mu, SEXP P, SEXP offset, SEXP wt,
      SEXP shape, SEXP rate, SEXP max_disp_perc, SEXP disp_lower,
      SEXP disp_upper, SEXP Gridtype, SEXP n_envopt, SEXP use_parallel,
      SEXP use_opencl, SEXP verbose, SEXP progbar) {
    return glmbayesCore_rIndepNormalGammaReg_fn()(
        n, y, x, mu, P, offset, wt, shape, rate, max_disp_perc, disp_lower,
        disp_upper, Gridtype, n_envopt, use_parallel, use_opencl, verbose,
        progbar);
  }

  /* Wraps glmbayes::progress::progress_bar. */
  static inline SEXP glmbayesCore_progress_bar(SEXP x, SEXP N, SEXP prefix) {
    return glmbayesCore_progress_bar_fn()(x, N, prefix);
  }

  /* Wraps glmbayes::progress::progress_bar_finish. */
  static inline SEXP glmbayesCore_progress_bar_finish(SEXP newline) {
    return glmbayesCore_progress_bar_finish_fn()(newline);
  }

#ifdef __cplusplus
}
#endif

#endif /* GLMBAYESCORE_CAPI_H */
