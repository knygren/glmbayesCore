# -------------------------------------------------------------------------
#  Rcpp interface wrappers for the two-block Gibbs engines and the
#  block ING (rIndepNormalGammaRegBlock*) envelope pipeline.
#  Same conventions as rcpp_wrappers.R: strictly positional, no validation.
# -------------------------------------------------------------------------


# =============================================================================
#  Tier 1b: Two-block batch (Block 1 / Block 2 piecewise; rGLMM_sweep)
#  Callers: two_block_batch_gibbs.R, build_mu_all.R
# =============================================================================

#' @noRd
#' @keywords internal
.two_block_build_mu_all_cpp <- function(x_hyper, fixef, re_names, group_levels) {
  .Call(`_glmbayesCore_two_block_build_mu_all_cpp_export`,
    x_hyper, fixef, re_names, group_levels
  )
}


#' @noRd
#' @keywords internal
.two_block_block2_one_chain_cpp <- function(
    b_i, fixef_rows, tau2_i, iters_i, x_hyper, group_levels,
    pfamily_list, ptypes, re_names
) {
  .Call(`_glmbayesCore_two_block_population_one_chain_cpp_export`,
    b_i, fixef_rows, tau2_i, iters_i, x_hyper, group_levels,
    pfamily_list, ptypes, re_names
  )
}


# =============================================================================
#  Independent-block separable-overbound variants (Ind)
#  See Appendix A of inst/BLOCK_ING_RINDEPNORMALGAMMA_REG.md for the theory.
# =============================================================================


#' Block ING independent-block pipeline: Centering → Build → DispersionBuildInd → SimInd
#' @noRd
#' @keywords internal
.rIndepNormalGammaRegBlockInd_cpp <- function(
    n,
    y,
    x,
    block,
    prior_list,
    prior_lists = NULL,
    offset,
    wt,
    p_re = -1L,
    n_rss_iter = 10L,
    Gridtype = 3L,
    n_envopt = -1L,
    RSS_ML = NA_real_,
    use_parallel = TRUE,
    use_opencl = FALSE,
    progbar = FALSE,
    verbose = FALSE,
    group_levels = character(0),
    re_names = character(0)
) {
  .Call(
    `_glmbayesCore_rIndepNormalGammaRegBlockInd_cpp_export`,
    n, y, x, block, prior_list, prior_lists,
    offset, wt, p_re, n_rss_iter, Gridtype, n_envopt, RSS_ML,
    use_parallel, use_opencl, progbar, verbose,
    group_levels, re_names
  )
}
