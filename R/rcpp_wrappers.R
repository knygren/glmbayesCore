# -------------------------------------------------------------------------
#  Rcpp Interface Wrappers for glmbayes
#
#  These functions provide the minimal, strictly positional R → C++ bridges
#  required by the package.  Each wrapper mirrors the exact argument order
#  expected by the corresponding C++ routine and performs no preprocessing,
#  validation, or postprocessing.  Their sole purpose is to ensure that
#  high‑level R code calls the correct compiled symbol with the correct
#  signature.
#
#  All wrappers are internal:
#    - They are not part of the public API.
#    - They exist only to guarantee stable, explicit R–C++ boundaries.
#    - They prevent accidental reliance on .Call() with named arguments,
#      which R ignores, and which can silently break when signatures change.
#
#  Any future C++ interface changes must be reflected here to maintain
#  positional consistency and avoid NULL → double coercion errors.
# -------------------------------------------------------------------------



#' @noRd
#' @keywords internal

.rnnorm_reg_std_cpp <- function(n, y, x, mu, P, alpha, wt,
                                f2, Envelope,
                                family, link,
                                progbar = 1L,
                                verbose = FALSE) {
  .Call(`_glmbayes_rNormalGLM_std_cpp_export`,
        n, y, x, mu, P, alpha, wt,
        f2, Envelope,
        family, link,
        progbar, verbose)
}


#' @noRd
#' @keywords internal

.rnorm_reg_cpp <- function(
    n, y, x, mu, P, offset, wt, dispersion,
    f2, f3, start,
    family = "gaussian",
    link = "identity",
    Gridtype = 2
) {
  .Call(
    "_glmbayes_rNormalReg_cpp_export",
    n, y, x, mu, P, offset, wt, dispersion,
    f2, f3, start,
    family, link, Gridtype
  )
}

#' @noRd
#' @keywords internal

.rindep_norm_gamma_reg_cpp <- function(n, y, x, mu, P, offset, wt, shape, rate, max_disp_perc, disp_lower, disp_upper, Gridtype, n_envopt, use_parallel, use_opencl, verbose, progbar) {
  .Call(`_glmbayes_rIndepNormalGammaReg_cpp_export`, n, y, x, mu, P, offset, wt, shape, rate, max_disp_perc, disp_lower, disp_upper, Gridtype, n_envopt, use_parallel, use_opencl, verbose, progbar)
}


##

#' @noRd
#' @keywords internal

.rindep_norm_gamma_reg_std_cpp <- function(n, y, x, mu, P, alpha, wt, f2, Envelope, gamma_list, UB_list, family, link, progbar, verbose) {
  .Call(`_glmbayes_rIndepNormalGammaReg_std_cpp_export`, n, y, x, mu, P, alpha, wt, f2, Envelope, gamma_list, UB_list, family, link, progbar, verbose)
}

.rindep_norm_gamma_reg_std_parallel_cpp <- function(n, y, x, mu, P, alpha, wt, f2, Envelope, gamma_list, UB_list, family, link, progbar, verbose) {
  .Call(`_glmbayes_rIndepNormalGammaReg_std_parallel_cpp_export`, n, y, x, mu, P, alpha, wt, f2, Envelope, gamma_list, UB_list, family, link, progbar, verbose)
}


.glmb_Standardize_Model_cpp <- function(y, x, P, bstar, A1) {
  .Call(`_glmbayes_glmb_Standardize_Model_cpp_export`, y, x, P, bstar, A1)
}


.rnnorm_reg_cpp <- function(n, y, x, mu, P, offset, wt, dispersion, f2, f3, start, family = "binomial", link = "logit", Gridtype = 2L, n_envopt = -1L, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE) {
  .Call(`_glmbayes_rNormalGLM_cpp_export`, n, y, x, mu, P, offset, wt, dispersion, f2, f3, start, family, link, Gridtype, n_envopt, use_parallel, use_opencl, verbose)
}






#' @noRd
#' @keywords internal



.EnvelopeSize <- function(a, G1, Gridtype, n, n_envopt, use_opencl, verbose) {
  .Call(`_glmbayes_EnvelopeSize_cpp_export`, a, G1, Gridtype, n, n_envopt, use_opencl, verbose)
}

#' @noRd
#' @keywords internal


.EnvelopeBuild_cpp<- function(bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, n_envopt, sortgrid, use_opencl, verbose) {
  .Call(`_glmbayes_EnvelopeBuild_cpp_export`, bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, n_envopt, sortgrid, use_opencl, verbose)
}


#' @noRd
#' @keywords internal

.EnvelopeBuild_Ind_Normal_Gamma <- function(bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, n_envopt, sortgrid, use_opencl, verbose) {
  .Call(`_glmbayes_EnvelopeBuild_Ind_Normal_Gamma_cpp_export`, bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, n_envopt, sortgrid, use_opencl, verbose)
}



#' @noRd
#' @keywords internal

.EnvelopeEval <- function(G4, y, x, mu, P, alpha, wt,
                          family, link,
                          use_opencl = FALSE,
                          verbose = FALSE) {
  .Call(`_glmbayes_EnvelopeEval_cpp_export`,
        G4, y, x, mu, P, alpha, wt,
        family, link,
        use_opencl, verbose)
}





#' @noRd
#' @keywords internal
.load_kernel_source_wrapper <- function(relative_path, package = "glmbayes") {
  .Call(`_glmbayes_load_kernel_source_wrapper_cpp_export`, relative_path, package)
}


#' @noRd
#' @keywords internal
.load_kernel_library_wrapper <- function(subdir, package = "glmbayes", verbose = FALSE) {
  .Call(`_glmbayes_load_kernel_library_wrapper_cpp_export`, subdir, package, verbose)
}


#' @noRd
#' @keywords internal
.get_opencl_core_count <- function() {
  .Call("_glmbayes_get_opencl_core_count_cpp_export")
}



#' Internal wrapper for RSS face evaluation
#'
#' @keywords internal
#' @noRd
.rss_face_at_disp <- function(dispersion, cache, cbars_j, y, x, alpha, wt) {
  .Call(`_glmbayes_rss_face_at_disp_cpp_export`,
        dispersion, cache, cbars_j, y, x, alpha, wt)
}



#' Internal wrapper for UB2 face evaluation
#'
#' @keywords internal
#' @noRd
.UB2 <- function(dispersion, cache, cbars_j, y, x, alpha, wt, rss_min_global) {
  .Call(`_glmbayes_UB2_cpp_export`,
        dispersion, cache, cbars_j, y, x, alpha, wt, rss_min_global)
}



#' @noRd
#' @keywords internal
.EnvelopeEval <- function(
    G4,
    y,
    x,
    mu,
    P,
    alpha,
    wt,
    family,
    link,
    use_opencl = FALSE,
    verbose = FALSE
) {
  .Call(
    "_glmbayes_EnvelopeEval_cpp_export",
    G4, y, x, mu, P, alpha, wt,
    family, link,
    use_opencl, verbose
  )
}


#' @noRd
#' @keywords internal
.EnvelopeDispersionBuild_cpp <- function(
    Env,
    Shape,
    Rate,
    P,
    y,
    x,
    alpha,
    n_obs,
    RSS_post,
    RSS_ML,
    mu,
    wt,
    max_disp_perc,
    disp_lower = NULL,
    disp_upper = NULL,
    verbose = FALSE,
    use_parallel = TRUE
) {
  .Call(
    "_glmbayes_EnvelopeDispersionBuild_cpp_export",
    Env,
    Shape,
    Rate,
    P,
    y,
    x,
    alpha,
    n_obs,
    RSS_post,
    RSS_ML,
    mu,
    wt,
    max_disp_perc,
    disp_lower,
    disp_upper,
    verbose,
    use_parallel
  )
}

#' @noRd
#' @keywords internal

.EnvelopeOrchestrator_cpp <- function(bstar2, A, y, x2, mu2, P2, alpha, wt, n, Gridtype, n_envopt, shape, rate, RSS_Post2, RSS_ML, max_disp_perc, disp_lower, disp_upper, use_parallel, use_opencl, verbose) {
  .Call(`_glmbayes_EnvelopeOrchestrator_cpp_export`, bstar2, A, y, x2, mu2, P2, alpha, wt, n, Gridtype, n_envopt, shape, rate, RSS_Post2, RSS_ML, max_disp_perc, disp_lower, disp_upper, use_parallel, use_opencl, verbose)
}


.Set_Grid_cpp <- function(GIndex, cbars, Lint) {
  .Call(`_glmbayes_EnvelopeSet_Grid_cpp_export`, GIndex, cbars, Lint)
}

.setlogP_cpp <- function(logP, NegLL, cbars, G3) {
  .Call(`_glmbayes_EnvelopeSet_LogP_cpp_export`, logP, NegLL, cbars, G3)
}



#' @noRd
#' @keywords internal
.has_opencl <- function() {
  .Call("_glmbayes_has_opencl_cpp_export")
}

#' @noRd
#' @keywords internal
.gpu_names <- function() {
  .Call("_glmbayes_gpu_names_cpp_export")
}



