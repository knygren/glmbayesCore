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
.rnnorm_reg_cpp <- function(
    n, y, x, mu, P, offset, wt, dispersion,
    f2, f3, start, family, link, Gridtype,
    n_envopt, use_parallel, use_opencl, verbose
) {
  .Call(
    "_glmbayes_rnnorm_reg_cpp",
    n, y, x, mu, P, offset, wt, dispersion,
    f2, f3, start, family, link, Gridtype,
    n_envopt, use_parallel, use_opencl, verbose
  )
}

.rnnorm_reg_std_cpp <- function(n, y, x, mu, P, alpha, wt,
                                f2, Envelope,
                                family, link,
                                progbar = 1L,
                                verbose = FALSE) {
  .Call(`_glmbayes_rnnorm_reg_std_cpp`,
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
    "_glmbayes_rnorm_reg_cpp",
    n, y, x, mu, P, offset, wt, dispersion,
    f2, f3, start,
    family, link, Gridtype
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
    "_glmbayes_EnvelopeDispersionBuild_cpp",
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
.EnvelopeSize <- function(
    a,
    G1,
    Gridtype,
    n,
    n_envopt,
    use_opencl,
    verbose
) {
  .Call(
    "_glmbayes_EnvelopeSize",
    a,
    G1,
    Gridtype,
    n,
    n_envopt,
    use_opencl,
    verbose
  )
}


#' @noRd
#' @keywords internal

.EnvelopeEval <- function(G4, y, x, mu, P, alpha, wt,
                          family, link,
                          use_opencl = FALSE,
                          verbose = FALSE) {
  .Call(`_glmbayes_EnvelopeEval`,
        G4, y, x, mu, P, alpha, wt,
        family, link,
        use_opencl, verbose)
}



#' @noRd
#' @keywords internal
.load_kernel_source_wrapper <- function(relative_path, package = "glmbayes") {
  .Call(`_glmbayes_load_kernel_source_wrapper`, relative_path, package)
}


#' @noRd
#' @keywords internal
.load_kernel_library_wrapper <- function(subdir, package = "glmbayes", verbose = FALSE) {
  .Call(`_glmbayes_load_kernel_library_wrapper`, subdir, package, verbose)
}


#' @noRd
#' @keywords internal
.get_opencl_core_count <- function() {
  .Call("_glmbayes_get_opencl_core_count")
}



#' Internal wrapper for RSS face evaluation
#'
#' @keywords internal
#' @noRd
.rss_face_at_disp <- function(dispersion, cache, cbars_j, y, x, alpha, wt) {
  .Call(`_glmbayes_rss_face_at_disp`,
        dispersion, cache, cbars_j, y, x, alpha, wt)
}



#' Internal wrapper for UB2 face evaluation
#'
#' @keywords internal
#' @noRd
.UB2 <- function(dispersion, cache, cbars_j, y, x, alpha, wt, rss_min_global) {
  .Call(`_glmbayes_UB2`,
        dispersion, cache, cbars_j, y, x, alpha, wt, rss_min_global)
}


#' @noRd
#' @keywords internal
.EnvelopeSize <- function(
    a,
    G1,
    Gridtype,
    n,
    n_envopt,
    use_opencl,
    verbose
) {
  .Call(
    "_glmbayes_EnvelopeSize",
    a,
    G1,
    Gridtype,
    n,
    n_envopt,
    use_opencl,
    verbose
  )
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
    "_glmbayes_EnvelopeEval",
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
    "_glmbayes_EnvelopeDispersionBuild_cpp",
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
.has_opencl <- function() {
  .Call("_glmbayes_has_opencl")
}

#' @noRd
#' @keywords internal
.gpu_names <- function() {
  .Call("_glmbayes_gpu_names")
}



