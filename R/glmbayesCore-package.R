#' @aliases glmbayesCore
#'
#' @title glmbayesCore: Core C++ Sampling Engine for glmbayes
#'
#' @description
#' Core C++ engine for envelope-based iid samplers, Gibbs building blocks, and
#' optional OpenCL acceleration. Developer backend for \pkg{glmbayes},
#' \pkg{lmebayes}, and related extensions. End users should install
#' \pkg{glmbayes} for formula-based modelling and S3 methods.
#'
#' @details
#' Low-level entry points include envelope construction
#' (\code{\link{EnvelopeBuild}}, \code{\link{EnvelopeOrchestrator}}),
#' registered simulation pipelines (\code{\link{rNormalGLM_std}},
#' \code{\link{rIndepNormalGammaReg_std}}), matrix-input samplers
#' (\code{\link{rglmb}}, \code{\link{rlmb}}), and OpenCL kernel loaders.
#' Formula interfaces \code{glmb} and \code{lmb} live in \pkg{glmbayes};
#' \pkg{glmbayesCore} supplies the sampling engine.
#'
#' @section OpenCL startup checks:
#' In interactive sessions, attaching the package with \code{library(glmbayesCore)}
#' may emit a short \code{\link{packageStartupMessage}} when \code{glmbayesCore_has_opencl()}
#' is \code{FALSE} but a GPU or OpenCL stack appears available on the host.
#' Set \code{options(glmbayes.quiet_opencl_startup = TRUE)} to suppress attach
#' notes (recommended for CI and \command{R CMD check}).
#'
#' @seealso \pkg{glmbayes} for the end-user modelling package.
#'
#' @references
#' \insertAllCited{}
#'
#' @import stats Rcpp
#' @importFrom Rcpp evalCpp
#' @importFrom MASS mvrnorm
#' @importFrom RcppParallel RcppParallelLibs
#' @import opencltools
#' @import nmathopencl
#' @useDynLib glmbayesCore, .registration = TRUE
"_PACKAGE"
