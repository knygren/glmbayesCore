#' @aliases glmbayes
#'
#' @title glmbayes: Bayesian Generalized Linear Models with iid Sampling
#'
#' @description
#' `glmbayes` provides independent and identically distributed (iid) samples for
#' Bayesian generalized linear models (GLMs), serving as a Bayesian analogue to
#' the base `glm()` function. Supported likelihood families include Gaussian,
#' Poisson, Binomial, and Gamma models with log-concave likelihoods.
#'
#' @details
#' The main user-facing interface is `glmb()`, which mirrors the structure of
#' `glm()` and supports prior specification through `pfamily` objects. Lower-level
#' functions such as `rglmb()` and `rGamma_reg()` provide direct access to the
#' underlying samplers and can be used in block Gibbs sampling or hierarchical
#' model implementations.
#'
#' For an introduction to the package, examples, and a complete set of vignettes,
#' see:
#'
#' - README: <https://github.com/knygren/glmbayes#readme>
#' - All vignettes: `browseVignettes("glmbayes")`
#'
#' The package includes extensive documentation on model fitting, prior
#' construction, diagnostics, and optional GPU acceleration using OpenCL.
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item GitHub: <https://github.com/knygren/glmbayes>
#'   \item R-Universe: <https://knygren.r-universe.dev/glmbayes>
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @author
#' Kjell Nygren
#'
#' @import stats Rcpp RcppArmadillo
#' @importFrom Rcpp evalCpp
#' @importFrom MASS mvrnorm
#' @importFrom Rdpack reprompt
#' @useDynLib glmbayes
"_PACKAGE"