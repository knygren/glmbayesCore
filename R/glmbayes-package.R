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
#' **Releases:** version **0.9.5** is on CRAN (`install.packages("glmbayes")`).
#' Source is available from GitHub; R-Universe (\url{https://knygren.r-universe.dev/glmbayes})
#' also builds binaries from that source.
#' Prebuilt CRAN and R-Universe binaries do not include OpenCL; for GPU support
#' with the CRAN release, install from source (see vignette *Chapter 12*).
#'
#' IID posterior simulation for non-Gaussian GLMs and several non-conjugate
#' linear-model setups uses the likelihood-subgradient envelope method of
#' \insertCite{Nygren2006}{glmbayes}. Introductory material and worked
#' examples are in \insertCite{glmbayesChapter00,glmbayesChapterA01}{glmbayes};
#' estimation and simulation background in
#' \insertCite{glmbayesChapterA02,glmbayesSimmethods,glmbayesChapterA08}{glmbayes};
#' prior derivations for \code{Prior_Setup()} in
#' \insertCite{glmbayesChapterA12}{glmbayes};
#' GPU/OpenCL topics in
#' \insertCite{glmbayesChapter12,glmbayesChapterA10}{glmbayes}.
#'
#' @example inst/examples/Ex_glmbayes-package.R
#'
#' @seealso
#' Main interfaces: \code{\link{glmb}}, \code{\link{lmb}},
#' \code{\link{rglmb}}, \code{\link{rlmb}}; low-level simulation API
#' \code{\link{simfuncs}}; envelope construction \code{\link{EnvelopeBuild}}.
#'
#' Useful links:
#' \itemize{
#'   \item CRAN: <https://CRAN.R-project.org/package=glmbayes>
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
#' @import stats Rcpp
#' @importFrom Rcpp evalCpp
#' @importFrom MASS mvrnorm
#' @importFrom Rdpack reprompt
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib glmbayes, .registration = TRUE
"_PACKAGE"