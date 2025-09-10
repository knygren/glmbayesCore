#' Builds Envelope function for simulation
#'
#' Builds an enveloping function for simulation using a grid and tangencies for
#' the posterior density. To construct this enveloping function, we follow
#' the approach in \insertCite{Nygren2006}{glmbayes}, which involves the
#' following steps when a maximally sized grid is constructed (if the prior
#' for some dimensions is relatively strong, this may not be needed):
#'
#' 1) For each dimension, a constant \code{omega_i} is found that depends on the
#'    corresponding diagonal element in the precision matrix.
#'
#' 2) Corresponding intervals \code{(thetastar_i - 0.5 * omega_i,
#'    thetastar_i + 0.5 * omega_i)} are constructed around the posterior mode
#'    \code{thetastar} for each dimension.
#'
#' 3) The mode as well as the points
#'    \code{thetastar_i - omega_i} and \code{thetastar_i + omega_i} are selected
#'    as the components of the points at which tangencies will be found.
#'
#' 4) A grid is constructed with all possible combinations of points and negative
#'    log-likelihood and gradient for the negative log-likelihood are evaluated.
#'
#' 5) The \code{\link{Set_Grid}} function is called to evaluate the log-density
#'    of each restricted multivariate normal by taking differences between its
#'    CDF at the lower and upper bounds.
#'
#' 6) The \code{\link{setlogP}} function is called to set sampling probabilities
#'    for each grid component, following Remark 6 in
#'    \insertCite{Nygren2006}{glmbayes}.
#'
#' Any constants needed by the sampling are added to a list and returned.
#'
#' @param bStar     Point at which envelope should be centered (typically posterior mode)
#' @param A         Diagonal precision matrix for the log-likelihood in standard form
#' @param y         A vector of observations of length \code{m}
#' @param x         A design matrix of dimension \code{m * p}
#' @param mu        A vector giving the prior means of the variables
#' @param P         Prior precision matrix of the variables (positive-definite)
#' @param alpha     Offset vector
#' @param wt        A vector of weights
#' @param family    Family for the envelope: binomial, quasibinomial, poisson, quasipoisson, or Gamma
#' @param link      Link function ("logit", "probit", "cloglog" for binomial; "log" for Poisson/Gamma)
#' @param Gridtype  Method to determine the number of subgradient densities in the grid
#' @param n         Number of draws from the posterior (used for grid sizing)
#' @param sortgrid  Logical; if TRUE, sort the envelope descending by component probability
#' @param use_opencl Logical; if TRUE, use OpenCL for gradient evaluations
#' @param verbose   Logical; if TRUE, print progress messages
#' @return A list with elements:
#'   \item{GridIndex}{Matrix indicating tail/center/line sampling per dimension}
#'   \item{thetabars}{Matrix of tangency points for each grid component}
#'   \item{cbars}{Matrix of negative log-likelihood gradients at each tangency}
#'   \item{logU}{Matrix of log CDF differences per dimension}
#'   \item{logrt}{Matrix of log right-tail probabilities}
#'   \item{loglt}{Matrix of log left-tail probabilities}
#'   \item{LLconst}{Vector of constants for each component used in the accept–reject step}
#'   \item{logP}{Matrix of log-probabilities for each grid component}
#'   \item{PLSD}{Vector of sampling probabilities for the grid components}
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @export
EnvelopeBuild <- function(
    bStar, A, y, x, mu, P, alpha, wt,
    family = "binomial", link = "logit",
    Gridtype = 2L, n = 1L, sortgrid = FALSE,
    use_opencl = FALSE, verbose = FALSE
) {
  if (family == "gaussian") {
    return(.EnvelopeBuild_Ind_Normal_Gamma(
      bStar, A, y, x, mu, P, alpha, wt,
      family = family, link = link,
      Gridtype = Gridtype, n = n, sortgrid
    ))
  }
  .EnvelopeBuild_cpp(
    bStar, A, y, x, mu, P, alpha, wt,
    family = family, link = link,
    Gridtype = Gridtype, n = n, sortgrid,
    use_opencl = use_opencl, verbose = verbose
  )
}