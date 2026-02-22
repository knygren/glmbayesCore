#' Envelope Construction Orchestrator for Bayesian Gaussian Regression
#'
#' @description
#' `EnvelopeOrchestrator()` provides a unified interface for constructing the
#' fixed‑dispersion and dispersion‑aware envelopes used in likelihood‑subgradient
#' simulation for Bayesian Gaussian regression with Normal–Gamma priors.
#'
#' This function coordinates:
#'
#' * fixed‑dispersion envelope construction via [`EnvelopeBuild()`],
#' * dispersion‑refined envelope construction via [`EnvelopeDispersionBuild()`],
#' * envelope sorting and reindexing via [`EnvelopeSort()`], and
#' * UB‑list alignment (reordered `lg_prob_factor` and `UB2min`).
#'
#' It is typically used inside higher‑level simulation routines such as
#' [`rindepNormalGamma_reg()`], but may also be called directly for
#' diagnostics, envelope visualization, or custom simulation workflows.
#'
#' @param bstar2 Numeric vector. Posterior mode of the standardized regression
#'   coefficients (from the standardized model).
#' @param A Numeric matrix. Posterior precision matrix (Hessian) at the mode.
#' @param y Numeric response vector of length \code{m}.
#' @param x2 Numeric matrix of standardized predictors (\code{m × p}).
#' @param mu2 Numeric vector. Standardized prior mean (typically a zero vector).
#' @param P2 Numeric matrix. Standardized prior precision component moved into
#'   the log‑likelihood.
#' @param alpha Numeric vector. Offset‑adjusted mean component.
#' @param wt Numeric vector of prior weights.
#' @param n Integer. Number of envelope grid points or simulation draws.
#' @param Gridtype Integer specifying the envelope grid construction method.
#' @param n_envopt Optional integer. Effective sample size passed to
#'   `EnvelopeOpt` during grid construction. Larger values encourage tighter
#'   envelopes.
#' @param shape Numeric. Shape parameter of the Gamma prior for the dispersion.
#' @param rate Numeric. Rate parameter of the Gamma prior for the dispersion.
#' @param RSS_Post2 Numeric. Posterior residual sum of squares used for
#'   dispersion anchoring.
#' @param RSS_ML Numeric. Maximum‑likelihood residual sum of squares.
#' @param max_disp_perc Numeric in \code{(0,1)}. Tail probability used to
#'   determine dispersion bounds when not explicitly supplied.
#' @param disp_lower Optional numeric. Lower bound for the dispersion
#'   (\eqn{\sigma^2}). If supplied, overrides quantile‑based bounds.
#' @param disp_upper Optional numeric. Upper bound for the dispersion
#'   (\eqn{\sigma^2}). Must be strictly greater than \code{disp_lower}.
#' @param use_parallel Logical. Whether to allow parallel computation inside
#'   [`EnvelopeDispersionBuild()`].
#' @param use_opencl Logical. Whether to allow OpenCL acceleration inside
#'   [`EnvelopeBuild()`].
#' @param verbose Logical. Whether to print detailed progress and timing
#'   messages.
#'
#' @return
#' A list with components:
#'
#' \describe{
#'   \item{\code{Env}}{The fully constructed and sorted envelope, including the
#'     PLSD component inserted by the dispersion‑aware refinement step.}
#'   \item{\code{gamma_list}}{Updated Gamma‑prior parameters for the dispersion
#'     (shape, rate, and dispersion bounds).}
#'   \item{\code{UB_list}}{Updated UB‑list including reordered
#'     \code{lg_prob_factor} and \code{UB2min}.}
#'   \item{\code{diagnostics}}{Diagnostic quantities returned by
#'     [`EnvelopeDispersionBuild()`], useful for debugging or envelope
#'     visualization.}
#'   \item{\code{low}}{Lower dispersion bound used.}
#'   \item{\code{upp}}{Upper dispersion bound used.}
#' }
#'
#' @details
#' `EnvelopeOrchestrator()` consolidates the envelope‑related steps that were
#' previously distributed across multiple R and C++ routines.  
#' It provides a stable, high‑level interface for envelope construction and
#' reduces the number of exported C++ functions required by the package.
#'
#' The function does **not** perform simulation.  
#' Simulation should be carried out afterward using either
#' `.rindep_norm_gamma_reg_std_cpp()` or
#' `.rindep_norm_gamma_reg_std_parallel_cpp()`, depending on the
#' \code{use_parallel} flag.
#'
#' @seealso
#' * [`EnvelopeBuild()`] – fixed‑dispersion envelope construction  
#' * [`EnvelopeDispersionBuild()`] – dispersion‑aware envelope refinement  
#' * [`EnvelopeSort()`] – envelope sorting and reindexing  
#' * [`rindepNormalGamma_reg()`] – full Normal–Gamma simulation routine  
#'
#' @examples
#' \dontrun{
#' env_out <- EnvelopeOrchestrator(
#'   bstar2 = bstar2,
#'   A = A,
#'   y = y,
#'   x2 = x2,
#'   mu2 = mu2,
#'   P2 = P2,
#'   alpha = alpha,
#'   wt = wt,
#'   n = 200,
#'   Gridtype = 2,
#'   n_envopt = NULL,
#'   shape = 2,
#'   rate = 1,
#'   RSS_Post2 = RSS_Post2,
#'   RSS_ML = RSS_ML,
#'   max_disp_perc = 0.99,
#'   disp_lower = NULL,
#'   disp_upper = NULL,
#'   use_parallel = TRUE,
#'   use_opencl = FALSE,
#'   verbose = TRUE
#' )
#'
#' Env3 <- env_out$Env
#' }
#'
#' @export

EnvelopeOrchestrator <- function(bstar2,
                                 A,
                                 y,
                                 x2,
                                 mu2,
                                 P2,
                                 alpha,
                                 wt,
                                 n,
                                 Gridtype,
                                 n_envopt,
                                 shape,
                                 rate,
                                 RSS_Post2,
                                 RSS_ML,
                                 max_disp_perc,
                                 disp_lower,
                                 disp_upper,
                                 use_parallel = TRUE,
                                 use_opencl  = FALSE,
                                 verbose     = FALSE) {
  
  # --- NEW: Call the C++ orchestrator directly ---
  out_cpp <- .EnvelopeOrchestrator_cpp(
    bstar2      = bstar2,
    A           = A,
    y           = y,
    x2          = x2,
    mu2         = as.matrix(mu2, ncol = 1),
    P2          = P2,
    alpha       = alpha,
    wt          = wt,
    n           = n,
    Gridtype    = Gridtype,
    n_envopt    = n_envopt,
    shape       = shape,
    rate        = rate,
    RSS_Post2   = RSS_Post2,
    RSS_ML      = RSS_ML,
    max_disp_perc = max_disp_perc,
    disp_lower  = disp_lower,
    disp_upper  = disp_upper,
    use_parallel = use_parallel,
    use_opencl   = use_opencl,
    verbose      = verbose
  )
  
  if (verbose) {
    cat("[EnvelopeOrchestrator] Using C++ orchestrator output\n")
  }
  
  return(out_cpp)
  
  

  
  } 
 


