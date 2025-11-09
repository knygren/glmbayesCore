#' Evaluate Negative Log-Likelihood and Gradients
#'
#' \code{EnvelopeEval()} evaluates the negative log-likelihood and gradients
#' at a grid of parameter values, optionally using OpenCL acceleration.
#'
#' The lower-level helpers \code{f2_f3_non_opencl()} and \code{f2_f3_opencl()}
#' expose the CPU and GPU kernels directly. The function \code{run_opencl_pilot()}
#' benchmarks OpenCL performance on a pilot subset of the grid to estimate
#' runtime before full evaluation.
#'
#' These functions implement the grid evaluation logic used in envelope
#' construction for rejection sampling. They make use of the theory described
#' in \insertCite{Nygren2006}{glmbayes} and the general implementation outlined
#' in \insertCite{glmbayesSimmethods}{glmbayes}.

#'
#' @param G4 Numeric matrix of parameter values (parameters × grid points).
#' @param y Numeric response vector.
#' @param x Numeric design matrix.
#' @param mu Numeric matrix of offsets or prior means.
#' @param P Numeric matrix representing the portion of the prior precision
#'   shifted into the likelihood.
#' @param alpha Numeric vector of prior shape parameters.
#' @param wt Numeric vector of weights.
#' @param family Character string; model family (e.g. \code{"gaussian"}).
#' @param link Character string; link function (e.g. \code{"identity"}).
#' @param use_opencl Logical; if \code{TRUE}, attempt OpenCL acceleration.
#' @param verbose Logical; if \code{TRUE}, print diagnostic output.
#' @param progbar Integer flag for progress bar control (internal use).
#' @param b Numeric matrix of parameter values (parameters × grid points).
#' @param threshold_sec Threshold seconds for run_opencl_pilot. If second exceeds this, a prompt for users is 
#' triggered allowing users to interrupt the run.
#' @details
#' The evaluation workflow has several layers:

#' **1. High-level dispatch (`EnvelopeEval`)**
#'
#' * `EnvelopeEval()` is the user-facing entry point. It accepts a grid of
#'   parameter values (`G4`) and the data (`y`, `x`, `mu`, `P`, `alpha`, `wt`).
#' * If the grid is large (≥ 14 columns), it first calls
#'   [`run_opencl_pilot()`] to benchmark OpenCL performance and optionally
#'   report estimated runtime.
#' * It then dispatches to either the CPU or GPU backend:
#'   - If `use_opencl = TRUE` and the family is not `"gaussian"`, it calls
#'     [`f2_f3_opencl()`].
#'   - Otherwise, it calls [`f2_f3_non_opencl()`].
#'
#' **2. CPU backend (`f2_f3_non_opencl`)**
#'
#' * This function evaluates the negative log-likelihood and gradients using
#'   standard CPU routines.
#' * It inspects the `family` and `link` arguments and routes to the correct
#'   pair of kernels (`f2_*` for the likelihood, `f3_*` for the gradient).
#' * For example:
#'   - `"binomial"` with `"logit"` calls `f2_binomial_logit()` and
#'     `f3_binomial_logit()`.
#'   - `"poisson"` calls `f2_poisson()` and `f3_poisson()`.
#'   - `"gaussian"` calls `f2_gaussian()` and `f3_gaussian()`.
#' * These kernels ultimately rely on the same C math routines that R itself
#'   uses (from the `nmath`/`rmath` libraries), ensuring numerical consistency
#'   with base R functions like `dnorm`, `dpois`, etc.
#'
#' **3. GPU backend (`f2_f3_opencl`)**
#'
#' * This function mirrors the CPU backend but executes the likelihood and
#'   gradient calculations on an OpenCL device (GPU or CPU).
#' * It flattens the input matrices/vectors and allocates output buffers.
#' * It then constructs a full OpenCL program by concatenating:
#'   - a generic OpenCL support header (`OPENCL.CL`),
#'   - OpenCL ports of R's `rmath`, `nmath`, and `dpq` libraries,
#'   - and the family/link-specific kernel source (e.g.
#'     `f2_f3_binomial_logit.cl`).
#' * The resulting program is compiled and passed to a kernel runner
#'   (`f2_f3_kernel_runner`) which executes the likelihood and gradient
#'   calculations in parallel on the device.
#' * This ensures that the GPU backend produces results consistent with the
#'   CPU backend, but can scale to much larger grids efficiently.
#'
#' **4. Pilot timing (`run_opencl_pilot`)**
#'
#' * This helper runs a small subset of the grid through the OpenCL backend
#'   to estimate runtime.
#' * It is used by `EnvelopeEval()` to inform users (when `verbose = TRUE`)
#'   whether OpenCL acceleration is likely to be beneficial.
#'
#' **5. Returned values**
#'
#' * All backends return a list with:
#'   - `NegLL`: numeric vector of negative log-likelihood values.
#'   - `cbars`: numeric matrix of gradients (parameters × grid points).
#' **6. Role of likelihood and gradients in sampling**
#'
#' * The outputs of `EnvelopeEval()` - the negative log-likelihood values
#'   (`NegLL`) and the gradient matrix (`cbars`) - are not endpoints in
#'   themselves. They form the *envelope* used in the rejection sampler
#'   implemented by internal functions such as
#'   `rnnorm_reg_std_cpp()` and `rnnorm_reg_std_cpp_parallel()`.
#' * These routines are called by `rnnorm_reg_cpp()`, which underlies the
#'   user-facing function `rNormal_reg()`. Together they implement
#'   envelope-based posterior sampling for GLMs with log-concave likelihoods
#'   and multivariate normal priors.
#'
#' **7. Simulation execution (accept/reject procedure)**
#'
#' The acceptance test is performed using
#' \deqn{
#'   \log(U_2) \leq
#'   \log f(y \mid \theta_i) -
#'   \Big(\log f(y \mid \bar{\theta}_{J(i)}) -
#'        c(\bar{\theta}_{J(i)})^T(\theta_i - \bar{\theta}_{J(i)})\Big) \leq 0
#' }
#'
#' **Connections between code and notation:**
#' * The arguments `G4` (in `EnvelopeEval`) and `b` (in `f2_f3_*`) both
#'   represent the grid of tangency points \eqn{\bar{\theta}_j}.
#' * The output `NegLL` corresponds to
#'   \eqn{-\log f(y \mid \bar{\theta}_{J(i)})}, i.e. the negative
#'   log-likelihood evaluated at each tangency point.
#' * The output `cbars` corresponds to the subgradient vectors
#'   \eqn{c(\bar{\theta}_{J(i)})}, which define the tangent hyperplanes
#'   used in the envelope construction.
#'
#' **Precomputation for efficiency:**
#' * Both `NegLL` and `cbars` are computed once during envelope construction,
#'   prior to the simulation stage.
#' * This means the sampler does not need to recompute likelihoods or
#'   gradients at every candidate draw - it simply reuses the stored values
#'   (`NegLL`, `cbars`, and `LLconst`) in the acceptance inequality.
#'
#' This design ensures that the envelope is tangent to the log-likelihood at
#' each \eqn{\bar{\theta}_j}, lies above it elsewhere, and that the
#' accept-reject procedure can run efficiently while still producing samples
#' from the true posterior \eqn{\pi(\theta \mid y)}.#'
#'
#' @references
#' \insertAllCited{}


#' @return
#' \describe{
#'   \item{EnvelopeEval}{List with components \code{NegLL} (numeric vector of
#'   negative log-likelihood values) and \code{cbars} (numeric matrix of gradients).}
#'   \item{f2_f3_non_opencl}{List with components \code{qf} (negative log-likelihood)
#'   and \code{grad} (gradients) from the CPU kernel.}
#'   \item{f2_f3_opencl}{List with components \code{qf} and \code{grad} from the
#'   OpenCL kernel.}
#'   \item{run_opencl_pilot}{Numeric scalar giving estimated runtime (seconds)
#'   for OpenCL evaluation on a pilot subset of the grid.}
#' }
#'
#' @rdname EnvelopeEval
#' @export
#' @usage EnvelopeEval(G4, y, x, mu, P, alpha, wt,
#'                     family, link,
#'                     use_opencl = FALSE, verbose = FALSE)
EnvelopeEval <- function(G4, y, x, mu, P, alpha, wt,
                         family, link,
                         use_opencl = FALSE,
                         verbose = FALSE) {
  if (!is.matrix(G4)) stop("G4 must be a numeric matrix")
  if (!is.numeric(y)) stop("y must be numeric")
  if (!is.matrix(x)) stop("x must be a numeric matrix")
  if (!is.matrix(mu)) stop("mu must be a numeric matrix")
  if (!is.matrix(P)) stop("P must be a numeric matrix")
  if (!is.numeric(alpha)) stop("alpha must be numeric")
  if (!is.numeric(wt)) stop("wt must be numeric")
  if (!is.character(family) || length(family) != 1L) stop("family must be a string")
  if (!is.character(link) || length(link) != 1L) stop("link must be a string")
  
  EnvelopeEval(G4, y, x, mu, P, alpha, wt,
               family, link,
               use_opencl, verbose)
}

#' @rdname EnvelopeEval
#' @export
#' @usage f2_f3_non_opencl(family, link, b, y, x, mu, P, alpha, wt, progbar)
f2_f3_non_opencl <- function(family, link, b, y, x, mu, P, alpha, wt, progbar = 0L) {
  f2_f3_non_opencl(family, link, b, y, x, mu, P, alpha, wt, progbar)
}

#' @rdname EnvelopeEval
#' @export
#' @usage f2_f3_opencl(family, link, b, y, x, mu, P, alpha, wt, progbar)
f2_f3_opencl <- function(family, link, b, y, x, mu, P, alpha, wt, progbar = 0L) {
  f2_f3_opencl(family, link, b, y, x, mu, P, alpha, wt, progbar)
}

#' @rdname EnvelopeEval
#' @export
#' @usage run_opencl_pilot(G4, y, x, mu, P, alpha, wt,
#'                         family, link,
#'                         use_opencl, verbose,
#'                         threshold_sec = 300)
run_opencl_pilot <- function(G4, y, x, mu, P, alpha, wt,
                             family, link,
                             use_opencl = FALSE,
                             verbose = FALSE,
                             threshold_sec = 300) {
  run_opencl_pilot(G4, y, x, mu, P, alpha, wt,
                   family, link, use_opencl, verbose, threshold_sec)
}