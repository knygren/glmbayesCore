## Pilot-chain cost trade-off for two-block GLMM sampling (analogous to
## EnvelopeOpt grid sizing: upfront setup vs per-draw inner sweeps).

#' Pilot / main chain cost calibration for two-block GLMM sampling
#'
#' @description
#' Tools for the pilot-vs-main trade-off in non-Gaussian two-block GLMM paths:
#' certificate the start distance when main chains begin at a pilot mean,
#' calibrate \code{m_convergence} from Nygren TV bounds, evaluate total
#' inner-sweep cost, and search for a cost-optimal \code{n_pilot}.
#'
#' @details
#' The functions form a bottom-up stack:
#' \enumerate{
#'   \item \code{\link{two_block_d0_pilot_start}}:
#'     \eqn{D_0 = \chi^2_{p,\alpha} / n_{\mathrm{pilot}}} for a pilot-mean start.
#'   \item \code{\link{two_block_m_convergence_for_pilot_start}}: main-stage
#'     \code{m_convergence} at that \eqn{D_0} (via \code{\link{two_block_l_for_tv}}).
#'   \item \code{\link{two_block_pilot_sampling_cost}}:
#'     \eqn{C = n_{\mathrm{pilot}} m_p + n\, m_{\mathrm{conv}}(n_{\mathrm{pilot}})}.
#'   \item \code{\link{two_block_optimize_pilot_cost}}: minimize \eqn{C} over
#'     \code{n_pilot} (advisory unless wired into \code{rGLMM_reg()} in
#'     \pkg{lmebayesCore}).
#' }
#' This is the Gibbs analogue of \code{\link[glmbayesCore]{EnvelopeOpt}}'s build-plus-sample cost.
#'
#' @param n_pilot Number of pilot replicate chains (positive integer).
#' @param p Dimension of the Block~2 hyper vector (positive integer).
#' @param pilot_start_tol One-sided chi-squared quantile in \code{(0, 1)}
#'   for the \eqn{D_0} certificate (default \code{0.95}).
#' @inheritParams two_block_tv_bound
#' @param tv_tol Target total-variation tolerance in \code{(0, 1)} passed to
#'   \code{\link{two_block_l_for_tv}} when calibrating \code{m_convergence}.
#' @param n Number of main-stage stored draws.
#' @param m_convergence_pilot Inner Gibbs sweeps per pilot stored draw
#'   (\eqn{m_p}).
#' @param n_pilot_min,n_pilot_max Integer bounds for the \code{n_pilot} search
#'   in \code{two_block_optimize_pilot_cost} (inclusive).  When
#'   \code{n_pilot_max} is \code{NULL}, defaults to \code{max(10000L, 10L * n)}.
#' @family simfuncs
#' @seealso \code{\link{two_block_rate}}, \code{\link{two_block_l_for_tv}},
#'   \code{\link[glmbayesCore]{EnvelopeOpt}}, \code{rGLMM_reg()} in \pkg{lmebayesCore}
#' @name two_block_optimize_pilot_cost
#' @aliases two_block_pilot_sampling_cost two_block_m_convergence_for_pilot_start
#'   two_block_d0_pilot_start
#' @keywords internal
NULL

#' @describeIn two_block_optimize_pilot_cost Squared standardized start distance
#'   from a pilot mean start: \eqn{D_0 = \chi^2_{p,\alpha} / n_{\mathrm{pilot}}}.
#' @return Non-negative scalar \eqn{D_0}.
two_block_d0_pilot_start <- function(n_pilot,
                                     p,
                                     pilot_start_tol = 0.95) {
  n_pilot <- as.integer(n_pilot[1L])
  p <- as.integer(p[1L])
  if (length(n_pilot) < 1L || !is.finite(n_pilot) || n_pilot < 1L) {
    stop("'n_pilot' must be a positive integer.", call. = FALSE)
  }
  if (length(p) < 1L || !is.finite(p) || p < 1L) {
    stop("'p' must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(pilot_start_tol) || length(pilot_start_tol) != 1L ||
      pilot_start_tol <= 0 || pilot_start_tol >= 1) {
    stop("'pilot_start_tol' must be a single value in (0, 1).", call. = FALSE)
  }
  stats::qchisq(pilot_start_tol, df = p) / n_pilot
}

#' @describeIn two_block_optimize_pilot_cost Main-stage inner sweeps for a
#'   pilot-mean chain start.  When \code{n_pilot} is large, \code{D_0} is small
#'   and the result approaches the mode-start minimum \code{m_min}.
#' @return List with \code{m_convergence} (integer sweeps per main draw),
#'   \code{m_min} (mode-start minimum), \code{delta_m}, \code{D_0}, and
#'   \code{n_pilot}.
two_block_m_convergence_for_pilot_start <- function(rate,
                                                    n_pilot,
                                                    tv_tol,
                                                    p,
                                                    pilot_start_tol = 0.95,
                                                    method = c("theorem3",
                                                               "corollary1")) {
  if (!inherits(rate, "two_block_rate")) {
    stop("'rate' must be a two_block_rate object.", call. = FALSE)
  }
  method <- match.arg(method)
  n_pilot <- as.integer(n_pilot[1L])
  if (length(n_pilot) < 1L || !is.finite(n_pilot) || n_pilot < 1L) {
    stop("'n_pilot' must be a positive integer.", call. = FALSE)
  }
  D0 <- two_block_d0_pilot_start(n_pilot, p, pilot_start_tol)
  m_min <- .two_block_cap_inner_sweeps(
    two_block_l_for_tv(rate, tv_tol, method = method, D0 = 0, warn = FALSE) + 1L
  )
  m_conv <- .two_block_cap_inner_sweeps(
    two_block_l_for_tv(rate, tv_tol, method = method, D0 = D0, warn = FALSE) + 1L
  )
  list(
    n_pilot         = n_pilot,
    D0              = D0,
    m_min           = m_min,
    m_convergence   = m_conv,
    delta_m         = m_conv - m_min,
    pilot_start_tol = pilot_start_tol
  )
}

#' @describeIn two_block_optimize_pilot_cost Total inner-sweep cost for a
#'   pilot / main sampling plan:
#'   \eqn{C = n_{\mathrm{pilot}} m_p + n\, m_{\mathrm{conv}}(n_{\mathrm{pilot}})}.
#' @return List with \code{total_cost}, \code{pilot_cost}, \code{main_cost},
#'   \code{m_convergence}, and the fields returned by
#'   \code{two_block_m_convergence_for_pilot_start}.
two_block_pilot_sampling_cost <- function(n,
                                          n_pilot,
                                          rate,
                                          tv_tol,
                                          m_convergence_pilot,
                                          p,
                                          pilot_start_tol = 0.95,
                                          method = c("theorem3", "corollary1")) {
  n <- as.integer(n[1L])
  m_convergence_pilot <- as.integer(m_convergence_pilot[1L])
  if (length(n) < 1L || !is.finite(n) || n < 1L) {
    stop("'n' must be a positive integer.", call. = FALSE)
  }
  if (length(m_convergence_pilot) < 1L || !is.finite(m_convergence_pilot) ||
      m_convergence_pilot < 1L) {
    stop("'m_convergence_pilot' must be a positive integer.", call. = FALSE)
  }
  conv <- two_block_m_convergence_for_pilot_start(
    rate            = rate,
    n_pilot         = n_pilot,
    tv_tol          = tv_tol,
    p               = p,
    pilot_start_tol = pilot_start_tol,
    method          = method
  )
  pilot_cost <- as.numeric(n_pilot) * as.numeric(m_convergence_pilot)
  main_cost  <- as.numeric(n) * as.numeric(conv$m_convergence)
  c(
    conv,
    list(
      n                     = n,
      m_convergence_pilot   = m_convergence_pilot,
      pilot_cost            = pilot_cost,
      main_cost             = main_cost,
      total_cost            = pilot_cost + main_cost
    )
  )
}

#' @describeIn two_block_optimize_pilot_cost Search over candidate
#'   \code{n_pilot} values for the minimum of
#'   \code{\link{two_block_pilot_sampling_cost}}.  The returned optimum is
#'   advisory only unless wired into the sampler.
#' @return List with \code{n_pilot_opt}, \code{m_convergence_opt},
#'   \code{total_cost_opt}, \code{cost_at_opt}, and \code{cost_curve} (data
#'   frame of evaluated \code{n_pilot} values near the optimum).
two_block_optimize_pilot_cost <- function(n,
                                          rate,
                                          tv_tol,
                                          m_convergence_pilot,
                                          p,
                                          pilot_start_tol = 0.95,
                                          method = c("theorem3", "corollary1"),
                                          n_pilot_min = 1L,
                                          n_pilot_max = NULL) {
  if (!inherits(rate, "two_block_rate")) {
    stop("'rate' must be a two_block_rate object.", call. = FALSE)
  }
  method <- match.arg(method)
  n <- as.integer(n[1L])
  n_pilot_min <- as.integer(n_pilot_min[1L])
  if (is.null(n_pilot_max)) {
    n_pilot_max <- max(10000L, 10L * n)
  } else {
    n_pilot_max <- as.integer(n_pilot_max[1L])
  }
  if (n_pilot_min < 1L || n_pilot_max < n_pilot_min) {
    stop("'n_pilot_min' and 'n_pilot_max' must satisfy 1 <= min <= max.",
         call. = FALSE)
  }

  cost_at <- function(n1) {
    two_block_pilot_sampling_cost(
      n                   = n,
      n_pilot             = n1,
      rate                = rate,
      tv_tol              = tv_tol,
      m_convergence_pilot = m_convergence_pilot,
      p                   = p,
      pilot_start_tol     = pilot_start_tol,
      method              = method
    )$total_cost
  }

  if (n_pilot_min == n_pilot_max) {
    n_opt <- n_pilot_min
  } else {
    opt_cont <- stats::optimize(
      function(x) cost_at(max(n_pilot_min, as.integer(round(x)))),
      lower = n_pilot_min,
      upper = n_pilot_max
    )
    n_cand <- unique(pmax(
      n_pilot_min,
      pmin(n_pilot_max, as.integer(round(opt_cont$minimum)) + (-1:1L))
    ))
    costs <- vapply(n_cand, cost_at, numeric(1L))
    n_opt <- n_cand[which.min(costs)]
  }

  cost_opt <- two_block_pilot_sampling_cost(
    n                   = n,
    n_pilot             = n_opt,
    rate                = rate,
    tv_tol              = tv_tol,
    m_convergence_pilot = m_convergence_pilot,
    p                   = p,
    pilot_start_tol     = pilot_start_tol,
    method              = method
  )

  curve_n <- unique(pmax(
    n_pilot_min,
    pmin(n_pilot_max, n_opt + (-5L:5L))
  ))
  cost_curve <- data.frame(
    n_pilot    = curve_n,
    total_cost = vapply(curve_n, cost_at, numeric(1L)),
    row.names  = NULL
  )

  list(
    n_pilot_opt       = n_opt,
    m_convergence_opt = cost_opt$m_convergence,
    total_cost_opt    = cost_opt$total_cost,
    cost_at_opt       = cost_opt,
    cost_curve        = cost_curve,
    n_pilot_min       = n_pilot_min,
    n_pilot_max       = n_pilot_max
  )
}

#' Whether a pilot stage will run from family and pilot arguments
#' @param any_non_normal When \code{TRUE} with \code{is_gaussian}, a pilot may
#'   still run (Gaussian LMM with ING Block~2 components).
#' @param random_measurement When \code{TRUE} with \code{is_gaussian}, a pilot
#'   may run for random measurement dispersion (dGamma / ING Block~1).
#' @noRd
.two_block_pilot_will_run <- function(
    is_gaussian,
    n_pilot_arg,
    gap_tol,
    tv_tol,
    any_non_normal = FALSE,
    random_measurement = FALSE
) {
  if (isTRUE(is_gaussian) && !isTRUE(any_non_normal) &&
      !isTRUE(random_measurement)) {
    return(FALSE)
  }
  if (!is.null(n_pilot_arg)) {
    return(as.integer(n_pilot_arg[1L]) > 0L)
  }
  if (!is.null(tv_tol)) {
    return(TRUE)
  }
  !is.null(.two_block_validate_gap_tol(gap_tol))
}

#' Resolve pilot chain count and main-stage \code{m_convergence}
#'
#' When \code{n_pilot} is \code{NULL} and \code{tv_tol} is set, uses
#' \code{\link{two_block_optimize_pilot_cost}}.  Otherwise falls back to
#' explicit \code{n_pilot}, legacy \code{gap_tol}, or no pilot.
#' @noRd
.two_block_resolve_pilot_plan <- function(
    is_gaussian,
    n,
    n_pilot_arg,
    gap_tol,
    tv_tol,
    m_convergence_user,
    m_convergence_pilot,
    rate,
    p_dim,
    m_min = NULL,
    pilot_start_tol = 0.95,
    n_pilot_max = NULL,
    any_non_normal = FALSE,
    random_measurement = FALSE) {
  gap_tol_validated <- .two_block_validate_gap_tol(gap_tol)
  n_pilot_gap_tol <- if (!is.null(gap_tol_validated)) {
    as.integer(ceiling((stats::qnorm(0.975) / gap_tol_validated)^2))
  } else {
    NULL
  }

  if (isTRUE(is_gaussian) && !isTRUE(any_non_normal) &&
      !isTRUE(random_measurement)) {
    m_conv <- if (!is.null(m_convergence_user)) {
      as.integer(m_convergence_user[1L])
    } else if (!is.null(m_min)) {
      m_min
    } else {
      10L
    }
    return(list(
      n_pilot          = 0L,
      n_pilot_source   = "gaussian",
      m_convergence    = m_conv,
      m_certificate    = NULL,
      pilot_cost_opt   = NULL,
      n_pilot_gap_tol  = n_pilot_gap_tol
    ))
  }

  pilot_cost_opt <- NULL
  n_pilot_source <- "none"

  if (!is.null(n_pilot_arg)) {
    n_pilot <- as.integer(n_pilot_arg[1L])
    if (n_pilot < 0L) {
      stop("'n_pilot' must be non-negative.", call. = FALSE)
    }
    n_pilot_source <- if (n_pilot > 0L) "explicit" else "none"
  } else if (!is.null(tv_tol) && !is.null(m_convergence_pilot)) {
    pilot_cost_opt <- two_block_optimize_pilot_cost(
      n                   = n,
      rate                = rate,
      tv_tol              = tv_tol,
      m_convergence_pilot = m_convergence_pilot,
      p                   = p_dim,
      pilot_start_tol     = pilot_start_tol,
      n_pilot_max         = n_pilot_max
    )
    n_pilot <- pilot_cost_opt$n_pilot_opt
    n_pilot_source <- "cost"
  } else if (!is.null(gap_tol_validated)) {
    n_pilot <- n_pilot_gap_tol
    n_pilot_source <- "gap_tol"
  } else {
    n_pilot <- 0L
    n_pilot_source <- "none"
  }

  run_pilot <- n_pilot > 0L
  m_certificate <- NULL
  if (run_pilot && !is.null(tv_tol)) {
    m_certificate <- two_block_m_convergence_for_pilot_start(
      rate            = rate,
      n_pilot         = n_pilot,
      tv_tol          = tv_tol,
      p               = p_dim,
      pilot_start_tol = pilot_start_tol
    )$m_convergence
  }

  if (!is.null(m_convergence_user)) {
    m_convergence <- as.integer(m_convergence_user[1L])
  } else if (run_pilot && !is.null(tv_tol)) {
    if (identical(n_pilot_source, "cost") && !is.null(pilot_cost_opt)) {
      m_convergence <- pilot_cost_opt$m_convergence_opt
    } else {
      m_convergence <- m_certificate
    }
  } else if (!is.null(m_min)) {
    m_convergence <- m_min
  } else {
    m_convergence <- 10L
  }

  if (!is.null(m_convergence_user)) {
    if (run_pilot && !is.null(m_certificate) && m_convergence < m_certificate) {
      warning(
        "rGLMM: m_convergence = ", m_convergence_user,
        " is below the pilot-start certificate ", m_certificate,
        " at n_pilot = ", n_pilot, "; using ", m_certificate, ".",
        call. = FALSE
      )
      m_convergence <- m_certificate
    } else if (!run_pilot && !is.null(m_min) && m_convergence < m_min) {
      warning(
        "rGLMM: m_convergence = ", m_convergence_user,
        " is below the mode-start minimum m_min = ", m_min,
        "; using m_min instead.",
        call. = FALSE
      )
      m_convergence <- m_min
    }
  }

  m_convergence <- .two_block_cap_inner_sweeps(m_convergence)
  if (!is.null(m_certificate)) {
    m_certificate <- .two_block_cap_inner_sweeps(m_certificate)
  }

  list(
    n_pilot          = n_pilot,
    n_pilot_source   = n_pilot_source,
    m_convergence    = m_convergence,
    m_certificate    = m_certificate,
    pilot_cost_opt   = pilot_cost_opt,
    n_pilot_gap_tol  = n_pilot_gap_tol
  )
}

#' Print resolved pilot / main sampling plan (before pilot stage)
#' @noRd
.two_block_print_pilot_plan <- function(pilot_plan,
                                      n,
                                      m_convergence_pilot,
                                      rate,
                                      tv_tol,
                                      p,
                                      pilot_start_tol = 0.95,
                                      verbose = TRUE) {
  if (!isTRUE(verbose) || is.null(pilot_plan) || pilot_plan$n_pilot < 1L) {
    return(invisible(NULL))
  }
  total_cost <- two_block_pilot_sampling_cost(
    n                   = n,
    n_pilot             = pilot_plan$n_pilot,
    rate                = rate,
    tv_tol              = tv_tol,
    m_convergence_pilot = m_convergence_pilot,
    p                   = p,
    pilot_start_tol     = pilot_start_tol
  )$total_cost
  cat(sprintf(
    paste0(
      "--- rGLMM: pilot / main plan [source = %s; n = %d, m_pilot = %d]:\n",
      "    n_pilot = %d => m_convergence = %d ",
      "(total inner sweeps = %d)\n"
    ),
    pilot_plan$n_pilot_source, n, m_convergence_pilot,
    pilot_plan$n_pilot, pilot_plan$m_convergence, total_cost
  ))
  flush.console()
  invisible(NULL)
}

