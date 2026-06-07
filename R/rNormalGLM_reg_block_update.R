#' One Gibbs block update via \code{block_rNormalGLM}
#'
#' Draw a single blockwise GLM posterior sample (\code{n = 1}) and return
#' latent-level draws (e.g. updated \code{theta}) for two-block Gibbs samplers
#' where each observation or group is its own block with scalar (or vector)
#' coefficients per block.
#'
#' This is a thin wrapper around \code{\link{block_rNormalGLM}} that always
#' uses \code{n = 1}. When \code{prior_lists} and \code{prior_list} are both
#' omitted, per-block \code{dNormal} priors are built from \code{mu_all} and
#' \code{sigma_theta_sq} (scalar intercept per block).
#'
#' @param mu_all Length-\code{k} vector of prior means per block (linear predictors
#'   from the population block). Required unless \code{prior_lists} or
#'   \code{prior_list} is supplied.
#' @param sigma_theta_sq Shared prior variance for scalar blocks when building
#'   \code{prior_lists} from \code{mu_all}.
#' @param seed Optional; passed to \code{\link{set.seed}} before sampling.
#' @param theta_coef_col Column index of \code{coefficients} to return as
#'   \code{theta} (default \code{1} for scalar intercept blocks).
#' @return A list with:
#'   \describe{
#'     \item{theta}{Vector of length \code{nrow(x)} (one draw per row/block).}
#'     \item{coefficients,coef.mode}{Matrices from \code{block_rNormalGLM}.}
#'     \item{block_rNormalGLM}{Full block sampler output.}
#'   }
#' @seealso \code{\link{block_rNormalGLM}}, \code{\link{normalize_block}}
#' @example inst/examples/Ex_block_rNormalGLM_update.R
#' @rdname block_simfuncs
#' @export
block_rNormalGLM_update <- function(mu_all,
                                    sigma_theta_sq = NULL,
                                    y,
                                    x,
                                    block,
                                    family = poisson(),
                                    prior_list = NULL,
                                    prior_lists = NULL,
                                    offset = NULL,
                                    weights = 1,
                                    Gridtype = 2L,
                                    n_envopt = 1L,
                                    use_parallel = TRUE,
                                    use_opencl = FALSE,
                                    verbose = FALSE,
                                    progbar = FALSE,
                                    seed = NULL,
                                    theta_coef_col = 1L) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(prior_lists) && is.null(prior_list)) {
    if (missing(mu_all)) {
      stop("Provide 'mu_all' (with 'sigma_theta_sq'), 'prior_lists', or 'prior_list'.",
           call. = FALSE)
    }
    if (is.null(sigma_theta_sq)) {
      stop("'sigma_theta_sq' is required when building per-block priors from 'mu_all'.",
           call. = FALSE)
    }
    prior_lists <- lapply(mu_all, function(m) {
      list(
        mu = m,
        Sigma = matrix(sigma_theta_sq, 1, 1),
        dispersion = 1,
        ddef = FALSE
      )
    })
  }

  out <- block_rNormalGLM(
    n = 1L,
    y = y,
    x = x,
    block = block,
    prior_list = prior_list,
    prior_lists = prior_lists,
    offset = offset,
    weights = weights,
    family = family,
    Gridtype = as.integer(Gridtype),
    n_envopt = n_envopt,
    use_parallel = use_parallel,
    use_opencl = use_opencl,
    verbose = verbose,
    progbar = progbar
  )

  theta <- as.vector(out$coefficients[, theta_coef_col, drop = TRUE])

  list(
    theta = theta,
    coefficients = out$coefficients,
    coef.mode = out$coef.mode,
    block_rNormalGLM = out
  )
}
