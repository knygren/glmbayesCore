#' @rdname pfamily_list
#' @order 2
#' @details
#' For a \code{\link{Prior_Setup_GLMM}} object, each group-effect
#' coefficient \eqn{k} uses parameters from
#' \code{object$pop.prior_list[[k]]}:
#' \itemize{
#'   \item \code{"dNormal"}: \code{dNormal(mu, Sigma, dispersion)} with
#'     group-effect variance \eqn{\tau^2_k} from
#'     \code{pop.prior_list[[k]]$dispersion} (A12 §3.3.4 marginal
#'     hyper-regression center for Gaussian models; classical VarCorr in
#'     \code{pop.dispersion.ref}).
#'   \item \code{"dIndependent_Normal_Gamma"}: the same \code{mu} and
#'     \code{Sigma}, plus a Gamma prior on precision \eqn{1/\tau^2_k}
#'     from \code{object$pop.ing_prior[[k]]} (calibrated by
#'     \code{\link{Prior_Setup_GLMM}} via \code{pop.dispersion.pwt} /
#'     \code{pop.dispersion.nprior}, mean-matched at that same marginal
#'     \eqn{\tau^2_k}).  With effective prior sample size
#'     \eqn{n_0} and \eqn{p_k} population predictors for component
#'     \eqn{k},
#'     \deqn{shape = (n_0 + 1 + p_k)/2, \qquad
#'           rate = \tau^2_k \, (n_0 + p_k - 1)/2.}
#'     Requires \eqn{n_0 \le J} (equivalently
#'     \code{pop.dispersion.pwt} \eqn{\le 0.5}). Truncation bounds
#'     \code{disp_lower} / \code{disp_upper} come from the limiting
#'     posterior window stored on \code{pop.ing_prior}; see
#'     \code{inst/ING_TRUNCATION_WINDOW.md}.
#' }
#'
#' This is separate from observation-level \eqn{\sigma^2} priors: those
#' use \code{\link{dGamma_list}} (per-group \code{dispformula}) or a
#' single \code{\link[glmbayesCore]{dGamma}} built from
#' \code{object$group.ing_prior} when dispersion is pooled.
#'
#' @export
#' @method pfamily_list Prior_Setup_GLMM
pfamily_list.Prior_Setup_GLMM <- function(object,
                                              ptypes = NULL,
                                              ...) {

  allowed <- c("dNormal", "dIndependent_Normal_Gamma")

  re_names <- names(object$pop.prior_list)
  p_re     <- length(re_names)

  if (is.null(ptypes)) {
    ptypes <- "dNormal"
  }

  if (is.list(ptypes)) {
    ok <- vapply(
      ptypes,
      function(p) is.character(p) && length(p) == 1L && !is.na(p),
      logical(1L)
    )
    if (!all(ok)) {
      stop("'ptypes' list elements must each be a single string.",
           call. = FALSE)
    }
    nms    <- names(ptypes)
    ptypes <- vapply(ptypes, identity, character(1L))
    names(ptypes) <- nms
  }
  if (!is.character(ptypes) || length(ptypes) < 1L || anyNA(ptypes)) {
    stop("'ptypes' must be a character vector or list of strings.",
         call. = FALSE)
  }
  bad <- setdiff(unique(ptypes), allowed)
  if (length(bad) > 0L) {
    stop(
      "Invalid 'ptypes' value(s): ", paste(bad, collapse = ", "),
      ". Allowed: ", paste(allowed, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (length(ptypes) == 1L) {
    ptypes <- stats::setNames(rep(unname(ptypes), p_re), re_names)
  } else {
    if (length(ptypes) != p_re) {
      stop(
        sprintf(
          "'ptypes' has length %d but the prior setup has %d random-effect component(s): %s.",
          length(ptypes), p_re, paste(re_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    if (!is.null(names(ptypes)) && any(nzchar(names(ptypes)))) {
      if (!setequal(names(ptypes), re_names)) {
        stop(
          "Names of 'ptypes' must match the random-effect coefficient names: ",
          paste(re_names, collapse = ", "), ".",
          call. = FALSE
        )
      }
      ptypes <- ptypes[re_names]
    } else {
      names(ptypes) <- re_names
    }
  }

  J   <- nlevels(object$design$group)
  npd <- object$pop.dispersion.nprior

  n_prior_for <- function(k) {
    if (!is.null(npd)) {
      return(unname(npd[[k]]))
    }
    w <- if (is.list(object$pop.pwt)) mean(object$pop.pwt[[k]]) else object$pop.pwt
    (w / (1 - w)) * J
  }

  out <- stats::setNames(vector("list", p_re), re_names)

  for (k in re_names) {
    pl    <- object$pop.prior_list[[k]]
    mu_k  <- pl$mu
    Sig_k <- pl$Sigma
    d_k   <- unname(pl$dispersion)
    p_k   <- length(mu_k)

    out[[k]] <- switch(
      ptypes[[k]],
      dNormal = dNormal(
        mu         = mu_k,
        Sigma      = Sig_k,
        dispersion = d_k
      ),
      dIndependent_Normal_Gamma = {
        n_prior_k <- n_prior_for(k)
        if (n_prior_k > J) {
          stop(
            "Component \"", k, "\": the dispersion prior has effective ",
            "prior sample size pop.dispersion.nprior = ", signif(n_prior_k, 4),
            ", but there are only J = ", J, " groups. ",
            "dIndependent_Normal_Gamma sampling requires ",
            "pop.dispersion.nprior <= J (pop.dispersion.pwt <= 0.5); lower ",
            "'pop.dispersion.pwt'/'pop.dispersion.nprior' in Prior_Setup_GLMM().",
            call. = FALSE
          )
        }
        ing_k <- object$pop.ing_prior[[k]]
        if (is.null(ing_k)) {
          shape_k <- (n_prior_k + 1) / 2 + p_k / 2
          rate_k  <- d_k * (n_prior_k + p_k - 1) / 2
          win_k <- .lmebayes_ing_limiting_posterior_window(d_k, J)
          ing_k <- list(
            shape      = shape_k,
            rate       = rate_k,
            disp_lower = win_k$disp_lower,
            disp_upper = win_k$disp_upper
          )
        }
        dIndependent_Normal_Gamma(
          mu         = mu_k,
          Sigma      = Sig_k,
          shape      = ing_k$shape,
          rate       = ing_k$rate,
          disp_lower = ing_k$disp_lower,
          disp_upper = ing_k$disp_upper
        )
      }
    )
  }

  attr(out, "ptypes") <- ptypes
  class(out) <- c("pfamily_list", "list")
  out
}
