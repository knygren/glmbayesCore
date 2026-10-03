#' Build a named list of dGamma measurement-dispersion priors
#'
#' Generic for constructing a named list of \code{\link{dGamma}} prior
#' objects from a prior-specification object.  For mixed models, one
#' \code{\link{dGamma}} per group level supplies observation-level
#' \eqn{\sigma^2_j} priors for \code{lmerb(..., group.dispersion = ...)}.
#'
#' @param object A prior-specification object.
#' @param ... Additional arguments passed to methods (e.g.
#'   \code{max_disp_perc_measurement}). For \code{print.dGamma_list()},
#'   passed to \code{print.pfamily}.
#' @return For \code{dGamma_list()}, an object of class
#'   \code{"dGamma_list"} (also inherits from \code{"list"}): a named
#'   list of \code{"pfamily"} objects (one per group).  Attributes such
#'   as \code{window_diagnostics} and \code{measurement_prior_group}
#'   are retained for downstream use but are not shown by
#'   \code{print.dGamma_list()}.  For \code{print.dGamma_list()},
#'   \code{x} invisibly.
#' @details
#' Where \code{\link{pfamily_list}} configures the \emph{between}-group
#' variance components \eqn{\tau^2_k}, \code{dGamma_list()} configures the
#' \emph{within}-group observation variance. It returns one
#' \code{\link{dGamma}} prior per group level, so each group
#' gets its own \eqn{\sigma^2_j} rather than sharing a pooled value.
#'
#' Each element places a Gamma prior on the \emph{precision}
#' \eqn{1/\sigma^2_j} (hence \code{Inv_Dispersion = TRUE}) and carries a
#' truncation window \code{disp_lower}/\code{disp_upper} calibrated from
#' that group's own data. Because groups differ in size and in residual
#' scatter, the windows differ across groups -- often by an order of
#' magnitude -- which is the point of the per-group treatment.
#'
#' @section Prerequisite: a per-group calibration:
#' \code{dGamma_list()} is not a transformation of an arbitrary prior
#' setup. It requires that \code{\link{Prior_Setup_GLMM}} was called with
#' a grouping \code{dispformula} (e.g. \code{dispformula = ~school_id}) so
#' that a per-group measurement-dispersion prior was actually computed.
#' Called on a \code{dispformula = ~1} setup it stops with a message
#' naming the fix rather than silently falling back to a pooled prior.
#'
#' @section What it costs:
#' Supplying the result as \code{dispersion_ranef} turns \eqn{J} residual
#' variances into parameters (see \dQuote{What the prior choices mean} in
#' \samp{rLMM_reg()}). The posterior is then no longer Gaussian, so
#' draws are no longer exact: they come from the final sweep of a
#' calibrated chain, and a pilot stage runs first. Because a value must be
#' drawn for every group at every sweep, per-group dispersion is
#' substantially more expensive than a fixed or pooled \eqn{\sigma^2}, and
#' it is worth the cost only when groups genuinely differ in residual
#' scatter.
#'
#' Posterior summaries for the per-group variances are produced by
#' \code{summary_sigma2()} in \pkg{lmebayes}, which reports the prior, the
#' truncation window, the number of observations per group, posterior
#' moments, and accept-reject candidates per draw.
#'
#' @seealso \code{\link{dGamma}}, \code{\link{Prior_Setup_GLMM}},
#'   \code{\link{pfamily_list}}, \samp{rLMM_reg()}
#' @example inst/examples/Ex_dGamma_list.R
#' @order 1
#' @export
dGamma_list <- function(object, ...) UseMethod("dGamma_list")

#' Resolve names/indices for a named list (dGamma_list / pfamily_list)
#' @noRd
.lmebayes_select_named_list_keys <- function(x, keys, arg = "groups",
                                             what = "name") {
  all_names <- names(x)
  if (is.null(all_names)) {
    all_names <- as.character(seq_along(x))
  }
  if (is.null(keys)) {
    return(all_names)
  }
  if (is.numeric(keys)) {
    if (any(is.na(keys)) || any(keys < 1) || any(keys > length(x)) ||
        any(keys != as.integer(keys))) {
      stop(
        "'", arg, "' numeric indices must be integers in 1:", length(x), ".",
        call. = FALSE
      )
    }
    return(all_names[as.integer(keys)])
  }
  if (!is.character(keys) || length(keys) < 1L || anyNA(keys)) {
    stop(
      "'", arg, "' must be NULL, a character vector of ", what, "s, ",
      "or integer indices.",
      call. = FALSE
    )
  }
  keys <- as.character(keys)
  missing <- setdiff(keys, all_names)
  if (length(missing) > 0L) {
    stop(
      "'", arg, "' contains unknown ", what, "(s): ",
      paste(missing, collapse = ", "),
      ".\nAvailable: ", paste(all_names, collapse = ", "), ".",
      call. = FALSE
    )
  }
  keys
}

#' @rdname dGamma_list
#' @method print dGamma_list
#' @param x An object of class \code{"dGamma_list"}.
#' @param groups For \code{print.dGamma_list()}: \code{NULL} (default;
#'   print all groups), a character vector of group names
#'   (\code{names(x)}), or integer indices into \code{x}. Each selected
#'   group's \code{dGamma()} is printed with
#'   \code{\link[glmbayesCore]{print.pfamily}}. List attributes
#'   (\code{window_diagnostics}, \code{group.dispersion.fit}, etc.) are
#'   omitted; use \code{attr(x, "...")} to inspect them.
#' @order 3
#' @export
print.dGamma_list <- function(x, groups = NULL, ...) {
  sel <- .lmebayes_select_named_list_keys(
    x, groups, arg = "groups", what = "group name"
  )
  for (nm in sel) {
    cat(sprintf("\n[[%s]]\n", nm))
    print(x[[nm]], ...)
  }
  invisible(x)
}
