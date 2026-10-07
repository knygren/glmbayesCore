#' Number of stored draws on a simfunc return
#' @noRd
.lmebayes_reg_n_draws <- function(x) {
  if (!is.null(x$n) && length(x$n) == 1L && is.finite(x$n) && x$n >= 0) {
    return(as.integer(x$n))
  }
  if (!is.null(x$popef) && length(x$popef) > 0L && !is.null(x$popef[[1L]])) {
    return(nrow(as.matrix(x$popef[[1L]])))
  }
  if (!is.null(x$groupef) && is.data.frame(x$groupef) && "draw" %in% names(x$groupef)) {
    return(length(unique(x$groupef$draw)))
  }
  0L
}

#' Resolve draw indices for print-only subsetting (object unchanged)
#' @noRd
.lmebayes_resolve_print_draws <- function(draws, n, arg = "draws") {
  if (is.null(draws)) {
    return(NULL)
  }
  if (n < 1L) {
    stop("Object has no stored draws to print.", call. = FALSE)
  }
  if (is.logical(draws)) {
    if (length(draws) != n) {
      stop(
        "Logical '", arg, "' must have length ", n,
        " (number of draws).",
        call. = FALSE
      )
    }
    draws <- which(draws)
  }
  if (!is.numeric(draws) || !length(draws) ||
      anyNA(draws) || any(draws != as.integer(draws))) {
    stop(
      "'", arg, "' must be integer draw indices (or a logical of length n), ",
      "e.g. 1:10.",
      call. = FALSE
    )
  }
  draws <- as.integer(draws)
  if (any(draws < 1L) || any(draws > n)) {
    stop("'", arg, "' indices must be in 1:", n, ".", call. = FALSE)
  }
  draws
}

#' Print simulated group coefficients from an LMM result
#'
#' S3 generic for distribution-sampler style display of \code{groupef} draws
#' (long data frame). Filters by grouping levels and optional RE columns.
#' Methods are registered for \code{rLMMNormal_reg},
#' \code{rLMMindepNormalGamma_reg}, and \code{rlmerb}; the \code{default}
#' method accepts any
#' list-like object with a usable \code{groupef} data frame.
#'
#' @param x An object with a \code{groupef} data frame (for example an
#'   \code{\link{rLMM_reg}} return).
#' @param groups \code{NULL} (default; all group levels), a character vector of
#'   levels in the grouping column, or integer indices into the unique levels
#'   present in \code{groupef}.
#' @param components \code{NULL} (all RE columns) or a character vector / integer
#'   indices selecting RE coefficient columns (excluding the draw and grouping
#'   columns).
#' @param draws \code{NULL} (all draws) or integer indices / logical of length
#'   \code{n} selecting which \code{draw} values to display (print-only).
#' @param digits Number of significant digits for printed numeric values.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' ## The groupef layout
#'
#' Population coefficients (\eqn{\gamma}) and group coefficients
#' (\eqn{\beta_j}) are stored in different shapes. \code{popef} is a
#' \emph{named list of matrices}, one per random-effect component, each
#' \code{n} draws by \eqn{q_k} level-2 predictors. \code{groupef} is a
#' single \strong{long-format data frame} with one row per
#' (draw, group level) pair:
#'
#' \tabular{ll}{
#'   \code{draw} \tab integer draw index, \code{1:n} \cr
#'   \emph{grouping column} \tab the group level, named after the grouping
#'     factor (\code{"Subject"}, \code{"state"}, ...) \cr
#'   \emph{one column per RE coefficient} \tab \code{"(Intercept)"} and any
#'     random slopes
#' }
#'
#' So a fit with \code{n = 1000} draws over 20 groups produces 20,000 rows.
#' Printing it whole is rarely what anyone wants, which is the reason this
#' function exists.
#'
#' ## Filtering
#'
#' The three filters are independent and compose; all of them are
#' \strong{print-only}, and \code{x} is never modified.
#'
#' \itemize{
#'   \item \code{draws} --- integer indices into \code{1:n}, or a logical mask
#'     of length \code{n}.  Filters rows by the \code{draw} column.
#'   \item \code{groups} --- level names as characters, or integer indices
#'     into the unique levels \emph{present in} \code{groupef} (which is the
#'     ordering the fit used, not necessarily \code{levels()} of the original
#'     factor).  Filters rows by the grouping column.
#'   \item \code{components} --- names or integer indices selecting which
#'     random-effect \emph{columns} to show.  The \code{draw} and grouping
#'     columns are always kept and are never counted in the indexing.
#' }
#'
#' \code{print_groupef(fit, draws = 1:5, groups = 1:3)} is the usual idiom:
#' the first few draws for the first few groups, enough to confirm the
#' sampler is producing sensible numbers without flooding the console.  An
#' empty result after filtering prints a placeholder rather than erroring.
#'
#' Numeric columns are formatted to \code{digits} significant figures for
#' display only; the stored values keep full precision.
#'
#' ## Dispatch
#'
#' The registered methods share the \code{default} implementation. They
#' exist so that the generic dispatches cleanly on the sampler return classes,
#' and so that class-specific layouts can be added later without changing
#' call sites. The \code{default} method works on any list-like object
#' carrying a usable \code{groupef} data frame.
#'
#' Point-estimate fits (\code{simulate = FALSE}) have no \code{groupef} at
#' all; use \code{groupef.mode} directly for those.
#'
#' @return No return value, called for side effects.
#' @seealso \code{\link{rLMM_reg}}
#' @name print_groupef
#' @export
print_groupef <- function(
    x,
    groups = NULL,
    components = NULL,
    draws = NULL,
    digits = max(3, getOption("digits") - 3),
    ...
) {
  UseMethod("print_groupef")
}

#' @rdname print_groupef
#' @export
print_groupef.default <- function(
    x,
    groups = NULL,
    components = NULL,
    draws = NULL,
    digits = max(3, getOption("digits") - 3),
    ...
) {
  gf <- x$groupef
  if (is.null(gf)) {
    stop("'x' has no 'groupef' draws to print.", call. = FALSE)
  }
  if (!is.data.frame(gf)) {
    stop("'x$groupef' must be a data frame.", call. = FALSE)
  }

  grp_col <- NULL
  if (!is.null(x$design$group_name) && nzchar(x$design$group_name) &&
      x$design$group_name %in% names(gf)) {
    grp_col <- x$design$group_name
  } else {
    skip <- c("draw", "draws")
    cand <- setdiff(names(gf), skip)
    if (!length(cand)) {
      stop(
        "Could not find a grouping column in 'x$groupef'.",
        call. = FALSE
      )
    }
    ## Prefer a non-numeric column; else first non-draw column.
    is_num <- vapply(gf[cand], is.numeric, logical(1L))
    grp_col <- if (any(!is_num)) cand[!is_num][1L] else cand[1L]
  }

  re_cols <- setdiff(names(gf), c("draw", grp_col))
  if (!length(re_cols)) {
    stop("'x$groupef' has no RE coefficient columns to print.", call. = FALSE)
  }

  if (!is.null(components)) {
    re_list <- stats::setNames(as.list(re_cols), re_cols)
    sel_re <- .lmebayes_select_named_list_keys(
      re_list, components, arg = "components", what = "component"
    )
  } else {
    sel_re <- re_cols
  }

  levels_all <- unique(as.character(gf[[grp_col]]))
  if (is.null(groups)) {
    sel_grp <- levels_all
  } else if (is.numeric(groups)) {
    if (any(is.na(groups)) || any(groups < 1) ||
        any(groups > length(levels_all)) ||
        any(groups != as.integer(groups))) {
      stop(
        "'groups' numeric indices must be integers in 1:",
        length(levels_all), ".",
        call. = FALSE
      )
    }
    sel_grp <- levels_all[as.integer(groups)]
  } else if (is.character(groups) || is.factor(groups)) {
    groups <- as.character(groups)
    missing <- setdiff(groups, levels_all)
    if (length(missing)) {
      stop(
        "'groups' contains unknown level(s): ",
        paste(missing, collapse = ", "),
        ".\nAvailable: ", paste(levels_all, collapse = ", "), ".",
        call. = FALSE
      )
    }
    sel_grp <- groups
  } else {
    stop(
      "'groups' must be NULL, a character vector of group levels, ",
      "or integer indices.",
      call. = FALSE
    )
  }

  keep <- as.character(gf[[grp_col]]) %in% sel_grp
  i <- .lmebayes_resolve_print_draws(draws, .lmebayes_reg_n_draws(x))
  if (!is.null(i) && "draw" %in% names(gf)) {
    keep <- keep & (as.integer(gf$draw) %in% i)
  }
  out_df <- gf[keep, c(intersect(c("draw", grp_col), names(gf)), sel_re),
               drop = FALSE]
  rownames(out_df) <- NULL

  cat("\nSimulated group coefficients:\n")
  if (!nrow(out_df)) {
    cat("  [no rows after filtering]\n\n")
    return(invisible(NULL))
  }

  ## Format numeric RE columns; keep draw / group labels as-is.
  show <- out_df
  for (nm in sel_re) {
    if (is.numeric(show[[nm]])) {
      show[[nm]] <- format(show[[nm]], digits = digits)
    }
  }
  print(show, row.names = FALSE, quote = FALSE)
  cat("\n")
  invisible(NULL)
}

#' @rdname print_groupef
#' @export
print_groupef.rLMMNormal_reg <- print_groupef.default

#' @rdname print_groupef
#' @export
print_groupef.rLMMindepNormalGamma_reg <- print_groupef.default

#' @rdname print_groupef
#' @export
print_groupef.rlmerb <- print_groupef.default
