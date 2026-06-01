#' Summarize and print \code{mlmb} fits
#'
#' @description
#' Methods for multi-response \code{\link{lmb}} fits (class \code{"mlmb"}).
#' \code{summary.mlmb} applies \code{\link{summary.glmb}} to each response;
#' printing follows \code{\link[stats]{summary.mlm}} with per-response sections.
#'
#' @param object An object of class \code{"mlmb"} from \code{\link{multi_lmb}}.
#' @param x An object of class \code{"summary.mlmb"}.
#' @param digits Number of significant digits for printing.
#' @param \ldots Passed to \code{\link{summary.glmb}} or print methods.
#' @return \code{summary.mlmb} returns a named list of \code{"summary.glmb"}
#'   objects with class \code{"summary.mlmb"}.
#' @seealso \code{\link{multi_lmb}}, \code{\link{print.mlmb}},
#'   \code{\link{summary.glmb}}, \code{\link[stats]{summary.mlm}}
#' @name summary.mlmb
#' @aliases summary.mlmb print.mlmb print.summary.mlmb
NULL

#' @rdname summary.mlmb
#' @export
#' @method summary mlmb
summary.mlmb <- function(object, ...) {
  res <- lapply(object, function(fit) {
    s <- summary(fit, ...)
    s$call <- fit$call
    s
  })
  names(res) <- names(object)
  attr(res, "mlmb_call") <- attr(object, "call")
  attr(res, "coef_means") <- .mlmb_coef_means_matrix(object)
  attr(res, "dic_table") <- .mlmb_dic_table(object)
  class(res) <- "summary.mlmb"
  res
}

#' @rdname summary.mlmb
#' @export
#' @method print mlmb
print.mlmb <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cl <- attr(x, "call")
  if (is.null(cl) && length(x) >= 1L) {
    cl <- x[[1L]]$call
  }
  cat("\nCall:\n")
  if (!is.null(cl)) {
    print(cl)
  } else {
    cat("  (not recorded)\n")
  }

  cm <- .mlmb_coef_means_matrix(x)
  if (!is.null(cm) && length(cm)) {
    cat("\nPosterior mean coefficients:\n")
    print.default(
      format(cm, digits = digits),
      print.gap = 2L,
      quote = FALSE
    )
  } else {
    cat("\nNo coefficients\n")
  }

  dic_tab <- .mlmb_dic_table(x)
  if (!is.null(dic_tab) && nrow(dic_tab) >= 1L) {
    cat("\nBayesian fit (per response; independent product model):\n")
    print.default(
      format(dic_tab, digits = digits),
      print.gap = 2L,
      quote = FALSE
    )
    cat(
      "Sum DIC:",
      format(sum(dic_tab[, "DIC"]), digits = digits),
      "  Sum pD:",
      format(sum(dic_tab[, "pD"]), digits = digits),
      "\n",
      sep = ""
    )
  }

  cat("\n")
  invisible(x)
}

#' @rdname summary.mlmb
#' @export
#' @method print summary.mlmb
print.summary.mlmb <- function(
    x,
    digits = max(3, getOption("digits") - 3),
    ...
) {
  cl <- attr(x, "mlmb_call")
  if (!is.null(cl)) {
    cat("\nCall:\n")
    print(cl)
  }

  cm <- attr(x, "coef_means")
  if (!is.null(cm) && length(cm)) {
    cat("\nPosterior mean coefficients:\n")
    print.default(
      format(cm, digits = digits),
      print.gap = 2L,
      quote = FALSE
    )
  }

  dic_tab <- attr(x, "dic_table")
  if (!is.null(dic_tab) && nrow(dic_tab) >= 1L) {
    cat("\nBayesian fit (per response):\n")
    print.default(
      format(dic_tab, digits = digits),
      print.gap = 2L,
      quote = FALSE
    )
    cat(
      "Sum DIC:",
      format(sum(dic_tab[, "DIC"]), digits = digits),
      "\n",
      sep = ""
    )
  }

  for (nm in names(x)) {
    cat("\nResponse", nm, ":\n")
    print(x[[nm]], digits = digits, ...)
  }
  invisible(x)
}
