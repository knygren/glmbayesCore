#' Summarizing mrglmb Objects
#'
#' These functions are \code{\link{methods}} for class \code{"mrglmb"} or
#' \code{"summary.mrglmb"} objects produced by
#' \code{\link{multi_rindepNormalGamma_reg}}, \code{\link{multi_rNormalGamma_reg}},
#'   \code{\link{multi_rNormal_reg}}, or \code{\link{multi_rlmb}}.
#'
#' @aliases
#' summary.mrglmb
#' print.summary.mrglmb
#' @param object An object of class \code{"mrglmb"}.
#' @param x An object of class \code{"summary.mrglmb"} for which printed
#'   output is desired.
#' @param digits The number of significant digits to use when printing.
#' @param \ldots Additional optional arguments passed to
#'   \code{\link{summary.rglmb}} or \code{\link{print.summary.rglmb}}.
#' @return \code{summary.mrglmb} returns a named list of
#'   \code{"summary.rglmb"} objects (one per response column), with class
#'   \code{"summary.mrglmb"}.  The names match \code{names(object)}.
#' @details
#' Mirrors the behavior of \code{\link[stats]{summary.mlm}}: each response
#' column is summarized independently using \code{\link{summary.rglmb}} and
#' printed with a \code{"Response <name> :"} header.
#' @seealso \code{\link{multi_rlmb}}, \code{\link{multi_rNormalGamma_reg}},
#'   \code{\link{multi_rNormal_reg}}, \code{\link{multi_rindepNormalGamma_reg}},
#'   \code{\link{summary.rglmb}}
#' @export
#' @method summary mrglmb
summary.mrglmb <- function(object, ...) {
  res <- lapply(object, summary, ...)
  names(res) <- names(object)
  class(res) <- "summary.mrglmb"
  res
}

#' @rdname summary.mrglmb
#' @export
#' @method print summary.mrglmb
print.summary.mrglmb <- function(x,
                                  digits = max(3, getOption("digits") - 3),
                                  ...) {
  for (nm in names(x)) {
    cat("\nResponse", nm, ":\n")
    print(x[[nm]], digits = digits, ...)
  }
  invisible(x)
}
