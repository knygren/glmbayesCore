#' Build a named list of pfamily objects from a prior-setup object
#'
#' @description
#' Generic for converting a prior-setup object into a named list of
#' \code{\link{pfamily}} objects (class \code{"pfamily_list"}) for block or
#' mixed-model samplers.
#'
#' @param object A prior-setup object. Dispatch is method-specific; see
#'   \strong{Methods} and linked help pages for required classes and arguments.
#' @param ptypes Prior-family name(s) per component; meaning depends on the
#'   active method (for row-block setup see \code{\link{Prior_SetupGroup}}).
#' @param ... Passed to methods (unused for \code{Prior_SetupGroup}).
#'
#' @return For \code{pfamily_list()}, an object of class \code{"pfamily_list"}
#'   (also inherits from \code{"list"}): a named list of \code{pfamily} objects.
#'   Attribute \code{"ptypes"} records the resolved prior-family name per
#'   component when the method sets it. For \code{print.pfamily_list()},
#'   \code{x} invisibly.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{Prior_SetupGroup}}{
#'     Row-block \code{\link{Prior_SetupGroup}} results: one \code{pfamily} per
#'     block. \code{ptypes}, allowed families, and examples are documented on
#'     \code{\link{Prior_SetupGroup}}.}
#' }
#'
#' @seealso \code{\link{Prior_SetupGroup}}, \code{\link{pfamily}},
#'   \code{\link{dNormal}}, \code{\link{dNormal_Gamma}},
#'   \code{\link{dIndependent_Normal_Gamma}}
#' @export
pfamily_list <- function(object, ptypes = NULL, ...) {
  UseMethod("pfamily_list")
}

#' @rdname pfamily_list
#' @method print pfamily_list
#' @param x An object of class \code{"pfamily_list"}.
#' @param components \code{NULL} (all blocks), block names, or integer indices.
#' @export
print.pfamily_list <- function(x, components = NULL, ...) {
  sel <- .pfamily_list_select_keys(
    x, components, arg = "components", what = "component"
  )
  for (nm in sel) {
    cat(sprintf("\n[[%s]]\n", nm))
    print(x[[nm]], ...)
  }
  invisible(x)
}

#' @noRd
.pfamily_list_select_keys <- function(x, keys, arg = "groups", what = "name") {
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
