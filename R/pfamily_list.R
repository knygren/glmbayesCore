#' Build a named list of pfamily objects
#'
#' Generic for constructing a named list of \code{\link{pfamily}} prior
#' objects from a prior-specification object.  Packages that define prior
#' setup containers (e.g. \code{Prior_Setup_lmebayes()} in \pkg{lmebayes})
#' provide methods that map each component of the container to a
#' \code{pfamily} constructor such as \code{\link{dNormal}} or
#' \code{\link{dIndependent_Normal_Gamma}}.
#'
#' @param object A prior-specification object.
#' @param ... Additional arguments passed to methods (e.g. \code{ptypes}).
#' @return A named list whose elements are objects of class
#'   \code{"pfamily"}.
#' @seealso \code{\link{pfamily}}, \code{\link{dNormal}},
#'   \code{\link{dIndependent_Normal_Gamma}}
#' @export
pfamily_list <- function(object, ...) UseMethod("pfamily_list")
