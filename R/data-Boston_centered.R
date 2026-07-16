#' Boston housing data with mean-centered predictors
#'
#' A copy of \code{\link[MASS]{Boston}} where all predictors (every column except
#' \code{medv}) have been mean-centered (subtract column means, no scaling).
#'
#' @format A data frame with 506 observations and 14 variables (same names as
#'   \code{\link[MASS]{Boston}}). See \code{?MASS::Boston} for variable
#'   descriptions.
#'
#' @source Derived from \code{MASS::Boston}. Original data described in Harrison
#'   and Rubinfeld (1978); see \code{?Boston} in \pkg{MASS}.
#'
#' @usage data("Boston_centered")
#'
#' @example inst/examples/Ex_Boston_centered_opencl.R
#'
#' @keywords datasets
#' @concept Bayesian linear regression
"Boston_centered"
