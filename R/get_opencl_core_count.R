#' Get the number of available OpenCL compute units
#'
#' Returns the number of compute units (cores) available on the default
#' OpenCL device. This can be useful for diagnostics or performance tuning.
#'
#' @return Integer count of compute units.
#' @export
get_opencl_core_count <- function() {
  .get_opencl_core_count_cpp()
}