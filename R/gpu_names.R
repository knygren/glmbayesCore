#' GPU Detection via nvidia-smi
#'
#' Calls the C++ implementation of `gpu_names()` to return GPU names
#' using `nvidia-smi --query-gpu=name --format=csv,noheader`.
#'
#' @title GPU Detection
#' @description Returns a character vector of GPU names detected by nvidia-smi.
#' @rdname GpuNames
#' @export
#' @usage gpu_names()
gpu_names <- function() {
  .Call(`_glmbayes_gpu_names`)
}