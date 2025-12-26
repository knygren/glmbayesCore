#' Verify OpenCL Runtime Availability
#'
#' This function attempts to compile and run a minimal C program that calls
#' \code{clGetPlatformIDs} to check whether a usable OpenCL runtime platform
#' is available on the system. It is designed to complement
#' \code{detect_compute_runtimes()}, which only checks for headers and libraries.
#'
#' The probe works by:
#' \enumerate{
#'   \item Writing a small C source file that includes \code{CL/cl.h}.
#'   \item Compiling it with \code{gcc} and linking against \code{-lOpenCL}.
#'   \item Executing the compiled binary to verify that at least one OpenCL
#'         platform is returned by \code{clGetPlatformIDs}.
#' }
#'
#' If compilation or execution fails, the function returns \code{FALSE}.
#' If compilation and execution succeed and at least one platform is detected,
#' the function returns \code{TRUE}.
#'
#' @return Logical. \code{TRUE} if a usable OpenCL runtime is detected,
#'         \code{FALSE} otherwise.
#'
#' @examples
#' \dontrun{
#' if (verify_opencl_runtime()) {
#'   message("OpenCL runtime is available.")
#' } else {
#'   message("No usable OpenCL runtime detected.")
#' }
#' }
#'
#' @export
verify_opencl_runtime <- function() {
  code <- '
  #define CL_TARGET_OPENCL_VERSION 300
  #include <CL/cl.h>
  int main() {
    cl_uint n = 0;
    cl_int status = clGetPlatformIDs(0, NULL, &n);
    if (status != CL_SUCCESS || n == 0) return 1;
    return 0;
  }'
  
  tf_c <- tempfile(fileext = ".c")
  tf_exe <- tempfile("ocltest")
  writeLines(code, tf_c)
  
  # Attempt to compile
  compile_status <- suppressWarnings(
    system2("gcc", c(tf_c, "-lOpenCL", "-o", tf_exe))
  )
  if (!is.numeric(compile_status) || compile_status != 0) {
    unlink(c(tf_c, tf_exe))
    return(FALSE)
  }
  
  # Attempt to run
  run_status <- suppressWarnings(system2(tf_exe))
  unlink(c(tf_c, tf_exe))
  
  return(is.numeric(run_status) && run_status == 0)
}