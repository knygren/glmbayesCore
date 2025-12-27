#' @title Verify OpenCL Runtime Availability
#'
#' @description
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
verify_opencl_runtime <- function(lib_dirs = NULL) {
  code <- '
  #define CL_TARGET_OPENCL_VERSION 300
  #include <CL/cl.h>
  int main() {
    cl_uint n = 0;
    cl_int status = clGetPlatformIDs(0, NULL, &n);
    if (status != CL_SUCCESS || n == 0) return 1;
    return 0;
  }'
  
  tf_c  <- tempfile(fileext = ".c")
  tf_exe <- tempfile("ocltest")
  writeLines(code, tf_c)
  
  # -----------------------------------------
  # 1. Locate actual libOpenCL.so* file
  # -----------------------------------------
  lib_path <- NULL
  
  if (!is.null(lib_dirs) && length(lib_dirs) > 0) {
    for (d in lib_dirs) {
      hits <- Sys.glob(file.path(d, "libOpenCL.so*"))
      if (length(hits) > 0) {
        lib_path <- hits[1]   # use first match
        break
      }
    }
  }
  
  # If still not found, try ldconfig
  if (is.null(lib_path)) {
    ld_hits <- try(
      system("ldconfig -p | grep -i opencl | awk '{print $NF}'", intern = TRUE),
      silent = TRUE
    )
    if (!inherits(ld_hits, "try-error") && length(ld_hits) > 0) {
      lib_path <- ld_hits[1]
    }
  }
  
  # If still not found, runtime probe cannot proceed
  if (is.null(lib_path)) {
    unlink(c(tf_c, tf_exe))
    return(FALSE)
  }
  
  lib_dir <- dirname(lib_path)
  
  # -----------------------------------------
  # 2. Attempt to compile using detected lib
  # -----------------------------------------
  compile_status <- suppressWarnings(
    system2("gcc", c(tf_c, paste0("-L", lib_dir), "-lOpenCL", "-o", tf_exe))
  )
  
  if (!is.numeric(compile_status) || compile_status != 0) {
    unlink(c(tf_c, tf_exe))
    return(FALSE)
  }
  
  # -----------------------------------------
  # 3. Attempt to run the test executable
  # -----------------------------------------
  run_status <- suppressWarnings(system2(tf_exe))
  unlink(c(tf_c, tf_exe))
  
  return(is.numeric(run_status) && run_status == 0)
}