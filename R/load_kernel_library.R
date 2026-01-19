#' Load OpenCL Kernel Source Files
#'
#' These functions provide a user-facing interface for loading OpenCL kernel
#' source files and kernel libraries from the package's `cl/` directory.
#' They call internal C++ routines that perform file lookup, dependency
#' resolution, and concatenation of kernel sources.
#'
#' OpenCL support is optional. If the package was built without OpenCL
#' (e.g., on systems lacking OpenCL headers or drivers), these functions
#' return a clear error message.
#'
#' @section OpenCL Availability:
#' Use \code{\link{has_opencl}} to check whether OpenCL support is available
#' in the current build of \pkg{glmbayes}.
#'
#'
#' @section How These Functions Assemble an OpenCL Program:
#'
#' The functions \code{load_kernel_source()} and \code{load_kernel_library()}
#' are the fundamental tools used by \pkg{glmbayes} to construct complete
#' OpenCL programs from modular components. OpenCL kernels in this package are
#' not stored as monolithic `.cl` files. Instead, they are built dynamically
#' by concatenating several layers of source code, each serving a distinct
#' purpose in the final GPU program.
#'
#' A typical OpenCL program used by \pkg{glmbayes} is assembled in the
#' following order:
#'
#' \enumerate{
#'
#'   \item \strong{Global configuration header}  
#'     The file \code{OPENCL.cl} defines global extensions, IEEE constants,
#'     helper macros, and device-side utilities. It plays a role analogous to
#'     a C/C++ header file and must appear at the very top of every combined
#'     kernel module.  
#'     It enables features such as double-precision arithmetic
#'     (\code{cl_khr_fp64}) and device-side debugging (\code{cl_khr_printf}).
#'
#'   \item \strong{Mathematical library modules}  
#'     Subdirectories such as \code{"rmath"}, \code{"dpq"}, and \code{"nmath"}
#'     contain collections of `.cl` files implementing mathematical functions
#'     used throughout the GLM likelihood and gradient computations.  
#'
#'     Each file may declare \code{@provides} and \code{@depends} tags.
#'     \code{load_kernel_library()} reads all files in a subdirectory,
#'     parses these annotations, performs a dependency-aware topological sort,
#'     and concatenates the files in an order that guarantees that upstream
#'     functions appear before downstream callers.  
#'
#'     This mechanism is conceptually similar to a sequence of
#'     \code{#include} statements in C/C++, but with automatic dependency
#'     resolution.
#'
#'   \item \strong{Model-specific helper functions}  
#'     Some kernels require additional device-side utilities that are not part
#'     of the shared libraries. These are typically loaded using
#'     \code{load_kernel_source()} and appended after the library modules.
#'
#'   \item \strong{Final kernel entry function}  
#'     The last component is the model-specific kernel that OpenCL will
#'     execute on the device. For example, in the function
#'     \code{f2_f3_opencl()}, the final kernel is selected based on the GLM
#'     family and link function (e.g., \code{"f2_f3_binomial_logit"}).  
#'
#'     This kernel must appear last in the combined program, after all helper
#'     functions and libraries have been defined, because OpenCL requires that
#'     all functions be defined before they are referenced.
#'
#' }
#'
#' The resulting program is a single, syntactically valid OpenCL source string
#' that is passed directly to the OpenCL compiler (e.g., via
#' \code{clBuildProgram}). The ordering performed by
#' \code{load_kernel_library()} is essential for successful compilation and
#' ensures that the GPU kernels used by \pkg{glmbayes} are reproducible,
#' modular, and maintainable.
#'
#' The function \code{f2_f3_opencl()} provides a concrete example of this
#' assembly process: it loads the global configuration header, then the
#' mathematical libraries, then the model-specific kernel file, and finally
#' concatenates these components into a complete OpenCL program that is sent
#' to the GPU for evaluation of the log-likelihood and gradient.
#' 
#' @param relative_path A file path inside the package's `cl/` directory.
#'   Used by \code{load_kernel_source()} to load a single `.cl` file.
#' @param subdir A subdirectory inside `cl/` containing a set of `.cl` files
#'   annotated with \code{@provides} and \code{@depends} tags. Used by
#'   \code{load_kernel_library()} to construct a dependency-resolved kernel
#'   library.
#' @param package Package name (default: \code{"glmbayes"}).
#' @param verbose Logical; print diagnostic information during dependency
#'   resolution (default: \code{FALSE}).
#'
#' @return
#' A character string containing the kernel source code or combined kernel
#' library.
#'
#' @examples
#' \dontrun{
#' if (has_opencl()) {
#'   src <- load_kernel_source("nmath/bd0.cl")
#'   lib <- load_kernel_library("nmath")
#' }
#' }
#'
#' @export
load_kernel_source <- function(relative_path, package = "glmbayes") {
  if (!has_opencl()) {
    stop("OpenCL support is not available in this build of glmbayes.")
  }
  .load_kernel_source_wrapper(relative_path, package)
}

#' @rdname load_kernel_source
#' @export
load_kernel_library <- function(subdir, package = "glmbayes", verbose = FALSE) {
  if (!has_opencl()) {
    stop("OpenCL support is not available in this build of glmbayes.")
  }
  .load_kernel_library_wrapper(subdir, package, verbose)
}
