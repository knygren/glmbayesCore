############################### Start of load_kernel_source example ####################
## Requires OpenCL in the opencltools build; kernels live in glmbayesCore inst/cl.

if (opencltools::has_opencl()) {
  src <- opencltools::load_kernel_source("nmath/bd0.cl", package = "glmbayesCore")
  lib <- opencltools::load_kernel_library("nmath", package = "glmbayesCore")
  cat("Loaded kernel source length:", nchar(src), "\n")
  cat("Loaded library length:", nchar(lib), "\n")
} else {
  message("OpenCL not enabled in this build of opencltools; skipping example.")
}

## End of load_kernel_source example
