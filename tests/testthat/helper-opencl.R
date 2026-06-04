skip_if_no_opencl <- function() {
  skip_if(!glmbayesCore_has_opencl(), "OpenCL not enabled in this build of glmbayesCore")
}
