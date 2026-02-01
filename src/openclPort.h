#ifndef OPENCLPORT_H
#define OPENCLPORT_H

#include <RcppArmadillo.h>
#include <string>
#include <vector>

#ifdef USE_OPENCL

// Ensure OpenCL types are available
#define CL_TARGET_OPENCL_VERSION 300
#include <CL/cl.h>
#include <string>
#endif 

#ifdef __linux__
#include <stdio.h>
#include <stdlib.h>
#endif

using namespace Rcpp;

// Dependencies:

// 1) OpenCL_helper.cpp
// 2) 



//
// -----------------------------------------------------------------------------
// openclPort: Public API for OpenCL kernel loading, device utilities,
//             and Rcpp → std::vector conversion helpers.
// -----------------------------------------------------------------------------
// Everything a user needs to write OpenCL-enabled wrappers lives here.
// -----------------------------------------------------------------------------
namespace openclPort {

// -------------------------------------------------------------------------
// Rcpp → std::vector conversion utilities
// -------------------------------------------------------------------------
std::vector<double> flattenMatrix(const Rcpp::NumericMatrix& mat);
std::vector<double> copyVector(const Rcpp::NumericVector& vec);

// -------------------------------------------------------------------------
// Device / OpenCL utilities
// -------------------------------------------------------------------------
Rcpp::CharacterVector gpu_names();

// Internal-only GPU detection (used by envelope scaling)
int detect_num_gpus_internal();


// -------------------------------------------------------------------------
// R-facing wrappers for kernel source loading
// -------------------------------------------------------------------------
std::string load_kernel_source_wrapper(
    std::string relative_path,
    std::string package = "glmbayes"
);

std::string load_kernel_library_wrapper(
    std::string subdir,
    std::string package = "glmbayes",
    bool verbose = false
);

// -------------------------------------------------------------------------
// Device / OpenCL utilities
// -------------------------------------------------------------------------

bool has_opencl();
int get_opencl_core_count();


// -------------------------------------------------------------------------
// Conditional declarations: only available when USE_OPENCL is defined
// -------------------------------------------------------------------------
#ifdef USE_OPENCL

// Ensure OpenCL types are available
#define CL_TARGET_OPENCL_VERSION 300
#include <CL/cl.h>
#include <string>


// Load a single .cl kernel file from inst/cl/<relative_path>
std::string load_kernel_source(
    const std::string& relative_path,
    const std::string& package = "glmbayes"
);

// Load and concatenate all .cl files in a subdirectory (inst/cl/<subdir>/)
std::string load_kernel_library(
    const std::string& subdir,
    const std::string& package = "glmbayes",
    bool verbose = false
);


struct OpenCLConfig {
  bool have_expm1;
  bool have_log1p;
  std::string buildOptions;
};

// Probe OpenCL device capabilities and construct build options
OpenCLConfig configureOpenCL(cl_context context,
                             cl_device_id device);

#endif // USE_OPENCL



} // namespace openclPort

#endif // OPENCLPORT_H


