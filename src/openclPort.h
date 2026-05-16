/**
 * @file openclPort.h
 * @brief Public OpenCL interface for glmbayes, including kernel loading,
 *        device discovery, capability probing, and Rcpp-to-std::vector
 *        conversion helpers.
 *
 * @namespace openclPort
 * @brief Lightweight OpenCL utility layer providing kernel management and
 *        device‑level information for optional GPU acceleration.
 *
 * @section ImplementedIn
 *   These declarations are implemented in:
 *     - OpenCL_helper.cpp
 *     - opencl_detect.cpp
 *     - kernel_loader.cpp
 *     - (optional) additional OpenCL backend files guarded by USE_OPENCL
 *
 * @section UsedBy
 *   These functions are consumed by:
 *     - Envelope construction routines (EnvelopeBuild, EnvelopeEval,
 *       EnvelopeDispersionBuild) when OpenCL acceleration is enabled
 *     - R wrappers that expose GPU availability and kernel loading to users
 *
 * @section Responsibilities
 *   Provides:
 *     - Rcpp → std::vector conversion utilities for kernel argument buffers
 *     - GPU/device enumeration and capability checks (gpu_names, has_opencl)
 *     - Kernel source and library loading from inst/cl/ directories
 *     - Conditional OpenCL configuration and build‑option generation
 *
 *   This module:
 *     - is optional and only active when compiled with USE_OPENCL,
 *     - isolates all OpenCL‑specific logic from the statistical code,
 *     - ensures safe fallback to CPU execution when no GPU is available.
 */


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

// Installed OpenCL source tree during migration (legacy copy under inst/cl_old/).
// All kernel loaders resolve paths via system.file(cl_inst_subdir(), ...).
// When inst/cl/ is repopulated, switch this to "cl".
inline const char* cl_inst_subdir() {
  return "cl_old";
}

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
// Kernel file loading (implemented in kernel_loader.cpp). Declared
// unconditionally so optional .cpp sources can compile without -DUSE_OPENCL;
// stubs return safe defaults when OpenCL is disabled at link time.
// -------------------------------------------------------------------------
std::string load_kernel_source(
    const std::string& relative_path,
    const std::string& package = "glmbayes"
);

std::string load_kernel_library(
    const std::string& subdir,
    const std::string& package = "glmbayes",
    bool verbose = false
);

std::string load_library_for_kernel(
    const std::string& kernel_relative_path,
    const std::string& library_subdir,
    const std::string& package = "glmbayes",
    const std::string& depends_tag = "depends_nmath"
);

bool kernel_all_depends_nmath_includes_qDiscrete_search(
    const std::string& kernel_relative_path,
    const std::string& package = "glmbayes"
);

std::string load_library_for_kernel_wrapper(
    std::string kernel_relative_path,
    std::string library_subdir,
    std::string package = "glmbayes",
    std::string depends_tag = "depends_nmath"
);

bool kernel_all_depends_nmath_includes_qDiscrete_search_wrapper(
    std::string kernel_relative_path,
    std::string package = "glmbayes"
);


// -------------------------------------------------------------------------
// Conditional declarations: only available when USE_OPENCL is defined
// -------------------------------------------------------------------------
#ifdef USE_OPENCL

// Ensure OpenCL types are available
#define CL_TARGET_OPENCL_VERSION 300
#include <CL/cl.h>
#include <string>


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


