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
 *     - R wrappers that expose GPU availability (glmbayesCore_has_opencl); kernel loading
 *       is provided by opencltools at the R level
 *
 * @section Responsibilities
 *   Provides:
 *     - Rcpp → std::vector conversion utilities for kernel argument buffers
 *     - GPU/device enumeration and capability checks (gpu_names, glmbayesCore_has_opencl)
 *     - Kernel source and library loading from inst/cl/ directories
 *       (generic loaders; full likelihood-subgradient program assembly is in
 *       glmbayes::opencl::load_likelihood_subgradient_program — see inst/cl/README.md)
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
#include "package_ns.h"

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

// OpenCL compute units via opencltools::get_opencl_core_count() (envelope scaling)
int opencl_core_count_for_scaling();

// -------------------------------------------------------------------------
// Device / OpenCL utilities
// -------------------------------------------------------------------------

bool glmbayesCore_has_opencl();


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
    const std::string& package = GLMBAYES_R_NS
);

// Load and concatenate all .cl files in a subdirectory (inst/cl/<subdir>/)
std::string load_kernel_library(
    const std::string& subdir,
    const std::string& package = GLMBAYES_R_NS,
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


