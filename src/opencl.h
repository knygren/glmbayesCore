#pragma once

#include <string>
#include <vector>
#include <Rcpp.h>

namespace glmbayes {
namespace opencl {

// -----------------------------------------------------------------------------
// Low-level OpenCL runner
// -----------------------------------------------------------------------------
//
// Executes the f2/f3 GLM kernel on the GPU.
// This is GLM-specific and belongs in glmbayes::opencl, not openclPort.
//
// kernel_source : full OpenCL program text (helpers + kernel)
// kernel_name   : e.g. "f2_f3_binomial_logit"
// l1, l2, m1    : dims (observations, coefficients, grid size)
// X_flat        : length = l1 * l2, column-major design matrix
// B_flat        : length = m1 * l2, row-major grid of β
// mu_flat       : length = l2, prior means
// P_flat        : length = l2 * l2, prior precision (row-major)
// alpha_flat    : length = l1, offsets
// y_flat        : length = l1, observed responses
// wt_flat       : length = l1, observation weights
// qf_flat       : OUT, length = m1, prior quadratic forms
// grad_flat     : OUT, length = m1 * l2, gradient w.r.t β
// progbar       : 0 = no progress, >0 = show progress bar
//
void f2_f3_kernel_runner(
    const std::string&            kernel_source,
    const char*                   kernel_name,
    int                           l1,
    int                           l2,
    int                           m1,
    const std::vector<double>&    X_flat,
    const std::vector<double>&    B_flat,
    const std::vector<double>&    mu_flat,
    const std::vector<double>&    P_flat,
    const std::vector<double>&    alpha_flat,
    const std::vector<double>&    y_flat,
    const std::vector<double>&    wt_flat,
    std::vector<double>&          qf_flat,
    std::vector<double>&          grad_flat,
    int                           progbar = 0
);


// -----------------------------------------------------------------------------
// High-level Rcpp wrapper
// -----------------------------------------------------------------------------
//
// Loads kernel sources, flattens R matrices, dispatches family/link,
// calls the runner, and reconstructs outputs for R.
//
// b      : l2 × m1 grid of β values
// y      : length l1 responses
// x      : l1 × l2 design matrix
// mu     : l2 × 1 prior mean
// P      : l2 × l2 prior precision
// alpha  : length l1 offsets
// wt     : length l1 weights
// progbar: 0 = no text bar, 1 = show progress
//
Rcpp::List f2_f3_opencl(
    std::string family,
    std::string link,
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar = 0
);

} // namespace opencl
} // namespace glmbayes