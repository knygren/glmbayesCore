#pragma once

#include <Rcpp.h>

// Rcpp‐exported wrapper that loads your .cl, flattens inputs, calls the runner,
// and reconstructs the xb/qf outputs for cpp.
//
// b      : l2 × m1 grid of β values
// y      : length l1 responses (unused in prep)
// x      : l1 × l2 design matrix
// mu     : l2 × 1 prior mean
// P      : l2 × l2 prior precision
// alpha  : length l1 offsets
// wt     : length l1 weights (unused here)
// progbar: 0 = no text bar, 1 = show progress


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
    int                  progbar=0
);

// Rcpp-exported wrapper that returns the number of detected OpenCL-compatible GPU cores.
// If OpenCL is not available, returns 1 as a safe fallback.
// Implementation lives in kernel_wrappers.cpp and wraps detect_num_gpus_internal().
int get_opencl_core_count();