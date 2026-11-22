// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

/**
 * @file Envelopefuncs_split.h
 * @brief EnvelopeBuild split around EnvelopeEval.
 *
 * @namespace glmbayes::env
 *
 * @section Parts
 *   1. EnvelopeGridSetup: tangent points, EnvelopeSize, expand.grid
 *                         (grid G3/G4, GIndex, Lint, E_draws).
 *   2. EnvelopeEval:      unchanged; NegLL and gradients at the grid points.
 *   3. EnvelopeFinalize:  EnvelopeSet_Grid, setlogP, PLSD, optional sort.
 *   Each part copies the matching lines of EnvelopeBuild, so
 *   setup + EnvelopeEval + finalize returns the same envelope. Keeping the
 *   evaluation separate allows the grids of many groups to be evaluated in
 *   one batch.
 *
 * @section ImplementedIn
 *   - EnvelopeBuild_split.cpp
 */

#ifndef GLMBAYES_ENVELOPEFUNCS_SPLIT_H
#define GLMBAYES_ENVELOPEFUNCS_SPLIT_H

#include "RcppArmadillo.h"

namespace glmbayes {

namespace env {

/// Grid produced by EnvelopeGridSetup and consumed by EnvelopeFinalize.
struct EnvGrid {
    Rcpp::NumericMatrix GIndex;  // l2 x l1 region codes
    Rcpp::NumericMatrix G3;      // l2 x l1 grid points
    Rcpp::NumericMatrix G4;      // l1 x l2 (EnvelopeEval input)
    Rcpp::NumericMatrix Lint1;   // 2 x l1 central-interval cutpoints
    Rcpp::NumericVector a_1;     // diag(A)
    double E_draws;
    int l1;
    int l2;
};

/// EnvelopeBuild's compile-time OpenCL fallback.
bool envelope_use_opencl(bool use_opencl, bool verbose);

EnvGrid EnvelopeGridSetup(Rcpp::NumericVector bStar,
                          Rcpp::NumericMatrix A,
                          int Gridtype,
                          int n,
                          int n_envopt,
                          bool use_opencl,
                          bool verbose);

/// x is only used for the diagnostic when the PLSD construction fails.
Rcpp::List EnvelopeFinalize(const EnvGrid& grid,
                            Rcpp::NumericVector NegLL,
                            const arma::mat& cbars_eval,
                            Rcpp::NumericMatrix x,
                            bool sortgrid,
                            bool verbose);

} // namespace env

} // namespace glmbayes

#endif
