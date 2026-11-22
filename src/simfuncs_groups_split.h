// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

/**
 * @file simfuncs_groups_split.h
 * @brief One group's rNormalGLM step split into mode, envelope and simulation.
 *
 * @namespace glmbayes::sim::group
 *
 * @section Phases
 *   1. group_mode:     shift the prior mean into the offset, optim (BFGS),
 *                      glmb_Standardize_Model.
 *   2. group_envelope: EnvelopeGridSetup, EnvelopeEval, EnvelopeFinalize
 *                      (EnvelopeBuild split around the evaluation).
 *   3. group_simulate: rNormalGLM_std (or _parallel), undo the
 *                      standardization, add the prior mean back.
 *   Phases 1 and 2 reproduce rNormalGLM step for step, so the mode and
 *   envelope equal rNormalGLM's for the same inputs.
 *
 * @section ImplementedIn
 *   - rNormalGLMGroups_split.cpp
 *
 * @section UsedBy
 *   - rNormalGLMGroups.cpp
 *   - glmbayesCore_ccallables.cpp (rNormalGLMGroupChains, C API)
 */

#ifndef GLMBAYES_SIM_GROUPS_SPLIT_H
#define GLMBAYES_SIM_GROUPS_SPLIT_H

#include "RcppArmadillo.h"

#include <string>

namespace glmbayes {

namespace sim {

namespace group {

/// One group's state after phase 1 (mode finding and standardization).
struct GroupStd {
    // Sliced inputs
    Rcpp::NumericVector y;
    Rcpp::NumericMatrix x;
    Rcpp::NumericVector offset;
    Rcpp::NumericVector wt;
    Rcpp::NumericVector mu;
    Rcpp::NumericMatrix P;
    double dispersion;

    // Shifted model (prior mean moved into alpha; wt2 = wt / dispersion)
    Rcpp::NumericVector wt2;
    Rcpp::NumericVector alpha;

    // optim output: b2a is the mode minus mu (l1 x 1), A1 the Hessian
    Rcpp::List fit;
    Rcpp::NumericMatrix b2a;
    Rcpp::NumericMatrix A1;

    // glmb_Standardize_Model output
    Rcpp::NumericVector bstar2;
    Rcpp::NumericMatrix A;
    Rcpp::NumericMatrix x2;
    Rcpp::NumericMatrix mu2;
    Rcpp::NumericMatrix P2;
    Rcpp::NumericMatrix L2Inv;
    Rcpp::NumericMatrix L3Inv;
};

GroupStd group_mode(
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector offset,
    Rcpp::NumericVector wt,
    Rcpp::NumericVector mu,
    Rcpp::NumericMatrix P,
    double dispersion,
    Rcpp::Function f2,
    Rcpp::Function f3,
    std::string family,
    bool verbose = false
);

/// Seconds spent in the three envelope parts (accumulated).
struct EnvTiming {
    double setup = 0.0;
    double eval = 0.0;
    double finalize = 0.0;
};

/// EnvelopeGridSetup, EnvelopeEval and EnvelopeFinalize in turn; same
/// result as EnvelopeBuild.
Rcpp::List group_envelope(
    const GroupStd& g,
    std::string family,
    std::string link,
    int Gridtype = 2,
    int n = 1,
    int n_envopt = -1,
    bool use_opencl = false,
    bool verbose = false,
    EnvTiming* timing = nullptr
);

/// Returns an rNormalGLM-shaped result list.
Rcpp::List group_simulate(
    const GroupStd& g,
    const Rcpp::List& envelope,
    int n,
    Rcpp::Function f2,
    std::string family,
    std::string link,
    bool use_parallel = true,
    bool verbose = false
);

/// One group drawn for every chain (one draw per chain). Column c of `mu`
/// is chain c's prior mean. Runs group_mode and EnvelopeGridSetup for all
/// chains, then EnvelopeEval for all chains (the step to batch across
/// chains), then EnvelopeFinalize and group_simulate per chain, releasing
/// each envelope after its draw. With `check`, each chain is also run
/// through rNormalGLM on cloned inputs and the call stops if the mode or
/// envelope differs.
/// Returns coefficients and coef.mode (n_chains x l1) and iters (n_chains).
Rcpp::List rNormalGLMGroupChains(
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector offset,
    Rcpp::NumericVector wt,
    double dispersion,
    Rcpp::NumericMatrix mu,
    Rcpp::NumericMatrix P,
    Rcpp::Function f2,
    Rcpp::Function f3,
    std::string family,
    std::string link,
    int Gridtype = 2,
    int n_envopt = -1,
    bool use_opencl = false,
    bool check = false,
    bool verbose = false
);

} // namespace group

} // namespace sim

} // namespace glmbayes

#endif
