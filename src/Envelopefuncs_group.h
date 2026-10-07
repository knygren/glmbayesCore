// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

/**
 * @file Envelopefuncs_group.h
 * @brief Group-factored independent Normal-Gamma envelopes.
 *
 * @namespace glmbayes::env::group
 *
 * @section ImplementedIn
 *   - block_rIndepNormalGammaReg.cpp
 *
 * @section UsedBy
 *   - export_wrappers.cpp
 */

#ifndef GLMBAYES_ENV_GROUP_H
#define GLMBAYES_ENV_GROUP_H

#include "RcppArmadillo.h"

using namespace Rcpp;

namespace glmbayes {

namespace env {

namespace group {

/// BlockEnvelopeCentering: pooled dispersion anchoring for block ING (Step 1).
List BlockEnvelopeCentering(
    NumericVector y,
    NumericMatrix x,
    SEXP block,
    SEXP prior_list,
    SEXP prior_lists,
    NumericVector offset,
    NumericVector wt,
    double shape,
    double rate,
    double max_disp_perc,
    Nullable<double> disp_lower,
    Nullable<double> disp_upper,
    int p_re,
    int n_rss_iter,
    bool verbose
);

/// BlockEnvelopeBuild: per-block beta envelopes + shared dispersion stub (Step 2).
List BlockEnvelopeBuild(
    const List& centering_out,
    NumericVector y,
    NumericMatrix x,
    SEXP block,
    SEXP prior_list,
    SEXP prior_lists,
    NumericVector offset,
    NumericVector wt,
    double max_disp_perc,
    Nullable<double> disp_lower,
    Nullable<double> disp_upper,
    int n,
    int Gridtype,
    Nullable<int> n_envopt,
    double RSS_ML,
    bool use_parallel,
    bool use_opencl,
    bool verbose
);

/// BlockEnvelopeDispersionBuild: per-block EDB + pooled sigma^2 constants (Step 3).
List BlockEnvelopeDispersionBuild(
    const List& build_out,
    const List& centering_out,
    NumericVector y,
    NumericMatrix x,
    SEXP block,
    NumericVector offset,
    NumericVector wt,
    double shape,
    double rate,
    double max_disp_perc,
    Nullable<double> disp_lower,
    Nullable<double> disp_upper,
    double RSS_ML,
    bool use_parallel,
    bool verbose
);

/// BlockEnvelopeSim: envelope-proposal draws per block (phase 1 auto-accept).
List BlockEnvelopeSim(
    const List& build_out,
    int n,
    bool progbar,
    bool verbose
);

/// Block ING pipeline orchestrator (Centering → Build → DispersionBuild → Sim).
List rIndepNormalGammaRegBlock(
    int n,
    NumericVector y,
    NumericMatrix x,
    SEXP block,
    SEXP prior_list,
    SEXP prior_lists,
    NumericVector offset,
    NumericVector wt,
    int p_re,
    int n_rss_iter,
    int Gridtype,
    Nullable<int> n_envopt,
    double RSS_ML,
    bool use_parallel,
    bool use_opencl,
    bool progbar,
    bool verbose,
    CharacterVector group_levels,
    CharacterVector re_names
);

/// Independent-block separable-overbound variant of BlockEnvelopeDispersionBuild.
/// Skips build_joint_product_face_slack; faces drawn from per-block PLSDs.
List BlockEnvelopeDispersionBuildInd(
    const List& build_out,
    const List& centering_out,
    NumericVector y,
    NumericMatrix x,
    SEXP block,
    NumericVector offset,
    NumericVector wt,
    double shape,
    double rate,
    double max_disp_perc,
    Nullable<double> disp_lower,
    Nullable<double> disp_upper,
    double RSS_ML,
    bool use_parallel,
    bool verbose
);

/// Independent-block variant of BlockEnvelopeSim.
/// Always draws faces independently from each block's own PLSD.
List BlockEnvelopeSimInd(
    const List& build_out,
    int n,
    bool progbar,
    bool verbose
);

/// Orchestrator for the independent-block pipeline
/// (Centering → Build → DispersionBuildInd → SimInd).
List rIndepNormalGammaRegBlockInd(
    int n,
    NumericVector y,
    NumericMatrix x,
    SEXP block,
    SEXP prior_list,
    SEXP prior_lists,
    NumericVector offset,
    NumericVector wt,
    int p_re,
    int n_rss_iter,
    int Gridtype,
    Nullable<int> n_envopt,
    double RSS_ML,
    bool use_parallel,
    bool use_opencl,
    bool progbar,
    bool verbose,
    CharacterVector group_levels,
    CharacterVector re_names
);

}  // namespace group

}  // namespace env

}  // namespace glmbayes

#endif
