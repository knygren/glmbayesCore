// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

/**
 * @file simfuncs_groups.h
 * @brief Row-group partition, prior layout, and groupwise simulation.
 *
 * @namespace glmbayes::sim::group
 * @brief Conditionally independent posterior draws over observation groups.
 *
 * @section ImplementedIn
 *   - group_utils.cpp
 *   - rNormalRegGroups.cpp
 *   - rNormalGLMGroups.cpp
 *
 * @section UsedBy
 *   - export_wrappers.cpp
 */

#ifndef GLMBAYES_SIM_GROUPS_H
#define GLMBAYES_SIM_GROUPS_H

#include "RcppArmadillo.h"

using namespace Rcpp;

namespace glmbayes {

namespace sim {

namespace group {

Rcpp::List normalize_group_cpp(SEXP group, int l2);

Rcpp::List normalize_prior_for_groups_cpp(
    SEXP prior_list_sexp,
    SEXP prior_lists_sexp,
    const List& group_info,
    int l1
);

Rcpp::List prior_payload_from_groups(const List& prior_block, int l1, int k);

Rcpp::List rNormalRegGroups(
    int n,
    NumericVector y,
    NumericMatrix x,
    NumericVector offset,
    NumericVector wt,
    NumericVector dispersion,
    NumericMatrix mu,
    List P_blocks,
    bool prior_by_block,
    List row_blocks,
    Function f2,
    Function f3,
    int Gridtype = 2
);

Rcpp::List rNormalGLMGroups(
    int n,
    NumericVector y,
    NumericMatrix x,
    NumericVector offset,
    NumericVector wt,
    NumericVector dispersion,
    NumericMatrix mu,
    List P_blocks,
    bool prior_by_block,
    List row_blocks,
    Function f2,
    Function f3,
    std::string family,
    std::string link,
    int Gridtype = 2,
    int n_envopt = -1,
    bool use_parallel = true,
    bool use_opencl = false,
    bool verbose = false
);

Rcpp::List group_rNormalReg_cpp_export(
    int n,
    const NumericVector& y,
    const NumericMatrix& x,
    SEXP group,
    SEXP prior_list,
    SEXP prior_lists,
    const NumericVector& offset,
    const NumericVector& wt,
    const Function& f2,
    const Function& f3,
    int Gridtype
);

Rcpp::List group_rNormalGLM_cpp_export(
    int n,
    const NumericVector& y,
    const NumericMatrix& x,
    SEXP group,
    SEXP prior_list,
    SEXP prior_lists,
    const NumericVector& offset,
    const NumericVector& wt,
    const Function& f2,
    const Function& f3,
    const std::string& family,
    const std::string& link,
    int Gridtype,
    int n_envopt,
    bool use_parallel,
    bool use_opencl,
    bool verbose
);

} // namespace group

} // namespace sim

} // namespace glmbayes

#endif
