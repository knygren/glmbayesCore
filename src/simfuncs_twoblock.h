// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

/**
 * @file simfuncs_twoblock.h
 * @brief Two-block Gibbs driver and the groups and population updates.
 *
 * @namespace glmbayes::sim::twoblock
 *
 * @section ImplementedIn
 *   - twoBlockGibbs.cpp
 *   - two_block_groups.cpp
 *   - two_block_population.cpp
 *
 * @section UsedBy
 *   - export_wrappers.cpp
 */

#ifndef GLMBAYES_SIM_TWOBLOCK_H
#define GLMBAYES_SIM_TWOBLOCK_H

#include "RcppArmadillo.h"

#include <vector>

using namespace Rcpp;

namespace glmbayes {

namespace sim {

namespace twoblock {

/// Two-block Gibbs driver (twoBlockGibbs.cpp): population priors are
/// pfamily objects (dNormal / dIndependent_Normal_Gamma); dispatch on the
/// pfamily type string. Returns dispersion_fixef_draws with the stored draws.
Rcpp::List two_block_rNormal_reg_cpp_export(
    int n,
    int m_convergence,
    const NumericVector& y,
    const NumericMatrix& x,
    SEXP block,
    const List& x_hyper,
    const List& prior_list_block1,
    SEXP dispersion_block1,
    SEXP ddef_block1,
    const List& pfamily_list,
    const List& fixef_start,
    const CharacterVector& group_levels,
    const std::string& family,
    const std::string& link,
    const Function& f2,
    const Function& f3,
    const Function& f2_gauss,
    const Function& f3_gauss,
    const NumericVector& offset,
    const NumericVector& wt,
    int Gridtype,
    int n_envopt,
    bool use_parallel,
    bool use_opencl,
    bool verbose,
    bool progbar
);


/// Cached \code{X_hyper} row lookups for \code{build_mu_all} (two_block_groups.cpp).
struct MuAllBuilder {
  int p_re;
  int J;
  std::vector<Rcpp::NumericMatrix> X;
  std::vector<std::vector<int>> row_idx;

  MuAllBuilder(
      const Rcpp::List& x_hyper,
      const Rcpp::CharacterVector& group_levels
  );
  Rcpp::NumericMatrix build(
      const std::vector<Rcpp::NumericVector>& fixef
  ) const;
};

/// Block~1: \code{mu_all} matrix (port of \code{build_mu_all()}).
Rcpp::NumericMatrix two_block_build_mu_all(
    const Rcpp::List& x_hyper,
    const Rcpp::List& fixef,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::CharacterVector& group_levels
);

/// Block~1: ING \code{P} refresh from \code{tau2} (two_block_groups.cpp).
Rcpp::List two_block_groups_prior_with_tau2(
    const Rcpp::List& base_prior,
    const Rcpp::NumericVector& tau2_vec,
    const Rcpp::CharacterVector& ptypes,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::NumericMatrix& mu_all
);

/// Block~1: mean envelope iters across groups (two_block_groups.cpp).
double two_block_groups_iters_mean(const Rcpp::List& block_out);

/// Extract chain-\code{i} fixef list (two_block_groups.cpp).
/// Port of \code{.two_block_batch_fixef_chain}.
Rcpp::List fixef_list_from_batch_chain(
    const Rcpp::List& batch_fixef,
    int chain_i,
    const Rcpp::CharacterVector& re_names
);

/// All-chains step A: \code{batch$tau2[chain_i, ]} (two_block_groups.cpp).
Rcpp::NumericVector batch_tau2_chain_row(
    const Rcpp::NumericMatrix& batch_tau2,
    int chain_i
);

/// All-chains step C: \code{batch$b[, , chain_i] <- b_draw} (two_block_groups.cpp).
void batch_b_assign_slice(
    Rcpp::NumericVector& b_store,
    int chain_i,
    const Rcpp::NumericMatrix& b_draw
);

/// \code{dim(b)[3]} for batch \code{b} (two_block_groups.cpp).
int batch_b_n_chains(const Rcpp::NumericVector& batch_b);

/// Restore \code{dim}/\code{dimnames} on batch \code{b} (two_block_groups.cpp).
void ensure_batch_b_dimnames(
    Rcpp::NumericVector& batch_b,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::CharacterVector& re_names,
    int n_chains
);

/// All-chains step D: \code{batch$iters_ranef[chain_i] += iters_mean}
/// (two_block_groups.cpp).
void batch_iters_ranef_add(
    Rcpp::NumericVector& iters_ranef,
    int chain_i,
    double iters_mean
);

/// Block~1: reorder coefficient rows to group_levels (two_block_groups.cpp).
Rcpp::NumericMatrix two_block_reorder_b_to_group_levels(
    Rcpp::NumericMatrix b_draw,
    SEXP block_ids,
    const Rcpp::CharacterVector& group_levels
);

/// Block~1 prep + draw for one chain (two_block_groups.cpp).
Rcpp::List two_block_groups_one_chain_impl(
    int chain_i,
    const Rcpp::List& batch_fixef,
    const Rcpp::NumericVector& tau2_i,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& Z,
    SEXP groups,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    const Rcpp::List& x_hyper,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::CharacterVector& ptypes,
    const Rcpp::List& block1_prior,
    bool is_gaussian,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::Function& f2_gauss,
    const Rcpp::Function& f3_gauss,
    const std::string& family,
    const std::string& link,
    int Gridtype,
    int n_envopt
);

/// Full per-chain Block~1 orchestration steps A→D (two_block_groups.cpp).
Rcpp::List two_block_groups_one_chain_orchestrate_impl(
    int chain_i,
    Rcpp::NumericVector b_store,
    Rcpp::NumericVector iters_ranef,
    const Rcpp::List& batch_fixef,
    const Rcpp::NumericMatrix& batch_tau2,
    const Rcpp::List& design,
    const Rcpp::List& block1_prior,
    SEXP family,
    const Rcpp::CharacterVector& ptypes,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::Function& f2_gauss,
    const Rcpp::Function& f3_gauss,
    bool use_cpp_tau2_row,
    bool use_cpp_b_slice,
    bool use_cpp_iters_ranef_add
);

/// Block~1 draw from R-built \code{mu_all} and \code{P} (two_block_groups.cpp).
Rcpp::List two_block_groups_one_chain_from_mu_P_impl(
    const Rcpp::NumericMatrix& mu_all,
    const Rcpp::NumericMatrix& P,
    SEXP dispersion,
    SEXP ddef,
    const Rcpp::List& design,
    SEXP family,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::Function& f2_gauss,
    const Rcpp::Function& f3_gauss
);


/// Population align (two_block_population.cpp): port of R two_block_align_b_to_xhyper().
Rcpp::NumericVector two_block_align_b_to_xhyper_cpp(
    Rcpp::NumericVector b_vec,
    Rcpp::NumericMatrix X_k,
    Rcpp::CharacterVector group_levels
);

/// Population one-chain update (two_block_population.cpp): align + rglmb per RE.
Rcpp::List two_block_population_one_chain_cpp_export(
    const Rcpp::NumericMatrix& b_i,
    const Rcpp::List& fixef_rows,
    const Rcpp::NumericVector& tau2_i,
    const Rcpp::NumericVector& iters_i,
    const Rcpp::List& x_hyper,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::List& pfamily_list,
    const Rcpp::CharacterVector& ptypes,
    const Rcpp::CharacterVector& re_names
);

}  // namespace twoblock

}  // namespace sim

}  // namespace glmbayes

#endif
