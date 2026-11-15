// C-callable API for packages that use `LinkingTo: glmbayesCore`.
// Declarations for consumers: inst/include/glmbayesCore/glmbayesCore_capi.h.

#include <R_ext/Rdynload.h>
#include <R_ext/Error.h>
#include "RcppArmadillo.h"
#include "simfuncs.h"
#include "simfuncs_groups.h"
#include "simfuncs_twoblock.h"
#include "progress_utils.h"

using namespace Rcpp;

namespace glmbayes {
namespace sim {
namespace twoblock {
Rcpp::List x_hyper_as_matrix(const Rcpp::List& x_hyper);
}
}
}

extern "C" {

int glmbayesCore_api_version(void) {
  return 1;
}

SEXP glmbayesCore_rNormalGLMGroups(SEXP n, SEXP y, SEXP x, SEXP offset,
                                   SEXP wt, SEXP dispersion, SEXP mu,
                                   SEXP P_blocks, SEXP prior_by_block,
                                   SEXP row_blocks, SEXP f2, SEXP f3,
                                   SEXP family, SEXP link, SEXP Gridtype,
                                   SEXP n_envopt, SEXP use_parallel,
                                   SEXP use_opencl, SEXP verbose) {
  try {
    return Rcpp::wrap(glmbayes::sim::group::rNormalGLMGroups(
        Rcpp::as<int>(n),
        Rcpp::as<NumericVector>(y),
        Rcpp::as<NumericMatrix>(x),
        Rcpp::as<NumericVector>(offset),
        Rcpp::as<NumericVector>(wt),
        Rcpp::as<NumericVector>(dispersion),
        Rcpp::as<NumericMatrix>(mu),
        Rcpp::as<List>(P_blocks),
        Rcpp::as<bool>(prior_by_block),
        Rcpp::as<List>(row_blocks),
        Rcpp::as<Function>(f2),
        Rcpp::as<Function>(f3),
        Rcpp::as<std::string>(family),
        Rcpp::as<std::string>(link),
        Rcpp::as<int>(Gridtype),
        Rcpp::as<int>(n_envopt),
        Rcpp::as<bool>(use_parallel),
        Rcpp::as<bool>(use_opencl),
        Rcpp::as<bool>(verbose)));
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  } catch (...) {
    Rf_error("glmbayesCore_rNormalGLMGroups failed");
  }
  return R_NilValue;
}

SEXP glmbayesCore_rNormalReg(SEXP n, SEXP y, SEXP x, SEXP mu, SEXP P,
                             SEXP offset, SEXP wt, SEXP dispersion, SEXP f2,
                             SEXP f3, SEXP start, SEXP family, SEXP link,
                             SEXP Gridtype) {
  try {
    return Rcpp::wrap(glmbayes::sim::rNormalReg(
        Rcpp::as<int>(n),
        Rcpp::as<NumericVector>(y),
        Rcpp::as<NumericMatrix>(x),
        Rcpp::as<NumericVector>(mu),
        Rcpp::as<NumericMatrix>(P),
        Rcpp::as<NumericVector>(offset),
        Rcpp::as<NumericVector>(wt),
        Rcpp::as<double>(dispersion),
        Rcpp::as<Function>(f2),
        Rcpp::as<Function>(f3),
        Rcpp::as<NumericVector>(start),
        Rcpp::as<std::string>(family),
        Rcpp::as<std::string>(link),
        Rcpp::as<int>(Gridtype)));
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  } catch (...) {
    Rf_error("glmbayesCore_rNormalReg failed");
  }
  return R_NilValue;
}

SEXP glmbayesCore_group_rNormalGLM(SEXP n, SEXP y, SEXP x, SEXP group,
                                   SEXP prior_list, SEXP prior_lists,
                                   SEXP offset, SEXP wt, SEXP f2, SEXP f3,
                                   SEXP family, SEXP link, SEXP Gridtype,
                                   SEXP n_envopt, SEXP use_parallel,
                                   SEXP use_opencl, SEXP verbose) {
  try {
    return Rcpp::wrap(glmbayes::sim::group::group_rNormalGLM_cpp_export(
        Rcpp::as<int>(n),
        Rcpp::as<NumericVector>(y),
        Rcpp::as<NumericMatrix>(x),
        group, prior_list, prior_lists,
        Rcpp::as<NumericVector>(offset),
        Rcpp::as<NumericVector>(wt),
        Rcpp::as<Function>(f2),
        Rcpp::as<Function>(f3),
        Rcpp::as<std::string>(family),
        Rcpp::as<std::string>(link),
        Rcpp::as<int>(Gridtype),
        Rcpp::as<int>(n_envopt),
        Rcpp::as<bool>(use_parallel),
        Rcpp::as<bool>(use_opencl),
        Rcpp::as<bool>(verbose)));
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  } catch (...) {
    Rf_error("glmbayesCore_group_rNormalGLM failed");
  }
  return R_NilValue;
}

SEXP glmbayesCore_group_rNormalReg(SEXP n, SEXP y, SEXP x, SEXP group,
                                   SEXP prior_list, SEXP prior_lists,
                                   SEXP offset, SEXP wt, SEXP f2, SEXP f3,
                                   SEXP Gridtype) {
  try {
    return Rcpp::wrap(glmbayes::sim::group::group_rNormalReg_cpp_export(
        Rcpp::as<int>(n),
        Rcpp::as<NumericVector>(y),
        Rcpp::as<NumericMatrix>(x),
        group, prior_list, prior_lists,
        Rcpp::as<NumericVector>(offset),
        Rcpp::as<NumericVector>(wt),
        Rcpp::as<Function>(f2),
        Rcpp::as<Function>(f3),
        Rcpp::as<int>(Gridtype)));
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  } catch (...) {
    Rf_error("glmbayesCore_group_rNormalReg failed");
  }
  return R_NilValue;
}

SEXP glmbayesCore_two_block_build_mu_all(SEXP x_hyper, SEXP fixef,
                                         SEXP re_names, SEXP group_levels) {
  try {
    return Rcpp::wrap(glmbayes::sim::twoblock::two_block_build_mu_all(
        Rcpp::as<List>(x_hyper),
        Rcpp::as<List>(fixef),
        Rcpp::as<CharacterVector>(re_names),
        Rcpp::as<CharacterVector>(group_levels)));
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  } catch (...) {
    Rf_error("glmbayesCore_two_block_build_mu_all failed");
  }
  return R_NilValue;
}

SEXP glmbayesCore_two_block_groups_prior_with_tau2(SEXP base_prior,
                                                   SEXP tau2_vec, SEXP ptypes,
                                                   SEXP re_names,
                                                   SEXP mu_all) {
  try {
    return Rcpp::wrap(glmbayes::sim::twoblock::two_block_groups_prior_with_tau2(
        Rcpp::as<List>(base_prior),
        Rcpp::as<NumericVector>(tau2_vec),
        Rcpp::as<CharacterVector>(ptypes),
        Rcpp::as<CharacterVector>(re_names),
        Rcpp::as<NumericMatrix>(mu_all)));
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  } catch (...) {
    Rf_error("glmbayesCore_two_block_groups_prior_with_tau2 failed");
  }
  return R_NilValue;
}

SEXP glmbayesCore_two_block_groups_one_chain_from_mu_P(
    SEXP mu_all, SEXP P, SEXP dispersion, SEXP ddef, SEXP design,
    SEXP family, SEXP re_names, SEXP group_levels, SEXP f2, SEXP f3,
    SEXP f2_gauss, SEXP f3_gauss) {
  try {
    return Rcpp::wrap(
        glmbayes::sim::twoblock::two_block_groups_one_chain_from_mu_P_impl(
            Rcpp::as<NumericMatrix>(mu_all),
            Rcpp::as<NumericMatrix>(P),
            dispersion, ddef,
            Rcpp::as<List>(design),
            family,
            Rcpp::as<CharacterVector>(re_names),
            Rcpp::as<CharacterVector>(group_levels),
            Rcpp::as<Function>(f2),
            Rcpp::as<Function>(f3),
            Rcpp::as<Function>(f2_gauss),
            Rcpp::as<Function>(f3_gauss)));
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  } catch (...) {
    Rf_error("glmbayesCore_two_block_groups_one_chain_from_mu_P failed");
  }
  return R_NilValue;
}

SEXP glmbayesCore_fixef_list_from_batch_chain(SEXP batch_fixef, SEXP chain_i,
                                              SEXP re_names) {
  try {
    return Rcpp::wrap(glmbayes::sim::twoblock::fixef_list_from_batch_chain(
        Rcpp::as<List>(batch_fixef),
        Rcpp::as<int>(chain_i),
        Rcpp::as<CharacterVector>(re_names)));
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  } catch (...) {
    Rf_error("glmbayesCore_fixef_list_from_batch_chain failed");
  }
  return R_NilValue;
}

SEXP glmbayesCore_batch_tau2_chain_row(SEXP batch_tau2, SEXP chain_i) {
  try {
    return Rcpp::wrap(glmbayes::sim::twoblock::batch_tau2_chain_row(
        Rcpp::as<NumericMatrix>(batch_tau2),
        Rcpp::as<int>(chain_i)));
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  } catch (...) {
    Rf_error("glmbayesCore_batch_tau2_chain_row failed");
  }
  return R_NilValue;
}

// b_store is modified in place; it must be a double array.
SEXP glmbayesCore_batch_b_assign_slice(SEXP b_store, SEXP chain_i,
                                       SEXP b_draw) {
  try {
    NumericVector b(b_store);
    glmbayes::sim::twoblock::batch_b_assign_slice(
        b, Rcpp::as<int>(chain_i), Rcpp::as<NumericMatrix>(b_draw));
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  } catch (...) {
    Rf_error("glmbayesCore_batch_b_assign_slice failed");
  }
  return R_NilValue;
}

SEXP glmbayesCore_batch_b_n_chains(SEXP batch_b) {
  try {
    return Rcpp::wrap(glmbayes::sim::twoblock::batch_b_n_chains(
        Rcpp::as<NumericVector>(batch_b)));
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  } catch (...) {
    Rf_error("glmbayesCore_batch_b_n_chains failed");
  }
  return R_NilValue;
}

// batch_b is modified in place; it must be a double array.
SEXP glmbayesCore_ensure_batch_b_dimnames(SEXP batch_b, SEXP group_levels,
                                          SEXP re_names, SEXP n_chains) {
  try {
    NumericVector b(batch_b);
    glmbayes::sim::twoblock::ensure_batch_b_dimnames(
        b,
        Rcpp::as<CharacterVector>(group_levels),
        Rcpp::as<CharacterVector>(re_names),
        Rcpp::as<int>(n_chains));
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  } catch (...) {
    Rf_error("glmbayesCore_ensure_batch_b_dimnames failed");
  }
  return R_NilValue;
}

// iters_ranef is modified in place; it must be a double vector.
SEXP glmbayesCore_batch_iters_ranef_add(SEXP iters_ranef, SEXP chain_i,
                                        SEXP iters_mean) {
  try {
    NumericVector it(iters_ranef);
    glmbayes::sim::twoblock::batch_iters_ranef_add(
        it, Rcpp::as<int>(chain_i), Rcpp::as<double>(iters_mean));
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  } catch (...) {
    Rf_error("glmbayesCore_batch_iters_ranef_add failed");
  }
  return R_NilValue;
}

SEXP glmbayesCore_x_hyper_as_matrix(SEXP x_hyper) {
  try {
    return Rcpp::wrap(glmbayes::sim::twoblock::x_hyper_as_matrix(
        Rcpp::as<List>(x_hyper)));
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  } catch (...) {
    Rf_error("glmbayesCore_x_hyper_as_matrix failed");
  }
  return R_NilValue;
}

SEXP glmbayesCore_rIndepNormalGammaReg(SEXP n, SEXP y, SEXP x, SEXP mu,
                                       SEXP P, SEXP offset, SEXP wt,
                                       SEXP shape, SEXP rate,
                                       SEXP max_disp_perc, SEXP disp_lower,
                                       SEXP disp_upper, SEXP Gridtype,
                                       SEXP n_envopt, SEXP use_parallel,
                                       SEXP use_opencl, SEXP verbose,
                                       SEXP progbar) {
  try {
    return Rcpp::wrap(glmbayes::sim::rIndepNormalGammaReg(
        Rcpp::as<int>(n),
        Rcpp::as<NumericVector>(y),
        Rcpp::as<NumericMatrix>(x),
        Rcpp::as<NumericVector>(mu),
        Rcpp::as<NumericMatrix>(P),
        Rcpp::as<NumericVector>(offset),
        Rcpp::as<NumericVector>(wt),
        Rcpp::as<double>(shape),
        Rcpp::as<double>(rate),
        Rcpp::as<double>(max_disp_perc),
        Rcpp::Nullable<NumericVector>(disp_lower),
        Rcpp::Nullable<NumericVector>(disp_upper),
        Rcpp::as<int>(Gridtype),
        Rcpp::as<int>(n_envopt),
        Rcpp::as<bool>(use_parallel),
        Rcpp::as<bool>(use_opencl),
        Rcpp::as<bool>(verbose),
        Rcpp::as<bool>(progbar)));
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  } catch (...) {
    Rf_error("glmbayesCore_rIndepNormalGammaReg failed");
  }
  return R_NilValue;
}

SEXP glmbayesCore_progress_bar(SEXP x, SEXP N, SEXP prefix) {
  try {
    glmbayes::progress::progress_bar(
        Rcpp::as<double>(x), Rcpp::as<double>(N),
        Rcpp::as<std::string>(prefix));
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  } catch (...) {
    Rf_error("glmbayesCore_progress_bar failed");
  }
  return R_NilValue;
}

SEXP glmbayesCore_progress_bar_finish(SEXP newline) {
  try {
    glmbayes::progress::progress_bar_finish(Rcpp::as<bool>(newline));
  } catch (const std::exception& e) {
    Rf_error("%s", e.what());
  } catch (...) {
    Rf_error("glmbayesCore_progress_bar_finish failed");
  }
  return R_NilValue;
}

} // extern "C"

// [[Rcpp::export]]
void register_glmbayesCore_ccallables_cpp_export() {
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_api_version",       (DL_FUNC) &glmbayesCore_api_version);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_rNormalGLMGroups",  (DL_FUNC) &glmbayesCore_rNormalGLMGroups);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_rNormalReg",        (DL_FUNC) &glmbayesCore_rNormalReg);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_group_rNormalGLM",  (DL_FUNC) &glmbayesCore_group_rNormalGLM);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_group_rNormalReg",  (DL_FUNC) &glmbayesCore_group_rNormalReg);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_two_block_build_mu_all", (DL_FUNC) &glmbayesCore_two_block_build_mu_all);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_two_block_groups_prior_with_tau2", (DL_FUNC) &glmbayesCore_two_block_groups_prior_with_tau2);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_two_block_groups_one_chain_from_mu_P", (DL_FUNC) &glmbayesCore_two_block_groups_one_chain_from_mu_P);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_fixef_list_from_batch_chain", (DL_FUNC) &glmbayesCore_fixef_list_from_batch_chain);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_batch_tau2_chain_row", (DL_FUNC) &glmbayesCore_batch_tau2_chain_row);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_batch_b_assign_slice", (DL_FUNC) &glmbayesCore_batch_b_assign_slice);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_batch_b_n_chains", (DL_FUNC) &glmbayesCore_batch_b_n_chains);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_ensure_batch_b_dimnames", (DL_FUNC) &glmbayesCore_ensure_batch_b_dimnames);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_batch_iters_ranef_add", (DL_FUNC) &glmbayesCore_batch_iters_ranef_add);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_x_hyper_as_matrix", (DL_FUNC) &glmbayesCore_x_hyper_as_matrix);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_rIndepNormalGammaReg", (DL_FUNC) &glmbayesCore_rIndepNormalGammaReg);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_progress_bar", (DL_FUNC) &glmbayesCore_progress_bar);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_progress_bar_finish", (DL_FUNC) &glmbayesCore_progress_bar_finish);
}
