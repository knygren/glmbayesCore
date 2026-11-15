// C-callable API for packages that use `LinkingTo: glmbayesCore`.
// Declarations for consumers: inst/include/glmbayesCore/glmbayesCore_capi.h.

#include <R_ext/Rdynload.h>
#include <R_ext/Error.h>
#include "RcppArmadillo.h"
#include "simfuncs.h"
#include "simfuncs_groups.h"

using namespace Rcpp;

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

} // extern "C"

// [[Rcpp::export]]
void register_glmbayesCore_ccallables_cpp_export() {
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_api_version",       (DL_FUNC) &glmbayesCore_api_version);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_rNormalGLMGroups",  (DL_FUNC) &glmbayesCore_rNormalGLMGroups);
  R_RegisterCCallable("glmbayesCore", "glmbayesCore_rNormalReg",        (DL_FUNC) &glmbayesCore_rNormalReg);
}
