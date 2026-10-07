#include "RcppArmadillo.h"
#include "Envelopefuncs.h"
#include "openclPort.h"
#include "simfuncs.h"
#include "simfuncs_groups.h"
#include "simfuncs_twoblock.h"
#include "Envelopefuncs_group.h"

using namespace openclPort;
using namespace glmbayes::env;
using namespace glmbayes::sim;


// -----------------------------------------------------------------------------
// Wrapper organization mirrors R/rcpp_wrappers.R:
//   Tier 1: Core Simulation   - Main sampling entry points
//   Tier 2: Envelope          - Build/eval; used by rNormalGLM
//   Tier 3: Indep NG std      - Split workflow samplers
//   Tier 4: Model Utilities   - Standardization
//   Tier 5: OpenCL/GPU        - Kernel loading, diagnostics
//   Phased out: rss_face_at_disp, UB2 (no R wrappers)
// -----------------------------------------------------------------------------


// =============================================================================
// Tier 1: Core Simulation
// Callers: rNormal_reg, rNormalGamma_reg, rindepNormalGamma_reg, rGamma_reg,
// User:    All users - primary paths via rglmb, rlmb, glmb, pfamily
// =============================================================================

// [[Rcpp::export]]
Rcpp::List rNormalGLM_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    double dispersion,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::NumericVector& start,
    const std::string& family = "binomial",
    const std::string& link   = "logit",
    int Gridtype = 2,
    int n_envopt = -1,
    bool use_parallel = true,
    bool use_opencl = false,
    bool verbose = false
) {
  return rNormalGLM(
    n, y, x, mu, P, offset, wt,
    dispersion,
    f2, f3, start,
    family, link, Gridtype,
    n_envopt, use_parallel, use_opencl, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List rNormalReg_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    double dispersion,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::NumericVector& start,
    const std::string& family = "gaussian",
    const std::string& link   = "identity",
    int Gridtype = 2
) {
  return rNormalReg(
    n, y, x, mu, P, offset, wt,
    dispersion, f2, f3, start,
    family, link, Gridtype
  );
}

// [[Rcpp::export]]
Rcpp::List rNormalGLMGroups_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    const Rcpp::NumericVector& dispersion,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::List& P_blocks,
    bool prior_by_block,
    const Rcpp::List& row_blocks,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const std::string& family = "binomial",
    const std::string& link   = "logit",
    int Gridtype = 2,
    int n_envopt = -1,
    bool use_parallel = true,
    bool use_opencl = false,
    bool verbose = false
) {
  return glmbayes::sim::group::rNormalGLMGroups(
    n, y, x, offset, wt,
    dispersion,
    mu, P_blocks, prior_by_block, row_blocks,
    f2, f3,
    family, link, Gridtype,
    n_envopt, use_parallel, use_opencl, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List rNormalRegGroups_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    const Rcpp::NumericVector& dispersion,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::List& P_blocks,
    bool prior_by_block,
    const Rcpp::List& row_blocks,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    int Gridtype = 2
) {
  return glmbayes::sim::group::rNormalRegGroups(
    n, y, x, offset, wt,
    dispersion,
    mu, P_blocks, prior_by_block, row_blocks,
    f2, f3,
    Gridtype
  );
}

// [[Rcpp::export]]
Rcpp::List group_rNormalReg_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    SEXP group,
    SEXP prior_list,
    SEXP prior_lists,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    int Gridtype = 2
) {
  return glmbayes::sim::group::group_rNormalReg_cpp_export(
    n, y, x, group, prior_list, prior_lists, offset, wt, f2, f3, Gridtype
  );
}

// [[Rcpp::export]]
Rcpp::List group_rNormalGLM_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    SEXP group,
    SEXP prior_list,
    SEXP prior_lists,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const std::string& family = "binomial",
    const std::string& link   = "logit",
    int Gridtype = 2,
    int n_envopt = -1,
    bool use_parallel = true,
    bool use_opencl = false,
    bool verbose = false
) {
  return glmbayes::sim::group::group_rNormalGLM_cpp_export(
    n, y, x, group, prior_list, prior_lists, offset, wt, f2, f3,
    family, link, Gridtype, n_envopt, use_parallel, use_opencl, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List rIndepNormalGammaReg_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    double shape,
    double rate,
    double max_disp_perc,
    Rcpp::Nullable<Rcpp::NumericVector> disp_lower,
    Rcpp::Nullable<Rcpp::NumericVector> disp_upper,
    int Gridtype,
    int n_envopt,
    bool use_parallel,
    bool use_opencl,
    bool verbose,
    bool progbar
) {
  glmbayes::env::check_disp_bounds_or_stop(
    disp_lower, disp_upper, "rIndepNormalGammaReg_cpp_export"
  );
  return rIndepNormalGammaReg(
    n,
    y,
    x,
    mu,
    P,
    offset,
    wt,
    shape,
    rate,
    max_disp_perc,
    disp_lower,
    disp_upper,
    Gridtype,
    n_envopt,
    use_parallel,
    use_opencl,
    verbose,
    progbar
  );
}

// [[Rcpp::export]]
Rcpp::List rNormalGammaReg_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    double shape,
    double rate,
    Rcpp::Nullable<double> max_disp_perc,
    Rcpp::Nullable<double> disp_lower,
    Rcpp::Nullable<double> disp_upper,
    bool verbose = false
) {
  return glmbayes::sim::rNormalGammaReg(
    n,
    y,
    x,
    mu,
    P,
    offset,
    wt,
    shape,
    rate,
    max_disp_perc,
    disp_lower,
    disp_upper,
    verbose
  );
}

// [[Rcpp::export]]
Rcpp::List rGammaGaussian_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& beta,
    const Rcpp::NumericVector& wt,
    const Rcpp::NumericVector& alpha,
    double shape,
    double rate,
    Rcpp::Nullable<double> disp_lower = R_NilValue,
    Rcpp::Nullable<double> disp_upper = R_NilValue,
    bool verbose = false
) {
  return glmbayes::sim::rGammaGaussian(
    n, y, x, beta, wt, alpha,
    shape, rate,
    disp_lower, disp_upper,
    verbose
  );
}

// [[Rcpp::export]]
Rcpp::List rGammaGamma_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& beta,
    const Rcpp::NumericVector& wt,
    const Rcpp::NumericVector& alpha,
    double shape,
    double rate,
    double max_disp_perc,
    Rcpp::Nullable<double> disp_lower = R_NilValue,
    Rcpp::Nullable<double> disp_upper = R_NilValue,
    bool verbose = false
) {
  return glmbayes::sim::rGammaGamma(
    n, y, x, beta, wt, alpha,
    shape, rate, max_disp_perc,
    disp_lower, disp_upper,
    verbose
  );
}


// =============================================================================
// Tier 2: Envelope & Standardization
// Callers: EnvelopeSize, EnvelopeBuild, EnvelopeEval, EnvelopeDispersionBuild,
//          EnvelopeOrchestrator, rNormalGLM_std; EnvelopeSet_* are internal
// User:    Advanced users - understanding algorithm, custom envelope workflows
// =============================================================================

// [[Rcpp::export]]
Rcpp::List rNormalGLM_std_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    const Rcpp::Function& f2,
    const Rcpp::List& Envelope,
    const Rcpp::CharacterVector& family,
    const Rcpp::CharacterVector& link,
    int progbar = 1,
    bool verbose = false
) {
  return rNormalGLM_std(
    n, y, x, mu, P, alpha, wt,
    f2, Envelope, family, link,
    progbar, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List EnvelopeSize_cpp_export(
    const arma::vec& a,
    const Rcpp::NumericMatrix& G1,
    int Gridtype,
    int n,
    int n_envopt,
    bool use_opencl,
    bool verbose
) {
  return glmbayes::env::EnvelopeSize(
    a, G1, Gridtype, n, n_envopt, use_opencl, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List EnvelopeBuild_cpp_export(
    Rcpp::NumericVector bStar,
    Rcpp::NumericMatrix A,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericMatrix mu,
    Rcpp::NumericMatrix P,
    Rcpp::NumericVector alpha,
    Rcpp::NumericVector wt,
    std::string family,
    std::string link,
    int Gridtype,
    int n,
    int n_envopt,
    bool sortgrid,
    bool use_opencl,
    bool verbose
) {
  return glmbayes::env::EnvelopeBuild(
    bStar, A, y, x, mu, P, alpha, wt,
    family, link, Gridtype, n, n_envopt,
    sortgrid, use_opencl, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List EnvelopeBuild_Ind_Normal_Gamma_cpp_export(
    const Rcpp::NumericVector& bStar,
    const Rcpp::NumericMatrix& A,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    const std::string& family = "binomial",
    const std::string& link   = "logit",
    int Gridtype              = 2,
    int n                     = 1,
    int n_envopt              = -1,
    bool sortgrid             = false,
    bool use_opencl           = false,
    bool verbose              = false
) {
  return EnvelopeBuild_Ind_Normal_Gamma(
    bStar, A, y, x, mu, P, alpha, wt,
    family, link,
    Gridtype, n, n_envopt,
    sortgrid, use_opencl, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List EnvelopeEval_cpp_export(
    const Rcpp::NumericMatrix& G4,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    const std::string& family,
    const std::string& link,
    bool use_opencl = false,
    bool verbose = false
) {
  return EnvelopeEval(
    G4, y, x, mu, P, alpha, wt,
    family, link,
    use_opencl, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List EnvelopeDispersionBuild_cpp_export(
    const Rcpp::List& Env,
    double Shape,
    double Rate,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& alpha,
    int n_obs,
    double RSS_post,
    double RSS_ML,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::NumericVector& wt,
    double max_disp_perc = 0.99,
    Rcpp::Nullable<double> disp_lower = R_NilValue,
    Rcpp::Nullable<double> disp_upper = R_NilValue,
    bool verbose = false,
    bool use_parallel = true
) {
  return EnvelopeDispersionBuild(
    Env,
    Shape,
    Rate,
    P,
    y,
    x,
    alpha,
    n_obs,
    RSS_post,
    RSS_ML,
    mu,
    wt,
    max_disp_perc,
    disp_lower,
    disp_upper,
    verbose,
    use_parallel
  );
}

// [[Rcpp::export]]
Rcpp::List EnvelopeOrchestrator_cpp_export(
    const Rcpp::NumericVector& bstar2,
    const Rcpp::NumericMatrix& A,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x2,
    const Rcpp::NumericMatrix& mu2,
    const Rcpp::NumericMatrix& P2,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    int n,
    int Gridtype,
    Rcpp::Nullable<int> n_envopt,
    double shape,
    double rate,
    double RSS_Post2,
    double RSS_ML,
    double max_disp_perc,
    Rcpp::Nullable<double> disp_lower,
    Rcpp::Nullable<double> disp_upper,
    bool use_parallel,
    bool use_opencl,
    bool verbose
) {
  return EnvelopeOrchestrator(
    bstar2,
    A,
    y,
    x2,
    mu2,
    P2,
    alpha,
    wt,
    n,
    Gridtype,
    n_envopt,
    shape,
    rate,
    RSS_Post2,
    RSS_ML,
    max_disp_perc,
    disp_lower,
    disp_upper,
    use_parallel,
    use_opencl,
    verbose
  );
}

// [[Rcpp::export]]
Rcpp::List EnvelopeCentering_cpp_export(
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    double shape,
    double rate,
    int Gridtype = 2,
    bool verbose = false
) {
  return glmbayes::env::EnvelopeCentering(
    y, x, mu, P, offset, wt,
    shape, rate,
    Gridtype, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List EnvelopeSet_Grid_cpp_export(
    const Rcpp::NumericMatrix& GIndex,
    const Rcpp::NumericMatrix& cbars,
    const Rcpp::NumericMatrix& Lint
) {
  return EnvelopeSet_Grid(
    GIndex,
    cbars,
    Lint
  );
}

// [[Rcpp::export]]
Rcpp::List EnvelopeSet_LogP_cpp_export(
    const Rcpp::NumericMatrix& logP,
    const Rcpp::NumericVector& NegLL,
    const Rcpp::NumericMatrix& cbars,
    const Rcpp::NumericMatrix& G3
) {
  return EnvelopeSet_LogP(
    logP,
    NegLL,
    cbars,
    G3
  );
}


// =============================================================================
// Tier 3: Standardized Samplers (Indep Normal-Gamma)
// Callers: C++ only (rIndepNormalGammaReg); R wrappers for custom split workflow
// User:    Advanced / developers - after EnvelopeOrchestrator, sample separately
// =============================================================================

// [[Rcpp::export]]
Rcpp::List rIndepNormalGammaReg_std_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    const Rcpp::Function& f2,
    const Rcpp::List& Envelope,
    const Rcpp::List& gamma_list,
    const Rcpp::List& UB_list,
    const Rcpp::CharacterVector& family,
    const Rcpp::CharacterVector& link,
    bool progbar = true,
    bool verbose = false
) {
  return rIndepNormalGammaReg_std(
    n, y, x, mu, P, alpha, wt,
    f2, Envelope, gamma_list, UB_list,
    family, link, progbar, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List rIndepNormalGammaReg_std_parallel_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& mu,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    const Rcpp::Function& f2,
    const Rcpp::List& Envelope,
    const Rcpp::List& gamma_list,
    const Rcpp::List& UB_list,
    const Rcpp::CharacterVector& family,
    const Rcpp::CharacterVector& link,
    bool progbar = true,
    bool verbose = false
) {
  return rIndepNormalGammaReg_std_parallel(
    n, y, x, mu, P, alpha, wt,
    f2, Envelope, gamma_list, UB_list,
    family, link, progbar, verbose
  );
}


// =============================================================================
// Tier 4: Model Utilities
// Callers: glmb_Standardize_Model
// User:    Advanced users - model preparation, standardization
// =============================================================================

// [[Rcpp::export]]
Rcpp::List glmb_Standardize_Model_cpp_export(
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericMatrix& P,
    const Rcpp::NumericMatrix& bstar,
    const Rcpp::NumericMatrix& A1
) {
  return glmb_Standardize_Model(
    y, x, P, bstar, A1
  );
}


// =============================================================================
// Tier 5: OpenCL / GPU
// Callers: glmbayesCore_has_opencl, gpu_names
// User:    Advanced users - GPU diagnostics for use_opencl
// Kernel loading: opencltools::load_kernel_source / load_kernel_library
// Core count:     opencltools::get_opencl_core_count
// =============================================================================

// [[Rcpp::export]]
bool glmbayesCore_has_opencl_cpp_export() {
  return glmbayesCore_has_opencl();
}

// [[Rcpp::export]]
Rcpp::CharacterVector gpu_names_cpp_export() {
  return gpu_names();
}


// =============================================================================
// Phased Out (no R wrappers; C++ exports commented out)
// - rss_face_at_disp, UB2: former RSS/UB2 minimization; active path uses
//   closed-form C++ bounds.
//
// To fully remove: delete this block, then (1) remove *.o from src/,
// (2) uninstall old glmbayes, (3) Rcpp::compileAttributes(),
// (4) devtools::document(), (5) devtools::install().
// =============================================================================

/*
// [[Rcpp::export]]
double rss_face_at_disp_cpp_export(
    double dispersion,
    const Rcpp::List& cache,
    const Rcpp::NumericVector& cbars_j,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt
) {
  return rss_face_at_disp(
    dispersion,
    cache,
    cbars_j,
    y,
    x,
    alpha,
    wt
  );
}

// [[Rcpp::export]]
double UB2_cpp_export(
    double dispersion,
    const Rcpp::List& cache,
    const Rcpp::NumericVector& cbars_j,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    double rss_min_global
) {
  return UB2(
    dispersion,
    cache,
    cbars_j,
    y,
    x,
    alpha,
    wt,
    rss_min_global
  );
}
*/


// =============================================================================
// Two-block Gibbs and block-ING entry points (lmerb engines).
// =============================================================================

// lmerb-unused: candidate for deletion (not on the lmerb path).
// [[Rcpp::export]]
Rcpp::List two_block_rNormal_reg_cpp_export(
    int n,
    int m_convergence,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    SEXP block,
    const Rcpp::List& x_hyper,
    const Rcpp::List& prior_list_block1,
    SEXP dispersion_block1,
    SEXP ddef_block1,
    const Rcpp::List& pfamily_list,
    const Rcpp::List& fixef_start,
    const Rcpp::CharacterVector& group_levels,
    const std::string& family,
    const std::string& link,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::Function& f2_gauss,
    const Rcpp::Function& f3_gauss,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    int Gridtype = 2,
    int n_envopt = 1,
    bool use_parallel = true,
    bool use_opencl = false,
    bool verbose = false,
    bool progbar = true
) {
  return glmbayes::sim::twoblock::two_block_rNormal_reg_cpp_export(
    n, m_convergence, y, x, block, x_hyper,
    prior_list_block1, dispersion_block1, ddef_block1,
    pfamily_list, fixef_start, group_levels,
    family, link, f2, f3, f2_gauss, f3_gauss,
    offset, wt, Gridtype, n_envopt,
    use_parallel, use_opencl, verbose, progbar
  );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix two_block_build_mu_all_cpp_export(
    const Rcpp::List& x_hyper,
    const Rcpp::List& fixef,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::CharacterVector& group_levels
) {
  return glmbayes::sim::twoblock::two_block_build_mu_all(
    x_hyper, fixef, re_names, group_levels
  );
}

// lmerb-unused: candidate for deletion (not on the lmerb path).
// [[Rcpp::export]]
Rcpp::List two_block_groups_prior_with_tau2_cpp_export(
    const Rcpp::List& base_prior,
    const Rcpp::NumericVector& tau2_vec,
    const Rcpp::CharacterVector& ptypes,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::NumericMatrix& mu_all
) {
  return glmbayes::sim::twoblock::two_block_groups_prior_with_tau2(
    base_prior, tau2_vec, ptypes, re_names, mu_all
  );
}

// lmerb-unused: candidate for deletion (not on the lmerb path).
// [[Rcpp::export]]
double two_block_groups_iters_mean_cpp_export(const Rcpp::List& block_out) {
  return glmbayes::sim::twoblock::two_block_groups_iters_mean(block_out);
}

// lmerb-unused: candidate for deletion (not on the lmerb path).
// [[Rcpp::export]]
Rcpp::List two_block_batch_fixef_chain_cpp_export(
    const Rcpp::List& batch_fixef,
    int chain_i,
    const Rcpp::CharacterVector& re_names
) {
  return glmbayes::sim::twoblock::fixef_list_from_batch_chain(
    batch_fixef, chain_i, re_names
  );
}

// lmerb-unused: candidate for deletion (not on the lmerb path).
// [[Rcpp::export]]
Rcpp::NumericVector two_block_batch_tau2_chain_row_cpp_export(
    const Rcpp::NumericMatrix& batch_tau2,
    int chain_i
) {
  return glmbayes::sim::twoblock::batch_tau2_chain_row(batch_tau2, chain_i);
}

// lmerb-unused: candidate for deletion (not on the lmerb path).
// [[Rcpp::export]]
Rcpp::NumericVector two_block_batch_b_assign_slice_cpp_export(
    Rcpp::NumericVector b_store,
    int chain_i,
    Rcpp::NumericMatrix b_draw
) {
  Rcpp::NumericVector out = Rcpp::clone(b_store);
  glmbayes::sim::twoblock::batch_b_assign_slice(out, chain_i, b_draw);
  return out;
}

// lmerb-unused: candidate for deletion (not on the lmerb path).
// [[Rcpp::export]]
Rcpp::NumericVector two_block_batch_iters_ranef_add_cpp_export(
    Rcpp::NumericVector iters_ranef,
    int chain_i,
    double iters_mean
) {
  Rcpp::NumericVector out = Rcpp::clone(iters_ranef);
  glmbayes::sim::twoblock::batch_iters_ranef_add(out, chain_i, iters_mean);
  return out;
}

// lmerb-unused: candidate for deletion (not on the lmerb path).
// [[Rcpp::export]]
Rcpp::List two_block_groups_one_chain_draw_cpp_export(
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
) {
  return glmbayes::sim::twoblock::two_block_groups_one_chain_impl(
    chain_i, batch_fixef, tau2_i, y, Z, groups, offset, wt,
    x_hyper, re_names, group_levels, ptypes, block1_prior,
    is_gaussian, f2, f3, f2_gauss, f3_gauss, family, link,
    Gridtype, n_envopt
  );
}

// lmerb-unused: candidate for deletion (not on the lmerb path).
// [[Rcpp::export]]
Rcpp::List two_block_groups_one_chain_cpp_export(
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
) {
  return glmbayes::sim::twoblock::two_block_groups_one_chain_orchestrate_impl(
    chain_i, b_store, iters_ranef, batch_fixef, batch_tau2, design,
    block1_prior, family, ptypes, re_names, group_levels,
    f2, f3, f2_gauss, f3_gauss,
    use_cpp_tau2_row, use_cpp_b_slice, use_cpp_iters_ranef_add
  );
}

// lmerb-unused: candidate for deletion (not on the lmerb path).
// [[Rcpp::export]]
Rcpp::List two_block_groups_one_chain_from_mu_P_cpp_export(
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
) {
  return glmbayes::sim::twoblock::two_block_groups_one_chain_from_mu_P_impl(
    mu_all, P, dispersion, ddef, design, family, re_names, group_levels,
    f2, f3, f2_gauss, f3_gauss
  );
}


// lmerb-unused: candidate for deletion (not on the lmerb path).
// [[Rcpp::export]]
Rcpp::NumericMatrix two_block_reorder_b_to_group_levels_cpp_export(
    Rcpp::NumericMatrix b_draw,
    SEXP block_ids,
    Rcpp::CharacterVector group_levels
) {
  return glmbayes::sim::twoblock::two_block_reorder_b_to_group_levels(
    b_draw, block_ids, group_levels
  );
}

// lmerb-unused export: candidate for deletion. The implementation stays;
// lmerb reaches it through another entry point.
// [[Rcpp::export]]
Rcpp::NumericVector two_block_align_b_to_xhyper_cpp_export(
    Rcpp::NumericVector b_vec,
    Rcpp::NumericMatrix X_k,
    Rcpp::CharacterVector group_levels
) {
  return glmbayes::sim::twoblock::two_block_align_b_to_xhyper_cpp(
    b_vec, X_k, group_levels
  );
}

// [[Rcpp::export]]
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
) {
  return glmbayes::sim::twoblock::two_block_population_one_chain_cpp_export(
    b_i, fixef_rows, tau2_i, iters_i, x_hyper, group_levels,
    pfamily_list, ptypes, re_names
  );
}

// lmerb-unused export: candidate for deletion. The implementation stays;
// lmerb reaches it through another entry point.
// [[Rcpp::export]]
Rcpp::List BlockEnvelopeCentering_cpp_export(
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    SEXP block,
    SEXP prior_list,
    SEXP prior_lists,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    double shape,
    double rate,
    double max_disp_perc,
    Rcpp::Nullable<double> disp_lower,
    Rcpp::Nullable<double> disp_upper,
    int p_re = -1,
    int n_rss_iter = 10,
    bool verbose = false
) {
  return glmbayes::env::group::BlockEnvelopeCentering(
    y, x, block, prior_list, prior_lists,
    offset, wt, shape, rate, max_disp_perc,
    disp_lower, disp_upper, p_re, n_rss_iter, verbose
  );
}

// lmerb-unused export: candidate for deletion. The implementation stays;
// lmerb reaches it through another entry point.
// [[Rcpp::export]]
Rcpp::List BlockEnvelopeBuild_cpp_export(
    const Rcpp::List& centering_out,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    SEXP block,
    SEXP prior_list,
    SEXP prior_lists,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    double max_disp_perc,
    Rcpp::Nullable<double> disp_lower,
    Rcpp::Nullable<double> disp_upper,
    int n = 1,
    int Gridtype = 3,
    int n_envopt = -1,
    double RSS_ML = NA_REAL,
    bool use_parallel = true,
    bool use_opencl = false,
    bool verbose = false
) {
  Rcpp::Nullable<int> n_envopt_nullable(
    n_envopt < 0 ? R_NilValue : Rcpp::wrap(n_envopt)
  );
  return glmbayes::env::group::BlockEnvelopeBuild(
    centering_out, y, x, block, prior_list, prior_lists,
    offset, wt, max_disp_perc, disp_lower, disp_upper,
    n, Gridtype, n_envopt_nullable, RSS_ML,
    use_parallel, use_opencl, verbose
  );
}

// lmerb-unused: candidate for deletion (not on the lmerb path).
// [[Rcpp::export]]
Rcpp::List BlockEnvelopeDispersionBuild_cpp_export(
    const Rcpp::List& build_out,
    const Rcpp::List& centering_out,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    SEXP block,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    double shape,
    double rate,
    double max_disp_perc,
    Rcpp::Nullable<double> disp_lower,
    Rcpp::Nullable<double> disp_upper,
    double RSS_ML = NA_REAL,
    bool use_parallel = true,
    bool verbose = false
) {
  return glmbayes::env::group::BlockEnvelopeDispersionBuild(
    build_out, centering_out, y, x, block,
    offset, wt, shape, rate, max_disp_perc,
    disp_lower, disp_upper, RSS_ML,
    use_parallel, verbose
  );
}

// lmerb-unused: candidate for deletion (not on the lmerb path).
// [[Rcpp::export]]
Rcpp::List BlockEnvelopeSim_cpp_export(
    const Rcpp::List& build_out,
    int n = 1,
    bool progbar = false,
    bool verbose = false
) {
  return glmbayes::env::group::BlockEnvelopeSim(
    build_out, n, progbar, verbose
  );
}

// lmerb-unused: candidate for deletion (not on the lmerb path).
// [[Rcpp::export]]
Rcpp::List rIndepNormalGammaRegBlock_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    SEXP block,
    SEXP prior_list,
    SEXP prior_lists,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    int p_re,
    int n_rss_iter,
    int Gridtype,
    int n_envopt,
    double RSS_ML,
    bool use_parallel,
    bool use_opencl,
    bool progbar,
    bool verbose,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::CharacterVector& re_names
) {
  Rcpp::Nullable<int> n_envopt_nullable(
    n_envopt < 0 ? R_NilValue : Rcpp::wrap(n_envopt)
  );
  return glmbayes::env::group::rIndepNormalGammaRegBlock(
    n, y, x, block, prior_list, prior_lists,
    offset, wt, p_re, n_rss_iter, Gridtype, n_envopt_nullable,
    RSS_ML, use_parallel, use_opencl, progbar, verbose,
    group_levels, re_names
  );
}

// lmerb-unused export: candidate for deletion. The implementation stays;
// lmerb reaches it through another entry point.
// [[Rcpp::export]]
Rcpp::List BlockEnvelopeDispersionBuildInd_cpp_export(
    const Rcpp::List& build_out,
    const Rcpp::List& centering_out,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    SEXP block,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    double shape,
    double rate,
    double max_disp_perc,
    Rcpp::Nullable<double> disp_lower,
    Rcpp::Nullable<double> disp_upper,
    double RSS_ML = NA_REAL,
    bool use_parallel = true,
    bool verbose = false
) {
  return glmbayes::env::group::BlockEnvelopeDispersionBuildInd(
    build_out, centering_out, y, x, block,
    offset, wt, shape, rate, max_disp_perc,
    disp_lower, disp_upper, RSS_ML,
    use_parallel, verbose
  );
}

// lmerb-unused export: candidate for deletion. The implementation stays;
// lmerb reaches it through another entry point.
// [[Rcpp::export]]
Rcpp::List BlockEnvelopeSimInd_cpp_export(
    const Rcpp::List& build_out,
    int n = 1,
    bool progbar = false,
    bool verbose = false
) {
  return glmbayes::env::group::BlockEnvelopeSimInd(
    build_out, n, progbar, verbose
  );
}

// [[Rcpp::export]]
Rcpp::List rIndepNormalGammaRegBlockInd_cpp_export(
    int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    SEXP block,
    SEXP prior_list,
    SEXP prior_lists,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    int p_re,
    int n_rss_iter,
    int Gridtype,
    int n_envopt,
    double RSS_ML,
    bool use_parallel,
    bool use_opencl,
    bool progbar,
    bool verbose,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::CharacterVector& re_names
) {
  Rcpp::Nullable<int> n_envopt_nullable(
    n_envopt < 0 ? R_NilValue : Rcpp::wrap(n_envopt)
  );
  return glmbayes::env::group::rIndepNormalGammaRegBlockInd(
    n, y, x, block, prior_list, prior_lists,
    offset, wt, p_re, n_rss_iter, Gridtype, n_envopt_nullable,
    RSS_ML, use_parallel, use_opencl, progbar, verbose,
    group_levels, re_names
  );
}
