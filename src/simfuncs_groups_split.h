// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

/**
 * @file simfuncs_groups_split.h
 * @brief One group's rNormalGLM step split into mode, envelope and simulation.
 *
 * @namespace glmbayes::sim::group
 *
 * @section Phases
 *   1. group_mode:     shift the prior mean into the offset, optim (BFGS)
 *                      or vmmin_mode, glmb_Standardize_Model.
 *   2. group_envelope: EnvelopeGridSetup, EnvelopeEval, EnvelopeFinalize
 *                      (EnvelopeBuild split around the evaluation).
 *   3. group_simulate: rNormalGLM_std (or _parallel), undo the
 *                      standardization, add the prior mean back.
 *   Phases 1 and 2 reproduce rNormalGLM step for step, so the mode and
 *   envelope equal rNormalGLM's for the same inputs.
 *
 * @section ImplementedIn
 *   - rNormalGLMGroups_split.cpp
 *   - group_mode_vmmin.cpp (vmmin_mode, mode_compare_update)
 *   - group_mode_hessian.cpp (exact_hessian)
 *
 * @section UsedBy
 *   - rNormalGLMGroups.cpp
 *   - glmbayesCore_ccallables.cpp (rNormalGLMGroupChains, C API)
 */

#ifndef GLMBAYES_SIM_GROUPS_SPLIT_H
#define GLMBAYES_SIM_GROUPS_SPLIT_H

#include "RcppArmadillo.h"

#include <string>
#include <vector>

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

/// Mode finder used by group_mode.
enum ModeMethod {
    MODE_OPTIM = 0,    // R optim (BFGS) with the R f2/f3 and optimhess
    MODE_VMMIN = 1,    // vmmin with the C++ f2/f3 and an optimhess replica
    MODE_COMPARE = 2   // optim kept downstream; both vmmin variants compared
};

/// Diagnostics accumulated over group_mode calls with MODE_COMPARE.
/// exact_mismatch counts calls where vmmin with the R f2/f3 differs from
/// optim in any bit of par, value or hessian. flagged counts calls where
/// vmmin with the C++ f2/f3 exceeds a tolerance (see group_mode_vmmin.cpp).
/// max_d: mode difference in posterior-SD units; max_df: |f difference| /
/// (1 + |f|); max_hess: relative Frobenius Hessian difference. Times in
/// seconds; fn/gr are evaluation counts.
struct ModeCompare {
    double calls = 0.0;
    double exact_mismatch = 0.0;
    double flagged = 0.0;
    double max_d = 0.0;
    double max_df = 0.0;
    double max_hess = 0.0;
    double t_optim = 0.0;
    double t_vmmin_r = 0.0;
    double t_vmmin_cpp = 0.0;
    double fn_optim = 0.0;
    double gr_optim = 0.0;
    double fn_vmmin = 0.0;
    double gr_vmmin = 0.0;
};

/// Starting point of vmmin in group_mode (vmmin path only).
enum WarmMethod {
    WARM_COLD = 0,     // par = 0 (the prior mean)
    WARM_START = 1,    // par = start - mu; cold rerun if that fails
    WARM_COMPARE = 2   // both; the cold result is kept downstream
};

/// Diagnostics accumulated over group_mode calls with a warm start.
/// calls and flagged are counted for WARM_COMPARE only (flagged uses the
/// MODE_COMPARE tolerances); fallback counts warm runs redone cold. max_d is
/// the warm-cold mode difference in the cold A1 metric, max_df the
/// |f difference| / (1 + |f_cold|), max_hess the relative Frobenius Hessian
/// difference. fn/gr are evaluation counts, times in seconds.
struct WarmCompare {
    double calls = 0.0;
    double flagged = 0.0;
    double fallback = 0.0;
    double max_d = 0.0;
    double max_df = 0.0;
    double max_hess = 0.0;
    double fn_cold = 0.0;
    double gr_cold = 0.0;
    double fn_warm = 0.0;
    double gr_warm = 0.0;
    double t_cold = 0.0;
    double t_warm = 0.0;
};

/// Hessian at the vmmin mode (vmmin path only).
enum HessMethod {
    HESS_NUMERIC = 0,  // optimhess replica (finite differences of f3)
    HESS_EXACT = 1,    // X' diag(w) X + P; numeric if not finite or not PD
    HESS_COMPARE = 2   // both; the numeric Hessian is kept downstream
};

/// Diagnostics accumulated over vmmin_mode calls with HESS_EXACT or
/// HESS_COMPARE. calls, flagged, max_hess and the times are counted for
/// HESS_COMPARE only; fallback counts exact Hessians that were not finite or
/// not positive definite (the numeric one is used instead). max_hess is the
/// relative Frobenius difference to the numeric Hessian; times in seconds.
struct HessCompare {
    double calls = 0.0;
    double flagged = 0.0;
    double fallback = 0.0;
    double max_hess = 0.0;
    double t_numeric = 0.0;
    double t_exact = 0.0;
};

/// Seconds in the parts of group_mode (accumulated); calls (group_mode
/// calls) and kernel_calls are counts. vmmin and hessian include their f2/f3 callbacks; the *_kernel
/// entries are the f2_f3_non_opencl calls inside them.
struct ModeTiming {
    double prep = 0.0;            // shift: alpha, clones, R lookups
    double optim = 0.0;           // R optim incl. optimhess
    double vmmin = 0.0;           // vmmin iterations
    double vmmin_kernel = 0.0;
    double hessian = 0.0;         // optimhess replica and / or exact Hessian
    double hessian_kernel = 0.0;
    double fit_list = 0.0;        // b2a, A1 and the optim-shaped fit list
    double compare = 0.0;         // MODE_COMPARE's extra vmmin runs
    double standardize = 0.0;     // glmb_Standardize_Model and unpacking
    double calls = 0.0;
    double kernel_calls = 0.0;
};

/// vmmin (R's C BFGS, as used by optim) from par = 0 (or `par0`, l1 values)
/// with optim's default settings, then the Hessian chosen by `hess_method`
/// (optimhess finite differences by default), on the shifted
/// model held in `g` (y, x, P, alpha, wt2; prior mean 0). `r_callbacks`
/// evaluates the R f2/f3 (same values as optim); otherwise the C++ kernels
/// for family/link. `fit` has the fields of an optim result.
struct VmminResult {
    Rcpp::NumericMatrix b2a;   // l1 x 1
    Rcpp::NumericMatrix A1;    // l1 x l1
    double value;
    int convergence;
    int fncount;
    int grcount;
    Rcpp::List fit;
};

/// Family/link pairs with C++ f2/f3 kernels (f2_f3_non_opencl). group_mode
/// uses optim for any other pair whatever the mode_method.
bool vmmin_supported(const std::string& family, const std::string& link);

VmminResult vmmin_mode(
    const GroupStd& g,
    Rcpp::Function f2,
    Rcpp::Function f3,
    const std::string& family,
    const std::string& link,
    bool r_callbacks,
    ModeTiming* timing = nullptr,
    const double* par0 = nullptr,
    int hess_method = HESS_NUMERIC,
    HessCompare* hc = nullptr
);

/// Exact Hessian X' diag(w) X + P of the shifted model's negative log
/// posterior at `par` (l1 values), symmetrized. w are the second derivatives
/// of the f2 kernels in eta = alpha + x par, times wt2 (observed information
/// for probit and cloglog). Returns false, leaving `H` unspecified, for an
/// unsupported family/link, a non-finite weight, or a Hessian whose Cholesky
/// factorization fails.
bool exact_hessian(
    const GroupStd& g,
    const std::vector<double>& par,
    const std::string& family,
    const std::string& link,
    arma::mat& H
);

/// Fold one WARM_COMPARE call into `wc` (`cold` and `warm` from vmmin_mode,
/// `t_cold`/`t_warm` their seconds).
void warm_compare_update(
    WarmCompare& wc,
    const VmminResult& cold,
    const VmminResult& warm,
    double t_cold,
    double t_warm
);

/// Fold one MODE_COMPARE call into `cmp`: `opt` is optim's result list,
/// `vr`/`vc` the vmmin results with R and C++ callbacks.
void mode_compare_update(
    ModeCompare& cmp,
    const GroupStd& g,
    const Rcpp::List& opt,
    const VmminResult& vr,
    const VmminResult& vc,
    Rcpp::Function f2
);

/// `link` is only needed for MODE_VMMIN / MODE_COMPARE. `cmp` receives the
/// MODE_COMPARE diagnostics, `mt` the time per part. With MODE_VMMIN,
/// `warm_start` (l1 values in the coordinates of `mu`, or nullptr for cold) and
/// `warm_method` choose vmmin's starting point; `wc` receives the
/// WarmCompare diagnostics; `hess_method` and `hc` are passed to
/// vmmin_mode. The optim path ignores them.
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
    bool verbose = false,
    const std::string& link = "",
    int mode_method = MODE_OPTIM,
    ModeCompare* cmp = nullptr,
    ModeTiming* mt = nullptr,
    const double* warm_start = nullptr,
    int warm_method = WARM_COLD,
    WarmCompare* wc = nullptr,
    int hess_method = HESS_NUMERIC,
    HessCompare* hc = nullptr
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
/// Returns coefficients and coef.mode (n_chains x l1), iters (n_chains), and
/// timing: seconds in mode, setup, eval, finalize, simulate (check excluded),
/// setup_timing: setup split into the GridSetupTiming fields,
/// mode_timing: mode split into the ModeTiming fields,
/// setup_compare: calls compared with setup_method SETUP_COMPARE (any
/// difference stops the call), and mode_compare: the ModeCompare fields (all zero unless mode_method is
/// MODE_COMPARE). `check` is ignored for MODE_VMMIN.
/// `start` is R_NilValue or an l1 x n_chains matrix of starting modes in
/// the coordinates of `mu` (a column with an NA starts cold); `warm_method`
/// is a WarmMethod and applies to MODE_VMMIN only. Also returned:
/// mode_hessian (l1*l1 x n_chains; column c is chain c's A1, column-major)
/// and warm_compare (the WarmCompare fields). `hess_method` is a HessMethod
/// (MODE_VMMIN only); hessian_compare returns the HessCompare fields.
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
    bool verbose = false,
    int mode_method = MODE_OPTIM,
    int setup_method = 0,
    SEXP start = R_NilValue,
    int warm_method = WARM_COLD,
    int hess_method = HESS_NUMERIC
);

} // namespace group

} // namespace sim

} // namespace glmbayes

#endif
