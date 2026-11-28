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
 *   - group_mode_newton.cpp (newton_minimize, newton_mode,
 *     newton_compare_update)
 *   - group_simulate_parallel.cpp (group_draw_in, group_draw_one,
 *     group_simulate_chains, sim_compare_update)
 *   - envelope_eval_batch.cpp (group_eval_method, eval_batch_pack,
 *     poisson_eval_batch)
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

namespace fam {
// As in famfuncs.h, which is not included here because of its
// namespace-scope using directives.
typedef double (*F2ValueFn)(int l1, int l2, const double* x, const double* y,
                            const double* alpha, const double* wt,
                            const double* mu, const double* P, const double* b,
                            double* work);
}  // namespace fam

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
    MODE_COMPARE = 2,  // optim kept downstream; both vmmin variants compared
    MODE_NEWTON = 3,   // Newton on f2_f3_non_opencl + exact_hessian
    MODE_NEWTON_COMPARE = 4,     // vmmin kept downstream; MODE_NEWTON compared
    MODE_NEWTON_F4 = 5,          // Newton on fam::f2_f3_f4_*
    MODE_NEWTON_F4_COMPARE = 6   // vmmin kept downstream; MODE_NEWTON_F4 compared
};

/// Evaluator behind newton_mode.
enum NewtonBackend {
    NEWTON_ENV = 0,    // env f2_f3_non_opencl for f and grad, exact_hessian for H
    NEWTON_FAM = 1     // fam::f2_f3_f4_* (allocation-free)
};

/// Diagnostics accumulated over group_mode calls with a Newton mode_method.
/// cold_retry counts warm-started Newton runs that failed and were redone
/// cold; fallback counts calls where the cold run also failed (vmmin is
/// used instead). The rest is counted for the compare modes only: calls,
/// flagged (the MODE_COMPARE tolerances in group_mode_vmmin.cpp, or a
/// failed Newton run), newton_lower (Newton reached the lower f2), max_d
/// (mode difference in the vmmin A1 metric), max_df (|f difference| /
/// (1 + |f_vmmin|), f2 from f2_f3_non_opencl at both points), max_hess
/// (relative Frobenius Hessian difference), the Newton iteration, value,
/// gradient and backtrack counts, and seconds in vmmin and Newton. With
/// MODE_NEWTON_F4_COMPARE, max_kf, max_kgrad and max_kH compare
/// fam::f2_f3_f4 at the Newton mode with f2_f3_non_opencl (relative value;
/// gradient in the H^-1 metric) and exact_hessian (relative Frobenius).
struct NewtonCompare {
    double calls = 0.0;
    double flagged = 0.0;
    double cold_retry = 0.0;
    double fallback = 0.0;
    double newton_lower = 0.0;
    double max_d = 0.0;
    double max_df = 0.0;
    double max_hess = 0.0;
    double max_kf = 0.0;
    double max_kgrad = 0.0;
    double max_kH = 0.0;
    double iters = 0.0;
    double fevals = 0.0;
    double gevals = 0.0;
    double backtracks = 0.0;
    double t_vmmin = 0.0;
    double t_newton = 0.0;
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
    double newton = 0.0;          // newton_mode (MODE_NEWTON, MODE_NEWTON_F4)
    double newton_kernel = 0.0;   // its evaluator calls
    double newton_iters = 0.0;
    double newton_backtracks = 0.0;
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

/// Evaluator for newton_minimize: f at b, plus the gradient (n values) and
/// Hessian (n x n, column-major) when `grad` / `H` are non-null. Returns
/// false when the evaluation failed.
struct NewtonEval {
    bool (*fn)(void* ctx, const double* b, double* f, double* grad, double* H);
    void* ctx;
};

struct NewtonOpts {
    double tol = 1e-10;       // stop when the Newton decrement g' H^-1 g / 2 <= tol
    int maxit = 50;
    int max_halvings = 40;    // backtracking steps per iteration
    double armijo = 1e-4;
};

struct NewtonStats {
    int iters = 0;
    int fevals = 0;           // all evaluator calls
    int gevals = 0;           // calls with gradient and Hessian
    int backtracks = 0;
};

enum NewtonStatus {
    NEWTON_OK = 0,
    NEWTON_MAXIT = 1,
    NEWTON_LINESEARCH = 2,
    NEWTON_CHOL = 3,
    NEWTON_EVAL = 4           // evaluation failed or not finite
};

/// newton_minimize workspace in doubles for n parameters.
inline int newton_work_size(int n) { return 3 * n + n * n; }

/// Damped Newton with Armijo backtracking from `b` (n values, overwritten
/// with the minimizer). On NEWTON_OK, `H` is the Hessian at the returned
/// `b` and `f_out` its value. Plain C style (no allocation, no exceptions;
/// `work` holds newton_work_size(n) doubles) so it can be ported to a GPU
/// work-item.
int newton_minimize(int n, const NewtonEval& ev, double* b, double* H,
                    double* f_out, double* work, const NewtonOpts& opts,
                    NewtonStats& stats);

/// Newton mode of the shifted model in `g`, from par = 0 (or `par0`), with
/// the evaluator chosen by `backend` (NewtonBackend). b2a, A1 and an
/// optim-shaped fit (counts are evaluator calls with and without the
/// gradient; convergence is the NewtonStatus) as in VmminResult. `work` is
/// resized as needed and can be reused across calls.
struct NewtonResult {
    Rcpp::NumericMatrix b2a;   // l1 x 1
    Rcpp::NumericMatrix A1;    // l1 x l1
    double value = 0.0;
    int status = NEWTON_EVAL;
    NewtonStats stats;
    Rcpp::List fit;
};

NewtonResult newton_mode(
    const GroupStd& g,
    const std::string& family,
    const std::string& link,
    const double* par0,
    int backend,
    std::vector<double>& work,
    ModeTiming* timing = nullptr
);

/// Fold one Newton compare call into `nc`: `vm` is the vmmin result kept
/// downstream, `nt` the Newton result, `t_vm`/`t_nt` their seconds.
/// `kernel_checks` adds the fam::f2_f3_f4 checks (max_kf, max_kgrad, max_kH).
void newton_compare_update(
    NewtonCompare& nc,
    const GroupStd& g,
    const std::string& family,
    const std::string& link,
    const VmminResult& vm,
    const NewtonResult& nt,
    double t_vm,
    double t_nt,
    bool kernel_checks
);

/// `link` is only needed for MODE_VMMIN / MODE_COMPARE. `cmp` receives the
/// MODE_COMPARE diagnostics, `mt` the time per part. With MODE_VMMIN,
/// `warm_start` (l1 values in the coordinates of `mu`, or nullptr for cold) and
/// `warm_method` choose vmmin's starting point; `wc` receives the
/// WarmCompare diagnostics; `hess_method` and `hc` are passed to
/// vmmin_mode. The optim path ignores them. The Newton modes start from the
/// same point as vmmin (warm when a warm start is given and `warm_method`
/// is not WARM_COLD); a failed warm run is redone cold, then vmmin is used.
/// `ncmp` receives the NewtonCompare diagnostics and `nwork` is Newton's
/// reusable workspace (a local one when nullptr).
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
    HessCompare* hc = nullptr,
    NewtonCompare* ncmp = nullptr,
    std::vector<double>* nwork = nullptr
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

/// Simulation path in rNormalGLMGroupChains, from
/// getOption("glmbayesCore.group_sim").
enum SimMethod {
    SIM_SERIAL = 0,    // group_simulate per chain (R's RNG)
    SIM_PARALLEL = 1,  // group_simulate_chains (thread-local RNG)
    SIM_COMPARE = 2    // serial draws kept; kernels, transform and timing compared
};

/// Raw views of one chain's standardized model and envelope for
/// group_draw_one. The pointers refer to the GroupStd and envelope List the
/// view was built from, which must outlive it. Matrices are column-major.
struct GroupDrawIn {
    int l1 = 0;               // coefficients
    int nobs = 0;             // observations
    int n_cells = 0;          // envelope cells
    const double* y = nullptr;
    const double* x2 = nullptr;      // nobs x l1
    const double* alpha = nullptr;
    const double* wt2 = nullptr;
    const double* mu2 = nullptr;     // l1
    const double* P2 = nullptr;      // l1 x l1
    const double* L2Inv = nullptr;   // l1 x l1
    const double* L3Inv = nullptr;   // l1 x l1
    const double* mu = nullptr;      // prior mean, l1
    const double* PLSD = nullptr;    // n_cells
    const double* LLconst = nullptr; // n_cells
    const double* loglt = nullptr;   // n_cells x l1
    const double* logrt = nullptr;   // n_cells x l1
    const double* cbars = nullptr;   // n_cells x l1
};

/// Views into `g` and `envelope` (an EnvelopeFinalize result). Uses the R
/// API, so call it outside parallel regions.
GroupDrawIn group_draw_in(const GroupStd& g, const Rcpp::List& envelope);

/// group_draw_one workspace in doubles.
int group_draw_work_size(int l1, int nobs);

/// fam::f2_value_select(family, link): nullptr for pairs without one.
glmbayes::fam::F2ValueFn group_f2_value(const std::string& family,
                                        const std::string& link);

/// SimMethod from getOption("glmbayesCore.group_sim"): "serial",
/// "parallel" (or unset) or "compare". Uses the R API.
int group_sim_method();

/// One accept-reject draw as in rNormalGLM_std (slice from PLSD, truncated
/// normals per coordinate, test LLconst + b'cbar - log U2 - f2(b) >= 0),
/// then L2Inv L3Inv b + mu into `coef_out` (l1). Draws from the thread-local
/// RNG (runif_safe), not R's. No R API, no allocation; `work` holds
/// group_draw_work_size doubles. Returns the number of candidates.
int group_draw_one(const GroupDrawIn& in, glmbayes::fam::F2ValueFn f2v,
                   double* work, double* coef_out);

/// group_draw_one for every chain in parallel (RcppParallel::parallelFor).
/// coef is n_chains x l1 (column-major), iters n_chains.
void group_simulate_chains(const std::vector<GroupDrawIn>& in,
                           glmbayes::fam::F2ValueFn f2v, double* coef,
                           double* iters);

/// Diagnostics accumulated with SIM_COMPARE: calls (chains), flagged
/// (max_kf or max_transform above 1e-10), max_kf (|f2_value - f2| /
/// (1 + |f2|) at bstar2 and the first and last tangent points), max_transform
/// (|L2Inv L3Inv bstar2 + mu - coef.mode| / (1 + |coef.mode|)), seconds in
/// the serial and parallel simulation, and candidates summed over draws.
struct SimCompare {
    double calls = 0.0;
    double flagged = 0.0;
    double max_kf = 0.0;
    double max_transform = 0.0;
    double t_serial = 0.0;
    double t_parallel = 0.0;
    double iters_serial = 0.0;
    double iters_parallel = 0.0;
};

/// Fold one chain into `sc` (kernel and transform checks).
void sim_compare_update(SimCompare& sc, const GroupStd& g,
                        const Rcpp::List& envelope, const GroupDrawIn& in,
                        glmbayes::fam::F2ValueFn f2v,
                        const std::string& family, const std::string& link);

/// Envelope evaluation path in rNormalGLMGroupChains, from
/// getOption("glmbayesCore.group_eval").
enum EvalMethod {
    EVAL_LEGACY = 0,   // EnvelopeEval per chain
    EVAL_BATCH = 1,    // poisson_eval_batch over every chain's grid points
    EVAL_COMPARE = 2   // legacy results kept; batch compared
};

/// EvalMethod from getOption("glmbayesCore.group_eval"): "legacy",
/// "batch" (or unset) or "compare". Uses the R API.
int group_eval_method();

/// Families poisson_eval_batch handles (poisson and quasipoisson, log link).
bool eval_batch_supported(const std::string& family);

/// Every chain's envelope evaluation for one group in flat, column-major
/// buffers. The group part (y, wt2, lg_const) is shared by all chains; the
/// chain part is stacked chain by chain: chain c's x2 starts at
/// c * nobs * l1, alpha at c * nobs, mu2 at c * l1, P2 at c * l1 * l1, and
/// its m[c] grid points (columns of G, l1 each) at column m_off[c]. Work
/// item t in [0, n_points) is grid point t - m_off[c] of the chain c with
/// m_off[c] <= t < m_off[c + 1].
struct EvalBatch {
    int nobs = 0;
    int l1 = 0;
    int n_chains = 0;
    int n_points = 0;
    std::vector<double> y;
    std::vector<double> wt2;
    double lg_const = 0.0;    // sum_i wt2[i] lgamma(y[i] + 1)
    std::vector<int> m;
    std::vector<int> m_off;   // n_chains + 1
    std::vector<double> x2;
    std::vector<double> alpha;
    std::vector<double> mu2;
    std::vector<double> P2;
    std::vector<double> G;
};

/// sum_i wt[i] lgamma(y[i] + 1): the y-only part of the Poisson negative
/// log-likelihood.
double poisson_eval_constant(const double* y, const double* wt, int nobs);

/// Pack st[c] (x2, alpha, mu2, P2) and G4[c] (l1 x m_c grid) for every
/// chain; y and wt2 come from st[0], which all chains of a call share. Uses
/// the R API, so call it outside parallel regions.
EvalBatch eval_batch_pack(const std::vector<GroupStd>& st,
                          const std::vector<Rcpp::NumericMatrix>& G4);

/// Negative log-posterior qf (n_points) and its gradient grad (l1 x
/// n_points, column t for work item t) at every grid point, with the values
/// of f2_f3_poisson. The data term is wt2 (exp(eta) - y eta) plus lg_const,
/// so each observation costs one exp. No R API, no allocation beyond an l1
/// workspace.
void poisson_eval_batch(const EvalBatch& b, double* qf, double* grad);

/// Diagnostics accumulated with EVAL_COMPARE: calls (chains), flagged
/// (max_qf or max_grad above 1e-10 for a chain), max_qf and max_grad
/// (|batch - legacy| / (1 + |legacy|)), and seconds in each path.
struct EvalCompare {
    double calls = 0.0;
    double flagged = 0.0;
    double max_qf = 0.0;
    double max_grad = 0.0;
    double t_legacy = 0.0;
    double t_batch = 0.0;
};

/// One group drawn for every chain (one draw per chain). Column c of `mu`
/// is chain c's prior mean. Runs group_mode and EnvelopeGridSetup for all
/// chains, then EnvelopeEval for all chains (the step to batch across
/// chains), then EnvelopeFinalize for all chains, then the draws.
/// getOption("glmbayesCore.group_sim") picks them (group_sim_method):
/// "serial" runs group_simulate per chain; "parallel" (default) runs
/// group_simulate_chains, whose draws come from the thread-local RNG
/// (runif_safe) rather than R's, so set.seed() does not fix them; "compare"
/// keeps the serial draws and fills sim_compare (SimCompare fields; all zero
/// otherwise). Pairs without an f2_value kernel always run serially. With
/// `check`, each chain is also run
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
/// (vmmin runs only); hessian_compare returns the HessCompare fields and
/// newton_compare the NewtonCompare fields (MODE_NEWTON* only). `check`
/// applies to MODE_OPTIM and MODE_COMPARE only.
/// getOption("glmbayesCore.group_eval") picks the evaluation
/// (group_eval_method): "batch" (default) runs poisson_eval_batch for
/// poisson and quasipoisson, "legacy" runs EnvelopeEval per chain, and
/// "compare" keeps EnvelopeEval and fills eval_compare (EvalCompare fields;
/// all zero otherwise). Other families and the OpenCL path always use
/// EnvelopeEval.
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
