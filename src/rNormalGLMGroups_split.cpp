// rNormalGLMGroups_split.cpp
// One group's rNormalGLM() step split into mode finding, envelope
// construction and simulation. Each phase mirrors the matching steps of
// rNormalGLM() (rNormalGLM.cpp), so the mode and envelope are identical
// for the same inputs.

#include "RcppArmadillo.h"
#include "Envelopefuncs.h"
#include "simfuncs.h"
#include "Envelopefuncs_split.h"
#include "simfuncs_groups_split.h"

#include <chrono>
#include <cmath>
#include <string>
#include <vector>

// Defined at global scope in rNormalGLM.cpp.
Rcpp::List rNormalGLM_std_parallel(
    int n,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericMatrix mu,
    Rcpp::NumericMatrix P,
    Rcpp::NumericVector alpha,
    Rcpp::NumericVector wt,
    Rcpp::Function f2,
    Rcpp::List Envelope,
    Rcpp::CharacterVector family,
    Rcpp::CharacterVector link,
    int progbar,
    bool verbose
);

namespace glmbayes {
namespace sim {
namespace group {

using Rcpp::CharacterVector;
using Rcpp::Function;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

namespace {
double seconds_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
      .count();
}

// Exact element-wise equality; NaN matches NaN.
bool same_numeric(SEXP a, SEXP b) {
  if (Rf_isNull(a) || Rf_isNull(b)) return Rf_isNull(a) && Rf_isNull(b);
  NumericVector va = Rcpp::as<NumericVector>(a);
  NumericVector vb = Rcpp::as<NumericVector>(b);
  if (va.size() != vb.size()) return false;
  for (R_xlen_t i = 0; i < va.size(); ++i) {
    const double p = va[i];
    const double q = vb[i];
    if (std::isnan(p) && std::isnan(q)) continue;
    if (p != q) return false;
  }
  return true;
}

// Stop when the split mode or envelope differs from rNormalGLM's.
void check_against_rNormalGLM(const List& ref, const NumericMatrix& mode,
                              const List& env, int chain) {
  if (!same_numeric(ref["coef.mode"], mode)) {
    Rcpp::stop("rNormalGLMGroupChains check: mode differs for chain %d",
               chain + 1);
  }
  const List ref_env = ref["Envelope"];
  const CharacterVector nm = ref_env.names();
  for (int f = 0; f < ref_env.size(); ++f) {
    const std::string field = Rcpp::as<std::string>(nm[f]);
    if (!env.containsElementNamed(field.c_str()) ||
        !same_numeric(ref_env[f], env[field])) {
      Rcpp::stop("rNormalGLMGroupChains check: envelope field '%s' differs "
                 "for chain %d", field, chain + 1);
    }
  }
}
}  // namespace

// rNormalGLM "Step 1" to "Step 3": shift, optim, standardize.
GroupStd group_mode(NumericVector y, NumericMatrix x, NumericVector offset,
                    NumericVector wt, NumericVector mu, NumericMatrix P,
                    double dispersion, Function f2, Function f3,
                    std::string family, bool verbose) {
  Rcpp::Function asMat("as.matrix");
  Rcpp::Function optfun("optim");
  Rcpp::Function tryfun("try");

  GroupStd g;
  g.y = y;
  g.x = x;
  g.offset = offset;
  g.wt = wt;
  g.mu = mu;
  g.P = P;

  const int l1 = x.ncol();
  const int l2 = x.nrow();

  g.dispersion =
      (family == "poisson" || family == "binomial") ? 1.0 : dispersion;
  g.wt2 = wt / g.dispersion;

  g.alpha = NumericVector(l2);
  NumericMatrix mu2a = asMat(mu);
  arma::mat x2(x.begin(), l2, l1, false);
  arma::vec alpha2(g.alpha.begin(), l2, false);
  arma::vec offset2b(offset.begin(), l2, false);
  arma::mat mu2(mu2a.begin(), mu2a.nrow(), mu2a.ncol(), false);
  alpha2 = x2 * mu2 + offset2b;

  NumericVector start = Rcpp::clone(mu);
  NumericVector parin = start - mu;
  NumericVector mu1 = mu - mu;

  Rcpp::Shield<SEXP> optSEXP(tryfun(
      optfun(Rcpp::_["par"] = parin, Rcpp::_["fn"] = f2, Rcpp::_["gr"] = f3,
             Rcpp::_["y"] = y, Rcpp::_["x"] = x, Rcpp::_["mu"] = mu1,
             Rcpp::_["P"] = P, Rcpp::_["alpha"] = g.alpha,
             Rcpp::_["wt"] = g.wt2, Rcpp::_["method"] = "BFGS",
             Rcpp::_["hessian"] = true),
      Rcpp::_["silent"] = true));
  if (Rf_inherits(optSEXP, "try-error")) {
    Rcpp::stop("Optimization failed in group_mode");
  }
  List opt(optSEXP);
  const int conv = Rcpp::as<int>(opt["convergence"]);
  if (conv != 0 && verbose) {
    Rcpp::Rcout << "optim() returned non-zero convergence code: " << conv
                << std::endl;
  }
  if (conv > 0) {
    Rcpp::stop("Posterior Optimization failed");
  }

  g.fit = opt;
  g.b2a = Rcpp::as<NumericMatrix>(asMat(opt[0]));
  g.A1 = Rcpp::as<NumericMatrix>(opt[5]);

  List sm = glmbayes::sim::glmb_Standardize_Model(y, x, P, g.b2a, g.A1);
  g.bstar2 = Rcpp::as<NumericVector>(sm[0]);
  g.A = Rcpp::as<NumericMatrix>(sm[1]);
  g.x2 = Rcpp::as<NumericMatrix>(sm[2]);
  g.mu2 = Rcpp::as<NumericMatrix>(sm[3]);
  g.P2 = Rcpp::as<NumericMatrix>(sm[4]);
  g.L2Inv = Rcpp::as<NumericMatrix>(sm[5]);
  g.L3Inv = Rcpp::as<NumericMatrix>(sm[6]);
  return g;
}

// rNormalGLM "Step 4", with EnvelopeBuild split around EnvelopeEval.
List group_envelope(const GroupStd& g, std::string family, std::string link,
                    int Gridtype, int n, int n_envopt, bool use_opencl,
                    bool verbose, EnvTiming* timing) {
  if (n_envopt < 0) n_envopt = n;
  const bool use_ocl = glmbayes::env::envelope_use_opencl(use_opencl, verbose);

  auto t0 = std::chrono::steady_clock::now();
  const glmbayes::env::EnvGrid grid = glmbayes::env::EnvelopeGridSetup(
      g.bstar2, g.A, Gridtype, n, n_envopt, use_ocl, verbose);
  if (timing) timing->setup += seconds_since(t0);

  t0 = std::chrono::steady_clock::now();
  List ev = glmbayes::env::EnvelopeEval(grid.G4, g.y, g.x2, g.mu2, g.P2,
                                        g.alpha, g.wt2, family, link, use_ocl,
                                        verbose);
  if (timing) timing->eval += seconds_since(t0);

  t0 = std::chrono::steady_clock::now();
  List env = glmbayes::env::EnvelopeFinalize(
      grid, Rcpp::as<NumericVector>(ev["NegLL"]),
      Rcpp::as<arma::mat>(ev["cbars"]), g.x2, n > 1, verbose);
  if (timing) timing->finalize += seconds_since(t0);
  return env;
}

// rNormalGLM "Step 5" and "Step 6".
List group_simulate(const GroupStd& g, const List& envelope, int n,
                    Function f2, std::string family, std::string link,
                    bool use_parallel, bool verbose) {
  CharacterVector family_cv(1);
  family_cv[0] = family;
  CharacterVector link_cv(1);
  link_cv[0] = link;
  const int progbar = 0;

  List sim = (use_parallel && n > 1)
      ? rNormalGLM_std_parallel(n, g.y, g.x2, g.mu2, g.P2, g.alpha, g.wt2, f2,
                                envelope, family_cv, link_cv, progbar, verbose)
      : glmbayes::sim::rNormalGLM_std(n, g.y, g.x2, g.mu2, g.P2, g.alpha,
                                      g.wt2, f2, envelope, family_cv, link_cv,
                                      progbar, verbose);

  const int l1 = g.x.ncol();
  NumericMatrix sim2 = sim[0];
  NumericMatrix outm(l1, n);
  arma::mat sim2b(sim2.begin(), sim2.nrow(), sim2.ncol(), false);
  arma::mat out2(outm.begin(), l1, n, false);
  arma::mat L2Inv(const_cast<double*>(g.L2Inv.begin()), g.L2Inv.nrow(),
                  g.L2Inv.ncol(), false);
  arma::mat L3Inv(const_cast<double*>(g.L3Inv.begin()), g.L3Inv.nrow(),
                  g.L3Inv.ncol(), false);
  out2 = L2Inv * L3Inv * arma::trans(sim2b);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < l1; ++j) {
      outm(j, i) = outm(j, i) + g.mu[j];
    }
  }

  NumericMatrix coef_mode(1, l1);
  for (int j = 0; j < l1; ++j) {
    coef_mode(0, j) = g.b2a[j] + g.mu[j];
  }

  return List::create(
      Rcpp::Named("coefficients") = Rcpp::transpose(outm),
      Rcpp::Named("coef.mode") = coef_mode,
      Rcpp::Named("dispersion") = g.dispersion,
      Rcpp::Named("Prior") = List::create(Rcpp::Named("mean") = g.mu,
                                          Rcpp::Named("Precision") = g.P),
      Rcpp::Named("offset") = g.offset, Rcpp::Named("prior.weights") = g.wt,
      Rcpp::Named("y") = g.y, Rcpp::Named("x") = g.x,
      Rcpp::Named("fit") = g.fit, Rcpp::Named("iters") = sim[1],
      Rcpp::Named("Envelope") = envelope);
}

List rNormalGLMGroupChains(NumericVector y, NumericMatrix x,
                           NumericVector offset, NumericVector wt,
                           double dispersion, NumericMatrix mu,
                           NumericMatrix P, Function f2, Function f3,
                           std::string family, std::string link, int Gridtype,
                           int n_envopt, bool use_opencl, bool check,
                           bool verbose) {
  const int l1 = x.ncol();
  const int l2 = x.nrow();
  const int nc = mu.ncol();
  if (mu.nrow() != l1) Rcpp::stop("nrow(mu) must equal ncol(x)");
  if (y.size() != l2 || offset.size() != l2 || wt.size() != l2) {
    Rcpp::stop("y, offset and wt must have nrow(x) elements");
  }
  if (P.nrow() != l1 || P.ncol() != l1) Rcpp::stop("P must be l1 x l1");

  const int n = 1;
  const int n_envopt_b = (n_envopt < 0) ? n : n_envopt;
  const bool use_ocl = glmbayes::env::envelope_use_opencl(use_opencl, verbose);

  // 1. Mode and grid setup for every chain.
  std::vector<GroupStd> st;
  std::vector<glmbayes::env::EnvGrid> grids;
  st.reserve(nc);
  grids.reserve(nc);
  for (int c = 0; c < nc; ++c) {
    NumericVector mu_c = mu(Rcpp::_, c);
    st.push_back(group_mode(y, x, offset, wt, mu_c, P, dispersion, f2, f3,
                            family, verbose && (c == 0)));
    grids.push_back(glmbayes::env::EnvelopeGridSetup(
        st[c].bstar2, st[c].A, Gridtype, n, n_envopt_b, use_ocl, false));
  }

  // 2. EnvelopeEval for every chain.
  std::vector<List> evals(nc);
  for (int c = 0; c < nc; ++c) {
    const GroupStd& g = st[c];
    evals[c] = glmbayes::env::EnvelopeEval(grids[c].G4, g.y, g.x2, g.mu2,
                                           g.P2, g.alpha, g.wt2, family, link,
                                           use_ocl, false);
  }

  // 3. Finalize, draw, release.
  NumericMatrix coef(nc, l1);
  NumericMatrix coef_mode(nc, l1);
  NumericVector iters(nc);
  for (int c = 0; c < nc; ++c) {
    List env = glmbayes::env::EnvelopeFinalize(
        grids[c], Rcpp::as<NumericVector>(evals[c]["NegLL"]),
        Rcpp::as<arma::mat>(evals[c]["cbars"]), st[c].x2, n > 1, false);
    const List sim = group_simulate(st[c], env, n, f2, family, link, false,
                                    false);
    const NumericMatrix cs = Rcpp::as<NumericMatrix>(sim["coefficients"]);
    const NumericMatrix cm = Rcpp::as<NumericMatrix>(sim["coef.mode"]);

    if (check) {
      NumericVector mu_c = Rcpp::clone(st[c].mu);
      const List ref = glmbayes::sim::rNormalGLM(
          n, Rcpp::clone(y), Rcpp::clone(x), mu_c, Rcpp::clone(P),
          Rcpp::clone(offset), Rcpp::clone(wt), dispersion, f2, f3,
          Rcpp::clone(mu_c), family, link, Gridtype, n_envopt_b, false,
          use_opencl, false);
      check_against_rNormalGLM(ref, cm, env, c);
    }

    for (int j = 0; j < l1; ++j) {
      coef(c, j) = cs(0, j);
      coef_mode(c, j) = cm(0, j);
    }
    iters[c] = Rcpp::as<NumericVector>(sim["iters"])[0];

    evals[c] = List();
    grids[c] = glmbayes::env::EnvGrid();
  }

  return List::create(Rcpp::Named("coefficients") = coef,
                      Rcpp::Named("coef.mode") = coef_mode,
                      Rcpp::Named("iters") = iters);
}

}  // namespace group
}  // namespace sim
}  // namespace glmbayes
