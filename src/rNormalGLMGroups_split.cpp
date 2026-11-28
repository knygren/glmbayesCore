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

#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
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
                    std::string family, bool verbose,
                    const std::string& link, int mode_method,
                    ModeCompare* cmp, ModeTiming* mt, const double* warm_start,
                    int warm_method, WarmCompare* wc, int hess_method,
                    HessCompare* hc, NewtonCompare* ncmp,
                    std::vector<double>* nwork) {
  auto t_part = std::chrono::steady_clock::now();

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
  arma::mat x2(x.begin(), l2, l1, false);
  arma::vec alpha2(g.alpha.begin(), l2, false);
  arma::vec offset2b(offset.begin(), l2, false);
  arma::vec mu2(mu.begin(), l1, false);
  alpha2 = x2 * mu2 + offset2b;

  if (mode_method != MODE_OPTIM && !vmmin_supported(family, link)) {
    mode_method = MODE_OPTIM;
  }
  if (mt) mt->prep += seconds_since(t_part);

  const bool newton_primary =
      mode_method == MODE_NEWTON || mode_method == MODE_NEWTON_F4;
  const bool newton_cmp = mode_method == MODE_NEWTON_COMPARE ||
                          mode_method == MODE_NEWTON_F4_COMPARE;
  const int backend = (mode_method == MODE_NEWTON_F4 ||
                       mode_method == MODE_NEWTON_F4_COMPARE)
                          ? NEWTON_FAM
                          : NEWTON_ENV;
  std::vector<double> par0;
  if (warm_start != nullptr && warm_method != WARM_COLD) {
    par0.resize(l1);
    for (int i = 0; i < l1; ++i) par0[i] = warm_start[i] - mu[i];
  }
  std::vector<double> nwork_local;
  std::vector<double>& nw = nwork ? *nwork : nwork_local;

  // Warm (when par0 is set) then cold; status != NEWTON_OK if both fail.
  auto run_newton = [&](ModeTiming* timing) {
    NewtonResult r;
    if (!par0.empty()) {
      r = newton_mode(g, family, link, par0.data(), backend, nw, timing);
      if (r.status == NEWTON_OK) return r;
      if (ncmp) ncmp->cold_retry += 1.0;
    }
    r = newton_mode(g, family, link, nullptr, backend, nw, timing);
    if (r.status != NEWTON_OK && ncmp) ncmp->fallback += 1.0;
    return r;
  };

  bool use_vmmin = mode_method == MODE_VMMIN || newton_cmp;
  if (newton_primary) {
    const NewtonResult r = run_newton(mt);
    if (r.status == NEWTON_OK) {
      g.fit = r.fit;
      g.b2a = r.b2a;
      g.A1 = r.A1;
    } else {
      use_vmmin = true;
    }
  }

  if (use_vmmin) {
    const auto t_vm0 = std::chrono::steady_clock::now();
    VmminResult v;
    if (par0.empty()) {
      v = vmmin_mode(g, f2, f3, family, link, false, mt, nullptr, hess_method,
                     hc);
    } else if (warm_method == WARM_COMPARE) {
      auto t0 = std::chrono::steady_clock::now();
      v = vmmin_mode(g, f2, f3, family, link, false, mt, nullptr, hess_method,
                     hc);
      const double t_cold = seconds_since(t0);
      t0 = std::chrono::steady_clock::now();
      bool ok = true;
      VmminResult w;
      try {
        w = vmmin_mode(g, f2, f3, family, link, false, nullptr, par0.data(),
                       hess_method, nullptr);
      } catch (std::exception&) {
        ok = false;
      }
      const double t_warm = seconds_since(t0);
      if (wc) {
        if (ok) {
          warm_compare_update(*wc, v, w, t_cold, t_warm);
        } else {
          wc->fallback += 1.0;
        }
      }
    } else {
      bool ok = true;
      try {
        v = vmmin_mode(g, f2, f3, family, link, false, mt, par0.data(),
                       hess_method, hc);
        ok = v.convergence == 0;
      } catch (std::exception&) {
        ok = false;
      }
      if (!ok) {
        if (wc) wc->fallback += 1.0;
        v = vmmin_mode(g, f2, f3, family, link, false, mt, nullptr,
                       hess_method, hc);
      }
    }
    if (v.convergence != 0 && verbose) {
      Rcpp::Rcout << "vmmin returned non-zero convergence code: "
                  << v.convergence << std::endl;
    }
    if (v.convergence > 0) {
      Rcpp::stop("Posterior Optimization failed");
    }
    g.fit = v.fit;
    g.b2a = v.b2a;
    g.A1 = v.A1;
    if (newton_cmp && ncmp) {
      const double t_vm = seconds_since(t_vm0);
      const auto t_cmp = std::chrono::steady_clock::now();
      const NewtonResult r = run_newton(nullptr);
      const double t_nt = seconds_since(t_cmp);
      newton_compare_update(*ncmp, g, family, link, v, r, t_vm, t_nt,
                            backend == NEWTON_FAM);
      if (mt) mt->compare += seconds_since(t_cmp);
    }
  } else if (!newton_primary) {
    auto t0 = std::chrono::steady_clock::now();
    Rcpp::Function asMat("as.matrix");
    Rcpp::Function optfun("optim");
    Rcpp::Function tryfun("try");
    NumericVector start = Rcpp::clone(mu);
    NumericVector parin = start - mu;
    NumericVector mu1 = mu - mu;
    Rcpp::Shield<SEXP> optSEXP(tryfun(
        optfun(Rcpp::_["par"] = parin, Rcpp::_["fn"] = f2,
               Rcpp::_["gr"] = f3, Rcpp::_["y"] = y, Rcpp::_["x"] = x,
               Rcpp::_["mu"] = mu1, Rcpp::_["P"] = P,
               Rcpp::_["alpha"] = g.alpha, Rcpp::_["wt"] = g.wt2,
               Rcpp::_["method"] = "BFGS", Rcpp::_["hessian"] = true),
        Rcpp::_["silent"] = true));
    if (Rf_inherits(optSEXP, "try-error")) {
      Rcpp::stop("Optimization failed in group_mode");
    }
    List opt(optSEXP);
    const double t_optim = seconds_since(t0);
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
    if (mt) mt->optim += seconds_since(t0);

    if (mode_method == MODE_COMPARE && cmp != nullptr) {
      const auto t_cmp = std::chrono::steady_clock::now();
      t0 = std::chrono::steady_clock::now();
      const VmminResult vr = vmmin_mode(g, f2, f3, family, link, true);
      const double t_r = seconds_since(t0);
      t0 = std::chrono::steady_clock::now();
      const VmminResult vc = vmmin_mode(g, f2, f3, family, link, false);
      const double t_c = seconds_since(t0);
      cmp->t_optim += t_optim;
      cmp->t_vmmin_r += t_r;
      cmp->t_vmmin_cpp += t_c;
      mode_compare_update(*cmp, g, opt, vr, vc, f2);
      if (mt) mt->compare += seconds_since(t_cmp);
    }
  }

  t_part = std::chrono::steady_clock::now();
  List sm = glmbayes::sim::glmb_Standardize_Model(y, x, P, g.b2a, g.A1);
  g.bstar2 = Rcpp::as<NumericVector>(sm[0]);
  g.A = Rcpp::as<NumericMatrix>(sm[1]);
  g.x2 = Rcpp::as<NumericMatrix>(sm[2]);
  g.mu2 = Rcpp::as<NumericMatrix>(sm[3]);
  g.P2 = Rcpp::as<NumericMatrix>(sm[4]);
  g.L2Inv = Rcpp::as<NumericMatrix>(sm[5]);
  g.L3Inv = Rcpp::as<NumericMatrix>(sm[6]);
  if (mt) {
    mt->standardize += seconds_since(t_part);
    mt->calls += 1.0;
  }
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
                           bool verbose, int mode_method, int setup_method,
                           SEXP start, int warm_method, int hess_method) {
  const int l1 = x.ncol();
  const int l2 = x.nrow();
  const int nc = mu.ncol();
  if (mu.nrow() != l1) Rcpp::stop("nrow(mu) must equal ncol(x)");
  NumericMatrix start_m;
  const bool has_start = !Rf_isNull(start);
  if (has_start) {
    start_m = Rcpp::as<NumericMatrix>(start);
    if (start_m.nrow() != l1 || start_m.ncol() != nc) {
      Rcpp::stop("start must be ncol(x) x ncol(mu)");
    }
  }
  if (y.size() != l2 || offset.size() != l2 || wt.size() != l2) {
    Rcpp::stop("y, offset and wt must have nrow(x) elements");
  }
  if (P.nrow() != l1 || P.ncol() != l1) Rcpp::stop("P must be l1 x l1");

  const int n = 1;
  const int n_envopt_b = (n_envopt < 0) ? n : n_envopt;
  const bool use_ocl = glmbayes::env::envelope_use_opencl(use_opencl, verbose);

  double t_mode = 0.0, t_setup = 0.0, t_eval = 0.0, t_fin = 0.0, t_sim = 0.0;
  ModeCompare cmp;
  ModeTiming mtm;
  WarmCompare wcm;
  HessCompare hcm;
  NewtonCompare ncm;
  std::vector<double> nwork;
  glmbayes::env::GridSetupTiming gst;
  // rNormalGLM uses optim, so its mode only matches when optim is kept.
  const bool do_check =
      check && (mode_method == MODE_OPTIM || mode_method == MODE_COMPARE);

  // 1. Mode and grid setup for every chain.
  std::vector<GroupStd> st;
  std::vector<glmbayes::env::EnvGrid> grids;
  st.reserve(nc);
  grids.reserve(nc);
  for (int c = 0; c < nc; ++c) {
    NumericVector mu_c = mu(Rcpp::_, c);
    const double* start_c = nullptr;
    if (has_start) {
      start_c = &start_m(0, c);
      for (int i = 0; i < l1; ++i) {
        if (ISNAN(start_c[i])) {
          start_c = nullptr;
          break;
        }
      }
    }
    auto t0 = std::chrono::steady_clock::now();
    st.push_back(group_mode(y, x, offset, wt, mu_c, P, dispersion, f2, f3,
                            family, verbose && (c == 0), link, mode_method,
                            &cmp, &mtm, start_c, warm_method, &wcm,
                            hess_method, &hcm, &ncm, &nwork));
    t_mode += seconds_since(t0);
    t0 = std::chrono::steady_clock::now();
    grids.push_back(glmbayes::env::EnvelopeGridSetup(
        st[c].bstar2, st[c].A, Gridtype, n, n_envopt_b, use_ocl, false,
        &gst, setup_method));
    t_setup += seconds_since(t0);
  }

  // 2. Envelope evaluation for every chain.
  std::vector<List> evals(nc);
  EvalCompare ecm;
  {
    int eval_method = group_eval_method();
    if (use_ocl || !eval_batch_supported(family)) eval_method = EVAL_LEGACY;

    if (eval_method != EVAL_BATCH) {
      const auto t0 = std::chrono::steady_clock::now();
      for (int c = 0; c < nc; ++c) {
        const GroupStd& g = st[c];
        evals[c] = glmbayes::env::EnvelopeEval(grids[c].G4, g.y, g.x2, g.mu2,
                                               g.P2, g.alpha, g.wt2, family,
                                               link, use_ocl, false);
      }
      const double t_l = seconds_since(t0);
      t_eval += t_l;
      ecm.t_legacy += (eval_method == EVAL_COMPARE) ? t_l : 0.0;
    }

    if (eval_method != EVAL_LEGACY) {
      const auto t0 = std::chrono::steady_clock::now();
      std::vector<NumericMatrix> G4(nc);
      for (int c = 0; c < nc; ++c) G4[c] = grids[c].G4;
      const EvalBatch eb = eval_batch_pack(st, G4);
      std::vector<double> qf(eb.n_points);
      std::vector<double> gr(static_cast<std::size_t>(eb.n_points) * l1);
      poisson_eval_batch(eb, qf.data(), gr.data());

      for (int c = 0; c < nc; ++c) {
        const int mc = eb.m[c];
        const int off = eb.m_off[c];
        if (eval_method == EVAL_BATCH) {
          NumericVector NegLL(mc);
          arma::mat cbars(mc, l1);
          for (int j = 0; j < mc; ++j) {
            NegLL[j] = qf[off + j];
            for (int k = 0; k < l1; ++k) {
              cbars(j, k) = gr[static_cast<std::size_t>(off + j) * l1 + k];
            }
          }
          evals[c] = List::create(Rcpp::Named("NegLL") = NegLL,
                                  Rcpp::Named("cbars") = cbars);
        } else {
          const NumericVector NegLL =
              Rcpp::as<NumericVector>(evals[c]["NegLL"]);
          const arma::mat cbars = Rcpp::as<arma::mat>(evals[c]["cbars"]);
          double dq = 0.0, dg = 0.0;
          for (int j = 0; j < mc; ++j) {
            dq = std::max(dq, std::abs(qf[off + j] - NegLL[j]) /
                                  (1.0 + std::abs(NegLL[j])));
            for (int k = 0; k < l1; ++k) {
              const double ref = cbars(j, k);
              const double v = gr[static_cast<std::size_t>(off + j) * l1 + k];
              dg = std::max(dg, std::abs(v - ref) / (1.0 + std::abs(ref)));
            }
          }
          ecm.calls += 1.0;
          if (dq > 1e-10 || dg > 1e-10) ecm.flagged += 1.0;
          ecm.max_qf = std::max(ecm.max_qf, dq);
          ecm.max_grad = std::max(ecm.max_grad, dg);
        }
      }
      const double t_b = seconds_since(t0);
      if (eval_method == EVAL_BATCH) {
        t_eval += t_b;
      } else {
        ecm.t_batch += t_b;
      }
    }
  }

  // 3. Finalize every chain, then draw (serial or parallel), then collect.
  NumericMatrix coef(nc, l1);
  NumericMatrix coef_mode(nc, l1);
  NumericMatrix mode_hessian(l1 * l1, nc);
  NumericVector iters(nc);
  for (int c = 0; c < nc; ++c) {
    const NumericMatrix& A1 = st[c].A1;
    for (int k = 0; k < l1 * l1; ++k) mode_hessian(k, c) = A1[k];
  }
  std::vector<List> envs(nc);
  for (int c = 0; c < nc; ++c) {
    const auto t0 = std::chrono::steady_clock::now();
    envs[c] = glmbayes::env::EnvelopeFinalize(
        grids[c], Rcpp::as<NumericVector>(evals[c]["NegLL"]),
        Rcpp::as<arma::mat>(evals[c]["cbars"]), st[c].x2, n > 1, false);
    t_fin += seconds_since(t0);
    evals[c] = List();
    grids[c] = glmbayes::env::EnvGrid();
  }

  int sim_method = group_sim_method();
  const glmbayes::fam::F2ValueFn f2v = group_f2_value(family, link);
  if (f2v == nullptr) sim_method = SIM_SERIAL;
  SimCompare scm;
  auto run_check = [&](int c) {
    NumericVector mu_c = Rcpp::clone(st[c].mu);
    const List ref = glmbayes::sim::rNormalGLM(
        n, Rcpp::clone(y), Rcpp::clone(x), mu_c, Rcpp::clone(P),
        Rcpp::clone(offset), Rcpp::clone(wt), dispersion, f2, f3,
        Rcpp::clone(mu_c), family, link, Gridtype, n_envopt_b, false,
        use_opencl, false);
    NumericMatrix cm(1, l1);
    for (int j = 0; j < l1; ++j) cm(0, j) = coef_mode(c, j);
    check_against_rNormalGLM(ref, cm, envs[c], c);
  };
  std::vector<GroupDrawIn> din;
  if (sim_method != SIM_SERIAL) {
    din.reserve(nc);
    for (int c = 0; c < nc; ++c) din.push_back(group_draw_in(st[c], envs[c]));
  }

  if (sim_method == SIM_PARALLEL) {
    const auto t0 = std::chrono::steady_clock::now();
    group_simulate_chains(din, f2v, coef.begin(), iters.begin());
    t_sim += seconds_since(t0);
    for (int c = 0; c < nc; ++c) {
      for (int j = 0; j < l1; ++j) {
        coef_mode(c, j) = st[c].b2a[j] + st[c].mu[j];
      }
      if (do_check) run_check(c);
    }
  } else {
    for (int c = 0; c < nc; ++c) {
      const auto t0 = std::chrono::steady_clock::now();
      const List sim = group_simulate(st[c], envs[c], n, f2, family, link,
                                      false, false);
      const double t_c = seconds_since(t0);
      t_sim += t_c;
      const NumericMatrix cs = Rcpp::as<NumericMatrix>(sim["coefficients"]);
      const NumericMatrix cm = Rcpp::as<NumericMatrix>(sim["coef.mode"]);
      for (int j = 0; j < l1; ++j) {
        coef(c, j) = cs(0, j);
        coef_mode(c, j) = cm(0, j);
      }
      iters[c] = Rcpp::as<NumericVector>(sim["iters"])[0];
      if (do_check) run_check(c);
      if (sim_method == SIM_COMPARE) {
        scm.t_serial += t_c;
        scm.iters_serial += iters[c];
        sim_compare_update(scm, st[c], envs[c], din[c], f2v, family, link);
      }
    }
    if (sim_method == SIM_COMPARE) {
      NumericMatrix coef_p(nc, l1);
      NumericVector iters_p(nc);
      const auto t0 = std::chrono::steady_clock::now();
      group_simulate_chains(din, f2v, coef_p.begin(), iters_p.begin());
      scm.t_parallel += seconds_since(t0);
      for (int c = 0; c < nc; ++c) scm.iters_parallel += iters_p[c];
    }
  }

  return List::create(
      Rcpp::Named("coefficients") = coef,
      Rcpp::Named("coef.mode") = coef_mode,
      Rcpp::Named("mode_hessian") = mode_hessian,
      Rcpp::Named("iters") = iters,
      Rcpp::Named("timing") = NumericVector::create(
          Rcpp::Named("mode") = t_mode, Rcpp::Named("setup") = t_setup,
          Rcpp::Named("eval") = t_eval, Rcpp::Named("finalize") = t_fin,
          Rcpp::Named("simulate") = t_sim),
      Rcpp::Named("setup_timing") = NumericVector::create(
          Rcpp::Named("arith") = gst.arith,
          Rcpp::Named("core_count") = gst.core_count,
          Rcpp::Named("envelope_opt") = gst.envelope_opt,
          Rcpp::Named("size_lists") = gst.size_lists,
          Rcpp::Named("gindex_expand") = gst.gindex_expand,
          Rcpp::Named("gindex_asmat") = gst.gindex_asmat,
          Rcpp::Named("g3_expand") = gst.g3_expand,
          Rcpp::Named("g3_asmat") = gst.g3_asmat,
          Rcpp::Named("g4") = gst.g4,
          Rcpp::Named("cpp_opt") = gst.cpp_opt,
          Rcpp::Named("cpp_grid") = gst.cpp_grid),
      Rcpp::Named("mode_timing") = NumericVector::create(
          Rcpp::Named("prep") = mtm.prep,
          Rcpp::Named("optim") = mtm.optim,
          Rcpp::Named("vmmin") = mtm.vmmin,
          Rcpp::Named("vmmin_kernel") = mtm.vmmin_kernel,
          Rcpp::Named("hessian") = mtm.hessian,
          Rcpp::Named("hessian_kernel") = mtm.hessian_kernel,
          Rcpp::Named("fit_list") = mtm.fit_list,
          Rcpp::Named("compare") = mtm.compare,
          Rcpp::Named("standardize") = mtm.standardize,
          Rcpp::Named("calls") = mtm.calls,
          Rcpp::Named("kernel_calls") = mtm.kernel_calls,
          Rcpp::Named("newton") = mtm.newton,
          Rcpp::Named("newton_kernel") = mtm.newton_kernel,
          Rcpp::Named("newton_iters") = mtm.newton_iters,
          Rcpp::Named("newton_backtracks") = mtm.newton_backtracks),
      Rcpp::Named("setup_compare") = NumericVector::create(
          Rcpp::Named("calls") = gst.compared),
      Rcpp::Named("mode_compare") = NumericVector::create(
          Rcpp::Named("calls") = cmp.calls,
          Rcpp::Named("exact_mismatch") = cmp.exact_mismatch,
          Rcpp::Named("flagged") = cmp.flagged,
          Rcpp::Named("max_d") = cmp.max_d,
          Rcpp::Named("max_df") = cmp.max_df,
          Rcpp::Named("max_hess") = cmp.max_hess,
          Rcpp::Named("t_optim") = cmp.t_optim,
          Rcpp::Named("t_vmmin_r") = cmp.t_vmmin_r,
          Rcpp::Named("t_vmmin_cpp") = cmp.t_vmmin_cpp,
          Rcpp::Named("fn_optim") = cmp.fn_optim,
          Rcpp::Named("gr_optim") = cmp.gr_optim,
          Rcpp::Named("fn_vmmin") = cmp.fn_vmmin,
          Rcpp::Named("gr_vmmin") = cmp.gr_vmmin),
      Rcpp::Named("warm_compare") = NumericVector::create(
          Rcpp::Named("calls") = wcm.calls,
          Rcpp::Named("flagged") = wcm.flagged,
          Rcpp::Named("fallback") = wcm.fallback,
          Rcpp::Named("max_d") = wcm.max_d,
          Rcpp::Named("max_df") = wcm.max_df,
          Rcpp::Named("max_hess") = wcm.max_hess,
          Rcpp::Named("fn_cold") = wcm.fn_cold,
          Rcpp::Named("gr_cold") = wcm.gr_cold,
          Rcpp::Named("fn_warm") = wcm.fn_warm,
          Rcpp::Named("gr_warm") = wcm.gr_warm,
          Rcpp::Named("t_cold") = wcm.t_cold,
          Rcpp::Named("t_warm") = wcm.t_warm),
      Rcpp::Named("hessian_compare") = NumericVector::create(
          Rcpp::Named("calls") = hcm.calls,
          Rcpp::Named("flagged") = hcm.flagged,
          Rcpp::Named("fallback") = hcm.fallback,
          Rcpp::Named("max_hess") = hcm.max_hess,
          Rcpp::Named("t_numeric") = hcm.t_numeric,
          Rcpp::Named("t_exact") = hcm.t_exact),
      Rcpp::Named("newton_compare") = NumericVector::create(
          Rcpp::Named("calls") = ncm.calls,
          Rcpp::Named("flagged") = ncm.flagged,
          Rcpp::Named("cold_retry") = ncm.cold_retry,
          Rcpp::Named("fallback") = ncm.fallback,
          Rcpp::Named("newton_lower") = ncm.newton_lower,
          Rcpp::Named("max_d") = ncm.max_d,
          Rcpp::Named("max_df") = ncm.max_df,
          Rcpp::Named("max_hess") = ncm.max_hess,
          Rcpp::Named("max_kf") = ncm.max_kf,
          Rcpp::Named("max_kgrad") = ncm.max_kgrad,
          Rcpp::Named("max_kH") = ncm.max_kH,
          Rcpp::Named("iters") = ncm.iters,
          Rcpp::Named("fevals") = ncm.fevals,
          Rcpp::Named("gevals") = ncm.gevals,
          Rcpp::Named("backtracks") = ncm.backtracks,
          Rcpp::Named("t_vmmin") = ncm.t_vmmin,
          Rcpp::Named("t_newton") = ncm.t_newton),
      Rcpp::Named("sim_compare") = NumericVector::create(
          Rcpp::Named("calls") = scm.calls,
          Rcpp::Named("flagged") = scm.flagged,
          Rcpp::Named("max_kf") = scm.max_kf,
          Rcpp::Named("max_transform") = scm.max_transform,
          Rcpp::Named("t_serial") = scm.t_serial,
          Rcpp::Named("t_parallel") = scm.t_parallel,
          Rcpp::Named("iters_serial") = scm.iters_serial,
          Rcpp::Named("iters_parallel") = scm.iters_parallel),
      Rcpp::Named("eval_compare") = NumericVector::create(
          Rcpp::Named("calls") = ecm.calls,
          Rcpp::Named("flagged") = ecm.flagged,
          Rcpp::Named("max_qf") = ecm.max_qf,
          Rcpp::Named("max_grad") = ecm.max_grad,
          Rcpp::Named("t_legacy") = ecm.t_legacy,
          Rcpp::Named("t_batch") = ecm.t_batch));
}

}  // namespace group
}  // namespace sim
}  // namespace glmbayes
