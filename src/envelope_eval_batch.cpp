// envelope_eval_batch.cpp
// Poisson envelope evaluation for every chain of one group at once, on flat
// buffers. The loop over work items (one per chain and grid point) has the
// shape of inst/cl/src/f2_f3_poisson.cl, so the buffers map onto a kernel
// with one work item per (chain, grid point).

#include "RcppArmadillo.h"
#include "simfuncs_groups_split.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace glmbayes {
namespace sim {
namespace group {

int group_eval_method() {
  SEXP opt = Rf_GetOption1(Rf_install("glmbayesCore.group_eval"));
  if (Rf_isNull(opt)) return EVAL_BATCH;
  if (!Rf_isString(opt) || Rf_length(opt) != 1) {
    Rcpp::stop("options(glmbayesCore.group_eval) must be \"legacy\", "
               "\"batch\" or \"compare\".");
  }
  const std::string s = CHAR(STRING_ELT(opt, 0));
  if (s == "legacy") return EVAL_LEGACY;
  if (s == "batch") return EVAL_BATCH;
  if (s == "compare") return EVAL_COMPARE;
  Rcpp::stop("options(glmbayesCore.group_eval) must be \"legacy\", "
             "\"batch\" or \"compare\".");
  return EVAL_LEGACY;
}

bool eval_batch_supported(const std::string& family) {
  return family == "poisson" || family == "quasipoisson";
}

double poisson_eval_constant(const double* y, const double* wt, int nobs) {
  double s = 0.0;
  for (int i = 0; i < nobs; ++i) s += wt[i] * std::lgamma(y[i] + 1.0);
  return s;
}

EvalBatch eval_batch_pack(const std::vector<GroupStd>& st,
                          const std::vector<Rcpp::NumericMatrix>& G4) {
  EvalBatch b;
  b.n_chains = static_cast<int>(st.size());
  if (b.n_chains == 0) return b;
  const GroupStd& g0 = st[0];
  b.nobs = g0.x2.nrow();
  b.l1 = g0.x2.ncol();
  const int nobs = b.nobs;
  const int l1 = b.l1;

  b.y.assign(g0.y.begin(), g0.y.end());
  b.wt2.assign(g0.wt2.begin(), g0.wt2.end());
  b.lg_const = poisson_eval_constant(b.y.data(), b.wt2.data(), nobs);

  b.m.resize(b.n_chains);
  b.m_off.resize(b.n_chains + 1);
  b.m_off[0] = 0;
  for (int c = 0; c < b.n_chains; ++c) {
    if (G4[c].nrow() != l1) Rcpp::stop("eval_batch_pack: nrow(G4) != l1");
    b.m[c] = G4[c].ncol();
    b.m_off[c + 1] = b.m_off[c] + b.m[c];
  }
  b.n_points = b.m_off[b.n_chains];

  b.x2.resize(static_cast<std::size_t>(b.n_chains) * nobs * l1);
  b.alpha.resize(static_cast<std::size_t>(b.n_chains) * nobs);
  b.mu2.resize(static_cast<std::size_t>(b.n_chains) * l1);
  b.P2.resize(static_cast<std::size_t>(b.n_chains) * l1 * l1);
  b.G.resize(static_cast<std::size_t>(b.n_points) * l1);
  for (int c = 0; c < b.n_chains; ++c) {
    const GroupStd& g = st[c];
    std::copy(g.x2.begin(), g.x2.end(),
              b.x2.begin() + static_cast<std::size_t>(c) * nobs * l1);
    std::copy(g.alpha.begin(), g.alpha.end(),
              b.alpha.begin() + static_cast<std::size_t>(c) * nobs);
    std::copy(g.mu2.begin(), g.mu2.end(),
              b.mu2.begin() + static_cast<std::size_t>(c) * l1);
    std::copy(g.P2.begin(), g.P2.end(),
              b.P2.begin() + static_cast<std::size_t>(c) * l1 * l1);
    std::copy(G4[c].begin(), G4[c].end(),
              b.G.begin() + static_cast<std::size_t>(b.m_off[c]) * l1);
  }
  return b;
}

void poisson_eval_batch(const EvalBatch& b, double* qf, double* grad) {
  const int nobs = b.nobs;
  const int l1 = b.l1;
  std::vector<double> d(l1);
  for (int c = 0; c < b.n_chains; ++c) {
    const double* x2 = b.x2.data() + static_cast<std::size_t>(c) * nobs * l1;
    const double* alpha = b.alpha.data() + static_cast<std::size_t>(c) * nobs;
    const double* mu2 = b.mu2.data() + static_cast<std::size_t>(c) * l1;
    const double* P2 = b.P2.data() + static_cast<std::size_t>(c) * l1 * l1;
    for (int t = b.m_off[c]; t < b.m_off[c + 1]; ++t) {
      const double* G = b.G.data() + static_cast<std::size_t>(t) * l1;
      double* g = grad + static_cast<std::size_t>(t) * l1;

      for (int k = 0; k < l1; ++k) d[k] = G[k] - mu2[k];
      double q = 0.0;
      for (int k = 0; k < l1; ++k) {
        double acc = 0.0;
        for (int e = 0; e < l1; ++e) acc += P2[k + l1 * e] * d[e];
        g[k] = acc;
        q += d[k] * acc;
      }
      q *= 0.5;

      for (int i = 0; i < nobs; ++i) {
        double eta = alpha[i];
        for (int k = 0; k < l1; ++k) eta += x2[i + nobs * k] * G[k];
        const double mi = std::exp(eta);
        const double w = b.wt2[i];
        q += w * (mi - b.y[i] * eta);
        const double r = (b.y[i] - mi) * w;
        for (int k = 0; k < l1; ++k) g[k] -= x2[i + nobs * k] * r;
      }
      qf[t] = q + b.lg_const;
    }
  }
}

}  // namespace group
}  // namespace sim
}  // namespace glmbayes
