// group_mode_hessian.cpp
// Exact Hessian of a group's shifted negative log posterior at the vmmin
// mode, X' diag(w) X + P, for the family/link pairs of vmmin_supported.
// The weights are the second derivatives in eta of the f2 kernels (whose
// derivatives are the f3 gradients); probit and cloglog use the observed
// information, which depends on y.

#include "RcppArmadillo.h"
#include "simfuncs_groups_split.h"

#include <cmath>
#include <string>
#include <vector>

namespace glmbayes {
namespace sim {
namespace group {

namespace {

// d2/deta2 of -[y log p + (1 - y) log(1 - p)] per unit weight, or NaN for an
// unsupported pair.
double eta_weight(double eta, double y, bool binom, bool pois, bool gam,
                  bool gauss, const std::string& link) {
  if (binom) {
    if (link == "logit") {
      double p, q;
      if (eta >= 0) {
        const double e = std::exp(-eta);
        p = 1.0 / (1.0 + e);
        q = e / (1.0 + e);
      } else {
        const double e = std::exp(eta);
        p = e / (1.0 + e);
        q = 1.0 / (1.0 + e);
      }
      return p * q;
    }
    if (link == "probit") {
      // lam1 = phi / Phi(eta), lam0 = phi / Phi(-eta), on the log scale so
      // that the tails do not underflow to 0 / 0.
      const double lphi = R::dnorm(eta, 0.0, 1.0, 1);
      const double lam1 = std::exp(lphi - R::pnorm(eta, 0.0, 1.0, 1, 1));
      const double lam0 = std::exp(lphi - R::pnorm(eta, 0.0, 1.0, 0, 1));
      return y * lam1 * (eta + lam1) + (1.0 - y) * lam0 * (lam0 - eta);
    }
    if (link == "cloglog") {
      // u = exp(eta); g = u / expm1(u) is -d/deta log p, with
      // d2/deta2 (-log p) = g (u + g - 1) and d2/deta2 (-log(1 - p)) = u.
      const double u = std::exp(eta);
      double gpart;
      if (u < 1e-4) {
        gpart = (1.0 - 0.5 * u) * (0.5 * u + u * u / 12.0);
      } else {
        const double em = std::expm1(u);
        const double g = R_FINITE(em) ? u / em : 0.0;
        gpart = g * (u + g - 1.0);
      }
      return y * gpart + (1.0 - y) * u;
    }
    return R_NaN;
  }
  if (pois) return std::exp(eta);
  if (gam) return y * std::exp(-eta);
  if (gauss) return 1.0;
  return R_NaN;
}

}  // namespace

bool exact_hessian(const GroupStd& g, const std::vector<double>& par,
                   const std::string& family, const std::string& link,
                   arma::mat& H) {
  const int l1 = g.x.ncol();
  const int l2 = g.x.nrow();
  const bool binom = family == "binomial" || family == "quasibinomial";
  const bool pois = family == "poisson" || family == "quasipoisson";
  const bool gam = family == "Gamma";
  const bool gauss = family == "gaussian";

  const arma::mat X(const_cast<double*>(g.x.begin()), l2, l1, false);
  const arma::mat P(const_cast<double*>(g.P.begin()), l1, l1, false);
  const arma::vec b(const_cast<double*>(par.data()), l1, false);
  const arma::vec alpha(const_cast<double*>(g.alpha.begin()), l2, false);
  const arma::vec eta = alpha + X * b;

  arma::vec w(l2);
  for (int j = 0; j < l2; ++j) {
    const double v =
        g.wt2[j] * eta_weight(eta[j], g.y[j], binom, pois, gam, gauss, link);
    if (!R_FINITE(v)) return false;
    // Every supported pair is log-concave, so a negative weight is rounding.
    w[j] = v < 0.0 ? 0.0 : v;
  }

  H = X.t() * (X.each_col() % w) + P;
  H = 0.5 * (H + H.t());
  arma::mat R;
  return arma::chol(R, H);
}

}  // namespace group
}  // namespace sim
}  // namespace glmbayes
