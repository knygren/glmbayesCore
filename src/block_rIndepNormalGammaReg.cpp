// block_rIndepNormalGammaReg.cpp
// Block envelope ING sampler: BlockEnvelopeCentering, BlockEnvelopeBuild,
// BlockEnvelopeSim, and rIndepNormalGammaRegBlock orchestrator.
// New code only — does not modify existing sampler implementations.

#include "RcppArmadillo.h"
#include "Envelopefuncs.h"
#include "Envelopefuncs_group.h"
#include "famfuncs.h"
#include "simfuncs.h"
#include "simfuncs_groups.h"
#include "package_ns.h"
#include "R_interface.h"
#include "rng_utils.h"
#include "progress_utils.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

Rcpp::NumericVector EnvBuildLinBound_cpp(
    Rcpp::NumericMatrix thetabars,
    Rcpp::NumericMatrix cbars,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericMatrix P,
    Rcpp::NumericVector alpha,
    double dispstar
);
Rcpp::NumericVector thetabar_const_cpp(
    Rcpp::NumericMatrix P,
    Rcpp::NumericMatrix cbars,
    Rcpp::NumericMatrix thetabars
);

namespace glmbayes {
namespace env {
namespace group {

namespace {

using Rcpp::CharacterVector;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::Nullable;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::RObject;
using Rcpp::Function;

using glmbayes::env::EnvelopeBuild;
using glmbayes::env::EnvelopeDispersionBuild;
using glmbayes::fam::Inv_f3_precompute_disp;

void add_ub2_bound_matrices_to_cache(
    List& cache,
    double low,
    double upp
) {
  if (cache.containsElementNamed("M_min") && cache.containsElementNamed("M_max")) {
    return;
  }
  arma::mat base_A = cache["base_A"];
  arma::mat Pmat = cache["Pmat"];
  Pmat = 0.5 * (Pmat + Pmat.t());

  arma::mat A_max = Pmat + base_A / low;
  A_max = 0.5 * (A_max + A_max.t());
  arma::mat Ainv_max = arma::inv_sympd(A_max);
  arma::mat M_min = Ainv_max.t() * base_A * Ainv_max;
  M_min = 0.5 * (M_min + M_min.t());

  arma::mat A_min = Pmat + base_A / upp;
  A_min = 0.5 * (A_min + A_min.t());
  arma::mat Ainv_min = arma::inv_sympd(A_min);
  arma::mat M_max = Ainv_min.t() * base_A * Ainv_min;
  M_max = 0.5 * (M_max + M_max.t());

  cache["A_max"] = A_max;
  cache["M_min"] = M_min;
  cache["A_min"] = A_min;
  cache["M_max"] = M_max;
}
using glmbayes::fam::Inv_f3_with_disp;
using glmbayes::fam::f2_gaussian;
using glmbayes::sim::glmb_Standardize_Model;
using glmbayes::progress::progress_bar;
using glmbayes::rng::runif_safe;
using glmbayes::rng::rnorm_ct;
using glmbayes::rng::rinvgamma_ct_safe;

#define BEB_DBG(verbose, msg) \
  do { \
    if (verbose) { \
      Rcpp::Rcout << msg << std::endl; \
    } \
  } while (0)

static inline double max_vec_local(const NumericVector& v) {
  double m = R_NegInf;
  for (int i = 0; i < v.size(); ++i) {
    if (v[i] > m) {
      m = v[i];
    }
  }
  return m;
}

static inline bool has_non_null(const List& pl, const char* name) {
  return pl.containsElementNamed(name) && !Rf_isNull(pl[name]);
}

static inline int resolve_n_envopt(const Nullable<int>& n_envopt) {
  if (!n_envopt.isSet() || n_envopt.isNull()) {
    return -1;
  }
  return Rcpp::as<int>(n_envopt.get());
}

static inline bool is_symmetric_mat(const NumericMatrix& M, double tol = 1e-10) {
  if (M.nrow() != M.ncol()) return false;
  for (int i = 0; i < M.nrow(); ++i) {
    for (int j = i + 1; j < M.ncol(); ++j) {
      if (std::fabs(M(i, j) - M(j, i)) > tol) return false;
    }
  }
  return true;
}

static void check_P_pd(const NumericMatrix& P, const char* what) {
  arma::mat P_arma(const_cast<double*>(P.begin()), P.nrow(), P.ncol(), false);
  try {
    arma::mat L = arma::chol(P_arma);
    (void)L;
  } catch (...) {
    Rcpp::stop("%s must be positive definite.", what);
  }
}

List prior_list_to_P_Sigma(List pl) {
  if (!has_non_null(pl, "mu")) {
    Rcpp::stop("prior_list must contain 'mu'.");
  }
  NumericVector mu = pl["mu"];
  List out = List::create(Rcpp::Named("mu") = mu);

  if (has_non_null(pl, "P")) {
    NumericMatrix P = pl["P"];
    if (!is_symmetric_mat(P)) {
      Rcpp::stop("prior precision matrix P must be symmetric.");
    }
    check_P_pd(P, "P");
    NumericMatrix Sigma = Rcpp::wrap(arma::inv(arma::mat(
      const_cast<double*>(P.begin()), P.nrow(), P.ncol(), false
    )));
    out["P"] = P;
    out["Sigma"] = Sigma;
    return out;
  }
  if (has_non_null(pl, "Sigma")) {
    NumericMatrix Sigma = pl["Sigma"];
    if (!is_symmetric_mat(Sigma)) {
      Rcpp::stop("prior covariance Sigma must be symmetric.");
    }
    arma::mat S(const_cast<double*>(Sigma.begin()), Sigma.nrow(), Sigma.ncol(), false);
    arma::mat P = arma::inv_sympd(S);
    out["Sigma"] = Sigma;
    out["P"] = NumericMatrix(Rcpp::wrap(0.5 * (P + P.t())));
    return out;
  }
  Rcpp::stop("prior_list must contain 'P' or 'Sigma'.");
}

List base_prior_block(List pl, int l1) {
  List ps = prior_list_to_P_Sigma(pl);
  NumericVector mu = ps["mu"];
  if (mu.size() != l1) {
    Rcpp::stop("length(mu) must equal ncol(x) (%d).", l1);
  }
  NumericMatrix P = ps["P"];
  if (P.nrow() != l1 || P.ncol() != l1) {
    Rcpp::stop("dim(P) or dim(Sigma) must be %d x %d.", l1, l1);
  }
  return List::create(
    Rcpp::Named("mu") = mu,
    Rcpp::Named("Sigma") = ps["Sigma"],
    Rcpp::Named("P") = P
  );
}

List rep_list_blocks(const List& one, int k) {
  List out(k);
  for (int j = 0; j < k; ++j) {
    out[j] = Rcpp::clone(one);
  }
  return out;
}

List normalize_prior_for_block_ing(
    SEXP prior_list_sexp,
    SEXP prior_lists_sexp,
    const List& group_info,
    int l1
) {
  const int k = group_info["k"];

  if (!Rf_isNull(prior_lists_sexp)) {
    List prior_lists(prior_lists_sexp);
    if (!Rcpp::is<List>(prior_lists)) {
      Rcpp::stop("'prior_lists' must be a list.");
    }
    if (prior_lists.size() == 1) {
      List one = base_prior_block(List(prior_lists[0]), l1);
      return rep_list_blocks(one, k);
    }
    if (prior_lists.size() != static_cast<R_xlen_t>(k)) {
      Rcpp::stop("'prior_lists' must have length 1 or k = %d.", k);
    }
    List out(k);
    for (int j = 0; j < k; ++j) {
      out[j] = base_prior_block(List(prior_lists[j]), l1);
    }
    return out;
  }

  if (Rf_isNull(prior_list_sexp)) {
    Rcpp::stop("Provide 'prior_list' or 'prior_lists'.");
  }

  List prior_list(prior_list_sexp);

  if (has_non_null(prior_list, "dispersion")) {
    Rcpp::stop(
      "BlockEnvelopeCentering requires unknown measurement dispersion; "
      "prior_list must not contain a fixed 'dispersion'."
    );
  }

  SEXP mu_sexp = has_non_null(prior_list, "mu") ? (SEXP)prior_list["mu"] : R_NilValue;
  if (Rcpp::is<NumericMatrix>(mu_sexp)) {
    NumericMatrix mu_mat(mu_sexp);
    if (mu_mat.nrow() != l1) {
      Rcpp::stop("nrow(prior_list$mu) must equal ncol(x) (%d).", l1);
    }
    if (mu_mat.ncol() == 1) {
      List pl_copy = Rcpp::clone(prior_list);
      pl_copy["mu"] = mu_mat(Rcpp::_, 0);
      List one = base_prior_block(pl_copy, l1);
      return rep_list_blocks(one, k);
    }
    if (mu_mat.ncol() != k) {
      Rcpp::stop("ncol(prior_list$mu) must equal number of blocks k = %d.", k);
    }

    List P_list;
    List Sigma_list;
    bool have_P_list = false;
    bool have_S_list = false;

    if (has_non_null(prior_list, "P")) {
      SEXP P_sexp = prior_list["P"];
      if (Rcpp::is<List>(P_sexp)) {
        P_list = List(P_sexp);
        have_P_list = P_list.size() == 1 || P_list.size() == static_cast<R_xlen_t>(k);
      } else if (Rcpp::is<NumericMatrix>(P_sexp)) {
        P_list = List::create(P_sexp);
        have_P_list = true;
      }
    }
    if (has_non_null(prior_list, "Sigma")) {
      SEXP S_sexp = prior_list["Sigma"];
      if (Rcpp::is<List>(S_sexp)) {
        Sigma_list = List(S_sexp);
        have_S_list = Sigma_list.size() == 1 || Sigma_list.size() == static_cast<R_xlen_t>(k);
      } else if (Rcpp::is<NumericMatrix>(S_sexp)) {
        Sigma_list = List::create(S_sexp);
        have_S_list = true;
      }
    }

    List out(k);
    for (int j = 0; j < k; ++j) {
      List pl_j = List::create(Rcpp::Named("mu") = mu_mat(Rcpp::_, j));
      if (have_P_list) {
        pl_j["P"] = P_list[std::min(j, static_cast<int>(P_list.size()) - 1)];
      } else if (have_S_list) {
        pl_j["Sigma"] = Sigma_list[std::min(j, static_cast<int>(Sigma_list.size()) - 1)];
      } else {
        Rcpp::stop("prior_list must contain 'P' or 'Sigma'.");
      }
      out[j] = base_prior_block(pl_j, l1);
    }
    return out;
  }

  List one = base_prior_block(prior_list, l1);
  return rep_list_blocks(one, k);
}

inline int check_index_1based(int idx, int n, const char* what) {
  if (idx < 1 || idx > n) {
    Rcpp::stop("%s index %d out of range [1, %d]", what, idx, n);
  }
  return idx - 1;
}

NumericVector slice_numeric(const NumericVector& v, const IntegerVector& rows) {
  const int m = rows.size();
  NumericVector out(m);
  for (int i = 0; i < m; ++i) {
    out[i] = v[check_index_1based(rows[i], v.size(), "row")];
  }
  return out;
}

NumericMatrix slice_matrix_rows(const NumericMatrix& x, const IntegerVector& rows) {
  const int m = rows.size();
  const int l1 = x.ncol();
  NumericMatrix out(m, l1);
  for (int i = 0; i < m; ++i) {
    const int r = check_index_1based(rows[i], x.nrow(), "row");
    for (int j = 0; j < l1; ++j) {
      out(i, j) = x(r, j);
    }
  }
  return out;
}

double sum_wt(const NumericVector& wt) {
  double s = 0.0;
  for (int i = 0; i < wt.size(); ++i) s += wt[i];
  return s;
}

struct BlockWlsInit {
  int rank;
  double rss_wls;
  double n_w_j;
  double dispersion_wls;
  bool identifiable;
};

BlockWlsInit block_wls_init(
    const NumericVector& y_j,
    const NumericMatrix& x_j,
    const NumericVector& offset_j,
    const NumericVector& wt_j
) {
  BlockWlsInit out;
  out.n_w_j = sum_wt(wt_j);

  const int n_j = y_j.size();
  NumericVector ystar(n_j);
  for (int i = 0; i < n_j; ++i) {
    ystar[i] = y_j[i] - offset_j[i];
  }

  Rcpp::Function lm_wfit("lm.wfit");
  List fit = lm_wfit(
    Rcpp::_["x"] = x_j,
    Rcpp::_["y"] = ystar,
    Rcpp::_["w"] = wt_j
  );

  NumericVector res = fit["residuals"];
  out.rss_wls = 0.0;
  for (int i = 0; i < res.size(); ++i) {
    out.rss_wls += res[i] * res[i];
  }
  out.rank = Rcpp::as<int>(fit["rank"]);
  out.identifiable = out.rank >= 1;
  const double denom = std::max(out.n_w_j - static_cast<double>(out.rank), 1.0);
  out.dispersion_wls = out.rss_wls / denom;
  return out;
}

List block_envelope_centering_one(
    const NumericVector& y_j,
    const NumericMatrix& x_j,
    const NumericVector& mu_j,
    const NumericMatrix& P_j,
    const NumericVector& offset_j,
    const NumericVector& wt_j,
    double dispersion2,
    int rank_j,
    bool identifiable
) {
  const int l1 = x_j.ncol();
  const int l2_j = x_j.nrow();
  const double n_w_j = sum_wt(wt_j);

  if (!identifiable || rank_j < 1) {
    return List::create(
      Rcpp::Named("RSS_post") = 0.0,
      Rcpp::Named("b_post_mean") = mu_j,
      Rcpp::Named("rank") = rank_j,
      Rcpp::Named("n_w_j") = n_w_j,
      Rcpp::Named("identifiable") = false
    );
  }

  const arma::mat X = Rcpp::as<arma::mat>(x_j);
  const arma::vec Y = Rcpp::as<arma::vec>(y_j);
  const arma::vec off = Rcpp::as<arma::vec>(offset_j);
  const arma::vec wv = Rcpp::as<arma::vec>(wt_j);
  const arma::vec mu_vec = Rcpp::as<arma::vec>(mu_j);
  const arma::mat P_arma = Rcpp::as<arma::mat>(P_j);

  arma::mat Xw(l2_j, l1);
  arma::vec yw(l2_j);
  for (int i = 0; i < l2_j; ++i) {
    const double sw = std::sqrt(wv[i]);
    Xw.row(i) = sw * X.row(i);
    yw[i] = (Y[i] - off[i]) * sw;
  }

  const arma::mat RA = arma::chol(P_arma);
  const arma::vec z_bot = RA * mu_vec;
  const arma::mat XtWX = X.t() * (arma::diagmat(wv) * X);

  const double s = 1.0 / std::sqrt(dispersion2);
  arma::mat W(l2_j + l1, l1);
  arma::vec z(l2_j + l1);
  W.rows(0, l2_j - 1) = s * Xw;
  W.rows(l2_j, l2_j + l1 - 1) = RA;
  z.rows(0, l2_j - 1) = s * yw;
  z.rows(l2_j, l2_j + l1 - 1) = z_bot;

  const arma::mat WtW = W.t() * W;
  const arma::mat IR = arma::inv(arma::trimatu(arma::chol(WtW)));
  const arma::mat Sigma = IR * arma::trans(IR);
  const arma::vec b2_fast = Sigma * (W.t() * z);

  const arma::vec r_fast = Y - X * b2_fast - off;
  const double rss_at_mean = arma::dot(wv, r_fast % r_fast);
  const double trace_term = arma::trace(XtWX * Sigma);
  const double RSS_post = rss_at_mean + trace_term;

  return List::create(
    Rcpp::Named("RSS_post") = RSS_post,
    Rcpp::Named("b_post_mean") = NumericVector(b2_fast.begin(), b2_fast.end()),
    Rcpp::Named("rank") = rank_j,
    Rcpp::Named("n_w_j") = n_w_j,
    Rcpp::Named("identifiable") = true
  );
}

Nullable<double> resolve_nullable_bound(
    Nullable<double> arg_val,
    const List& prior_list,
    const char* name
) {
  if (arg_val.isSet() && arg_val.isNotNull()) {
    return arg_val;
  }
  if (has_non_null(prior_list, name)) {
    NumericVector v = prior_list[name];
    if (v.size() >= 1) {
      return Nullable<double>(Rcpp::wrap(static_cast<double>(v[0])));
    }
  }
  return R_NilValue;
}

List block_envelope_build_one(
    const NumericVector& y_j,
    const NumericMatrix& x_j,
    const NumericVector& mu_j,
    const NumericMatrix& P_j,
    const NumericVector& offset_j,
    const NumericVector& wt_j,
    const NumericVector& b_post_mean,
    double dispersion2,
    bool identifiable,
    const std::string& block_id,
    Function& optim,
    Function& f2,
    Function& f3,
    int Gridtype,
    int n,
    int n_envopt_val,
    bool use_opencl,
    bool verbose
) {
  const int l1 = x_j.ncol();
  BEB_DBG(verbose, "[BEB 4.0] block_envelope_build_one id=" << block_id
          << " identifiable=" << identifiable
          << " l1=" << l1 << " n_j=" << y_j.size());

  if (!identifiable) {
    BEB_DBG(verbose, "[BEB 4.1] prior-only stub (not identifiable)");
    return List::create(
      Rcpp::Named("block_envelope") = R_NilValue,
      Rcpp::Named("block_standardization") = List::create(
        Rcpp::Named("mu") = mu_j,
        Rcpp::Named("P") = P_j,
        Rcpp::Named("prior_only") = true
      ),
      Rcpp::Named("block_id") = block_id,
      Rcpp::Named("identifiable") = false
    );
  }

  const int n_j = y_j.size();
  NumericVector wt2(n_j);
  for (int i = 0; i < n_j; ++i) {
    wt2[i] = wt_j[i] / dispersion2;
  }

  const arma::mat X = Rcpp::as<arma::mat>(x_j);
  const arma::vec alpha_vec =
    X * Rcpp::as<arma::vec>(mu_j) + Rcpp::as<arma::vec>(offset_j);
  const NumericVector alpha = Rcpp::wrap(alpha_vec);

  NumericVector mu2(l1);
  std::fill(mu2.begin(), mu2.end(), 0.0);

  NumericVector parin(l1);
  for (int i = 0; i < l1; ++i) {
    parin[i] = b_post_mean[i];
  }

  BEB_DBG(verbose, "[BEB 4.2] before optim block_id=" << block_id);
  List opt_out = optim(
    Rcpp::_["par"] = parin,
    Rcpp::_["fn"] = f2,
    Rcpp::_["gr"] = f3,
    Rcpp::_["y"] = y_j,
    Rcpp::_["x"] = x_j,
    Rcpp::_["mu"] = mu2,
    Rcpp::_["P"] = P_j,
    Rcpp::_["alpha"] = alpha,
    Rcpp::_["wt"] = wt2,
    Rcpp::_["method"] = "BFGS",
    Rcpp::_["hessian"] = true
  );

  NumericVector bstar = opt_out["par"];
  NumericMatrix A1 = opt_out["hessian"];
  BEB_DBG(verbose, "[BEB 4.3] after optim block_id=" << block_id
          << " bstar_len=" << bstar.size());

  NumericMatrix bstar_mat(l1, 1);
  for (int i = 0; i < l1; ++i) {
    bstar_mat(i, 0) = bstar[i];
  }

  NumericMatrix x2_mat = x_j;
  NumericMatrix P2_mat = P_j;

  BEB_DBG(verbose, "[BEB 4.4] before glmb_Standardize_Model block_id=" << block_id);
  List Standard_Mod = glmb_Standardize_Model(
    y_j, x2_mat, P2_mat, bstar_mat, A1
  );
  BEB_DBG(verbose, "[BEB 4.5] after glmb_Standardize_Model block_id=" << block_id);

  NumericVector bstar2 = Standard_Mod["bstar2"];
  NumericMatrix A = Standard_Mod["A"];
  NumericMatrix x2_std = Standard_Mod["x2"];
  NumericMatrix mu2_std = Standard_Mod["mu2"];
  NumericMatrix P2_std = Standard_Mod["P2"];

  BEB_DBG(verbose, "[BEB 4.6] before EnvelopeBuild block_id=" << block_id
          << " Gridtype=" << Gridtype << " n=" << n);
  List Env2 = EnvelopeBuild(
    bstar2,
    A,
    y_j,
    x2_std,
    mu2_std,
    P2_std,
    alpha,
    wt2,
    "gaussian",
    "identity",
    Gridtype,
    n,
    n_envopt_val,
    false,
    use_opencl,
    verbose
  );
  BEB_DBG(verbose, "[BEB 4.7] after EnvelopeBuild block_id=" << block_id);

  List block_standardization = List::create(
    Rcpp::Named("bstar") = bstar,
    Rcpp::Named("bstar2") = bstar2,
    Rcpp::Named("A") = A,
    Rcpp::Named("x2") = x2_std,
    Rcpp::Named("P2") = P2_std,
    Rcpp::Named("mu2") = mu2_std,
    Rcpp::Named("L2Inv") = Standard_Mod["L2Inv"],
    Rcpp::Named("L3Inv") = Standard_Mod["L3Inv"],
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("y") = y_j,
    Rcpp::Named("wt") = wt_j,
    Rcpp::Named("mu") = mu_j,
    Rcpp::Named("P") = P_j,
    Rcpp::Named("prior_only") = false
  );

  return List::create(
    Rcpp::Named("block_envelope") = Env2,
    Rcpp::Named("block_standardization") = block_standardization,
    Rcpp::Named("block_id") = block_id,
    Rcpp::Named("identifiable") = true
  );
}

void compute_sigma2_bounds_cpp(
    double shape2,
    double rate3,
    double max_disp_perc,
    Nullable<double> disp_lower,
    Nullable<double> disp_upper,
    double* low_out,
    double* upp_out
) {
  if (!disp_lower.isUsable() || !disp_upper.isUsable()) {
    Function qgamma("qgamma");
    NumericVector q_low = qgamma(
      Rcpp::_["p"] = max_disp_perc,
      Rcpp::_["shape"] = shape2,
      Rcpp::_["rate"] = rate3
    );
    NumericVector q_upp = qgamma(
      Rcpp::_["p"] = 1.0 - max_disp_perc,
      Rcpp::_["shape"] = shape2,
      Rcpp::_["rate"] = rate3
    );
    *low_out = 1.0 / static_cast<double>(q_low[0]);
    *upp_out = 1.0 / static_cast<double>(q_upp[0]);
  } else {
    *low_out = Rcpp::as<double>(disp_lower.get());
    *upp_out = Rcpp::as<double>(disp_upper.get());
    if (!R_finite(*low_out) || !R_finite(*upp_out)) {
      Rcpp::stop("disp_lower/disp_upper must be finite.");
    }
    if (*low_out <= 0.0 || *upp_out <= 0.0) {
      Rcpp::stop("disp_lower/disp_upper must be positive.");
    }
    if (*upp_out <= *low_out) {
      Rcpp::stop("disp_upper must be strictly greater than disp_lower.");
    }
  }
}

List block_envelope_dispersion_one(
    List Env_j,
    double shape,
    double rate,
    const NumericMatrix& P2_j,
    const NumericVector& y_j,
    const NumericMatrix& x2_j,
    const NumericMatrix& mu2_j,
    const NumericVector& alpha_j,
    const NumericVector& wt_j,
    double RSS_post_j,
    double RSS_ML_j,
    double max_disp_perc,
    Nullable<double> disp_lower,
    Nullable<double> disp_upper,
    bool use_parallel,
    bool verbose,
    const std::string& block_id
) {
  const int n_obs = y_j.size();

  List edb = EnvelopeDispersionBuild(
    Env_j,
    shape,
    rate,
    P2_j,
    y_j,
    x2_j,
    alpha_j,
    n_obs,
    RSS_post_j,
    RSS_ML_j,
    mu2_j,
    wt_j,
    max_disp_perc,
    disp_lower,
    disp_upper,
    verbose,
    use_parallel
  );

  List Env_out = edb["Env_out"];
  List gamma_list_j = edb["gamma_list"];
  List UB_list_j = edb["UB_list"];
  List diagnostics_j = edb["diagnostics"];

  NumericMatrix cbars = Env_out["cbars"];
  List cache_j = Inv_f3_precompute_disp(
    cbars, y_j, x2_j, mu2_j, P2_j, alpha_j, wt_j
  );
  const double low_j = Rcpp::as<double>(gamma_list_j["disp_lower"]);
  const double upp_j = Rcpp::as<double>(gamma_list_j["disp_upper"]);
  add_ub2_bound_matrices_to_cache(cache_j, low_j, upp_j);

  const double rss_min_global_j = Rcpp::as<double>(UB_list_j["RSS_Min"]);
  double rss_ml_j = Rcpp::as<double>(UB_list_j["RSS_ML"]);
  if (!R_finite(RSS_ML_j)) {
    rss_ml_j = Rcpp::as<double>(UB_list_j["RSS_ML"]);
  } else {
    rss_ml_j = RSS_ML_j;
  }

  if (verbose) {
    Rcpp::Rcout << "BlockEnvelopeDispersionBuild block id=" << block_id
                << " gs=" << cbars.nrow() << std::endl;
  }

  return List::create(
    Rcpp::Named("Env_out") = Env_out,
    Rcpp::Named("gamma_list") = gamma_list_j,
    Rcpp::Named("UB_list") = UB_list_j,
    Rcpp::Named("diagnostics") = diagnostics_j,
    Rcpp::Named("cache") = cache_j,
    Rcpp::Named("lg_prob_factor") = UB_list_j["lg_prob_factor"],
    Rcpp::Named("UB2min") = UB_list_j["UB2min"],
    Rcpp::Named("rss_min_global") = rss_min_global_j,
    Rcpp::Named("RSS_ML") = rss_ml_j,
    Rcpp::Named("block_id") = block_id
  );
}

struct BlockFaceGeom {
  NumericVector slope;
  NumericVector upp_apprx;
  NumericVector low_apprx;
  NumericVector logP1;
  NumericVector ub2_min;
  int gs = 0;
};

// Per-block arrays hoisted once before the O(prod gs_t) product-face loop in
// build_joint_product_face_slack (avoids repeated SEXP/List indexing and
// recomputing ||cbar||^2 for every occurrence of the same block face).
struct BlockJointSlackCache {
  NumericVector ub2_at_low;
  NumericVector ub2_at_upp;
  NumericVector logP1;
  NumericVector norm2;
  int gs = 0;
};

NumericVector cbars_row_norm2_sq(const NumericMatrix& cbars) {
  const int gs = cbars.nrow();
  const int p = cbars.ncol();
  NumericVector out(gs);
  for (int f = 0; f < gs; ++f) {
    double s = 0.0;
    for (int c = 0; c < p; ++c) {
      const double cjk = cbars(f, c);
      s += cjk * cjk;
    }
    out[f] = s;
  }
  return out;
}

// max_upp and mean_slope are taken as pre-reduced scalars rather than being
// recomputed here from the full length-gs_total product-face arrays. This is
// not an approximation: for a product face flat = (j1,...,jk), joint_upp_apprx
// and joint_slope are sums of independent per-block terms (no cross terms,
// see build_joint_face_product_geometry), so
//   max over product faces of a separable sum = sum of per-block maxes, and
//   mean over product faces of a separable sum = sum of per-block means.
// The caller therefore reduces each block's own (small) upp_apprx/slope
// vectors once and sums the results, which is exactly equal to scanning the
// full O(prod gs_t) joint array but costs only O(sum gs_t).
List ub3_geometry_from_joint_faces(
    double max_upp,
    double mean_slope,
    const NumericVector& joint_low_apprx,
    const NumericVector& joint_upp_apprx,
    double n_w_global,
    double low,
    double upp
) {
  const double lmc2_max =
    (n_w_global / 2.0) * (std::log(upp) - std::log(low)) / (upp - low);
  if (mean_slope > lmc2_max) {
    Rcpp::Rcout
      << "[BlockEnvelopeDispersionBuild] UB3A mean slope (" << mean_slope
      << ") exceeds UB3B-compatible maximum (" << lmc2_max << ").\n"
      << "          Capping global slope to preserve lm_log2 <= n_w/2 "
      << "and ensure shape3 > 0.\n";
  }
  const double m_New_LL_Slope = std::min(mean_slope, lmc2_max);
  const double max_low_mean = max_upp - m_New_LL_Slope * (upp - low);
  const double new_slope = (max_upp - max_low_mean) / (upp - low);
  const double new_int = max_low_mean - new_slope * low;
  const double d2_star = (upp - low) / (std::log(upp / low));

  return List::create(
    Rcpp::Named("thetabar_const_low_apprx") = joint_low_apprx,
    Rcpp::Named("thetabar_const_upp_apprx") = joint_upp_apprx,
    Rcpp::Named("max_low") = max_upp - mean_slope * (upp - low),
    Rcpp::Named("max_upp") = max_upp,
    Rcpp::Named("new_slope") = new_slope,
    Rcpp::Named("new_int") = new_int,
    Rcpp::Named("d2_star") = d2_star
  );
}

List gamma_ub_from_joint_geometry(
    const List& geom,
    double shape2_global,
    double rate_prior,
    double RSS_Min_global,
    double RSS_ML_global,
    double low,
    double upp
) {
  const double max_low = Rcpp::as<double>(geom["max_low"]);
  const double max_upp = Rcpp::as<double>(geom["max_upp"]);
  const double new_slope = Rcpp::as<double>(geom["new_slope"]);
  const double new_int = Rcpp::as<double>(geom["new_int"]);
  const double d2_star = Rcpp::as<double>(geom["d2_star"]);

  const double lm_log2 = new_slope * d2_star;
  const double lm_log1 =
    new_int + new_slope * d2_star - new_slope * std::log(d2_star);
  const double shape3 = shape2_global - lm_log2;
  const double rate2 = rate_prior + RSS_Min_global / 2.0;

  if (!R_finite(lm_log1) || !R_finite(lm_log2)) {
    Rcpp::stop(
      "BlockEnvelopeDispersionBuild: joint lm_log1/lm_log2 non-finite."
    );
  }
  if (!R_finite(shape3) || shape3 <= 0.0) {
    Rcpp::stop(
      "BlockEnvelopeDispersionBuild: joint shape3 <= 0; invalid tilted IG."
    );
  }
  if (!R_finite(rate2) || rate2 <= 0.0) {
    Rcpp::stop(
      "BlockEnvelopeDispersionBuild: joint rate2 <= 0; invalid tilted IG."
    );
  }

  List gamma_list = List::create(
    Rcpp::Named("shape3") = shape3,
    Rcpp::Named("rate2") = rate2,
    Rcpp::Named("disp_lower") = low,
    Rcpp::Named("disp_upper") = upp
  );

  List UB_list = List::create(
    Rcpp::Named("RSS_ML") = RSS_ML_global,
    Rcpp::Named("RSS_Min") = RSS_Min_global,
    Rcpp::Named("max_New_LL_UB") = max_upp,
    Rcpp::Named("max_LL_log_disp") = lm_log1 + lm_log2 * std::log(upp),
    Rcpp::Named("lm_log1") = lm_log1,
    Rcpp::Named("lm_log2") = lm_log2,
    Rcpp::Named("lmc1") = new_int,
    Rcpp::Named("lmc2") = new_slope
  );

  return List::create(
    Rcpp::Named("gamma_list") = gamma_list,
    Rcpp::Named("UB_list") = UB_list,
    Rcpp::Named("max_low") = max_low,
    Rcpp::Named("max_upp") = max_upp
  );
}

// O(sum gs_t) scalar-only reduction of global UB3 geometry for the Ind path.
// Exact by the same additivity argument as ub3_geometry_from_joint_faces:
//   max(sum_t f_t(j_t)) = sum_t max_t f_t  and  mean = sum_t mean_t
// for separable sums over product faces.  The Ind path does not need the
// full joint_upp_apprx / joint_low_apprx arrays, so we pass empty vectors.
List build_global_dispersion_constants_scalar(
    int n_identifiable,
    const List& single_gamma_list,
    const List& single_ub_list,
    const IntegerVector& identifiable_idx,
    const List& block_envelopes,
    const List& block_standardization,
    double shape2_global,
    double rate_prior,
    double rate3_global,
    double n_w_global,
    double RSS_Min_global,
    double RSS_ML_global,
    double low,
    double upp,
    double RSS_post_global
) {
  if (n_identifiable == 1) {
    return List::create(
      Rcpp::Named("gamma_list") = single_gamma_list,
      Rcpp::Named("UB_list") = single_ub_list,
      Rcpp::Named("source") = "single_block_edb"
    );
  }

  const int n_blocks = identifiable_idx.size();
  const double d1_star = rate3_global / (shape2_global - 1.0);

  double max_upp = 0.0;
  double mean_slope = 0.0;
  for (int t = 0; t < n_blocks; ++t) {
    const int j = identifiable_idx[t];
    List Env = block_envelopes[j];
    List std_j = block_standardization[j];
    NumericMatrix cbars = Env["cbars"];
    NumericMatrix thetabars = Env["thetabars"];
    NumericMatrix x2 = std_j["x2"];
    NumericMatrix P2 = std_j["P2"];
    NumericVector alpha = std_j["alpha"];
    NumericVector y = std_j["y"];

    NumericVector slope = EnvBuildLinBound_cpp(
      thetabars, cbars, y, x2, P2, alpha, d1_star
    );
    NumericVector const_base = thetabar_const_cpp(P2, cbars, thetabars);
    const int gs = slope.size();
    double max_upp_t = R_NegInf;
    for (int f = 0; f < gs; ++f) {
      double v = const_base[f] + (upp - d1_star) * slope[f];
      if (v > max_upp_t) max_upp_t = v;
    }
    max_upp += max_upp_t;
    mean_slope += static_cast<double>(Rcpp::mean(slope));
  }

  List geom = ub3_geometry_from_joint_faces(
    max_upp, mean_slope,
    NumericVector(0), NumericVector(0),
    n_w_global, low, upp
  );

  List out = gamma_ub_from_joint_geometry(
    geom, shape2_global, rate_prior,
    RSS_Min_global, RSS_ML_global, low, upp
  );

  const double max_low = max_upp - mean_slope * (upp - low);
  return List::create(
    Rcpp::Named("gamma_list") = out["gamma_list"],
    Rcpp::Named("UB_list") = out["UB_list"],
    Rcpp::Named("source") = "joint_face_scalar_edb",
    Rcpp::Named("RSS_post") = RSS_post_global,
    Rcpp::Named("prob_max_low") = max_low
  );
}

List build_joint_face_product_geometry(
    const IntegerVector& identifiable_idx,
    const List& block_envelopes,
    const List& block_standardization,
    const List& block_dispersion,
    double shape2_global,
    double rate3_global,
    double n_w_global,
    double low,
    double upp
) {
  const int n_blocks = identifiable_idx.size();
  const double d1_star = rate3_global / (shape2_global - 1.0);

  std::vector<BlockFaceGeom> block_geom(static_cast<size_t>(n_blocks));
  int gs_total = 1;
  for (int t = 0; t < n_blocks; ++t) {
    const int j = identifiable_idx[t];
    List Env = block_envelopes[j];
    List std_j = block_standardization[j];
    List disp_j = block_dispersion[j];

    NumericMatrix cbars = Env["cbars"];
    NumericMatrix thetabars = Env["thetabars"];
    NumericVector logP1 = Env["logP"];
    NumericVector ub2_min = disp_j["UB2min"];
    NumericMatrix x2 = std_j["x2"];
    NumericMatrix P2 = std_j["P2"];
    NumericVector alpha = std_j["alpha"];
    NumericVector y = std_j["y"];

    BlockFaceGeom bg;
    bg.gs = cbars.nrow();
    bg.slope = EnvBuildLinBound_cpp(
      thetabars, cbars, y, x2, P2, alpha, d1_star
    );
    NumericVector const_base = thetabar_const_cpp(P2, cbars, thetabars);
    bg.upp_apprx = NumericVector(bg.gs);
    bg.low_apprx = NumericVector(bg.gs);
    for (int f = 0; f < bg.gs; ++f) {
      bg.upp_apprx[f] = const_base[f] + (upp - d1_star) * bg.slope[f];
      bg.low_apprx[f] = const_base[f] + (low - d1_star) * bg.slope[f];
    }
    bg.logP1 = logP1;
    bg.ub2_min = ub2_min;
    block_geom[static_cast<size_t>(t)] = bg;
    gs_total *= bg.gs;
  }

  // max_upp / mean_slope: O(sum gs_t) additive shortcut (see comment on
  // ub3_geometry_from_joint_faces). Avoids ever materializing a length-
  // gs_total slope array just to reduce it, which is the only thing
  // joint_slope was used for.
  double max_upp = 0.0;
  double mean_slope = 0.0;
  for (int t = 0; t < n_blocks; ++t) {
    const BlockFaceGeom& bg = block_geom[static_cast<size_t>(t)];
    max_upp += max_vec_local(bg.upp_apprx);
    mean_slope += static_cast<double>(Rcpp::mean(bg.slope));
  }

  NumericVector joint_upp_apprx(gs_total);
  NumericVector joint_low_apprx(gs_total);
  std::vector<int> face_idx(static_cast<size_t>(n_blocks), 0);

  for (int flat = 0; flat < gs_total; ++flat) {
    double upp_sum = 0.0;
    double low_sum = 0.0;
    for (int t = 0; t < n_blocks; ++t) {
      const BlockFaceGeom& bg = block_geom[static_cast<size_t>(t)];
      const int f = face_idx[static_cast<size_t>(t)];
      upp_sum += bg.upp_apprx[f];
      low_sum += bg.low_apprx[f];
    }
    joint_upp_apprx[flat] = upp_sum;
    joint_low_apprx[flat] = low_sum;

    for (int t = n_blocks - 1; t >= 0; --t) {
      ++face_idx[static_cast<size_t>(t)];
      if (face_idx[static_cast<size_t>(t)] < block_geom[static_cast<size_t>(t)].gs) {
        break;
      }
      face_idx[static_cast<size_t>(t)] = 0;
    }
  }

  List geom = ub3_geometry_from_joint_faces(
    max_upp, mean_slope, joint_low_apprx, joint_upp_apprx, n_w_global, low, upp
  );
  geom["joint_upp_apprx"] = joint_upp_apprx;
  geom["joint_low_apprx"] = joint_low_apprx;
  return geom;
}

List build_global_dispersion_constants(
    int n_identifiable,
    const List& single_gamma_list,
    const List& single_ub_list,
    const IntegerVector& identifiable_idx,
    const List& block_envelopes,
    const List& block_standardization,
    const List& block_dispersion,
    double shape2_global,
    double rate_prior,
    double rate3_global,
    double n_w_global,
    double RSS_Min_global,
    double RSS_ML_global,
    double low,
    double upp,
    double RSS_post_global
) {
  if (n_identifiable == 1) {
    return List::create(
      Rcpp::Named("gamma_list") = single_gamma_list,
      Rcpp::Named("UB_list") = single_ub_list,
      Rcpp::Named("source") = "single_block_edb"
    );
  }

  List geom = build_joint_face_product_geometry(
    identifiable_idx,
    block_envelopes,
    block_standardization,
    block_dispersion,
    shape2_global,
    rate3_global,
    n_w_global,
    low,
    upp
  );

  List out = gamma_ub_from_joint_geometry(
    geom,
    shape2_global,
    rate_prior,
    RSS_Min_global,
    RSS_ML_global,
    low,
    upp
  );

  return List::create(
    Rcpp::Named("gamma_list") = out["gamma_list"],
    Rcpp::Named("UB_list") = out["UB_list"],
    Rcpp::Named("source") = "joint_face_product_edb",
    Rcpp::Named("RSS_post") = RSS_post_global,
    Rcpp::Named("prob_max_low") = geom["max_low"],
    Rcpp::Named("joint_upp_apprx") = geom["joint_upp_apprx"],
    Rcpp::Named("joint_low_apprx") = geom["joint_low_apprx"]
  );
}

// Per-block upp/low extrapolations at global d1_star + joint prob_factor anchors
// (max_upp, max_low) via additive shortcuts — no product-face loop.
List compute_block_face_apprx_and_prob_anchors(
    const IntegerVector& identifiable_idx,
    const List& block_envelopes,
    const List& block_standardization,
    const List& block_dispersion,
    double shape2_global,
    double rate3_global,
    double low,
    double upp
) {
  const int n_blocks = identifiable_idx.size();
  const double d1_star = rate3_global / (shape2_global - 1.0);
  double max_upp_prob = 0.0;
  double mean_slope_sum = 0.0;
  List block_apprx(n_blocks);

  for (int t = 0; t < n_blocks; ++t) {
    const int j = identifiable_idx[t];
    List Env = block_envelopes[j];
    List std_j = block_standardization[j];
    NumericMatrix cbars = Env["cbars"];
    NumericMatrix thetabars = Env["thetabars"];
    NumericMatrix x2 = std_j["x2"];
    NumericMatrix P2 = std_j["P2"];
    NumericVector alpha = std_j["alpha"];
    NumericVector y = std_j["y"];
    NumericVector wt = std_j["wt"];
    List disp_j = block_dispersion[j];
    List cache = disp_j["cache"];
    const double rss_min = Rcpp::as<double>(disp_j["rss_min_global"]);

    NumericVector slope = EnvBuildLinBound_cpp(
      thetabars, cbars, y, x2, P2, alpha, d1_star
    );
    NumericVector const_base = thetabar_const_cpp(P2, cbars, thetabars);
    NumericVector upp_apprx(slope.size());
    NumericVector low_apprx(slope.size());
    NumericVector ub2_at_low(slope.size());
    NumericVector ub2_at_upp(slope.size());
    for (int f = 0; f < slope.size(); ++f) {
      upp_apprx[f] = const_base[f] + (upp - d1_star) * slope[f];
      low_apprx[f] = const_base[f] + (low - d1_star) * slope[f];
      NumericVector cbars_j = cbars(f, _);
      ub2_at_low[f] = UB2(low, cache, cbars_j, y, x2, alpha, wt, rss_min);
      ub2_at_upp[f] = UB2(upp, cache, cbars_j, y, x2, alpha, wt, rss_min);
    }

    max_upp_prob += max_vec_local(upp_apprx);
    mean_slope_sum += static_cast<double>(Rcpp::mean(slope));

    block_apprx[t] = List::create(
      Rcpp::Named("block_index") = j,
      Rcpp::Named("upp_apprx") = upp_apprx,
      Rcpp::Named("low_apprx") = low_apprx,
      Rcpp::Named("ub2_at_low") = ub2_at_low,
      Rcpp::Named("ub2_at_upp") = ub2_at_upp
    );
  }

  const double max_low_prob = max_upp_prob - mean_slope_sum * (upp - low);

  return List::create(
    Rcpp::Named("block_apprx") = block_apprx,
    Rcpp::Named("prob_max_upp") = max_upp_prob,
    Rcpp::Named("prob_max_low") = max_low_prob
  );
}

void patch_block_dispersion_apprx(
    List& block_dispersion,
    const List& apprx_out
) {
  List block_apprx = apprx_out["block_apprx"];
  for (int t = 0; t < block_apprx.size(); ++t) {
    List one = block_apprx[t];
    const int j = Rcpp::as<int>(one["block_index"]);
    List bd = Rcpp::as<List>(block_dispersion[j]);
    bd = Rcpp::clone(bd);
    bd["upp_apprx"] = one["upp_apprx"];
    bd["low_apprx"] = one["low_apprx"];
    bd["ub2_at_low"] = one["ub2_at_low"];
    bd["ub2_at_upp"] = one["ub2_at_upp"];
    block_dispersion[j] = bd;
  }
}

int product_face_flat_index(
    const std::vector<int>& J_draw,
    const IntegerVector& identifiable_idx,
    const IntegerVector& gs_per_block
) {
  const int n_blocks = identifiable_idx.size();
  int flat = 0;
  int stride = 1;
  for (int t = n_blocks - 1; t >= 0; --t) {
    const int block_j = identifiable_idx[t];
    const int J = J_draw[static_cast<size_t>(block_j)];
    if (J < 0) {
      Rcpp::stop("product_face_flat_index: invalid face index for identifiable block.");
    }
    if (J >= gs_per_block[t]) {
      Rcpp::stop("product_face_flat_index: face index out of range for block.");
    }
    flat += J * stride;
    stride *= gs_per_block[t];
  }
  return flat;
}

void decode_product_face_flat_index(
    int flat,
    const IntegerVector& gs_per_block,
    std::vector<int>& face_idx_out
) {
  const int n_blocks = gs_per_block.size();
  face_idx_out.assign(static_cast<size_t>(n_blocks), 0);
  int rem = flat;
  for (int t = n_blocks - 1; t >= 0; --t) {
    const int gs = gs_per_block[t];
    if (gs <= 0) {
      Rcpp::stop("decode_product_face_flat_index: invalid gs_per_block entry.");
    }
    face_idx_out[static_cast<size_t>(t)] = rem % gs;
    rem /= gs;
  }
  if (rem != 0) {
    Rcpp::stop("decode_product_face_flat_index: flat index out of range.");
  }
}

List build_joint_product_face_slack(
    const IntegerVector& identifiable_idx,
    const List& block_envelopes,
    const List& block_dispersion,
    const NumericVector& joint_upp_apprx,
    const NumericVector& joint_low_apprx,
    double prob_max_upp,
    double prob_max_low
) {
  const int n_blocks = identifiable_idx.size();
  const int gs_total = joint_upp_apprx.size();

  std::vector<BlockJointSlackCache> block_cache(static_cast<size_t>(n_blocks));
  for (int t = 0; t < n_blocks; ++t) {
    const int block_j = identifiable_idx[t];
    List bd = block_dispersion[block_j];
    List Env = block_envelopes[block_j];
    NumericMatrix cbars = Env["cbars"];
    BlockJointSlackCache bc;
    bc.ub2_at_low = bd["ub2_at_low"];
    bc.ub2_at_upp = bd["ub2_at_upp"];
    bc.logP1 = Env["logP"];
    bc.norm2 = cbars_row_norm2_sq(cbars);
    bc.gs = bc.ub2_at_low.size();
    if (bc.ub2_at_upp.size() != bc.gs ||
        bc.logP1.size() != bc.gs ||
        bc.norm2.size() != bc.gs) {
      Rcpp::stop(
        "build_joint_product_face_slack: per-block face table length mismatch."
      );
    }
    block_cache[static_cast<size_t>(t)] = bc;
  }

  NumericVector joint_lg_prob_factor(gs_total);
  NumericVector joint_ub2min_product(gs_total);
  NumericVector joint_logw(gs_total);
  std::vector<int> face_idx(static_cast<size_t>(n_blocks), 0);

  for (int flat = 0; flat < gs_total; ++flat) {
    const double pf_upp = joint_upp_apprx[flat] - prob_max_upp;
    const double pf_low = joint_low_apprx[flat] - prob_max_low;
    joint_lg_prob_factor[flat] = (pf_upp > pf_low ? pf_upp : pf_low);

    double ub2_low_sum = 0.0;
    double ub2_upp_sum = 0.0;
    double logp_sum = 0.0;
    double norm2_sum = 0.0;
    for (int t = 0; t < n_blocks; ++t) {
      const BlockJointSlackCache& bc = block_cache[static_cast<size_t>(t)];
      const int f = face_idx[static_cast<size_t>(t)];
      ub2_low_sum += bc.ub2_at_low[f];
      ub2_upp_sum += bc.ub2_at_upp[f];
      logp_sum += bc.logP1[f];
      norm2_sum += bc.norm2[f];
    }
    joint_ub2min_product[flat] =
      (ub2_low_sum <= ub2_upp_sum ? ub2_low_sum : ub2_upp_sum);
    joint_logw[flat] = logp_sum + 0.5 * norm2_sum +
      (joint_lg_prob_factor[flat] - joint_ub2min_product[flat]);

    for (int t = n_blocks - 1; t >= 0; --t) {
      ++face_idx[static_cast<size_t>(t)];
      if (face_idx[static_cast<size_t>(t)] < block_cache[static_cast<size_t>(t)].gs) {
        break;
      }
      face_idx[static_cast<size_t>(t)] = 0;
    }
  }

  NumericVector joint_PLSD(gs_total);
  double max_logw = joint_logw[0];
  for (int flat = 1; flat < gs_total; ++flat) {
    if (joint_logw[flat] > max_logw) {
      max_logw = joint_logw[flat];
    }
  }
  double sumP = 0.0;
  for (int flat = 0; flat < gs_total; ++flat) {
    joint_PLSD[flat] = std::exp(joint_logw[flat] - max_logw);
    sumP += joint_PLSD[flat];
  }
  if (sumP <= 0.0 || !R_finite(sumP)) {
    Rcpp::stop("build_joint_product_face_slack: joint_PLSD normalization failed.");
  }
  for (int flat = 0; flat < gs_total; ++flat) {
    joint_PLSD[flat] /= sumP;
  }

  return List::create(
    Rcpp::Named("joint_lg_prob_factor") = joint_lg_prob_factor,
    Rcpp::Named("joint_ub2min_product") = joint_ub2min_product,
    Rcpp::Named("joint_PLSD") = joint_PLSD,
    Rcpp::Named("n_product_faces") = gs_total
  );
}

double joint_lg_prob_factor_at_draw(
    const NumericVector& joint_lg_prob_factor,
    const std::vector<int>& J_draw,
    const IntegerVector& identifiable_idx,
    const IntegerVector& gs_per_block
) {
  const int flat = product_face_flat_index(J_draw, identifiable_idx, gs_per_block);
  if (flat < 0 || flat >= joint_lg_prob_factor.size()) {
    Rcpp::stop("joint_lg_prob_factor_at_draw: product face index out of range.");
  }
  return joint_lg_prob_factor[flat];
}

double joint_ub2min_at_draw(
    const NumericVector& joint_ub2min_product,
    const std::vector<int>& J_draw,
    const IntegerVector& identifiable_idx,
    const IntegerVector& gs_per_block
) {
  const int flat = product_face_flat_index(J_draw, identifiable_idx, gs_per_block);
  if (flat < 0 || flat >= joint_ub2min_product.size()) {
    Rcpp::stop("joint_ub2min_at_draw: product face index out of range.");
  }
  return joint_ub2min_product[flat];
}

// g1_j(d) = -0.5 * theta_j(d)^T P theta_j(d) + c_j^T theta_j(d)  (matches rIndepNormalGammaReg_std)
double g1_face_at_disp(
    double dispersion,
    int j,
    const List& cache,
    const arma::mat& P2,
    const NumericMatrix& cbars
) {
  const int p = cbars.ncol();
  NumericMatrix cbars_small(1, p);
  for (int k = 0; k < p; ++k) {
    cbars_small(0, k) = cbars(j, k);
  }
  NumericMatrix cbars_small_t = Rcpp::transpose(cbars_small);
  arma::mat theta_mat = Inv_f3_with_disp(cache, dispersion, cbars_small_t);
  arma::rowvec theta_row = theta_mat.row(0);
  arma::vec theta = theta_row.t();
  arma::vec c_j(p);
  for (int k = 0; k < p; ++k) {
    c_j(k) = cbars(j, k);
  }
  return arma::as_scalar(-0.5 * theta.t() * P2 * theta + c_j.t() * theta);
}

// Cumulative sum of a (normalized) probability vector, precomputed once so
// that repeated draws from the same PLSD can use binary search instead of a
// linear scan (see draw_face_index_from_cdf below).
std::vector<double> build_face_cdf(const NumericVector& p) {
  std::vector<double> cdf(static_cast<size_t>(p.size()));
  double running = 0.0;
  for (int i = 0; i < p.size(); ++i) {
    running += p[i];
    cdf[static_cast<size_t>(i)] = running;
  }
  return cdf;
}

// O(log gs) equivalent of draw_face_index(PLSD) via binary search on a
// precomputed cumulative sum. Same semantics (smallest J with U <= CDF[J]),
// but avoids re-scanning the full face array on every resample-until-accept
// attempt, which matters because draw_face_index is called once per attempt
// (not once per build) and gs can be O(prod gs_t) for k > 1 blocks. Clamps to
// the last index if floating-point error leaves U fractionally above
// cdf.back() (~1), which is strictly safer than the unbounded linear scan.
int draw_face_index_from_cdf(const std::vector<double>& cdf) {
  const double U = runif_safe();
  const auto it = std::lower_bound(cdf.begin(), cdf.end(), U);
  if (it == cdf.end()) {
    return static_cast<int>(cdf.size()) - 1;
  }
  return static_cast<int>(it - cdf.begin());
}

NumericVector draw_beta_std_face(
    const List& Env,
    int J,
    int l1
) {
  NumericMatrix loglt = Env["loglt"];
  NumericMatrix logrt = Env["logrt"];
  NumericMatrix cbars = Env["cbars"];
  NumericVector beta_std(l1);
  for (int j = 0; j < l1; ++j) {
    beta_std[j] = rnorm_ct(logrt(J, j), loglt(J, j), -cbars(J, j), 1.0);
  }
  return beta_std;
}

void block_ar_accumulate_one(
    int J,
    const NumericVector& beta_std,
    double sigma2,
    const List& Env,
    const List& block_std,
    const List& block_disp,
    double* LL_sum,
    double* UB1_sum,
    double* quad_sum,
    double* ub3a_block_sum,
    bool include_block_lg_in_ub3a
) {
  NumericMatrix x2 = block_std["x2"];
  NumericMatrix mu2 = block_std["mu2"];
  NumericMatrix P2 = block_std["P2"];
  NumericVector alpha = block_std["alpha"];
  NumericVector y = block_std["y"];
  NumericVector wt = block_std["wt"];
  List cache = block_disp["cache"];
  NumericVector lg_prob_factor = block_disp["lg_prob_factor"];
  NumericMatrix cbars = Env["cbars"];

  const int l1 = beta_std.size();
  NumericVector wt2 = wt / sigma2;

  NumericMatrix out_mat(1, l1);
  for (int t = 0; t < l1; ++t) {
    out_mat(0, t) = beta_std[t];
  }
  const double LL_Test =
    -f2_gaussian(Rcpp::transpose(out_mat), y, x2, mu2, P2, alpha, wt2)[0];
  *LL_sum += LL_Test;

  NumericMatrix cbars_small(1, l1);
  for (int t = 0; t < l1; ++t) {
    cbars_small(0, t) = cbars(J, t);
  }
  arma::mat theta2 = Inv_f3_with_disp(cache, sigma2, Rcpp::transpose(cbars_small));
  NumericMatrix thetabars_new = Rcpp::wrap(theta2);
  const double LL_New2 =
    -f2_gaussian(Rcpp::transpose(thetabars_new), y, x2, mu2, P2, alpha, wt2)[0];

  arma::vec betadiff(l1);
  arma::vec cbars_j(l1);
  for (int t = 0; t < l1; ++t) {
    betadiff(t) = beta_std[t] - thetabars_new(0, t);
    cbars_j(t) = cbars(J, t);
  }
  *UB1_sum += LL_New2 - arma::as_scalar(cbars_j.t() * betadiff);

  NumericVector cbars_j_nv = cbars(J, _);
  *quad_sum += rss_face_at_disp(sigma2, cache, cbars_j_nv, y, x2, alpha, wt);

  arma::mat P2_arma(P2.begin(), l1, l1, false);
  const double g1j = g1_face_at_disp(sigma2, J, cache, P2_arma, cbars);
  if (include_block_lg_in_ub3a) {
    *ub3a_block_sum += lg_prob_factor[J] - g1j;
  } else {
    *ub3a_block_sum -= g1j;
  }
}

bool block_ar_check_signs(
    double test1,
    double UB1,
    double UB2,
    double UB3A,
    double UB3B,
    double test,
    double quad_sum,
    double RSS_Min,
    double UB2_raw,
    double UB2min_sum,
    bool verbose
) {
  bool bad = false;
  std::ostringstream msg;

  const double tol1 = 1e-9 * std::max(1.0, std::abs(UB1));
  if (test1 > tol1) {
    bad = true;
    msg << "Sign violation: test1 = " << test1 << " > 0\n";
  }
  if (UB2 < 0.0) {
    const double ratio = std::abs(UB2) / std::max(std::abs(test), 1e-15);
    if (ratio >= 1e-2) {
      bad = true;
      msg << "Sign violation: UB2 = " << UB2 << " < 0\n";
    } else if (ratio >= 1e-4 && verbose) {
      Rcpp::Rcout << "BlockEnvelopeSim warning [UB2]: UB2=" << UB2
                  << " ratio=" << ratio << "\n";
    }
  }
  if (UB3A < 0.0) {
    bad = true;
    msg << "Sign violation: UB3A = " << UB3A << " < 0\n";
  }
  if (UB3B < 0.0) {
    bad = true;
    msg << "Sign violation: UB3B = " << UB3B << " < 0\n";
  }
  if (bad && verbose) {
    Rcpp::Rcout << "BlockEnvelopeSim sign check:\n" << msg.str()
                << " quad_sum=" << quad_sum
                << " RSS_Min=" << RSS_Min
                << " UB2_raw=" << UB2_raw
                << " UB2min_sum=" << UB2min_sum
                << " test=" << test << "\n";
  }
  return bad;
}

namespace {
constexpr double kBlockArTimeLimitSec = 5.0;

inline double block_ar_elapsed_seconds(
    const std::chrono::steady_clock::time_point& t0
) {
  return std::chrono::duration<double>(
    std::chrono::steady_clock::now() - t0
  ).count();
}

inline void block_ar_warn_time_limit_accept(
    const char* fn_name,
    int draw_index,
    double iters,
    double elapsed_sec,
    double limit_sec
) {
  Rcpp::warning(
    "%s: draw %d accept-reject exceeded %.0f s (elapsed %.2f s); "
    "accepting last proposal after %.0f attempts.",
    fn_name,
    draw_index + 1,
    limit_sec,
    elapsed_sec,
    iters
  );
}
}  // namespace

// Chapter A05 steps 1–2: draw face J ~ PLSD, then beta_j ~ TN(loglt, logrt, -cbars).
NumericVector envelope_draw_beta_std_one(const List& Env, int l1) {
  NumericVector PLSD = Env["PLSD"];
  NumericMatrix cbars = Env["cbars"];
  NumericMatrix loglt = Env["loglt"];
  NumericMatrix logrt = Env["logrt"];

  double U = runif_safe();
  int J = 0;
  int a2 = 0;
  while (a2 == 0) {
    if (U <= PLSD[J]) {
      a2 = 1;
    } else {
      U -= PLSD[J];
      ++J;
    }
  }

  NumericVector beta_std(l1);
  for (int j = 0; j < l1; ++j) {
    beta_std[j] = rnorm_ct(logrt(J, j), loglt(J, j), -cbars(J, j), 1.0);
  }
  return beta_std;
}

NumericVector unstandardize_beta_one(
    const NumericVector& beta_std,
    const NumericMatrix& L2Inv,
    const NumericMatrix& L3Inv,
    const NumericVector& mu
) {
  const arma::vec b = Rcpp::as<arma::vec>(beta_std);
  const arma::mat L2 = Rcpp::as<arma::mat>(L2Inv);
  const arma::mat L3 = Rcpp::as<arma::mat>(L3Inv);
  const arma::vec mu_v = Rcpp::as<arma::vec>(mu);
  const arma::vec out = L2 * L3 * b + mu_v;
  return Rcpp::wrap(out);
}

List block_envelope_sim_one(
    const List& block_envelope,
    const List& block_standardization,
    const std::string& block_id,
    int l1,
    int n
) {
  const bool prior_only = Rcpp::as<bool>(block_standardization["prior_only"]);
  NumericVector mu_j = block_standardization["mu"];
  NumericMatrix beta_draw(l1, n);
  NumericVector iters_out(n);

  if (prior_only || Rf_isNull(block_envelope)) {
    for (int i = 0; i < n; ++i) {
      beta_draw(Rcpp::_, i) = mu_j;
      iters_out[i] = 1.0;
    }
    return List::create(
      Rcpp::Named("block_id") = block_id,
      Rcpp::Named("identifiable") = false,
      Rcpp::Named("beta") = beta_draw,
      Rcpp::Named("iters_out") = iters_out,
      Rcpp::Named("prior_only") = true
    );
  }

  NumericMatrix L2Inv = block_standardization["L2Inv"];
  NumericMatrix L3Inv = block_standardization["L3Inv"];
  List Env = block_envelope;

  for (int i = 0; i < n; ++i) {
    NumericVector beta_std = envelope_draw_beta_std_one(Env, l1);
    beta_draw(Rcpp::_, i) = unstandardize_beta_one(beta_std, L2Inv, L3Inv, mu_j);
    iters_out[i] = 1.0;
  }

  return List::create(
    Rcpp::Named("block_id") = block_id,
    Rcpp::Named("identifiable") = true,
    Rcpp::Named("beta") = beta_draw,
    Rcpp::Named("iters_out") = iters_out,
    Rcpp::Named("prior_only") = false
  );
}

}  // anonymous namespace

List BlockEnvelopeCentering(
    NumericVector y,
    NumericMatrix x,
    SEXP block,
    SEXP prior_list_sexp,
    SEXP prior_lists_sexp,
    NumericVector offset,
    NumericVector wt,
    double shape,
    double rate,
    double max_disp_perc,
    Nullable<double> disp_lower,
    Nullable<double> disp_upper,
    int p_re,
    int n_rss_iter,
    bool verbose
) {
  const int l2 = y.size();
  const int l1 = x.ncol();
  if (x.nrow() != l2) {
    Rcpp::stop("nrow(x) must equal length(y).");
  }
  if (n_rss_iter < 1) {
    Rcpp::stop("'n_rss_iter' must be at least 1.");
  }

  if (offset.size() == 1) offset = Rcpp::rep(offset[0], l2);
  if (wt.size() == 1) wt = Rcpp::rep(wt[0], l2);
  if (offset.size() != l2) {
    Rcpp::stop("length(offset) must be 1 or length(y).");
  }
  if (wt.size() != l2) {
    Rcpp::stop("length(wt) must be 1 or length(y).");
  }

  List group_info = glmbayes::sim::group::normalize_group_cpp(block, l2);
  const int k = group_info["k"];
  CharacterVector ids = group_info["ids"];
  List row_blocks = group_info["rows"];

  List prior_block = normalize_prior_for_block_ing(
    prior_list_sexp, prior_lists_sexp, group_info, l1
  );

  List prior_list = Rf_isNull(prior_list_sexp) ? List() : List(prior_list_sexp);
  if (max_disp_perc <= 0.0 || max_disp_perc >= 1.0) {
    if (has_non_null(prior_list, "max_disp_perc")) {
      max_disp_perc = Rcpp::as<double>(prior_list["max_disp_perc"]);
    } else {
      max_disp_perc = 0.99;
    }
  }

  disp_lower = resolve_nullable_bound(disp_lower, prior_list, "disp_lower");
  disp_upper = resolve_nullable_bound(disp_upper, prior_list, "disp_upper");

  if (p_re < 1) {
    p_re = l1;
  }

  const double n_w = sum_wt(wt);
  const double n_prior_implied = 2.0 * shape - 1.0 - static_cast<double>(p_re);
  if (n_prior_implied > n_w) {
    Rcpp::stop(
      "dIndependent_Normal_Gamma prior implies n_prior = %g effective prior "
      "observations, but the data supply only n_w = sum(weights) = %g. The "
      "dispersion envelope requires n_prior <= n_w (prior weight pwt <= 0.5); "
      "weaken the prior (smaller shape) or supply more data.",
      n_prior_implied, n_w
    );
  }

  std::vector<BlockWlsInit> wls_init(k);
  double rss_wls_sum = 0.0;
  double denom_sum = 0.0;
  for (int j = 0; j < k; ++j) {
    IntegerVector rows = row_blocks[j];
    List pb = prior_block[j];
    NumericVector y_j = slice_numeric(y, rows);
    NumericMatrix x_j = slice_matrix_rows(x, rows);
    NumericVector offset_j = slice_numeric(offset, rows);
    NumericVector wt_j = slice_numeric(wt, rows);

    wls_init[j] = block_wls_init(y_j, x_j, offset_j, wt_j);
    rss_wls_sum += wls_init[j].rss_wls;
    denom_sum += std::max(wls_init[j].n_w_j - static_cast<double>(wls_init[j].rank), 1.0);

    if (verbose) {
      Rcpp::Rcout << "BlockEnvelopeCentering block " << (j + 1)
                  << " rank=" << wls_init[j].rank
                  << " n_obs=" << rows.size() << std::endl;
    }
  }

  double dispersion2 = rss_wls_sum / denom_sum;
  double RSS_post_pooled = NA_REAL;

  List blocks(k);
  for (int iter = 0; iter < n_rss_iter; ++iter) {
    RSS_post_pooled = 0.0;
    for (int j = 0; j < k; ++j) {
      IntegerVector rows = row_blocks[j];
      List pb = prior_block[j];
      NumericVector mu_j = pb["mu"];
      NumericMatrix P_j = pb["P"];
      NumericVector y_j = slice_numeric(y, rows);
      NumericMatrix x_j = slice_matrix_rows(x, rows);
      NumericVector offset_j = slice_numeric(offset, rows);
      NumericVector wt_j = slice_numeric(wt, rows);

      List one = block_envelope_centering_one(
        y_j, x_j, mu_j, P_j, offset_j, wt_j,
        dispersion2,
        wls_init[j].rank,
        wls_init[j].identifiable
      );

      const double RSS_j = Rcpp::as<double>(one["RSS_post"]);
      RSS_post_pooled += RSS_j;

      if (iter == n_rss_iter - 1) {
        blocks[j] = List::create(
          Rcpp::Named("id") = Rcpp::wrap(Rcpp::as<std::string>(ids[j])),
          Rcpp::Named("n_obs") = rows.size(),
          Rcpp::Named("n_w_j") = one["n_w_j"],
          Rcpp::Named("rank") = one["rank"],
          Rcpp::Named("identifiable") = one["identifiable"],
          Rcpp::Named("RSS_post") = RSS_j,
          Rcpp::Named("mu") = mu_j,
          Rcpp::Named("b_post_mean") = one["b_post_mean"],
          Rcpp::Named("dispersion_wls") = wls_init[j].dispersion_wls
        );
      }
    }

    const double shape2 = shape + n_w / 2.0;
    const double rate2 = rate + RSS_post_pooled / 2.0;
    dispersion2 = rate2 / (shape2 - 1.0);
  }

  const double shape2 = shape + n_w / 2.0;
  const double rate3 = rate + RSS_post_pooled / 2.0;

  RObject disp_lower_out = disp_lower.isUsable()
    ? RObject(Rcpp::wrap(Rcpp::as<double>(disp_lower.get())))
    : RObject(R_NilValue);
  RObject disp_upper_out = disp_upper.isUsable()
    ? RObject(Rcpp::wrap(Rcpp::as<double>(disp_upper.get())))
    : RObject(R_NilValue);

  return List::create(
    Rcpp::Named("dispersion") = dispersion2,
    Rcpp::Named("RSS_post") = RSS_post_pooled,
    Rcpp::Named("n_w") = n_w,
    Rcpp::Named("p_re") = p_re,
    Rcpp::Named("l1") = l1,
    Rcpp::Named("l2") = l2,
    Rcpp::Named("k") = k,
    Rcpp::Named("shape2") = shape2,
    Rcpp::Named("rate3") = rate3,
    Rcpp::Named("max_disp_perc") = max_disp_perc,
    Rcpp::Named("disp_lower") = disp_lower_out,
    Rcpp::Named("disp_upper") = disp_upper_out,
    Rcpp::Named("group_info") = group_info,
    Rcpp::Named("blocks") = blocks,
    Rcpp::Named("prior_lists") = prior_block
  );
}

List BlockEnvelopeBuild(
    const List& centering_out,
    NumericVector y,
    NumericMatrix x,
    SEXP block,
    SEXP prior_list_sexp,
    SEXP prior_lists_sexp,
    NumericVector offset,
    NumericVector wt,
    double max_disp_perc,
    Nullable<double> disp_lower,
    Nullable<double> disp_upper,
    int n,
    int Gridtype,
    Nullable<int> n_envopt,
    double RSS_ML,
    bool use_parallel,
    bool use_opencl,
    bool verbose
) {
  BEB_DBG(verbose, "[BEB 1.0] BlockEnvelopeBuild enter");

  if (!centering_out.containsElementNamed("dispersion") ||
      !centering_out.containsElementNamed("blocks") ||
      !centering_out.containsElementNamed("group_info") ||
      !centering_out.containsElementNamed("k")) {
    Rcpp::stop("'centering_out' must be a full BlockEnvelopeCentering return list.");
  }

  const double dispersion2 = Rcpp::as<double>(centering_out["dispersion"]);
  const double RSS_post = Rcpp::as<double>(centering_out["RSS_post"]);
  const double shape2 = Rcpp::as<double>(centering_out["shape2"]);
  const double rate3 = Rcpp::as<double>(centering_out["rate3"]);
  const double n_w = Rcpp::as<double>(centering_out["n_w"]);
  const int p_re = Rcpp::as<int>(centering_out["p_re"]);
  const int l1 = Rcpp::as<int>(centering_out["l1"]);
  const int l2 = Rcpp::as<int>(centering_out["l2"]);
  const int k = Rcpp::as<int>(centering_out["k"]);
  BEB_DBG(verbose, "[BEB 1.1] parsed centering_out k=" << k << " l1=" << l1
          << " l2=" << l2 << " dispersion=" << dispersion2);

  if (y.size() != l2) {
    Rcpp::stop("length(y) must match centering_out$l2.");
  }
  if (x.nrow() != l2 || x.ncol() != l1) {
    Rcpp::stop("dim(x) must match centering_out$l2 x l1.");
  }

  if (max_disp_perc <= 0.0 || max_disp_perc >= 1.0) {
    if (centering_out.containsElementNamed("max_disp_perc")) {
      max_disp_perc = Rcpp::as<double>(centering_out["max_disp_perc"]);
    } else {
      max_disp_perc = 0.99;
    }
  }

  if (!disp_lower.isUsable() && centering_out.containsElementNamed("disp_lower")) {
    SEXP dl = centering_out["disp_lower"];
    if (!Rf_isNull(dl)) {
      disp_lower = Nullable<double>(Rcpp::wrap(Rcpp::as<double>(dl)));
    }
  }
  if (!disp_upper.isUsable() && centering_out.containsElementNamed("disp_upper")) {
    SEXP du = centering_out["disp_upper"];
    if (!Rf_isNull(du)) {
      disp_upper = Nullable<double>(Rcpp::wrap(Rcpp::as<double>(du)));
    }
  }

  List group_info = centering_out["group_info"];
  List row_blocks = group_info["rows"];
  List blocks_centering = centering_out["blocks"];

  List group_info_data = glmbayes::sim::group::normalize_group_cpp(block, l2);
  if (Rcpp::as<int>(group_info_data["k"]) != k) {
    Rcpp::stop("'block' partition k must match centering_out$k.");
  }

  List prior_block = normalize_prior_for_block_ing(
    prior_list_sexp, prior_lists_sexp, group_info_data, l1
  );
  BEB_DBG(verbose, "[BEB 1.2] partition + prior expand done");

  if (offset.size() == 1) offset = Rcpp::rep(offset[0], l2);
  if (wt.size() == 1) wt = Rcpp::rep(wt[0], l2);
  if (offset.size() != l2) {
    Rcpp::stop("length(offset) must be 1 or length(y).");
  }
  if (wt.size() != l2) {
    Rcpp::stop("length(wt) must be 1 or length(y).");
  }

  const int n_envopt_val = resolve_n_envopt(n_envopt);

  BEB_DBG(verbose, "[BEB 1.3] before optim() lookup");
  Function optim("optim");
  BEB_DBG(verbose, "[BEB 1.4] before r_glmbfamfunc()");
  Function gaussian("gaussian");
  Function glmbfamfunc = glmbayes_R::r_glmbfamfunc();
  BEB_DBG(verbose, "[BEB 1.5] before glmbfamfunc(gaussian())");
  List famfunc = glmbfamfunc(gaussian());
  BEB_DBG(verbose, "[BEB 1.6] after glmbfamfunc(gaussian())");
  Function f2 = famfunc["f2"];
  Function f3 = famfunc["f3"];
  BEB_DBG(verbose, "[BEB 1.7] f2/f3 extracted");

  List block_envelopes(k);
  List block_standardization(k);
  int n_identifiable = 0;

  for (int j = 0; j < k; ++j) {
    IntegerVector rows = row_blocks[j];
    List pb = prior_block[j];
    List bc = blocks_centering[j];
    NumericVector mu_j = pb["mu"];
    NumericMatrix P_j = pb["P"];
    NumericVector y_j = slice_numeric(y, rows);
    NumericMatrix x_j = slice_matrix_rows(x, rows);
    NumericVector offset_j = slice_numeric(offset, rows);
    NumericVector wt_j = slice_numeric(wt, rows);
    NumericVector b_post_mean = bc["b_post_mean"];
    const bool identifiable = Rcpp::as<bool>(bc["identifiable"]);
    const std::string block_id = Rcpp::as<std::string>(bc["id"]);

    BEB_DBG(verbose, "[BEB 2.0] block j=" << (j + 1) << "/" << k
            << " id=" << block_id
            << " identifiable=" << identifiable
            << " n_obs=" << rows.size());

    BEB_DBG(verbose, "[BEB 2.1] before block_envelope_build_one id=" << block_id);
    List one = block_envelope_build_one(
      y_j, x_j, mu_j, P_j, offset_j, wt_j,
      b_post_mean, dispersion2, identifiable, block_id,
      optim, f2, f3, Gridtype, n, n_envopt_val,
      use_opencl, verbose
    );
    BEB_DBG(verbose, "[BEB 2.2] after block_envelope_build_one id=" << block_id);

    if (identifiable) {
      List env_j = Rcpp::as<List>(one["block_envelope"]);
      env_j["block_id"] = block_id;
      block_envelopes[j] = env_j;
      ++n_identifiable;
    } else {
      block_envelopes[j] = R_NilValue;
    }
    block_standardization[j] = one["block_standardization"];

    if (verbose) {
      Rcpp::Rcout << "BlockEnvelopeBuild block " << (j + 1)
                  << " id=" << block_id
                  << " identifiable=" << identifiable << std::endl;
    }
  }

  BEB_DBG(verbose, "[BEB 3.0] all blocks done n_identifiable=" << n_identifiable);

  RObject disp_lower_out = disp_lower.isUsable()
    ? RObject(Rcpp::wrap(Rcpp::as<double>(disp_lower.get())))
    : RObject(R_NilValue);
  RObject disp_upper_out = disp_upper.isUsable()
    ? RObject(Rcpp::wrap(Rcpp::as<double>(disp_upper.get())))
    : RObject(R_NilValue);

  List dispersion_envelope = List::create(
    Rcpp::Named("gamma_list") = List::create(
      Rcpp::Named("shape2") = shape2,
      Rcpp::Named("rate2") = rate3,
      Rcpp::Named("disp_lower") = disp_lower_out,
      Rcpp::Named("disp_upper") = disp_upper_out
    ),
    Rcpp::Named("UB_list") = List::create(
      Rcpp::Named("status") = "stub_v1"
    ),
    Rcpp::Named("low") = disp_lower_out,
    Rcpp::Named("upp") = disp_upper_out,
    Rcpp::Named("RSS_post") = RSS_post,
    Rcpp::Named("RSS_ML") = RSS_ML,
    Rcpp::Named("status") = "stub_pending_BlockEnvelopeDispersionBuild"
  );

  List meta = List::create(
    Rcpp::Named("k") = k,
    Rcpp::Named("l1") = l1,
    Rcpp::Named("l2") = l2,
    Rcpp::Named("n_w") = n_w,
    Rcpp::Named("p_re") = p_re,
    Rcpp::Named("n") = n,
    Rcpp::Named("Gridtype") = Gridtype,
    Rcpp::Named("d1_star") = dispersion2,
    Rcpp::Named("dispersion") = dispersion2,
    Rcpp::Named("RSS_post") = RSS_post,
    Rcpp::Named("shape2") = shape2,
    Rcpp::Named("rate3") = rate3,
    Rcpp::Named("group_info") = group_info,
    Rcpp::Named("n_identifiable") = n_identifiable,
    Rcpp::Named("use_parallel") = use_parallel
  );

  BEB_DBG(verbose, "[BEB 3.1] BlockEnvelopeBuild return");
  return List::create(
    Rcpp::Named("block_envelopes") = block_envelopes,
    Rcpp::Named("dispersion_envelope") = dispersion_envelope,
    Rcpp::Named("block_standardization") = block_standardization,
    Rcpp::Named("meta") = meta
  );
}

// lmerb-unused: candidate for deletion (not on the lmerb path).
List BlockEnvelopeDispersionBuild(
    const List& build_out,
    const List& centering_out,
    NumericVector y,
    NumericMatrix x,
    SEXP block,
    NumericVector offset,
    NumericVector wt,
    double shape,
    double rate,
    double max_disp_perc,
    Nullable<double> disp_lower,
    Nullable<double> disp_upper,
    double RSS_ML,
    bool use_parallel,
    bool verbose
) {
  if (!build_out.containsElementNamed("block_envelopes") ||
      !build_out.containsElementNamed("block_standardization") ||
      !build_out.containsElementNamed("meta")) {
    Rcpp::stop("'build_out' must be a full BlockEnvelopeBuild return list.");
  }
  if (!centering_out.containsElementNamed("dispersion") ||
      !centering_out.containsElementNamed("blocks") ||
      !centering_out.containsElementNamed("group_info") ||
      !centering_out.containsElementNamed("RSS_post") ||
      !centering_out.containsElementNamed("shape2") ||
      !centering_out.containsElementNamed("rate3")) {
    Rcpp::stop("'centering_out' must be a full BlockEnvelopeCentering return list.");
  }

  List block_envelopes = Rcpp::as<List>(build_out["block_envelopes"]);
  block_envelopes = Rcpp::clone(block_envelopes);
  List block_standardization = build_out["block_standardization"];
  List meta = build_out["meta"];

  const double RSS_post_global = Rcpp::as<double>(centering_out["RSS_post"]);
  const double shape2_global = Rcpp::as<double>(centering_out["shape2"]);
  const double rate3_global = Rcpp::as<double>(centering_out["rate3"]);
  const double n_w_global = Rcpp::as<double>(centering_out["n_w"]);
  const int l1 = Rcpp::as<int>(centering_out["l1"]);
  const int l2 = Rcpp::as<int>(centering_out["l2"]);
  const int k = Rcpp::as<int>(centering_out["k"]);

  if (y.size() != l2) {
    Rcpp::stop("length(y) must match centering_out$l2.");
  }
  if (x.nrow() != l2 || x.ncol() != l1) {
    Rcpp::stop("dim(x) must match centering_out$l2 x l1.");
  }

  if (max_disp_perc <= 0.0 || max_disp_perc >= 1.0) {
    if (centering_out.containsElementNamed("max_disp_perc")) {
      max_disp_perc = Rcpp::as<double>(centering_out["max_disp_perc"]);
    } else {
      max_disp_perc = 0.99;
    }
  }

  if (!disp_lower.isUsable() && centering_out.containsElementNamed("disp_lower")) {
    SEXP dl = centering_out["disp_lower"];
    if (!Rf_isNull(dl)) {
      disp_lower = Nullable<double>(Rcpp::wrap(Rcpp::as<double>(dl)));
    }
  }
  if (!disp_upper.isUsable() && centering_out.containsElementNamed("disp_upper")) {
    SEXP du = centering_out["disp_upper"];
    if (!Rf_isNull(du)) {
      disp_upper = Nullable<double>(Rcpp::wrap(Rcpp::as<double>(du)));
    }
  }

  List group_info = centering_out["group_info"];
  List row_blocks = group_info["rows"];
  List blocks_centering = centering_out["blocks"];

  List group_info_data = glmbayes::sim::group::normalize_group_cpp(block, l2);
  if (Rcpp::as<int>(group_info_data["k"]) != k) {
    Rcpp::stop("'block' partition k must match centering_out$k.");
  }

  if (offset.size() == 1) offset = Rcpp::rep(offset[0], l2);
  if (wt.size() == 1) wt = Rcpp::rep(wt[0], l2);
  if (offset.size() != l2) {
    Rcpp::stop("length(offset) must be 1 or length(y).");
  }
  if (wt.size() != l2) {
    Rcpp::stop("length(wt) must be 1 or length(y).");
  }

  double low = 0.0;
  double upp = 0.0;
  compute_sigma2_bounds_cpp(
    shape2_global, rate3_global, max_disp_perc,
    disp_lower, disp_upper, &low, &upp
  );

  List block_dispersion(k);
  IntegerVector identifiable_idx;
  double RSS_Min_global = 0.0;
  double RSS_ML_global = 0.0;
  bool rss_ml_provided = R_finite(RSS_ML);
  if (rss_ml_provided) {
    RSS_ML_global = RSS_ML;
  }

  List single_gamma_list;
  List single_ub_list;

  for (int j = 0; j < k; ++j) {
    List bc = blocks_centering[j];
    const bool identifiable = Rcpp::as<bool>(bc["identifiable"]);
    const std::string block_id = Rcpp::as<std::string>(bc["id"]);

    if (!identifiable) {
      block_dispersion[j] = R_NilValue;
      continue;
    }

    IntegerVector rows = row_blocks[j];
    List std_j = block_standardization[j];
    List Env_j = block_envelopes[j];

    NumericVector y_j = slice_numeric(y, rows);
    NumericMatrix x2_j = std_j["x2"];
    NumericMatrix P2_j = std_j["P2"];
    NumericMatrix mu2_j = std_j["mu2"];
    NumericVector alpha_j = std_j["alpha"];
    NumericVector wt_j = slice_numeric(wt, rows);
    const double RSS_post_j = Rcpp::as<double>(bc["RSS_post"]);

    double RSS_ML_j = NA_REAL;
    if (rss_ml_provided) {
      RSS_ML_j = RSS_ML;
    }

    List one = block_envelope_dispersion_one(
      Env_j,
      shape,
      rate,
      P2_j,
      y_j,
      x2_j,
      mu2_j,
      alpha_j,
      wt_j,
      RSS_post_j,
      RSS_ML_j,
      max_disp_perc,
      disp_lower,
      disp_upper,
      use_parallel,
      verbose,
      block_id
    );

    block_envelopes[j] = one["Env_out"];
    List ub_j = one["UB_list"];
    List gamma_j = one["gamma_list"];
    block_dispersion[j] = List::create(
      Rcpp::Named("cache") = one["cache"],
      Rcpp::Named("lg_prob_factor") = one["lg_prob_factor"],
      Rcpp::Named("UB2min") = one["UB2min"],
      Rcpp::Named("UB_list") = ub_j,
      Rcpp::Named("gamma_list") = gamma_j,
      Rcpp::Named("rss_min_global") = one["rss_min_global"],
      Rcpp::Named("RSS_ML") = one["RSS_ML"],
      Rcpp::Named("diagnostics") = one["diagnostics"],
      Rcpp::Named("block_id") = block_id,
      Rcpp::Named("identifiable") = true
    );

    RSS_Min_global += Rcpp::as<double>(one["rss_min_global"]);
    if (!rss_ml_provided) {
      RSS_ML_global += Rcpp::as<double>(one["RSS_ML"]);
    }
    if (identifiable_idx.size() == 0) {
      single_gamma_list = gamma_j;
      single_ub_list = ub_j;
    }

    identifiable_idx.push_back(j);
  }

  if (identifiable_idx.size() < 1) {
    Rcpp::stop("BlockEnvelopeDispersionBuild: no identifiable blocks.");
  }

  List global_constants = build_global_dispersion_constants(
    identifiable_idx.size(),
    single_gamma_list,
    single_ub_list,
    identifiable_idx,
    block_envelopes,
    block_standardization,
    block_dispersion,
    shape2_global,
    rate,
    rate3_global,
    n_w_global,
    RSS_Min_global,
    RSS_ML_global,
    low,
    upp,
    RSS_post_global
  );

  List gamma_list = global_constants["gamma_list"];
  List UB_list_global = global_constants["UB_list"];

  double prob_max_upp = Rcpp::as<double>(UB_list_global["max_New_LL_UB"]);
  double prob_max_low = NA_REAL;
  List joint_slack;
  if (identifiable_idx.size() > 1) {
    List apprx = compute_block_face_apprx_and_prob_anchors(
      identifiable_idx,
      block_envelopes,
      block_standardization,
      block_dispersion,
      shape2_global,
      rate3_global,
      low,
      upp
    );
    patch_block_dispersion_apprx(block_dispersion, apprx);
    if (global_constants.containsElementNamed("prob_max_low")) {
      prob_max_low = Rcpp::as<double>(global_constants["prob_max_low"]);
    } else {
      prob_max_low = Rcpp::as<double>(apprx["prob_max_low"]);
    }
    joint_slack = build_joint_product_face_slack(
      identifiable_idx,
      block_envelopes,
      block_dispersion,
      global_constants["joint_upp_apprx"],
      global_constants["joint_low_apprx"],
      prob_max_upp,
      prob_max_low
    );
  }

  RObject disp_lower_out = disp_lower.isUsable()
    ? RObject(Rcpp::wrap(Rcpp::as<double>(disp_lower.get())))
    : RObject(R_NilValue);
  RObject disp_upper_out = disp_upper.isUsable()
    ? RObject(Rcpp::wrap(Rcpp::as<double>(disp_upper.get())))
    : RObject(R_NilValue);

  IntegerVector gs_per_block(identifiable_idx.size());
  for (int t = 0; t < identifiable_idx.size(); ++t) {
    List env_t = block_envelopes[identifiable_idx[t]];
    NumericMatrix cbars_t = env_t["cbars"];
    gs_per_block[t] = cbars_t.nrow();
  }

  List dispersion_envelope = List::create(
    Rcpp::Named("gamma_list") = gamma_list,
    Rcpp::Named("UB_list") = UB_list_global,
    Rcpp::Named("block_dispersion") = block_dispersion,
    Rcpp::Named("cross_face_meta") = List::create(
      Rcpp::Named("gs_per_block") = gs_per_block,
      Rcpp::Named("identifiable_idx") = identifiable_idx,
      Rcpp::Named("n_identifiable") = identifiable_idx.size(),
      Rcpp::Named("aggregation") = global_constants["source"],
      Rcpp::Named("ub3a_lg_mode") = (identifiable_idx.size() > 1)
        ? "joint_product_face_lookup_v2"
        : "single_block_lg",
      Rcpp::Named("ub2_mode") = (identifiable_idx.size() > 1)
        ? "joint_product_face_ub2_rss_v2"
        : "single_block_ub2min",
      Rcpp::Named("face_draw_mode") = (identifiable_idx.size() > 1)
        ? "joint_product_plsd_v1"
        : "per_block_plsd"
    ),
    Rcpp::Named("prob_max_upp") = prob_max_upp,
    Rcpp::Named("prob_max_low") = prob_max_low,
    Rcpp::Named("low") = low,
    Rcpp::Named("upp") = upp,
    Rcpp::Named("RSS_post") = RSS_post_global,
    Rcpp::Named("RSS_ML") = RSS_ML_global,
    Rcpp::Named("RSS_Min") = RSS_Min_global,
    Rcpp::Named("status") = "v1"
  );

  if (identifiable_idx.size() > 1) {
    dispersion_envelope["joint_lg_prob_factor"] =
      joint_slack["joint_lg_prob_factor"];
    dispersion_envelope["joint_ub2min_product"] =
      joint_slack["joint_ub2min_product"];
    dispersion_envelope["joint_PLSD"] = joint_slack["joint_PLSD"];
    dispersion_envelope["n_product_faces"] = joint_slack["n_product_faces"];
  }

  List meta_out = Rcpp::clone(meta);
  meta_out["dispersion_envelope_status"] = "v1";

  List build_out_patched = List::create(
    Rcpp::Named("block_envelopes") = block_envelopes,
    Rcpp::Named("dispersion_envelope") = dispersion_envelope,
    Rcpp::Named("block_standardization") = block_standardization,
    Rcpp::Named("meta") = meta_out
  );

  if (verbose) {
    Rcpp::Rcout << "BlockEnvelopeDispersionBuild: k=" << k
                << " n_identifiable=" << identifiable_idx.size()
                << " RSS_Min=" << RSS_Min_global
                << " low=" << low << " upp=" << upp << std::endl;
  }

  return List::create(
    Rcpp::Named("build_out") = build_out_patched,
    Rcpp::Named("dispersion_envelope") = dispersion_envelope,
    Rcpp::Named("low") = low,
    Rcpp::Named("upp") = upp,
    Rcpp::Named("diagnostics") = List::create(
      Rcpp::Named("RSS_Min_global") = RSS_Min_global,
      Rcpp::Named("RSS_ML_global") = RSS_ML_global,
      Rcpp::Named("n_identifiable") = identifiable_idx.size(),
      Rcpp::Named("aggregation") = global_constants["source"]
    )
  );
}

// lmerb-unused: candidate for deletion (not on the lmerb path).
List BlockEnvelopeSim(
    const List& build_out,
    int n,
    bool progbar,
    bool verbose
) {
  if (n < 1) {
    Rcpp::stop("'n' must be at least 1.");
  }
  if (!build_out.containsElementNamed("block_envelopes") ||
      !build_out.containsElementNamed("block_standardization") ||
      !build_out.containsElementNamed("meta")) {
    Rcpp::stop("'build_out' must be a full BlockEnvelopeBuild return list.");
  }

  List block_envelopes = build_out["block_envelopes"];
  List block_standardization = build_out["block_standardization"];
  List meta = build_out["meta"];
  if (!meta.containsElementNamed("group_info") ||
      !meta.containsElementNamed("l1") ||
      !meta.containsElementNamed("dispersion")) {
    Rcpp::stop("'build_out$meta' must contain group_info, l1, and dispersion.");
  }

  List group_info = meta["group_info"];
  const int k = Rcpp::as<int>(group_info["k"]);
  const int l1 = Rcpp::as<int>(meta["l1"]);
  const double dispersion_anchor = Rcpp::as<double>(meta["dispersion"]);
  CharacterVector ids = group_info["ids"];

  if (block_envelopes.size() != k ||
      block_standardization.size() != k) {
    Rcpp::stop("block_envelopes and block_standardization must have length k.");
  }

  const bool has_disp_env =
    build_out.containsElementNamed("dispersion_envelope") &&
    !Rf_isNull(build_out["dispersion_envelope"]);
  std::string disp_status = "missing";
  if (has_disp_env) {
    List dispersion_envelope = build_out["dispersion_envelope"];
    if (dispersion_envelope.containsElementNamed("status") &&
        !Rf_isNull(dispersion_envelope["status"])) {
      disp_status = Rcpp::as<std::string>(dispersion_envelope["status"]);
    }
  }
  const bool use_ar_compute = (disp_status == "v1");

  List group_results(k);
  NumericVector dispersion_out(n);
  NumericVector iters_out(n);

  if (!use_ar_compute) {
    std::fill(dispersion_out.begin(), dispersion_out.end(), dispersion_anchor);
    std::fill(iters_out.begin(), iters_out.end(), 1.0);
    for (int j = 0; j < k; ++j) {
      if (progbar && verbose) {
        Rcpp::Rcout << "BlockEnvelopeSim block " << (j + 1) << "/" << k << std::endl;
      }
      const std::string block_id = Rcpp::as<std::string>(ids[j]);
      SEXP env_j = block_envelopes[j];
      List env_list = Rf_isNull(env_j) ? List() : List(env_j);
      List std_j = block_standardization[j];
      group_results[j] = block_envelope_sim_one(
        env_list, std_j, block_id, l1, n
      );
    }
    List meta_out = Rcpp::clone(meta);
    meta_out["accept_mode"] = "auto_v1";
    meta_out["n"] = n;
    return List::create(
      Rcpp::Named("group_results") = group_results,
      Rcpp::Named("dispersion") = dispersion_out,
      Rcpp::Named("iters_out") = iters_out,
      Rcpp::Named("group_info") = group_info,
      Rcpp::Named("meta") = meta_out
    );
  }

  List dispersion_envelope = build_out["dispersion_envelope"];
  List gamma_list = dispersion_envelope["gamma_list"];
  List UB_global = dispersion_envelope["UB_list"];
  List block_dispersion = dispersion_envelope["block_dispersion"];

  const double shape3 = Rcpp::as<double>(gamma_list["shape3"]);
  const double rate2 = Rcpp::as<double>(gamma_list["rate2"]);
  const double disp_lower = Rcpp::as<double>(gamma_list["disp_lower"]);
  const double disp_upper = Rcpp::as<double>(gamma_list["disp_upper"]);
  const double RSS_Min_G = Rcpp::as<double>(dispersion_envelope["RSS_Min"]);
  const double max_New_LL_UB = Rcpp::as<double>(UB_global["max_New_LL_UB"]);
  const double max_LL_log_disp = Rcpp::as<double>(UB_global["max_LL_log_disp"]);
  const double lm_log1 = Rcpp::as<double>(UB_global["lm_log1"]);
  const double lm_log2 = Rcpp::as<double>(UB_global["lm_log2"]);
  const double lmc1 = Rcpp::as<double>(UB_global["lmc1"]);
  const double lmc2 = Rcpp::as<double>(UB_global["lmc2"]);

  List cross_face_meta = dispersion_envelope["cross_face_meta"];
  const int n_identifiable = Rcpp::as<int>(cross_face_meta["n_identifiable"]);
  const bool use_joint_lg_ub3a = (n_identifiable > 1);
  IntegerVector identifiable_idx;
  IntegerVector gs_per_block;
  NumericVector joint_lg_prob_factor;
  NumericVector joint_ub2min_product;
  NumericVector joint_PLSD;
  if (use_joint_lg_ub3a) {
    identifiable_idx = cross_face_meta["identifiable_idx"];
    gs_per_block = cross_face_meta["gs_per_block"];
    if (!dispersion_envelope.containsElementNamed("joint_lg_prob_factor") ||
        Rf_isNull(dispersion_envelope["joint_lg_prob_factor"]) ||
        !dispersion_envelope.containsElementNamed("joint_ub2min_product") ||
        Rf_isNull(dispersion_envelope["joint_ub2min_product"]) ||
        !dispersion_envelope.containsElementNamed("joint_PLSD") ||
        Rf_isNull(dispersion_envelope["joint_PLSD"])) {
      Rcpp::stop(
        "BlockEnvelopeSim: joint product-face tables missing "
        "(joint_lg_prob_factor / joint_ub2min_product / joint_PLSD)."
      );
    }
    joint_lg_prob_factor = dispersion_envelope["joint_lg_prob_factor"];
    joint_ub2min_product = dispersion_envelope["joint_ub2min_product"];
    joint_PLSD = dispersion_envelope["joint_PLSD"];
  }

  // Precompute face-draw CDFs once per sim call (joint_PLSD / each block's own
  // PLSD are fixed for the whole resample-until-accept loop below), so every
  // attempt does an O(log gs) binary-search draw instead of an O(gs) linear
  // scan. gs can be O(prod gs_t) for k > 1 blocks, so this matters most there.
  std::vector<double> joint_plsd_cdf;
  if (use_joint_lg_ub3a) {
    joint_plsd_cdf = build_face_cdf(joint_PLSD);
  }
  std::vector<std::vector<double>> block_plsd_cdf(k);

  for (int j = 0; j < k; ++j) {
    List std_j = block_standardization[j];
    const bool prior_only = Rcpp::as<bool>(std_j["prior_only"]);
    const std::string block_id = Rcpp::as<std::string>(ids[j]);

    if (prior_only || Rf_isNull(block_envelopes[j])) {
      NumericVector mu_j = std_j["mu"];
      NumericMatrix beta_block(l1, n);
      NumericVector block_iters(n);
      for (int i = 0; i < n; ++i) {
        beta_block(Rcpp::_, i) = mu_j;
        block_iters[i] = 1.0;
      }
      group_results[j] = List::create(
        Rcpp::Named("block_id") = block_id,
        Rcpp::Named("identifiable") = false,
        Rcpp::Named("beta") = beta_block,
        Rcpp::Named("iters_out") = block_iters,
        Rcpp::Named("prior_only") = true
      );
    } else {
      if (!use_joint_lg_ub3a) {
        List Env = block_envelopes[j];
        block_plsd_cdf[j] = build_face_cdf(Env["PLSD"]);
      }
      group_results[j] = List::create(
        Rcpp::Named("block_id") = block_id,
        Rcpp::Named("identifiable") = true,
        Rcpp::Named("beta") = NumericMatrix(l1, n),
        Rcpp::Named("iters_out") = NumericVector(n),
        Rcpp::Named("prior_only") = false,
        Rcpp::Named("face_J_last") = NA_INTEGER
      );
    }
  }

  std::vector<int> J_draw(k, 0);
  std::vector<NumericVector> beta_std_draw(k);
  std::vector<NumericVector> beta_orig_draw(k);

  for (int i = 0; i < n; ++i) {
    Rcpp::checkUserInterrupt();
    if (progbar && n > 1) {
      progress_bar(static_cast<double>(i), static_cast<double>(n - 1));
      if (i == n - 1) {
        Rcpp::Rcout << "" << std::endl;
      }
    }

    iters_out[i] = 1.0;
    int accept = 0;
    double sigma2_accept = NA_REAL;
    const auto ar_t0 = std::chrono::steady_clock::now();

    while (accept == 0) {
      if (use_joint_lg_ub3a) {
        const int flat = draw_face_index_from_cdf(joint_plsd_cdf);
        std::vector<int> faces_identifiable(static_cast<size_t>(n_identifiable));
        decode_product_face_flat_index(flat, gs_per_block, faces_identifiable);
        for (int t = 0; t < n_identifiable; ++t) {
          const int j = identifiable_idx[t];
          List std_j = block_standardization[j];
          List Env = block_envelopes[j];
          J_draw[j] = faces_identifiable[static_cast<size_t>(t)];
          beta_std_draw[j] = draw_beta_std_face(Env, J_draw[j], l1);
          beta_orig_draw[j] = unstandardize_beta_one(
            beta_std_draw[j],
            std_j["L2Inv"],
            std_j["L3Inv"],
            std_j["mu"]
          );
        }
        for (int j = 0; j < k; ++j) {
          List std_j = block_standardization[j];
          if (Rcpp::as<bool>(std_j["prior_only"]) || Rf_isNull(block_envelopes[j])) {
            beta_orig_draw[j] = std_j["mu"];
            J_draw[j] = -1;
          }
        }
      } else {
        for (int j = 0; j < k; ++j) {
          List std_j = block_standardization[j];
          if (Rcpp::as<bool>(std_j["prior_only"]) || Rf_isNull(block_envelopes[j])) {
            beta_orig_draw[j] = std_j["mu"];
            J_draw[j] = -1;
            continue;
          }
          List Env = block_envelopes[j];
          J_draw[j] = draw_face_index_from_cdf(block_plsd_cdf[j]);
          beta_std_draw[j] = draw_beta_std_face(Env, J_draw[j], l1);
          beta_orig_draw[j] = unstandardize_beta_one(
            beta_std_draw[j],
            std_j["L2Inv"],
            std_j["L3Inv"],
            std_j["mu"]
          );
        }
      }

      const double sigma2 =
        rinvgamma_ct_safe(shape3, rate2, disp_upper, disp_lower);
      if (!R_finite(sigma2)) {
        Rcpp::stop("BlockEnvelopeSim: non-finite sigma2 draw.");
      }

      double LL_total = 0.0;
      double UB1_total = 0.0;
      double quad_sum = 0.0;
      double ub3a_block_sum = 0.0;
      double UB2min_used = 0.0;

      for (int j = 0; j < k; ++j) {
        List std_j = block_standardization[j];
        if (Rcpp::as<bool>(std_j["prior_only"]) || Rf_isNull(block_envelopes[j])) {
          continue;
        }
        List Env = block_envelopes[j];
        List block_disp_j = block_dispersion[j];
        block_ar_accumulate_one(
          J_draw[j], beta_std_draw[j], sigma2, Env, std_j, block_disp_j,
          &LL_total, &UB1_total, &quad_sum, &ub3a_block_sum,
          !use_joint_lg_ub3a
        );
        if (!use_joint_lg_ub3a) {
          NumericVector UB2min_j = block_disp_j["UB2min"];
          UB2min_used += UB2min_j[J_draw[j]];
        }
      }

      if (use_joint_lg_ub3a) {
        ub3a_block_sum += joint_lg_prob_factor_at_draw(
          joint_lg_prob_factor, J_draw, identifiable_idx, gs_per_block
        );
        UB2min_used = joint_ub2min_at_draw(
          joint_ub2min_product, J_draw, identifiable_idx, gs_per_block
        );
      }

      const double UB2_raw = 0.5 * (1.0 / sigma2) * (quad_sum - RSS_Min_G);
      const double UB2 = UB2_raw - UB2min_used;
      const double UB3A = ub3a_block_sum + lmc1 + lmc2 * sigma2;
      const double New_LL_log_disp = lm_log1 + lm_log2 * std::log(sigma2);
      const double UB3B =
        (max_New_LL_UB - max_LL_log_disp + New_LL_log_disp) -
        (lmc1 + lmc2 * sigma2);

      const double test1 = LL_total - UB1_total;
      double test = test1 - (UB2 + UB3A + UB3B);

      block_ar_check_signs(
        test1, UB1_total, UB2, UB3A, UB3B, test,
        quad_sum, RSS_Min_G, UB2_raw, UB2min_used, verbose
      );

      const double U2 = runif_safe();
      test -= std::log(U2);

      if (test >= 0.0) {
        accept = 1;
        sigma2_accept = sigma2;

        if (verbose && i == 0) {
          Rcpp::Rcout << "BlockEnvelopeSim resample_until_accept draw 0: sigma2="
                      << sigma2 << " test=" << test
                      << " iters=" << iters_out[i]
                      << " LL=" << LL_total << " UB1=" << UB1_total
                      << " UB2=" << UB2 << " UB3A=" << UB3A
                      << " UB3B=" << UB3B << std::endl;
        }
      } else {
        iters_out[i] = iters_out[i] + 1.0;
        const double ar_elapsed = block_ar_elapsed_seconds(ar_t0);
        if (ar_elapsed > kBlockArTimeLimitSec) {
          accept = 1;
          sigma2_accept = sigma2;
          block_ar_warn_time_limit_accept(
            "BlockEnvelopeSim",
            i,
            iters_out[i],
            ar_elapsed,
            kBlockArTimeLimitSec
          );
        }
      }
    }

    dispersion_out[i] = sigma2_accept;

    for (int j = 0; j < k; ++j) {
      List std_j = block_standardization[j];
      if (Rcpp::as<bool>(std_j["prior_only"]) || Rf_isNull(block_envelopes[j])) {
        List res_j = group_results[j];
        NumericMatrix beta_block = res_j["beta"];
        beta_block(Rcpp::_, i) = beta_orig_draw[j];
        res_j["beta"] = beta_block;
        NumericVector block_iters = res_j["iters_out"];
        block_iters[i] = iters_out[i];
        res_j["iters_out"] = block_iters;
        group_results[j] = res_j;
        continue;
      }
      List res_j = group_results[j];
      NumericMatrix beta_block = res_j["beta"];
      beta_block(Rcpp::_, i) = beta_orig_draw[j];
      res_j["beta"] = beta_block;
      NumericVector block_iters = res_j["iters_out"];
      block_iters[i] = iters_out[i];
      res_j["iters_out"] = block_iters;
      res_j["face_J_last"] = J_draw[j];
      group_results[j] = res_j;
    }
  }

  List meta_out = Rcpp::clone(meta);
  meta_out["accept_mode"] = use_joint_lg_ub3a
    ? "resample_until_accept_joint_product_slack_v2"
    : "resample_until_accept_v1";
  meta_out["face_draw_mode"] = use_joint_lg_ub3a
    ? "joint_product_plsd_v1"
    : "per_block_plsd";
  meta_out["n"] = n;

  return List::create(
    Rcpp::Named("group_results") = group_results,
    Rcpp::Named("dispersion") = dispersion_out,
    Rcpp::Named("iters_out") = iters_out,
    Rcpp::Named("group_info") = group_info,
    Rcpp::Named("meta") = meta_out
  );
}

static void prior_hyperparams_from_list(
    SEXP prior_list_sexp,
    double* shape,
    double* rate,
    double* max_disp_perc,
    Nullable<double>* disp_lower,
    Nullable<double>* disp_upper
) {
  List prior_list(prior_list_sexp);
  if (!has_non_null(prior_list, "shape") || !has_non_null(prior_list, "rate")) {
    Rcpp::stop("'prior_list' must contain 'shape' and 'rate'.");
  }
  *shape = Rcpp::as<double>(prior_list["shape"]);
  *rate = Rcpp::as<double>(prior_list["rate"]);
  if (has_non_null(prior_list, "max_disp_perc")) {
    *max_disp_perc = Rcpp::as<double>(prior_list["max_disp_perc"]);
  } else {
    *max_disp_perc = 0.99;
  }
  *disp_lower = resolve_nullable_bound(Nullable<double>(), prior_list, "disp_lower");
  *disp_upper = resolve_nullable_bound(Nullable<double>(), prior_list, "disp_upper");
}

NumericMatrix map_block_sim_to_b_draw(
    const List& sim,
    const CharacterVector& group_levels,
    int p_re,
    int draw_i
) {
  List group_results = sim["group_results"];
  const int k = group_results.size();
  if (k < 1) {
    Rcpp::stop("BlockEnvelopeSim returned no group_results.");
  }

  std::unordered_map<std::string, NumericVector> by_id;
  by_id.reserve(static_cast<std::size_t>(k));
  for (int j = 0; j < k; ++j) {
    List br = group_results[j];
    const std::string block_id = Rcpp::as<std::string>(br["block_id"]);
    NumericMatrix beta = br["beta"];
    if (draw_i < 0 || draw_i >= beta.ncol()) {
      Rcpp::stop("map_block_sim_to_b_draw: draw index out of range.");
    }
    NumericVector v = beta(Rcpp::_, draw_i);
    if (v.size() != p_re) {
      Rcpp::stop(
        "BlockEnvelopeSim beta length (%d) must equal p_re (%d).",
        v.size(), p_re
      );
    }
    by_id[block_id] = v;
  }

  if (group_levels.size() < 1) {
    NumericMatrix b_draw(k, p_re);
    for (int j = 0; j < k; ++j) {
      List br = group_results[j];
      const std::string block_id = Rcpp::as<std::string>(br["block_id"]);
      b_draw(j, Rcpp::_) = by_id.at(block_id);
    }
    return b_draw;
  }

  NumericMatrix b_draw(group_levels.size(), p_re);
  for (int i = 0; i < group_levels.size(); ++i) {
    const std::string gl = Rcpp::as<std::string>(group_levels[i]);
    const auto it = by_id.find(gl);
    if (it == by_id.end()) {
      Rcpp::stop(
        "BlockEnvelopeSim block ids do not cover all group levels."
      );
    }
    b_draw(i, Rcpp::_) = it->second;
  }
  return b_draw;
}

NumericMatrix block_beta_to_out_k1(
    const List& sim,
    int p_re,
    int n
) {
  List group_results = sim["group_results"];
  if (group_results.size() != 1) {
    Rcpp::stop("'out' is only defined when k = 1 identifiable block layout.");
  }
  List br = group_results[0];
  NumericMatrix beta = br["beta"];
  if (beta.nrow() != p_re || beta.ncol() != n) {
    Rcpp::stop("k = 1 beta matrix must be p_re x n.");
  }
  NumericMatrix out(p_re, n);
  for (int j = 0; j < p_re; ++j) {
    for (int i = 0; i < n; ++i) {
      out(j, i) = beta(j, i);
    }
  }
  return out;
}

// Orchestrator: BlockEnvelopeCentering → Build → DispersionBuild → Sim
// (same sequence as .two_block_block1_envelope_draw_one_chain() in R).
// lmerb-unused: candidate for deletion (not on the lmerb path).
List rIndepNormalGammaRegBlock(
    int n,
    NumericVector y,
    NumericMatrix x,
    SEXP block,
    SEXP prior_list_sexp,
    SEXP prior_lists_sexp,
    NumericVector offset,
    NumericVector wt,
    int p_re,
    int n_rss_iter,
    int Gridtype,
    Nullable<int> n_envopt,
    double RSS_ML,
    bool use_parallel,
    bool use_opencl,
    bool progbar,
    bool verbose,
    CharacterVector group_levels,
    CharacterVector re_names
) {
  if (n < 1) {
    Rcpp::stop("'n' must be at least 1.");
  }
  if (n_rss_iter < 1) {
    Rcpp::stop("'n_rss_iter' must be at least 1.");
  }

  double shape = NA_REAL;
  double rate = NA_REAL;
  double max_disp_perc = 0.99;
  Nullable<double> disp_lower;
  Nullable<double> disp_upper;
  prior_hyperparams_from_list(
    prior_list_sexp, &shape, &rate, &max_disp_perc, &disp_lower, &disp_upper
  );

  List center = BlockEnvelopeCentering(
    y, x, block, prior_list_sexp, prior_lists_sexp,
    offset, wt, shape, rate, max_disp_perc,
    disp_lower, disp_upper, p_re, n_rss_iter, verbose
  );

  List build = BlockEnvelopeBuild(
    center, y, x, block, prior_list_sexp, prior_lists_sexp,
    offset, wt, max_disp_perc, disp_lower, disp_upper,
    n, Gridtype, n_envopt, RSS_ML, use_parallel, use_opencl, verbose
  );

  List disp = BlockEnvelopeDispersionBuild(
    build, center, y, x, block, offset, wt,
    shape, rate, max_disp_perc, disp_lower, disp_upper,
    RSS_ML, use_parallel, verbose
  );
  build = disp["build_out"];

  List sim = BlockEnvelopeSim(build, n, progbar, verbose);

  const int p_re_eff = Rcpp::as<int>(center["p_re"]);
  const int k = Rcpp::as<int>(center["k"]);
  NumericVector dispersion = sim["dispersion"];
  NumericVector iters_out = sim["iters_out"];

  NumericMatrix b_draw = map_block_sim_to_b_draw(
    sim, group_levels, p_re_eff, 0
  );
  if (re_names.size() == p_re_eff) {
    colnames(b_draw) = re_names;
  } else if (re_names.size() > 0 && re_names.size() != p_re_eff) {
    Rcpp::stop("'re_names' must have length 0 or p_re.");
  }
  if (group_levels.size() == b_draw.nrow()) {
    rownames(b_draw) = group_levels;
  }

  NumericVector weight_out(n);
  std::fill(weight_out.begin(), weight_out.end(), 1.0);

  const double low = Rcpp::as<double>(disp["low"]);
  const double upp = Rcpp::as<double>(disp["upp"]);

  List ret = List::create(
    Rcpp::Named("b") = b_draw,
    Rcpp::Named("dispersion_ranef") = dispersion[0],
    Rcpp::Named("iters_mean") = iters_out[0],
    Rcpp::Named("disp_out") = dispersion,
    Rcpp::Named("iters_out") = iters_out,
    Rcpp::Named("weight_out") = weight_out,
    Rcpp::Named("low") = low,
    Rcpp::Named("upp") = upp,
    Rcpp::Named("k") = k,
    Rcpp::Named("p_re") = p_re_eff,
    Rcpp::Named("sim") = sim,
    Rcpp::Named("centering_out") = center,
    Rcpp::Named("build_out") = build
  );

  if (k == 1) {
    List blocks = center["blocks"];
    List block0 = blocks[0];
    NumericVector betastar = block0["b_post_mean"];
    ret["out"] = block_beta_to_out_k1(sim, p_re_eff, n);
    ret["betastar"] = betastar;
  }

  return ret;
}

// ============================================================
// Independent-block separable-overbound sampler (Ind variants)
// Appendix A of BLOCK_ING_RINDEPNORMALGAMMA_REG.md.
//
// These functions implement the valid independent proposal where
//   lg_ind(j1,j2)      = sum_t max(pf_upp_t, pf_low_t)   [>= lg_joint]
//   ub2min_ind(j1,j2)  = sum_t min(ub2_low_t, ub2_upp_t) [<= ub2min_joint]
// so PLSD_ind(j1,j2) = PLSD_1(j1) * PLSD_2(j2) (faces drawn independently).
//
// The existing joint-product functions (BlockEnvelopeDispersionBuild,
// BlockEnvelopeSim, rIndepNormalGammaRegBlock) are NOT modified.
// ============================================================

List BlockEnvelopeDispersionBuildInd(
    const List& build_out,
    const List& centering_out,
    NumericVector y,
    NumericMatrix x,
    SEXP block,
    NumericVector offset,
    NumericVector wt,
    double shape,
    double rate,
    double max_disp_perc,
    Nullable<double> disp_lower,
    Nullable<double> disp_upper,
    double RSS_ML,
    bool use_parallel,
    bool verbose
) {
  if (!build_out.containsElementNamed("block_envelopes") ||
      !build_out.containsElementNamed("block_standardization") ||
      !build_out.containsElementNamed("meta")) {
    Rcpp::stop("'build_out' must be a full BlockEnvelopeBuild return list.");
  }
  if (!centering_out.containsElementNamed("dispersion") ||
      !centering_out.containsElementNamed("blocks") ||
      !centering_out.containsElementNamed("group_info") ||
      !centering_out.containsElementNamed("RSS_post") ||
      !centering_out.containsElementNamed("shape2") ||
      !centering_out.containsElementNamed("rate3")) {
    Rcpp::stop("'centering_out' must be a full BlockEnvelopeCentering return list.");
  }

  List block_envelopes = Rcpp::as<List>(build_out["block_envelopes"]);
  block_envelopes = Rcpp::clone(block_envelopes);
  List block_standardization = build_out["block_standardization"];
  List meta = build_out["meta"];

  const double RSS_post_global = Rcpp::as<double>(centering_out["RSS_post"]);
  const double shape2_global = Rcpp::as<double>(centering_out["shape2"]);
  const double rate3_global = Rcpp::as<double>(centering_out["rate3"]);
  const double n_w_global = Rcpp::as<double>(centering_out["n_w"]);
  const int l1 = Rcpp::as<int>(centering_out["l1"]);
  const int l2 = Rcpp::as<int>(centering_out["l2"]);
  const int k = Rcpp::as<int>(centering_out["k"]);

  if (y.size() != l2) {
    Rcpp::stop("length(y) must match centering_out$l2.");
  }
  if (x.nrow() != l2 || x.ncol() != l1) {
    Rcpp::stop("dim(x) must match centering_out$l2 x l1.");
  }

  if (max_disp_perc <= 0.0 || max_disp_perc >= 1.0) {
    if (centering_out.containsElementNamed("max_disp_perc")) {
      max_disp_perc = Rcpp::as<double>(centering_out["max_disp_perc"]);
    } else {
      max_disp_perc = 0.99;
    }
  }

  if (!disp_lower.isUsable() && centering_out.containsElementNamed("disp_lower")) {
    SEXP dl = centering_out["disp_lower"];
    if (!Rf_isNull(dl)) {
      disp_lower = Nullable<double>(Rcpp::wrap(Rcpp::as<double>(dl)));
    }
  }
  if (!disp_upper.isUsable() && centering_out.containsElementNamed("disp_upper")) {
    SEXP du = centering_out["disp_upper"];
    if (!Rf_isNull(du)) {
      disp_upper = Nullable<double>(Rcpp::wrap(Rcpp::as<double>(du)));
    }
  }

  List group_info = centering_out["group_info"];
  List row_blocks = group_info["rows"];
  List blocks_centering = centering_out["blocks"];

  List group_info_data = glmbayes::sim::group::normalize_group_cpp(block, l2);
  if (Rcpp::as<int>(group_info_data["k"]) != k) {
    Rcpp::stop("'block' partition k must match centering_out$k.");
  }

  if (offset.size() == 1) offset = Rcpp::rep(offset[0], l2);
  if (wt.size() == 1) wt = Rcpp::rep(wt[0], l2);
  if (offset.size() != l2) {
    Rcpp::stop("length(offset) must be 1 or length(y).");
  }
  if (wt.size() != l2) {
    Rcpp::stop("length(wt) must be 1 or length(y).");
  }

  double low = 0.0;
  double upp = 0.0;
  compute_sigma2_bounds_cpp(
    shape2_global, rate3_global, max_disp_perc,
    disp_lower, disp_upper, &low, &upp
  );

  List block_dispersion(k);
  IntegerVector identifiable_idx;
  double RSS_Min_global = 0.0;
  double RSS_ML_global = 0.0;
  bool rss_ml_provided = R_finite(RSS_ML);
  if (rss_ml_provided) {
    RSS_ML_global = RSS_ML;
  }

  List single_gamma_list;
  List single_ub_list;

  for (int j = 0; j < k; ++j) {
    List bc = blocks_centering[j];
    const bool identifiable = Rcpp::as<bool>(bc["identifiable"]);
    const std::string block_id = Rcpp::as<std::string>(bc["id"]);

    if (!identifiable) {
      block_dispersion[j] = R_NilValue;
      continue;
    }

    IntegerVector rows = row_blocks[j];
    List std_j = block_standardization[j];
    List Env_j = block_envelopes[j];

    NumericVector y_j = slice_numeric(y, rows);
    NumericMatrix x2_j = std_j["x2"];
    NumericMatrix P2_j = std_j["P2"];
    NumericMatrix mu2_j = std_j["mu2"];
    NumericVector alpha_j = std_j["alpha"];
    NumericVector wt_j = slice_numeric(wt, rows);
    const double RSS_post_j = Rcpp::as<double>(bc["RSS_post"]);

    double RSS_ML_j = NA_REAL;
    if (rss_ml_provided) {
      RSS_ML_j = RSS_ML;
    }

    List one = block_envelope_dispersion_one(
      Env_j,
      shape,
      rate,
      P2_j,
      y_j,
      x2_j,
      mu2_j,
      alpha_j,
      wt_j,
      RSS_post_j,
      RSS_ML_j,
      max_disp_perc,
      disp_lower,
      disp_upper,
      use_parallel,
      verbose,
      block_id
    );

    block_envelopes[j] = one["Env_out"];
    List ub_j = one["UB_list"];
    List gamma_j = one["gamma_list"];
    block_dispersion[j] = List::create(
      Rcpp::Named("cache") = one["cache"],
      Rcpp::Named("lg_prob_factor") = one["lg_prob_factor"],
      Rcpp::Named("UB2min") = one["UB2min"],
      Rcpp::Named("UB_list") = ub_j,
      Rcpp::Named("gamma_list") = gamma_j,
      Rcpp::Named("rss_min_global") = one["rss_min_global"],
      Rcpp::Named("RSS_ML") = one["RSS_ML"],
      Rcpp::Named("diagnostics") = one["diagnostics"],
      Rcpp::Named("block_id") = block_id,
      Rcpp::Named("identifiable") = true
    );

    RSS_Min_global += Rcpp::as<double>(one["rss_min_global"]);
    if (!rss_ml_provided) {
      RSS_ML_global += Rcpp::as<double>(one["RSS_ML"]);
    }
    if (identifiable_idx.size() == 0) {
      single_gamma_list = gamma_j;
      single_ub_list = ub_j;
    }

    identifiable_idx.push_back(j);
  }

  if (identifiable_idx.size() < 1) {
    Rcpp::stop("BlockEnvelopeDispersionBuildInd: no identifiable blocks.");
  }

  // Ind path: use O(sum gs_t) scalar shortcut — avoids materializing the
  // O(prod gs_t) product-face arrays that overflow for large k.
  List global_constants = build_global_dispersion_constants_scalar(
    identifiable_idx.size(),
    single_gamma_list,
    single_ub_list,
    identifiable_idx,
    block_envelopes,
    block_standardization,
    shape2_global,
    rate,
    rate3_global,
    n_w_global,
    RSS_Min_global,
    RSS_ML_global,
    low,
    upp,
    RSS_post_global
  );

  List gamma_list = global_constants["gamma_list"];
  List UB_list_global = global_constants["UB_list"];

  double prob_max_upp = Rcpp::as<double>(UB_list_global["max_New_LL_UB"]);
  double prob_max_low = NA_REAL;

  // Independent variant: when k > 1 we still need the apprx patch so that
  // per-block ub2_at_low / ub2_at_upp are in the block_dispersion entries
  // (they are used by block_ar_accumulate_one at draw time).  We do NOT call
  // build_joint_product_face_slack — no O(prod gs_t) enumeration required.
  if (identifiable_idx.size() > 1) {
    List apprx = compute_block_face_apprx_and_prob_anchors(
      identifiable_idx,
      block_envelopes,
      block_standardization,
      block_dispersion,
      shape2_global,
      rate3_global,
      low,
      upp
    );
    patch_block_dispersion_apprx(block_dispersion, apprx);
    if (global_constants.containsElementNamed("prob_max_low")) {
      prob_max_low = Rcpp::as<double>(global_constants["prob_max_low"]);
    } else {
      prob_max_low = Rcpp::as<double>(apprx["prob_max_low"]);
    }
  }

  RObject disp_lower_out = disp_lower.isUsable()
    ? RObject(Rcpp::wrap(Rcpp::as<double>(disp_lower.get())))
    : RObject(R_NilValue);
  RObject disp_upper_out = disp_upper.isUsable()
    ? RObject(Rcpp::wrap(Rcpp::as<double>(disp_upper.get())))
    : RObject(R_NilValue);

  IntegerVector gs_per_block(identifiable_idx.size());
  for (int t = 0; t < identifiable_idx.size(); ++t) {
    List env_t = block_envelopes[identifiable_idx[t]];
    NumericMatrix cbars_t = env_t["cbars"];
    gs_per_block[t] = cbars_t.nrow();
  }

  List dispersion_envelope = List::create(
    Rcpp::Named("gamma_list") = gamma_list,
    Rcpp::Named("UB_list") = UB_list_global,
    Rcpp::Named("block_dispersion") = block_dispersion,
    Rcpp::Named("cross_face_meta") = List::create(
      Rcpp::Named("gs_per_block") = gs_per_block,
      Rcpp::Named("identifiable_idx") = identifiable_idx,
      Rcpp::Named("n_identifiable") = identifiable_idx.size(),
      Rcpp::Named("aggregation") = global_constants["source"],
      Rcpp::Named("ub3a_lg_mode") = (identifiable_idx.size() > 1)
        ? "per_block_lg_sum_v1"
        : "single_block_lg",
      Rcpp::Named("ub2_mode") = (identifiable_idx.size() > 1)
        ? "per_block_ub2min_sum_v1"
        : "single_block_ub2min",
      Rcpp::Named("face_draw_mode") = (identifiable_idx.size() > 1)
        ? "per_block_plsd_ind_v1"
        : "per_block_plsd"
    ),
    Rcpp::Named("prob_max_upp") = prob_max_upp,
    Rcpp::Named("prob_max_low") = prob_max_low,
    Rcpp::Named("low") = low,
    Rcpp::Named("upp") = upp,
    Rcpp::Named("RSS_post") = RSS_post_global,
    Rcpp::Named("RSS_ML") = RSS_ML_global,
    Rcpp::Named("RSS_Min") = RSS_Min_global,
    Rcpp::Named("status") = "v1"
  );

  List meta_out = Rcpp::clone(meta);
  meta_out["dispersion_envelope_status"] = "v1";

  List build_out_patched = List::create(
    Rcpp::Named("block_envelopes") = block_envelopes,
    Rcpp::Named("dispersion_envelope") = dispersion_envelope,
    Rcpp::Named("block_standardization") = block_standardization,
    Rcpp::Named("meta") = meta_out
  );

  if (verbose) {
    Rcpp::Rcout << "BlockEnvelopeDispersionBuildInd: k=" << k
                << " n_identifiable=" << identifiable_idx.size()
                << " RSS_Min=" << RSS_Min_global
                << " low=" << low << " upp=" << upp
                << " (independent per-block proposal)" << std::endl;
  }

  return List::create(
    Rcpp::Named("build_out") = build_out_patched,
    Rcpp::Named("dispersion_envelope") = dispersion_envelope,
    Rcpp::Named("low") = low,
    Rcpp::Named("upp") = upp,
    Rcpp::Named("diagnostics") = List::create(
      Rcpp::Named("RSS_Min_global") = RSS_Min_global,
      Rcpp::Named("RSS_ML_global") = RSS_ML_global,
      Rcpp::Named("n_identifiable") = identifiable_idx.size(),
      Rcpp::Named("aggregation") = global_constants["source"]
    )
  );
}

// Independent-block simulator. Identical to BlockEnvelopeSim except that
// use_joint_lg_ub3a is always false: faces are drawn independently from each
// block's own PLSD and per-block lg_prob_factor / UB2min are summed, using
// the separable overbounds proved valid in Appendix A of the README.
List BlockEnvelopeSimInd(
    const List& build_out,
    int n,
    bool progbar,
    bool verbose
) {
  if (n < 1) {
    Rcpp::stop("'n' must be at least 1.");
  }
  if (!build_out.containsElementNamed("block_envelopes") ||
      !build_out.containsElementNamed("block_standardization") ||
      !build_out.containsElementNamed("meta")) {
    Rcpp::stop("'build_out' must be a full BlockEnvelopeBuild return list.");
  }

  List block_envelopes = build_out["block_envelopes"];
  List block_standardization = build_out["block_standardization"];
  List meta = build_out["meta"];
  if (!meta.containsElementNamed("group_info") ||
      !meta.containsElementNamed("l1") ||
      !meta.containsElementNamed("dispersion")) {
    Rcpp::stop("'build_out$meta' must contain group_info, l1, and dispersion.");
  }

  List group_info = meta["group_info"];
  const int k = Rcpp::as<int>(group_info["k"]);
  const int l1 = Rcpp::as<int>(meta["l1"]);
  const double dispersion_anchor = Rcpp::as<double>(meta["dispersion"]);
  CharacterVector ids = group_info["ids"];

  if (block_envelopes.size() != k ||
      block_standardization.size() != k) {
    Rcpp::stop("block_envelopes and block_standardization must have length k.");
  }

  const bool has_disp_env =
    build_out.containsElementNamed("dispersion_envelope") &&
    !Rf_isNull(build_out["dispersion_envelope"]);
  std::string disp_status = "missing";
  if (has_disp_env) {
    List dispersion_envelope = build_out["dispersion_envelope"];
    if (dispersion_envelope.containsElementNamed("status") &&
        !Rf_isNull(dispersion_envelope["status"])) {
      disp_status = Rcpp::as<std::string>(dispersion_envelope["status"]);
    }
  }
  const bool use_ar_compute = (disp_status == "v1");

  List group_results(k);
  NumericVector dispersion_out(n);
  NumericVector iters_out(n);

  if (!use_ar_compute) {
    std::fill(dispersion_out.begin(), dispersion_out.end(), dispersion_anchor);
    std::fill(iters_out.begin(), iters_out.end(), 1.0);
    for (int j = 0; j < k; ++j) {
      if (progbar && verbose) {
        Rcpp::Rcout << "BlockEnvelopeSimInd block " << (j + 1) << "/" << k << std::endl;
      }
      const std::string block_id = Rcpp::as<std::string>(ids[j]);
      SEXP env_j = block_envelopes[j];
      List env_list = Rf_isNull(env_j) ? List() : List(env_j);
      List std_j = block_standardization[j];
      group_results[j] = block_envelope_sim_one(
        env_list, std_j, block_id, l1, n
      );
    }
    List meta_out = Rcpp::clone(meta);
    meta_out["accept_mode"] = "auto_v1";
    meta_out["n"] = n;
    return List::create(
      Rcpp::Named("group_results") = group_results,
      Rcpp::Named("dispersion") = dispersion_out,
      Rcpp::Named("iters_out") = iters_out,
      Rcpp::Named("group_info") = group_info,
      Rcpp::Named("meta") = meta_out
    );
  }

  List dispersion_envelope = build_out["dispersion_envelope"];
  List gamma_list = dispersion_envelope["gamma_list"];
  List UB_global = dispersion_envelope["UB_list"];
  List block_dispersion = dispersion_envelope["block_dispersion"];

  const double shape3 = Rcpp::as<double>(gamma_list["shape3"]);
  const double rate2 = Rcpp::as<double>(gamma_list["rate2"]);
  const double disp_lower = Rcpp::as<double>(gamma_list["disp_lower"]);
  const double disp_upper = Rcpp::as<double>(gamma_list["disp_upper"]);
  const double RSS_Min_G = Rcpp::as<double>(dispersion_envelope["RSS_Min"]);
  const double max_New_LL_UB = Rcpp::as<double>(UB_global["max_New_LL_UB"]);
  const double max_LL_log_disp = Rcpp::as<double>(UB_global["max_LL_log_disp"]);
  const double lm_log1 = Rcpp::as<double>(UB_global["lm_log1"]);
  const double lm_log2 = Rcpp::as<double>(UB_global["lm_log2"]);
  const double lmc1 = Rcpp::as<double>(UB_global["lmc1"]);
  const double lmc2 = Rcpp::as<double>(UB_global["lmc2"]);

  // Independent variant: always draw from each block's own PLSD.
  // No joint product-face tables are read or required.
  std::vector<std::vector<double>> block_plsd_cdf(k);

  for (int j = 0; j < k; ++j) {
    List std_j = block_standardization[j];
    const bool prior_only = Rcpp::as<bool>(std_j["prior_only"]);
    const std::string block_id = Rcpp::as<std::string>(ids[j]);

    if (prior_only || Rf_isNull(block_envelopes[j])) {
      NumericVector mu_j = std_j["mu"];
      NumericMatrix beta_block(l1, n);
      NumericVector block_iters(n);
      for (int i = 0; i < n; ++i) {
        beta_block(Rcpp::_, i) = mu_j;
        block_iters[i] = 1.0;
      }
      group_results[j] = List::create(
        Rcpp::Named("block_id") = block_id,
        Rcpp::Named("identifiable") = false,
        Rcpp::Named("beta") = beta_block,
        Rcpp::Named("iters_out") = block_iters,
        Rcpp::Named("prior_only") = true
      );
    } else {
      List Env = block_envelopes[j];
      block_plsd_cdf[j] = build_face_cdf(Env["PLSD"]);
      group_results[j] = List::create(
        Rcpp::Named("block_id") = block_id,
        Rcpp::Named("identifiable") = true,
        Rcpp::Named("beta") = NumericMatrix(l1, n),
        Rcpp::Named("iters_out") = NumericVector(n),
        Rcpp::Named("prior_only") = false,
        Rcpp::Named("face_J_last") = NA_INTEGER
      );
    }
  }

  std::vector<int> J_draw(k, 0);
  std::vector<NumericVector> beta_std_draw(k);
  std::vector<NumericVector> beta_orig_draw(k);

  for (int i = 0; i < n; ++i) {
    Rcpp::checkUserInterrupt();
    if (progbar && n > 1) {
      progress_bar(static_cast<double>(i), static_cast<double>(n - 1));
      if (i == n - 1) {
        Rcpp::Rcout << "" << std::endl;
      }
    }

    iters_out[i] = 1.0;
    int accept = 0;
    double sigma2_accept = NA_REAL;
    const auto ar_t0 = std::chrono::steady_clock::now();

    while (accept == 0) {
      // Draw faces independently from each block's own PLSD.
      for (int j = 0; j < k; ++j) {
        List std_j = block_standardization[j];
        if (Rcpp::as<bool>(std_j["prior_only"]) || Rf_isNull(block_envelopes[j])) {
          beta_orig_draw[j] = std_j["mu"];
          J_draw[j] = -1;
          continue;
        }
        List Env = block_envelopes[j];
        J_draw[j] = draw_face_index_from_cdf(block_plsd_cdf[j]);
        beta_std_draw[j] = draw_beta_std_face(Env, J_draw[j], l1);
        if (!is_true(all(is_finite(beta_std_draw[j])))) {
          // A non-finite draw here means rnorm_ct() produced NaN/Inf for
          // this face's truncated-normal proposal (most likely log-space
          // cancellation for an extreme/degenerate [loglt, logrt] window).
          // This is a real numerical failure, not something to paper over
          // by resampling: report the offending face so it can be
          // diagnosed, rather than silently retrying or letting NaN
          // propagate into the downstream Cholesky solves.
          const std::string block_id_j =
            Env.containsElementNamed("block_id")
              ? Rcpp::as<std::string>(Env["block_id"])
              : std::string("<unknown>");
          NumericMatrix loglt_j = Env["loglt"];
          NumericMatrix logrt_j = Env["logrt"];
          NumericMatrix cbars_j = Env["cbars"];
          std::ostringstream diag;
          diag << "BlockEnvelopeSimInd: non-finite beta_std draw for block id="
               << block_id_j << ", draw i=" << (i + 1) << ", face J=" << J_draw[j]
               << " (l1=" << l1 << "). beta_std=[";
          for (int t = 0; t < l1; ++t) {
            if (t > 0) diag << ", ";
            diag << beta_std_draw[j][t];
          }
          diag << "], loglt(J,.)=[";
          for (int t = 0; t < l1; ++t) {
            if (t > 0) diag << ", ";
            diag << loglt_j(J_draw[j], t);
          }
          diag << "], logrt(J,.)=[";
          for (int t = 0; t < l1; ++t) {
            if (t > 0) diag << ", ";
            diag << logrt_j(J_draw[j], t);
          }
          diag << "], cbars(J,.)=[";
          for (int t = 0; t < l1; ++t) {
            if (t > 0) diag << ", ";
            diag << cbars_j(J_draw[j], t);
          }
          diag << "].";
          Rcpp::stop(diag.str());
        }
        beta_orig_draw[j] = unstandardize_beta_one(
          beta_std_draw[j],
          std_j["L2Inv"],
          std_j["L3Inv"],
          std_j["mu"]
        );
      }

      const double sigma2 =
        rinvgamma_ct_safe(shape3, rate2, disp_upper, disp_lower);
      if (!R_finite(sigma2)) {
        Rcpp::stop("BlockEnvelopeSimInd: non-finite sigma2 draw.");
      }

      double LL_total = 0.0;
      double UB1_total = 0.0;
      double quad_sum = 0.0;
      double ub3a_block_sum = 0.0;
      double UB2min_used = 0.0;

      // Accumulate per-block: lg_prob_factor and UB2min are summed directly
      // (separable overbound — Lemmas A1 and A2 in the README appendix).
      for (int j = 0; j < k; ++j) {
        List std_j = block_standardization[j];
        if (Rcpp::as<bool>(std_j["prior_only"]) || Rf_isNull(block_envelopes[j])) {
          continue;
        }
        List Env = block_envelopes[j];
        List block_disp_j = block_dispersion[j];
        block_ar_accumulate_one(
          J_draw[j], beta_std_draw[j], sigma2, Env, std_j, block_disp_j,
          &LL_total, &UB1_total, &quad_sum, &ub3a_block_sum,
          true  // always accumulate per-block lg_prob_factor
        );
        NumericVector UB2min_j = block_disp_j["UB2min"];
        UB2min_used += UB2min_j[J_draw[j]];
      }

      const double UB2_raw = 0.5 * (1.0 / sigma2) * (quad_sum - RSS_Min_G);
      const double UB2 = UB2_raw - UB2min_used;
      const double UB3A = ub3a_block_sum + lmc1 + lmc2 * sigma2;
      const double New_LL_log_disp = lm_log1 + lm_log2 * std::log(sigma2);
      const double UB3B =
        (max_New_LL_UB - max_LL_log_disp + New_LL_log_disp) -
        (lmc1 + lmc2 * sigma2);

      const double test1 = LL_total - UB1_total;
      double test = test1 - (UB2 + UB3A + UB3B);

      block_ar_check_signs(
        test1, UB1_total, UB2, UB3A, UB3B, test,
        quad_sum, RSS_Min_G, UB2_raw, UB2min_used, verbose
      );

      const double U2 = runif_safe();
      test -= std::log(U2);

      if (test >= 0.0) {
        accept = 1;
        sigma2_accept = sigma2;

        if (verbose && i == 0) {
          Rcpp::Rcout << "BlockEnvelopeSimInd resample_until_accept draw 0: sigma2="
                      << sigma2 << " test=" << test
                      << " iters=" << iters_out[i]
                      << " LL=" << LL_total << " UB1=" << UB1_total
                      << " UB2=" << UB2 << " UB3A=" << UB3A
                      << " UB3B=" << UB3B << std::endl;
        }
      } else {
        iters_out[i] = iters_out[i] + 1.0;
        const double ar_elapsed = block_ar_elapsed_seconds(ar_t0);
        if (ar_elapsed > kBlockArTimeLimitSec) {
          accept = 1;
          sigma2_accept = sigma2;
          block_ar_warn_time_limit_accept(
            "BlockEnvelopeSimInd",
            i,
            iters_out[i],
            ar_elapsed,
            kBlockArTimeLimitSec
          );
        }
      }
    }

    dispersion_out[i] = sigma2_accept;

    for (int j = 0; j < k; ++j) {
      List std_j = block_standardization[j];
      if (Rcpp::as<bool>(std_j["prior_only"]) || Rf_isNull(block_envelopes[j])) {
        List res_j = group_results[j];
        NumericMatrix beta_block = res_j["beta"];
        beta_block(Rcpp::_, i) = beta_orig_draw[j];
        res_j["beta"] = beta_block;
        NumericVector block_iters = res_j["iters_out"];
        block_iters[i] = iters_out[i];
        res_j["iters_out"] = block_iters;
        group_results[j] = res_j;
        continue;
      }
      List res_j = group_results[j];
      NumericMatrix beta_block = res_j["beta"];
      beta_block(Rcpp::_, i) = beta_orig_draw[j];
      res_j["beta"] = beta_block;
      NumericVector block_iters = res_j["iters_out"];
      block_iters[i] = iters_out[i];
      res_j["iters_out"] = block_iters;
      res_j["face_J_last"] = J_draw[j];
      group_results[j] = res_j;
    }
  }

  List meta_out = Rcpp::clone(meta);
  meta_out["accept_mode"] = "resample_until_accept_ind_v1";
  meta_out["face_draw_mode"] = "per_block_plsd_ind_v1";
  meta_out["n"] = n;

  return List::create(
    Rcpp::Named("group_results") = group_results,
    Rcpp::Named("dispersion") = dispersion_out,
    Rcpp::Named("iters_out") = iters_out,
    Rcpp::Named("group_info") = group_info,
    Rcpp::Named("meta") = meta_out
  );
}

// Orchestrator for the independent-block variant:
//   BlockEnvelopeCentering → BlockEnvelopeBuild
//   → BlockEnvelopeDispersionBuildInd → BlockEnvelopeSimInd
List rIndepNormalGammaRegBlockInd(
    int n,
    NumericVector y,
    NumericMatrix x,
    SEXP block,
    SEXP prior_list_sexp,
    SEXP prior_lists_sexp,
    NumericVector offset,
    NumericVector wt,
    int p_re,
    int n_rss_iter,
    int Gridtype,
    Nullable<int> n_envopt,
    double RSS_ML,
    bool use_parallel,
    bool use_opencl,
    bool progbar,
    bool verbose,
    CharacterVector group_levels,
    CharacterVector re_names
) {
  if (n < 1) {
    Rcpp::stop("'n' must be at least 1.");
  }
  if (n_rss_iter < 1) {
    Rcpp::stop("'n_rss_iter' must be at least 1.");
  }

  double shape = NA_REAL;
  double rate = NA_REAL;
  double max_disp_perc = 0.99;
  Nullable<double> disp_lower;
  Nullable<double> disp_upper;
  prior_hyperparams_from_list(
    prior_list_sexp, &shape, &rate, &max_disp_perc, &disp_lower, &disp_upper
  );

  List center = BlockEnvelopeCentering(
    y, x, block, prior_list_sexp, prior_lists_sexp,
    offset, wt, shape, rate, max_disp_perc,
    disp_lower, disp_upper, p_re, n_rss_iter, verbose
  );

  List build = BlockEnvelopeBuild(
    center, y, x, block, prior_list_sexp, prior_lists_sexp,
    offset, wt, max_disp_perc, disp_lower, disp_upper,
    n, Gridtype, n_envopt, RSS_ML, use_parallel, use_opencl, verbose
  );

  List disp = BlockEnvelopeDispersionBuildInd(
    build, center, y, x, block, offset, wt,
    shape, rate, max_disp_perc, disp_lower, disp_upper,
    RSS_ML, use_parallel, verbose
  );
  build = disp["build_out"];

  List sim = BlockEnvelopeSimInd(build, n, progbar, verbose);

  const int p_re_eff = Rcpp::as<int>(center["p_re"]);
  const int k = Rcpp::as<int>(center["k"]);
  NumericVector dispersion = sim["dispersion"];
  NumericVector iters_out = sim["iters_out"];

  NumericMatrix b_draw = map_block_sim_to_b_draw(
    sim, group_levels, p_re_eff, 0
  );
  if (re_names.size() == p_re_eff) {
    colnames(b_draw) = re_names;
  } else if (re_names.size() > 0 && re_names.size() != p_re_eff) {
    Rcpp::stop("'re_names' must have length 0 or p_re.");
  }
  if (group_levels.size() == b_draw.nrow()) {
    rownames(b_draw) = group_levels;
  }

  NumericVector weight_out(n);
  std::fill(weight_out.begin(), weight_out.end(), 1.0);

  const double low = Rcpp::as<double>(disp["low"]);
  const double upp = Rcpp::as<double>(disp["upp"]);

  List ret = List::create(
    Rcpp::Named("b") = b_draw,
    Rcpp::Named("dispersion_ranef") = dispersion[0],
    Rcpp::Named("iters_mean") = iters_out[0],
    Rcpp::Named("disp_out") = dispersion,
    Rcpp::Named("iters_out") = iters_out,
    Rcpp::Named("weight_out") = weight_out,
    Rcpp::Named("low") = low,
    Rcpp::Named("upp") = upp,
    Rcpp::Named("k") = k,
    Rcpp::Named("p_re") = p_re_eff,
    Rcpp::Named("sim") = sim,
    Rcpp::Named("centering_out") = center,
    Rcpp::Named("build_out") = build
  );

  if (k == 1) {
    List blocks = center["blocks"];
    List block0 = blocks[0];
    NumericVector betastar = block0["b_post_mean"];
    ret["out"] = block_beta_to_out_k1(sim, p_re_eff, n);
    ret["betastar"] = betastar;
  }

  return ret;
}

}  // namespace group
}  // namespace env
}  // namespace glmbayes
