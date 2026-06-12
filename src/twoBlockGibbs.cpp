// twoBlockGibbs.cpp
// C++ port of the two_block_rNormal_reg() Gibbs loop (port-only).
//
// Mirrors the R loop in R/two_block_rNormal_reg.R verbatim:
//   for (i in 1:n) {
//     fixef <- fixef_start                          # "replicate" sampling
//     for (m in 1:m_convergence) {
//       mu_all  <- .two_block_mu_all(fixef, x_hyper, re_names, group_levels)
//       Block 1 <- block_rNormalReg() / block_rNormalGLM()   (n = 1)
//       Block 2 <- multi_rNormal_reg() list-x branch (rNormal_reg gaussian
//                  per RE component, n = 1)
//     }
//   }
// Block 1 goes through the existing block_rNormalReg_cpp_export /
// block_rNormalGLM_cpp_export (per-iteration prior normalization unchanged).
// Block 2 replicates rNormal_reg()'s gaussian-branch prep per component per
// iteration (Sigma -> P via inv_sympd == chol2inv(chol()), PD checks), then
// calls the rNormalReg core.  f2/f3 R closures and the optim call inside
// rNormalGLM are untouched.

#include "RcppArmadillo.h"
#include "simfuncs.h"
#include "progress_utils.h"

#include <string>
#include <vector>

namespace glmbayes {
namespace sim {

namespace {

using Rcpp::CharacterVector;
using Rcpp::Function;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

// Mirror R's is.null(): a present-but-NULL element counts as absent.
inline bool has_non_null(const List& pl, const char* name) {
  return pl.containsElementNamed(name) && !Rf_isNull(pl[name]);
}

inline bool is_symmetric_mat(const NumericMatrix& M, double tol = 1e-8) {
  if (M.nrow() != M.ncol()) return false;
  for (int i = 0; i < M.nrow(); ++i) {
    for (int j = i + 1; j < M.ncol(); ++j) {
      if (std::fabs(M(i, j) - M(j, i)) > tol) return false;
    }
  }
  return true;
}

// R .check_symmetric_pd() / rNormal_reg PD check: eigen, tol 1e-6.
inline void check_pd(const NumericMatrix& M, const char* label) {
  arma::mat M2(const_cast<double*>(M.begin()), M.nrow(), M.ncol(), false);
  arma::vec ev = arma::eig_sym(M2);
  const double tol = 1e-6;
  if (ev.min() < -tol * std::fabs(ev.max())) {
    Rcpp::stop("'%s' is not positive definite.", label);
  }
}

// Port of .two_block_mu_all(): mu_all[i, j] = sum(X_k[row_j, ] * gamma_k),
// where row_j is the positional row j when X_k has no rownames, or the row
// named group_levels[j] otherwise.  Row lookups are resolved once up front
// (x_hyper and group_levels are constant across iterations).
struct MuAllBuilder {
  int p_re;
  int J;
  std::vector<NumericMatrix> X;          // per-component J x q_k design
  std::vector<std::vector<int>> row_idx; // per-component 0-based row for group j

  MuAllBuilder(const List& x_hyper, const CharacterVector& group_levels) {
    p_re = x_hyper.size();
    J = group_levels.size();
    X.reserve(p_re);
    row_idx.resize(p_re);
    for (int i = 0; i < p_re; ++i) {
      NumericMatrix X_k = Rcpp::as<NumericMatrix>(x_hyper[i]);
      X.push_back(X_k);

      CharacterVector rn;
      SEXP dn = X_k.attr("dimnames");
      if (!Rf_isNull(dn)) {
        List dnl(dn);
        if (dnl.size() >= 1 && !Rf_isNull(dnl[0])) {
          rn = CharacterVector(dnl[0]);
        }
      }

      std::vector<int>& idx = row_idx[i];
      idx.resize(J);
      if (rn.size() == 0) {
        if (X_k.nrow() != J) {
          Rcpp::stop(
            "nrow(x_hyper[[%d]]) must equal length(group_levels).", i + 1
          );
        }
        for (int j = 0; j < J; ++j) idx[j] = j;
      } else {
        for (int j = 0; j < J; ++j) {
          const std::string lev = Rcpp::as<std::string>(group_levels[j]);
          int found = -1;
          for (int r = 0; r < rn.size(); ++r) {
            if (!CharacterVector::is_na(rn[r]) &&
                Rcpp::as<std::string>(rn[r]) == lev) {
              found = r;
              break;
            }
          }
          if (found < 0) {
            Rcpp::stop(
              "x_hyper[[%d]] has no row named \"%s\" (group_levels[%d]).",
              i + 1, lev.c_str(), j + 1
            );
          }
          idx[j] = found;
        }
      }
    }
  }

  NumericMatrix build(const std::vector<NumericVector>& fixef) const {
    NumericMatrix mu_all(p_re, J);
    for (int i = 0; i < p_re; ++i) {
      const NumericMatrix& X_k = X[i];
      const NumericVector& gamma_k = fixef[i];
      if (gamma_k.size() != X_k.ncol()) {
        Rcpp::stop(
          "length(fixef[[%d]]) (%d) must equal ncol(x_hyper[[%d]]) (%d).",
          i + 1, gamma_k.size(), i + 1, X_k.ncol()
        );
      }
      const std::vector<int>& idx = row_idx[i];
      for (int j = 0; j < J; ++j) {
        const int r = idx[j];
        double s = 0.0;
        for (int c = 0; c < X_k.ncol(); ++c) {
          s += X_k(r, c) * gamma_k[c];
        }
        mu_all(i, j) = s;
      }
    }
    return mu_all;
  }
};

// Port of .two_block_block1_prior_list(): list(mu, dispersion, ddef, P/Sigma).
// dispersion/ddef are forwarded as-is (possibly NULL); the block exports
// treat present-but-NULL as absent, matching R semantics.
List block1_prior_list(
    const NumericMatrix& mu_all,
    const List& prior_list_block1,
    SEXP dispersion_block1,
    SEXP ddef_block1
) {
  List out = List::create(
    Rcpp::Named("mu") = mu_all,
    Rcpp::Named("dispersion") = dispersion_block1,
    Rcpp::Named("ddef") = ddef_block1
  );
  if (has_non_null(prior_list_block1, "P")) {
    out["P"] = prior_list_block1["P"];
  }
  if (has_non_null(prior_list_block1, "Sigma")) {
    out["Sigma"] = prior_list_block1["Sigma"];
  }
  if (!has_non_null(out, "P") && !has_non_null(out, "Sigma")) {
    Rcpp::stop("prior_list_block1 must contain 'P' or 'Sigma'.");
  }
  return out;
}

// Port of .validate_normal_prior_list() + rNormal_reg() gaussian-branch prior
// prep for one Block 2 component: mu, P (from chol2inv(chol(Sigma)) when P is
// absent; inv_sympd is the same LAPACK path), required dispersion.
struct Block2Prior {
  NumericVector mu;
  NumericMatrix P;
  double dispersion;
};

Block2Prior block2_prior_prep(const List& pl, int j1 /*1-based*/, int p) {
  if (!has_non_null(pl, "mu")) {
    Rcpp::stop("prior_list[[%d]] must contain 'mu'.", j1);
  }
  if (!has_non_null(pl, "Sigma") && !has_non_null(pl, "P")) {
    Rcpp::stop("prior_list[[%d]] must contain 'Sigma' or 'P'.", j1);
  }

  NumericVector mu = Rcpp::as<NumericVector>(pl["mu"]);
  if (mu.size() != p) {
    Rcpp::stop("prior_list[[%d]]$mu must have length ncol(x) = %d.", j1, p);
  }

  NumericMatrix P;
  if (has_non_null(pl, "Sigma")) {
    NumericMatrix S = Rcpp::as<NumericMatrix>(pl["Sigma"]);
    if (S.nrow() != p || S.ncol() != p) {
      Rcpp::stop("prior_list[[%d]]$Sigma must be %d x %d.", j1, p, p);
    }
    if (!is_symmetric_mat(S)) {
      Rcpp::stop("prior_list[[%d]]$Sigma must be symmetric.", j1);
    }
    check_pd(S, "Sigma");
    if (!has_non_null(pl, "P")) {
      // rNormal_reg(): P <- 0.5 * (chol2inv(chol(Sigma)) + t(...))
      arma::mat S2(const_cast<double*>(S.begin()), p, p, false);
      arma::mat Pinv = arma::inv_sympd(S2);
      P = NumericMatrix(Rcpp::wrap(0.5 * (Pinv + Pinv.t())));
    }
  }
  if (has_non_null(pl, "P")) {
    P = Rcpp::as<NumericMatrix>(pl["P"]);
    if (P.nrow() != p || P.ncol() != p) {
      Rcpp::stop("prior_list[[%d]]$P must be %d x %d.", j1, p, p);
    }
  }
  if (!is_symmetric_mat(P)) {
    Rcpp::stop("matrix P must be symmetric");
  }
  check_pd(P, "P");

  // rNormal_reg(): gaussian requires an explicit dispersion (ddef rules).
  bool ddef;
  if (pl.containsElementNamed("ddef")) {
    SEXP dd = pl["ddef"];
    ddef = Rf_isLogical(dd) && Rf_length(dd) >= 1 &&
           LOGICAL(dd)[0] == TRUE;
  } else {
    ddef = !has_non_null(pl, "dispersion");
  }
  if (ddef || !has_non_null(pl, "dispersion")) {
    Rcpp::stop(
      "For gaussian() models, dNormal() requires an explicit dispersion "
      "(prior_list[[%d]]). Omitted or NULL dispersion is not allowed", j1
    );
  }
  double dispersion = Rcpp::as<NumericVector>(pl["dispersion"])[0];

  Block2Prior out;
  out.mu = mu;
  out.P = P;
  out.dispersion = dispersion;
  return out;
}

} // anonymous namespace

List two_block_rNormal_reg_cpp_export(
    int n,
    int m_convergence,
    const NumericVector& y,
    const NumericMatrix& x,
    SEXP block,
    const List& x_hyper,
    const List& prior_list_block1,
    SEXP dispersion_block1,
    SEXP ddef_block1,
    const List& prior_list_block2,
    const List& fixef_start,
    const CharacterVector& group_levels,
    const std::string& family,
    const std::string& link,
    const Function& f2,
    const Function& f3,
    const Function& f2_gauss,
    const Function& f3_gauss,
    const NumericVector& offset,
    const NumericVector& wt,
    int Gridtype,
    int n_envopt,
    bool use_parallel,
    bool use_opencl,
    bool verbose,
    bool progbar
) {
  if (n < 1) {
    Rcpp::stop("'n' must be at least 1.");
  }
  if (m_convergence < 1) {
    Rcpp::stop("'m_convergence' must be at least 1.");
  }
  const int p_re = x.ncol();
  const int J = group_levels.size();
  if (x_hyper.size() != p_re) {
    Rcpp::stop("length(x_hyper) must equal ncol(x) = %d.", p_re);
  }
  if (prior_list_block2.size() != p_re) {
    Rcpp::stop("length(prior_list_block2) must equal ncol(x) = %d.", p_re);
  }
  if (fixef_start.size() != p_re) {
    Rcpp::stop("length(fixef_start) must equal ncol(x) = %d.", p_re);
  }

  const bool is_gaussian = (family == "gaussian");

  MuAllBuilder mu_builder(x_hyper, group_levels);

  // fixef state (component order fixed by the R wrapper)
  std::vector<NumericVector> fixef_start_v(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_start_v[j] = Rcpp::as<NumericVector>(fixef_start[j]);
  }
  std::vector<NumericVector> fixef = fixef_start_v;

  // Storage: fixef draws per component (n x q_k); b draws as J x p_re x n.
  List fixef_draws(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_draws[j] = NumericMatrix(n, fixef_start_v[j].size());
  }
  NumericVector b_arr(Rcpp::Dimension(J, p_re, n));

  NumericMatrix mu_all;
  NumericMatrix b_i;
  CharacterVector group_ids;
  bool have_ids = false;

  for (int i = 0; i < n; ++i) {
    Rcpp::checkUserInterrupt();
    if (progbar) {
      glmbayes::progress::progress_bar(
        static_cast<double>(i + 1), static_cast<double>(n)
      );
    }

    fixef = fixef_start_v;

    for (int m = 0; m < m_convergence; ++m) {

      mu_all = mu_builder.build(fixef);
      List pl1 = block1_prior_list(
        mu_all, prior_list_block1, dispersion_block1, ddef_block1
      );

      List block_i;
      if (is_gaussian) {
        // R: block_rNormalReg(n = 1, ...) with default Gridtype = 2L.
        block_i = block_rNormalReg_cpp_export(
          1, y, x, block, pl1, R_NilValue, offset, wt, f2, f3, 2
        );
      } else {
        block_i = block_rNormalGLM_cpp_export(
          1, y, x, block, pl1, R_NilValue, offset, wt, f2, f3,
          family, link, Gridtype, n_envopt,
          use_parallel, use_opencl, verbose
        );
      }

      b_i = Rcpp::as<NumericMatrix>(block_i["coefficients"]);
      if (b_i.nrow() != J || b_i.ncol() != p_re) {
        Rcpp::stop(
          "Block 1 returned a %d x %d coefficient matrix; expected %d x %d.",
          b_i.nrow(), b_i.ncol(), J, p_re
        );
      }
      if (!have_ids) {
        List bi_info = block_i["block_info"];
        group_ids = Rcpp::as<CharacterVector>(bi_info["ids"]);
        have_ids = true;
      }

      // Block 2: multi_rNormal_reg() list-x branch, n = 1 per component.
      for (int j = 0; j < p_re; ++j) {
        const NumericMatrix& X_j = mu_builder.X[j];
        if (X_j.nrow() != b_i.nrow()) {
          Rcpp::stop(
            "nrow(x_hyper[[%d]]) (%d) must equal nrow(y) (%d).",
            j + 1, X_j.nrow(), b_i.nrow()
          );
        }
        Block2Prior pr = block2_prior_prep(
          List(prior_list_block2[j]), j + 1, X_j.ncol()
        );
        NumericVector y_j = b_i(Rcpp::_, j);
        NumericVector offset_j(X_j.nrow(), 0.0);
        NumericVector wt_j(X_j.nrow(), 1.0);

        List out_j = rNormalReg(
          1, y_j, X_j, pr.mu, pr.P, offset_j, wt_j,
          pr.dispersion, f2_gauss, f3_gauss, pr.mu,
          "gaussian", "identity", 2
        );
        NumericMatrix coef_j = Rcpp::as<NumericMatrix>(out_j["coefficients"]);
        fixef[j] = NumericVector(coef_j(0, Rcpp::_));
      }
    }

    for (int j = 0; j < p_re; ++j) {
      NumericMatrix fd = fixef_draws[j];
      const NumericVector& fj = fixef[j];
      if (fj.size() != fd.ncol()) {
        Rcpp::stop(
          "Block 2 draw for component %d has length %d; expected %d.",
          j + 1, fj.size(), fd.ncol()
        );
      }
      fd(i, Rcpp::_) = fj;
    }
    for (int j = 0; j < p_re; ++j) {
      for (int g = 0; g < J; ++g) {
        b_arr[g + J * (j + p_re * i)] = b_i(g, j);
      }
    }
  }

  List fixef_last(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_last[j] = fixef[j];
  }

  return List::create(
    Rcpp::Named("fixef_draws") = fixef_draws,
    Rcpp::Named("b_draws") = b_arr,
    Rcpp::Named("fixef_last") = fixef_last,
    Rcpp::Named("b_last") = b_i,
    Rcpp::Named("mu_all_last") = mu_all,
    Rcpp::Named("group_ids") = group_ids
  );
}

// ===========================================================================
// v2: pfamily-based Block 2 contract (development track).
//
// Differences vs the v1 driver above:
//   * Block 2 priors arrive as glmbayesCore pfamily objects (dNormal or
//     dIndependent_Normal_Gamma); the type string in pf$pfamily selects the
//     Block 2 update (native analogue of rglmb()'s simfun dispatch).
//   * Block 2 priors are parsed once before the loop (no RNG is consumed in
//     the prep, so the dNormal draw stream is identical to v1).
//   * A tau2 working vector tracks the Block 2 dispersion per RE component;
//     for ING components it will be updated by the joint (gamma_k, tau2_k)
//     draw and fed back into the Block 1 prior precision (next milestone).
//     The values at each stored draw are returned as dispersion_fixef_draws.
// ===========================================================================

namespace {

struct Block2PriorV2 {
  NumericVector mu;
  NumericMatrix P;
  double dispersion;       // dNormal: fixed tau2_k; ING: initial tau2_k
  bool is_ing;
  double shape;
  double rate;
  double max_disp_perc;
  SEXP disp_lower;         // NumericVector or R_NilValue (rglmb semantics)
  SEXP disp_upper;
};

// Shared mu / Sigma-or-P parsing (same checks and LAPACK path as
// block2_prior_prep above; duplicated so the v1 code stays byte-identical
// until v2 is promoted).
void block2_mu_P_prep_v2(
    const List& pl, int j1 /*1-based*/, int p,
    NumericVector& mu_out, NumericMatrix& P_out
) {
  if (!has_non_null(pl, "mu")) {
    Rcpp::stop("pfamily_list[[%d]]$prior_list must contain 'mu'.", j1);
  }
  if (!has_non_null(pl, "Sigma") && !has_non_null(pl, "P")) {
    Rcpp::stop(
      "pfamily_list[[%d]]$prior_list must contain 'Sigma' or 'P'.", j1
    );
  }

  NumericVector mu = Rcpp::as<NumericVector>(pl["mu"]);
  if (mu.size() != p) {
    Rcpp::stop(
      "pfamily_list[[%d]]$prior_list$mu must have length ncol(x) = %d.",
      j1, p
    );
  }

  NumericMatrix P;
  if (has_non_null(pl, "Sigma")) {
    NumericMatrix S = Rcpp::as<NumericMatrix>(pl["Sigma"]);
    if (S.nrow() != p || S.ncol() != p) {
      Rcpp::stop(
        "pfamily_list[[%d]]$prior_list$Sigma must be %d x %d.", j1, p, p
      );
    }
    if (!is_symmetric_mat(S)) {
      Rcpp::stop(
        "pfamily_list[[%d]]$prior_list$Sigma must be symmetric.", j1
      );
    }
    check_pd(S, "Sigma");
    if (!has_non_null(pl, "P")) {
      arma::mat S2(const_cast<double*>(S.begin()), p, p, false);
      arma::mat Pinv = arma::inv_sympd(S2);
      P = NumericMatrix(Rcpp::wrap(0.5 * (Pinv + Pinv.t())));
    }
  }
  if (has_non_null(pl, "P")) {
    P = Rcpp::as<NumericMatrix>(pl["P"]);
    if (P.nrow() != p || P.ncol() != p) {
      Rcpp::stop(
        "pfamily_list[[%d]]$prior_list$P must be %d x %d.", j1, p, p
      );
    }
  }
  if (!is_symmetric_mat(P)) {
    Rcpp::stop("matrix P must be symmetric");
  }
  check_pd(P, "P");

  mu_out = mu;
  P_out = P;
}

Block2PriorV2 block2_prior_prep_v2(const List& pf, int j1 /*1-based*/, int p) {
  if (!has_non_null(pf, "pfamily") || !has_non_null(pf, "prior_list")) {
    Rcpp::stop(
      "pfamily_list[[%d]] must be a pfamily object (e.g. dNormal() or "
      "dIndependent_Normal_Gamma()) with 'pfamily' and 'prior_list' fields.",
      j1
    );
  }
  const std::string ptype = Rcpp::as<std::string>(pf["pfamily"]);
  List pl(pf["prior_list"]);

  Block2PriorV2 out;
  out.is_ing = false;
  out.shape = NA_REAL;
  out.rate = NA_REAL;
  out.max_disp_perc = 0.99;
  out.disp_lower = R_NilValue;
  out.disp_upper = R_NilValue;

  block2_mu_P_prep_v2(pl, j1, p, out.mu, out.P);

  if (ptype == "dNormal") {
    // Same dispersion requirement as the v1 gaussian Block 2 path.
    bool ddef;
    if (pl.containsElementNamed("ddef")) {
      SEXP dd = pl["ddef"];
      ddef = Rf_isLogical(dd) && Rf_length(dd) >= 1 && LOGICAL(dd)[0] == TRUE;
    } else {
      ddef = !has_non_null(pl, "dispersion");
    }
    if (ddef || !has_non_null(pl, "dispersion")) {
      Rcpp::stop(
        "For gaussian() Block 2, dNormal() requires an explicit dispersion "
        "(pfamily_list[[%d]]).", j1
      );
    }
    out.dispersion = Rcpp::as<NumericVector>(pl["dispersion"])[0];
  } else if (ptype == "dIndependent_Normal_Gamma") {
    out.is_ing = true;
    if (!has_non_null(pl, "shape") || !has_non_null(pl, "rate")) {
      Rcpp::stop(
        "pfamily_list[[%d]]: dIndependent_Normal_Gamma requires 'shape' "
        "and 'rate'.", j1
      );
    }
    out.shape = Rcpp::as<NumericVector>(pl["shape"])[0];
    out.rate = Rcpp::as<NumericVector>(pl["rate"])[0];
    if (!(out.shape > 0.0) || !(out.rate > 0.0)) {
      Rcpp::stop(
        "pfamily_list[[%d]]: 'shape' and 'rate' must be positive.", j1
      );
    }
    if (has_non_null(pl, "max_disp_perc")) {
      out.max_disp_perc = Rcpp::as<NumericVector>(pl["max_disp_perc"])[0];
    }
    if (has_non_null(pl, "disp_lower")) {
      out.disp_lower = pl["disp_lower"];
    }
    if (has_non_null(pl, "disp_upper")) {
      out.disp_upper = pl["disp_upper"];
    }
    // Initial tau2_k for the first Block 1 sweep of each replicate chain:
    // the conservative disp_lower plug-in (consistent with the lmebayes
    // calibration); refined by the joint Block 2 draw within the sweep.
    out.dispersion = has_non_null(pl, "disp_lower")
      ? Rcpp::as<NumericVector>(pl["disp_lower"])[0]
      : NA_REAL;
  } else {
    Rcpp::stop(
      "pfamily_list[[%d]]: unsupported pfamily '%s' (allowed: dNormal, "
      "dIndependent_Normal_Gamma).", j1, ptype.c_str()
    );
  }

  return out;
}

} // anonymous namespace

List two_block_rNormal_reg_v2_cpp_export(
    int n,
    int m_convergence,
    const NumericVector& y,
    const NumericMatrix& x,
    SEXP block,
    const List& x_hyper,
    const List& prior_list_block1,
    SEXP dispersion_block1,
    SEXP ddef_block1,
    const List& pfamily_list,
    const List& fixef_start,
    const CharacterVector& group_levels,
    const std::string& family,
    const std::string& link,
    const Function& f2,
    const Function& f3,
    const Function& f2_gauss,
    const Function& f3_gauss,
    const NumericVector& offset,
    const NumericVector& wt,
    int Gridtype,
    int n_envopt,
    bool use_parallel,
    bool use_opencl,
    bool verbose,
    bool progbar
) {
  if (n < 1) {
    Rcpp::stop("'n' must be at least 1.");
  }
  if (m_convergence < 1) {
    Rcpp::stop("'m_convergence' must be at least 1.");
  }
  const int p_re = x.ncol();
  const int J = group_levels.size();
  if (x_hyper.size() != p_re) {
    Rcpp::stop("length(x_hyper) must equal ncol(x) = %d.", p_re);
  }
  if (pfamily_list.size() != p_re) {
    Rcpp::stop("length(pfamily_list) must equal ncol(x) = %d.", p_re);
  }
  if (fixef_start.size() != p_re) {
    Rcpp::stop("length(fixef_start) must equal ncol(x) = %d.", p_re);
  }

  const bool is_gaussian = (family == "gaussian");

  MuAllBuilder mu_builder(x_hyper, group_levels);

  // Block 2 priors parsed once (constant across draws; no RNG consumed).
  std::vector<Block2PriorV2> pr2(p_re);
  bool any_ing = false;
  for (int j = 0; j < p_re; ++j) {
    const NumericMatrix& X_j = mu_builder.X[j];
    pr2[j] = block2_prior_prep_v2(List(pfamily_list[j]), j + 1, X_j.ncol());
    if (pr2[j].is_ing) any_ing = true;
  }

  // fixef state (component order fixed by the R wrapper)
  std::vector<NumericVector> fixef_start_v(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_start_v[j] = Rcpp::as<NumericVector>(fixef_start[j]);
  }
  std::vector<NumericVector> fixef = fixef_start_v;

  // tau2 working vector: dNormal components stay fixed; ING components are
  // re-drawn jointly with gamma_k each sweep (and reset per stored draw).
  std::vector<double> tau2_start(p_re);
  for (int j = 0; j < p_re; ++j) {
    tau2_start[j] = pr2[j].dispersion;
  }
  std::vector<double> tau2 = tau2_start;

  // Base Block 1 prior precision (p_re x p_re), needed only when ING
  // components feed the current tau2 back into Block 1 each sweep.
  // ING rows/cols are overridden with diag(1/tau2_j); dNormal rows keep
  // their fixed entries from prior_list_block1.
  NumericMatrix base_P1;
  if (any_ing) {
    if (has_non_null(prior_list_block1, "P")) {
      base_P1 = Rcpp::as<NumericMatrix>(prior_list_block1["P"]);
    } else if (has_non_null(prior_list_block1, "Sigma")) {
      NumericMatrix S1 = Rcpp::as<NumericMatrix>(prior_list_block1["Sigma"]);
      arma::mat S1a(const_cast<double*>(S1.begin()),
                    S1.nrow(), S1.ncol(), false);
      arma::mat P1inv = arma::inv_sympd(S1a);
      base_P1 = NumericMatrix(Rcpp::wrap(0.5 * (P1inv + P1inv.t())));
    } else {
      Rcpp::stop("prior_list_block1 must contain 'P' or 'Sigma'.");
    }
    if (base_P1.nrow() != p_re || base_P1.ncol() != p_re) {
      Rcpp::stop(
        "prior_list_block1 P/Sigma must be %d x %d.", p_re, p_re
      );
    }
  }

  // Storage: fixef draws per component (n x q_k); b draws as J x p_re x n;
  // tau2 values at each stored draw (n x p_re); total Block 2 candidates
  // generated per stored draw (n x p_re; summed over the m_convergence
  // inner sweeps -- ING components count envelope candidates until
  // acceptance, dNormal components count 1 per conjugate draw).
  List fixef_draws(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_draws[j] = NumericMatrix(n, fixef_start_v[j].size());
  }
  NumericVector b_arr(Rcpp::Dimension(J, p_re, n));
  NumericMatrix disp_draws(n, p_re);
  NumericMatrix iters_draws(n, p_re);

  NumericMatrix mu_all;
  NumericMatrix b_i;
  CharacterVector group_ids;
  bool have_ids = false;

  for (int i = 0; i < n; ++i) {
    Rcpp::checkUserInterrupt();
    if (progbar) {
      glmbayes::progress::progress_bar(
        static_cast<double>(i + 1), static_cast<double>(n)
      );
    }

    fixef = fixef_start_v;
    tau2 = tau2_start;

    for (int m = 0; m < m_convergence; ++m) {

      mu_all = mu_builder.build(fixef);
      List pl1;
      if (any_ing) {
        // Refresh the Block 1 prior precision from the current tau2: ING
        // rows/cols become diag(1/tau2_j); dNormal rows keep base entries.
        NumericMatrix P1 = Rcpp::clone(base_P1);
        for (int j = 0; j < p_re; ++j) {
          if (!pr2[j].is_ing) continue;
          for (int c = 0; c < p_re; ++c) {
            P1(j, c) = 0.0;
            P1(c, j) = 0.0;
          }
          P1(j, j) = 1.0 / tau2[j];
        }
        List pl1_base = List::create(
          Rcpp::Named("P") = P1
        );
        pl1 = block1_prior_list(
          mu_all, pl1_base, dispersion_block1, ddef_block1
        );
      } else {
        pl1 = block1_prior_list(
          mu_all, prior_list_block1, dispersion_block1, ddef_block1
        );
      }

      List block_i;
      if (is_gaussian) {
        block_i = block_rNormalReg_cpp_export(
          1, y, x, block, pl1, R_NilValue, offset, wt, f2, f3, 2
        );
      } else {
        block_i = block_rNormalGLM_cpp_export(
          1, y, x, block, pl1, R_NilValue, offset, wt, f2, f3,
          family, link, Gridtype, n_envopt,
          use_parallel, use_opencl, verbose
        );
      }

      b_i = Rcpp::as<NumericMatrix>(block_i["coefficients"]);
      if (b_i.nrow() != J || b_i.ncol() != p_re) {
        Rcpp::stop(
          "Block 1 returned a %d x %d coefficient matrix; expected %d x %d.",
          b_i.nrow(), b_i.ncol(), J, p_re
        );
      }
      if (!have_ids) {
        List bi_info = block_i["block_info"];
        group_ids = Rcpp::as<CharacterVector>(bi_info["ids"]);
        have_ids = true;
      }

      // Block 2: per-component dispatch on the pfamily type.
      for (int j = 0; j < p_re; ++j) {
        const NumericMatrix& X_j = mu_builder.X[j];
        if (X_j.nrow() != b_i.nrow()) {
          Rcpp::stop(
            "nrow(x_hyper[[%d]]) (%d) must equal nrow(y) (%d).",
            j + 1, X_j.nrow(), b_i.nrow()
          );
        }
        NumericVector y_j = b_i(Rcpp::_, j);
        NumericVector offset_j(X_j.nrow(), 0.0);
        NumericVector wt_j(X_j.nrow(), 1.0);

        if (!pr2[j].is_ing) {
          // dNormal: conjugate draw at fixed dispersion (identical to v1).
          List out_j = rNormalReg(
            1, y_j, X_j, pr2[j].mu, pr2[j].P, offset_j, wt_j,
            pr2[j].dispersion, f2_gauss, f3_gauss, pr2[j].mu,
            "gaussian", "identity", 2
          );
          NumericMatrix coef_j =
            Rcpp::as<NumericMatrix>(out_j["coefficients"]);
          fixef[j] = NumericVector(coef_j(0, Rcpp::_));
          iters_draws(i, j) += 1.0;
        } else {
          // ING: joint (gamma_k, tau2_k) draw via the likelihood-subgradient
          // envelope sampler (same path as rglmb with a
          // dIndependent_Normal_Gamma pfamily; rIndepNormalGammaReg, n = 1).
          List out_j = rIndepNormalGammaReg(
            1, y_j, X_j, pr2[j].mu, pr2[j].P, offset_j, wt_j,
            pr2[j].shape, pr2[j].rate, pr2[j].max_disp_perc,
            Rcpp::Nullable<NumericVector>(pr2[j].disp_lower),
            Rcpp::Nullable<NumericVector>(pr2[j].disp_upper),
            2 /*Gridtype*/, 1 /*n_envopt*/,
            false /*use_parallel: n = 1 is serial anyway*/,
            false /*use_opencl*/, false /*verbose*/, false /*progbar*/
          );
          // "out" is p x n (coefficients in columns, original scale).
          NumericMatrix coef_j = Rcpp::as<NumericMatrix>(out_j["out"]);
          NumericVector gamma_j(coef_j.nrow());
          for (int c = 0; c < coef_j.nrow(); ++c) {
            gamma_j[c] = coef_j(c, 0);
          }
          fixef[j] = gamma_j;
          NumericVector disp_j = Rcpp::as<NumericVector>(out_j["disp_out"]);
          tau2[j] = disp_j[0];
          // Candidates generated by the envelope accept-reject sampler for
          // this accepted joint draw (>= 1; reading it consumes no RNG).
          NumericVector iters_j =
            Rcpp::as<NumericVector>(out_j["iters_out"]);
          iters_draws(i, j) += iters_j[0];
        }
      }
    }

    for (int j = 0; j < p_re; ++j) {
      NumericMatrix fd = fixef_draws[j];
      const NumericVector& fj = fixef[j];
      if (fj.size() != fd.ncol()) {
        Rcpp::stop(
          "Block 2 draw for component %d has length %d; expected %d.",
          j + 1, fj.size(), fd.ncol()
        );
      }
      fd(i, Rcpp::_) = fj;
    }
    for (int j = 0; j < p_re; ++j) {
      for (int g = 0; g < J; ++g) {
        b_arr[g + J * (j + p_re * i)] = b_i(g, j);
      }
      disp_draws(i, j) = tau2[j];
    }
  }

  List fixef_last(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_last[j] = fixef[j];
  }

  return List::create(
    Rcpp::Named("fixef_draws") = fixef_draws,
    Rcpp::Named("b_draws") = b_arr,
    Rcpp::Named("fixef_last") = fixef_last,
    Rcpp::Named("b_last") = b_i,
    Rcpp::Named("mu_all_last") = mu_all,
    Rcpp::Named("group_ids") = group_ids,
    Rcpp::Named("dispersion_fixef_draws") = disp_draws,
    Rcpp::Named("iters_fixef_draws") = iters_draws,
    Rcpp::Named("any_ing") = any_ing
  );
}

} // namespace sim
} // namespace glmbayes
