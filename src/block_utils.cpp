// block_utils.cpp
// C++ block partition + prior normalization for block_rNormalReg.

#include "RcppArmadillo.h"
#include "simfuncs.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace glmbayes {
namespace sim {

namespace {

using Rcpp::CharacterVector;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::RObject;

bool approx_equal_vec(const NumericVector& a, const NumericVector& b, double tol = 1e-8) {
  if (a.size() != b.size()) return false;
  for (int i = 0; i < a.size(); ++i) {
    if (std::fabs(a[i] - b[i]) > tol * (1.0 + std::fabs(a[i]) + std::fabs(b[i]))) {
      return false;
    }
  }
  return true;
}

bool approx_equal_mat(const NumericMatrix& a, const NumericMatrix& b, double tol = 1e-8) {
  if (a.nrow() != b.nrow() || a.ncol() != b.ncol()) return false;
  for (int i = 0; i < a.nrow(); ++i) {
    for (int j = 0; j < a.ncol(); ++j) {
      const double da = a(i, j);
      const double db = b(i, j);
      if (std::fabs(da - db) > tol * (1.0 + std::fabs(da) + std::fabs(db))) {
        return false;
      }
    }
  }
  return true;
}

void check_P_pd(const NumericMatrix& P, const char* label) {
  arma::mat P2(const_cast<double*>(P.begin()), P.nrow(), P.ncol(), false);
  arma::vec ev = arma::eig_sym(P2);
  const double tol = 1e-6;
  if (ev.min() < -tol * std::fabs(ev.max())) {
    Rcpp::stop("'%s' is not positive definite.", label);
  }
}

bool is_symmetric_mat(const NumericMatrix& M, double tol = 1e-8) {
  if (M.nrow() != M.ncol()) return false;
  for (int i = 0; i < M.nrow(); ++i) {
    for (int j = i + 1; j < M.ncol(); ++j) {
      if (std::fabs(M(i, j) - M(j, i)) > tol) return false;
    }
  }
  return true;
}

List prior_list_to_P_Sigma(List pl) {
  if (!pl.containsElementNamed("mu")) {
    Rcpp::stop("prior_list must contain 'mu'.");
  }
  NumericVector mu = pl["mu"];
  List out = List::create(Rcpp::Named("mu") = mu);

  if (pl.containsElementNamed("P")) {
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
  if (pl.containsElementNamed("Sigma")) {
    NumericMatrix Sigma = pl["Sigma"];
    if (!is_symmetric_mat(Sigma)) {
      Rcpp::stop("prior covariance Sigma must be symmetric.");
    }
    arma::mat S(const_cast<double*>(Sigma.begin()), Sigma.nrow(), Sigma.ncol(), false);
    arma::mat R = arma::chol(S);
    arma::mat Pinv = arma::inv(R);
    arma::mat P = Pinv.t() * Pinv;
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
  List out = List::create(
    Rcpp::Named("mu") = mu,
    Rcpp::Named("Sigma") = ps["Sigma"],
    Rcpp::Named("P") = P
  );
  if (pl.containsElementNamed("dispersion")) {
    out["dispersion"] = pl["dispersion"];
  }
  if (pl.containsElementNamed("ddef")) {
    out["ddef"] = pl["ddef"];
  }
  return out;
}

List rep_list_blocks(const List& one, int k) {
  List out(k);
  for (int j = 0; j < k; ++j) {
    out[j] = Rcpp::clone(one);
  }
  return out;
}

List normalize_prior_for_blocks_cpp(
    SEXP prior_list_sexp,
    SEXP prior_lists_sexp,
    const List& block_info,
    int l1
) {
  const int k = block_info["k"];
  CharacterVector ids = block_info["ids"];

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

  if (prior_list.containsElementNamed("blocks")) {
    List bl = prior_list["blocks"];
    if (bl.size() != static_cast<R_xlen_t>(k)) {
      Rcpp::stop("'prior_list$blocks' must have length k = %d.", k);
    }
    CharacterVector bl_names = bl.names();
    if (bl_names.size() == bl.size() && !ids.isNULL()) {
      bool all_named = true;
      for (int j = 0; j < bl_names.size(); ++j) {
        if (CharacterVector::is_na(bl_names[j]) || as<std::string>(bl_names[j]).empty()) {
          all_named = false;
          break;
        }
      }
      if (all_named) {
        bool all_match = true;
        for (int j = 0; j < k; ++j) {
          if (!CharacterVector::is_na(ids[j]) &&
              bl.containsElementNamed(as<std::string>(ids[j]).c_str())) {
            continue;
          }
          all_match = false;
          break;
        }
        if (all_match) {
          List out(k);
          for (int j = 0; j < k; ++j) {
            out[j] = base_prior_block(List(bl[as<std::string>(ids[j])]), l1);
          }
          return out;
        }
      }
    }
    List out(k);
    for (int j = 0; j < k; ++j) {
      out[j] = base_prior_block(List(bl[j]), l1);
    }
    return out;
  }

  SEXP mu_sexp = prior_list["mu"];
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

    if (prior_list.containsElementNamed("P")) {
      SEXP P_sexp = prior_list["P"];
      if (Rcpp::is<List>(P_sexp)) {
        P_list = List(P_sexp);
        have_P_list = P_list.size() == 1 || P_list.size() == static_cast<R_xlen_t>(k);
      } else if (Rcpp::is<NumericMatrix>(P_sexp)) {
        P_list = List::create(P_sexp);
        have_P_list = true;
      }
    }
    if (prior_list.containsElementNamed("Sigma")) {
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
      if (prior_list.containsElementNamed("dispersion")) {
        NumericVector disp = prior_list["dispersion"];
        if (disp.size() == k) pl_j["dispersion"] = disp[j];
        else pl_j["dispersion"] = disp[0];
      }
      if (prior_list.containsElementNamed("ddef")) {
        pl_j["ddef"] = prior_list["ddef"];
      }
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

List prior_payload_from_blocks(const List& prior_block, int l1, int k) {
  if (prior_block.size() != static_cast<R_xlen_t>(k)) {
    Rcpp::stop("prior_block must have length k = %d.", k);
  }

  List pb1 = prior_block[0];
  NumericVector disp_v(k);
  for (int j = 0; j < k; ++j) {
    List pb = prior_block[j];
    if (pb.containsElementNamed("dispersion")) {
      NumericVector d = pb["dispersion"];
      disp_v[j] = d[0];
    } else {
      disp_v[j] = 1.0;
    }
  }

  bool prior_by_block = false;
  for (int j = 1; j < k; ++j) {
    List pb = prior_block[j];
    List pb0 = pb1;
    if (!approx_equal_vec(NumericVector(pb["mu"]), NumericVector(pb0["mu"])) ||
        !approx_equal_mat(NumericMatrix(pb["P"]), NumericMatrix(pb0["P"])) ||
        std::fabs(disp_v[j] - disp_v[0]) > 1e-12) {
      prior_by_block = true;
      break;
    }
  }

  if (!prior_by_block) {
    NumericVector mu1 = pb1["mu"];
    NumericMatrix mu_mat(l1, 1);
    for (int i = 0; i < l1; ++i) mu_mat(i, 0) = mu1[i];
    return List::create(
      Rcpp::Named("mu") = mu_mat,
      Rcpp::Named("P_blocks") = List::create(pb1["P"]),
      Rcpp::Named("dispersion") = NumericVector::create(disp_v[0]),
      Rcpp::Named("prior_by_block") = false
    );
  }

  NumericMatrix mu_mat(l1, k);
  List P_blocks(k);
  for (int j = 0; j < k; ++j) {
    List pb = prior_block[j];
    NumericVector muj = pb["mu"];
    for (int i = 0; i < l1; ++i) mu_mat(i, j) = muj[i];
    P_blocks[j] = pb["P"];
  }

  return List::create(
    Rcpp::Named("mu") = mu_mat,
    Rcpp::Named("P_blocks") = P_blocks,
    Rcpp::Named("dispersion") = disp_v,
    Rcpp::Named("prior_by_block") = true
  );
}

List split_factor_rows(const IntegerVector& block, int l2) {
  CharacterVector lev = block.attr("levels");
  const int k = lev.size();
  List rows(k);
  IntegerVector counts(k, 0);
  for (int i = 0; i < l2; ++i) {
    counts[block[i] - 1]++;
  }
  for (int j = 0; j < k; ++j) {
    rows[j] = IntegerVector(counts[j]);
  }
  IntegerVector pos(k, 0);
  for (int i = 0; i < l2; ++i) {
    const int j = block[i] - 1;
    IntegerVector r = rows[j];
    r[pos[j]++] = i + 1;
    rows[j] = r;
  }
  IntegerVector l2_blocks(k);
  IntegerVector starts(k);
  int cum = 0;
  for (int j = 0; j < k; ++j) {
    l2_blocks[j] = counts[j];
    starts[j] = cum + 1;
    cum += counts[j];
  }
  return List::create(
    Rcpp::Named("k") = k,
    Rcpp::Named("ids") = lev,
    Rcpp::Named("l2_blocks") = l2_blocks,
    Rcpp::Named("starts") = starts,
    Rcpp::Named("rows") = rows
  );
}

} // anonymous namespace

List normalize_block_cpp(SEXP block_sexp, int l2) {
  if (l2 < 1) {
    Rcpp::stop("'l2' must be a positive integer (length of y).");
  }

  if (Rf_isNewList(block_sexp)) {
    List block_list(block_sexp);
    if (block_list.size() < 1) {
      Rcpp::stop("'block' list must have at least one element.");
    }
    const int k = block_list.size();
    List rows(k);
    std::vector<int> all_idx;
    all_idx.reserve(l2);
    for (int j = 0; j < k; ++j) {
      IntegerVector idx = block_list[j];
      std::vector<int> seen;
      IntegerVector rv;
      for (int i = 0; i < idx.size(); ++i) {
        const int r = idx[i];
        if (Rcpp::IntegerVector::is_na(r) || r < 1 || r > l2) {
          Rcpp::stop("Row indices in 'block' must be integers in 1:l2.");
        }
        if (std::find(seen.begin(), seen.end(), r) != seen.end()) continue;
        seen.push_back(r);
        rv.push_back(r);
        all_idx.push_back(r);
      }
      rows[j] = rv;
    }
    std::sort(all_idx.begin(), all_idx.end());
    if (static_cast<int>(all_idx.size()) != l2) {
      Rcpp::stop("Row indices in 'block' list must be disjoint and cover 1:l2.");
    }
    for (int i = 0; i < l2; ++i) {
      if (all_idx[i] != i + 1) {
        Rcpp::stop("Row indices in 'block' list must cover exactly 1:l2.");
      }
    }
    CharacterVector ids = block_list.names();
    bool replace_ids = ids.size() != block_list.size();
    if (!replace_ids) {
      for (int j = 0; j < ids.size(); ++j) {
        if (CharacterVector::is_na(ids[j]) || as<std::string>(ids[j]).empty()) {
          replace_ids = true;
          break;
        }
      }
    }
    if (replace_ids) {
      ids = CharacterVector(k);
      for (int j = 0; j < k; ++j) {
        ids[j] = ("block" + std::to_string(j + 1)).c_str();
      }
    }
    IntegerVector l2_blocks(k);
    IntegerVector starts(k);
    int cum = 0;
    for (int j = 0; j < k; ++j) {
      IntegerVector r = rows[j];
      l2_blocks[j] = r.size();
      starts[j] = cum + 1;
      cum += r.size();
    }
    return List::create(
      Rcpp::Named("k") = k,
      Rcpp::Named("ids") = ids,
      Rcpp::Named("l2_blocks") = l2_blocks,
      Rcpp::Named("starts") = starts,
      Rcpp::Named("rows") = rows
    );
  }

  if (Rf_isFactor(block_sexp)) {
    IntegerVector block(block_sexp);
    if (block.size() != l2) {
      Rcpp::stop("'block' factor length must equal l2.");
    }
    return split_factor_rows(block, l2);
  }

  IntegerVector block;
  if (Rcpp::is<IntegerVector>(block_sexp)) {
    block = IntegerVector(block_sexp);
  } else if (Rcpp::is<NumericVector>(block_sexp)) {
    block = Rcpp::as<IntegerVector>(block_sexp);
  } else {
    Rcpp::stop("'block' must be a factor, integer vector, or list of row indices.");
  }

  if (block.size() == l2) {
    Rcpp::Environment base_env = Rcpp::Environment::base_env();
    Rcpp::Function factor_fun = base_env["factor"];
    SEXP fblock = factor_fun(block);
    return split_factor_rows(IntegerVector(fblock), l2);
  }

  if (block.size() >= 1 && block.size() < l2) {
    int sum = 0;
    for (int i = 0; i < block.size(); ++i) {
      if (block[i] < 1) Rcpp::stop("l2_blocks counts must be positive.");
      sum += block[i];
    }
    if (sum != l2) {
      Rcpp::stop("sum(block) must equal length(y) (%d).", l2);
    }
    const int k = block.size();
    List rows(k);
    IntegerVector l2_blocks(k);
    IntegerVector starts(k);
    int start = 1;
    CharacterVector ids(k);
    for (int j = 0; j < k; ++j) {
      const int len = block[j];
      IntegerVector r(len);
      for (int t = 0; t < len; ++t) r[t] = start + t;
      rows[j] = r;
      l2_blocks[j] = len;
      starts[j] = start;
      ids[j] = ("block" + std::to_string(j + 1)).c_str();
      start += len;
    }
    return List::create(
      Rcpp::Named("k") = k,
      Rcpp::Named("ids") = ids,
      Rcpp::Named("l2_blocks") = l2_blocks,
      Rcpp::Named("starts") = starts,
      Rcpp::Named("rows") = rows
    );
  }

  Rcpp::stop(
    "'block' must be a factor or integer vector of length l2, a list of row indices, ",
    "or an integer vector of l2_blocks counts."
  );
}

List block_rNormalReg_cpp_export(
    int n,
    const NumericVector& y,
    const NumericMatrix& x,
    SEXP block,
    SEXP prior_list,
    SEXP prior_lists,
    const NumericVector& offset,
    const NumericVector& wt,
    const Function& f2,
    const Function& f3,
    int Gridtype
) {
  if (n < 1) {
    Rcpp::stop("'n' must be at least 1.");
  }
  const int l2 = y.size();
  const int l1 = x.ncol();
  if (x.nrow() != l2) {
    Rcpp::stop("nrow(x) must equal length(y).");
  }

  NumericVector offset2 = offset;
  NumericVector wt2 = wt;
  if (offset2.size() == 1) offset2 = Rcpp::rep(offset2[0], l2);
  if (wt2.size() == 1) wt2 = Rcpp::rep(wt2[0], l2);
  if (offset2.size() != l2) {
    Rcpp::stop("length(offset) must be 1 or length(y).");
  }
  if (wt2.size() != l2) {
    Rcpp::stop("length(weights) must be 1 or length(y).");
  }

  List block_info = normalize_block_cpp(block, l2);
  const int k = block_info["k"];

  List prior_norm = normalize_prior_for_blocks_cpp(
    prior_list, prior_lists, block_info, l1
  );
  List prior_block = prior_norm;

  for (int j = 0; j < k; ++j) {
    List pb = prior_block[j];
    if (!pb.containsElementNamed("dispersion")) {
      Rcpp::stop(
        "prior_list must contain 'dispersion' for block_rNormalReg. Block %d has no dispersion.",
        j + 1
      );
    }
  }

  List prior_cpp = prior_payload_from_blocks(prior_block, l1, k);

  List row_blocks = block_info["rows"];
  NumericMatrix mu = prior_cpp["mu"];
  List P_blocks = prior_cpp["P_blocks"];
  NumericVector dispersion = prior_cpp["dispersion"];
  bool prior_by_block = prior_cpp["prior_by_block"];

  List cpp_out = rNormalRegBlocks(
    n, y, x, offset2, wt2, dispersion,
    mu, P_blocks, prior_by_block, row_blocks, f2, f3, Gridtype
  );

  return List::create(
    Rcpp::Named("coefficients") = cpp_out["coefficients"],
    Rcpp::Named("coef.mode") = cpp_out["coef.mode"],
    Rcpp::Named("dispersion") = cpp_out["dispersion"],
    Rcpp::Named("block_results") = cpp_out["block_results"],
    Rcpp::Named("block_info") = block_info,
    Rcpp::Named("prior_lists") = prior_block,
    Rcpp::Named("n") = n,
    Rcpp::Named("k") = k,
    Rcpp::Named("l1") = l1,
    Rcpp::Named("l2") = l2,
    Rcpp::Named("y") = y,
    Rcpp::Named("x") = x,
    Rcpp::Named("offset") = offset2,
    Rcpp::Named("prior.weights") = wt2
  );
}

} // namespace sim
} // namespace glmbayes
