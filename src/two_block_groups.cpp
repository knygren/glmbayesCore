// two_block_groups.cpp
// Block 1 batch helpers for rGLMM_sweep (R-oracle semantics).
// Migrated incrementally from the R driver in R/two_block_batch_gibbs.R.
//
// R .Call names for these helpers no longer match lmebayesCore.
// See inst/TWO_BLOCK_CPP_CALLS.md. Argument name block1_prior is unchanged.

#include "simfuncs.h"
#include "simfuncs_twoblock.h"
#include "simfuncs_groups.h"
#include "progress_utils.h"
#include "package_ns.h"

#include <string>
#include <vector>

namespace glmbayes {
namespace sim {
namespace twoblock {

namespace {

void set_matrix_dimnames(
    Rcpp::NumericMatrix& M,
    const Rcpp::CharacterVector& rownames,
    SEXP colnames_src
) {
  if (rownames.size() == 0 && Rf_isNull(colnames_src)) {
    return;
  }
  Rcpp::List dn(2);
  if (rownames.size() > 0) {
    dn[0] = rownames;
  } else {
    dn[0] = R_NilValue;
  }
  if (Rf_isNull(colnames_src)) {
    dn[1] = R_NilValue;
  } else {
    dn[1] = colnames_src;
  }
  M.attr("dimnames") = dn;
}

std::vector<Rcpp::NumericVector> fixef_vec_from_list(
    const Rcpp::List& fixef,
    const Rcpp::CharacterVector& re_names
) {
  const int p_re = re_names.size();
  std::vector<Rcpp::NumericVector> out(static_cast<size_t>(p_re));
  for (int i = 0; i < p_re; ++i) {
    const std::string nm = Rcpp::as<std::string>(re_names[i]);
    if (!fixef.containsElementNamed(nm.c_str())) {
      Rcpp::stop("fixef must contain element '%s'.", nm.c_str());
    }
    out[static_cast<size_t>(i)] =
      Rcpp::as<Rcpp::NumericVector>(fixef[nm.c_str()]);
  }
  return out;
}

/// Forward optional list element (mirror R \code{$} on absent name → \code{NULL}).
SEXP optional_list_elt(const Rcpp::List& pl, const char* name) {
  return pl.containsElementNamed(name) ? pl[name] : R_NilValue;
}

/// Mirror \code{rNormalGLM_reg_group()} / \code{rNormal_reg_group()} dimnames on
/// \code{coefficients} (\code{colnames(x)}, \code{rownames} from
/// \code{group_info$ids}).
void set_block_draw_coefficient_dimnames(
    Rcpp::NumericMatrix& coef_draw,
    const Rcpp::NumericMatrix& x,
    const Rcpp::List& group_info
) {
  Rcpp::List dn(2);
  SEXP old_dn = coef_draw.attr("dimnames");
  if (!Rf_isNull(old_dn)) {
    Rcpp::List old_dnl(old_dn);
    if (old_dnl.size() >= 1 && !Rf_isNull(old_dnl[0])) {
      dn[0] = old_dnl[0];
    }
    if (old_dnl.size() >= 2 && !Rf_isNull(old_dnl[1])) {
      dn[1] = old_dnl[1];
    }
  }

  SEXP x_dn = x.attr("dimnames");
  if (!Rf_isNull(x_dn)) {
    Rcpp::List x_dnl(x_dn);
    if (x_dnl.size() >= 2 && !Rf_isNull(x_dnl[1])) {
      dn[1] = x_dnl[1];
    }
  }
  if (group_info.containsElementNamed("ids")) {
    SEXP ids = group_info["ids"];
    if (!Rf_isNull(ids)) {
      dn[0] = ids;
    }
  }
  if (!Rf_isNull(dn[0]) || !Rf_isNull(dn[1])) {
    coef_draw.attr("dimnames") = dn;
  }
}

/// Row ids for reorder (default \code{block_ids = rownames(b_draw)} in R).
SEXP matrix_rownames_sexp(const Rcpp::NumericMatrix& M) {
  SEXP dn = M.attr("dimnames");
  if (Rf_isNull(dn)) {
    return R_NilValue;
  }
  Rcpp::List dnl(dn);
  if (dnl.size() < 1 || Rf_isNull(dnl[0])) {
    return R_NilValue;
  }
  return dnl[0];
}

} // namespace

/// Port of \code{.two_block_batch_fixef_chain(batch, i)}.
// lmerb-unused: candidate for deletion (not on the lmerb path).
Rcpp::List fixef_list_from_batch_chain(
    const Rcpp::List& batch_fixef,
    int chain_i,
    const Rcpp::CharacterVector& re_names
) {
  if (chain_i < 1) {
    Rcpp::stop("'chain_i' must be at least 1.");
  }
  const int row_i = chain_i - 1;
  const int p_re = re_names.size();
  Rcpp::List out(p_re);
  out.attr("names") = re_names;
  for (int k = 0; k < p_re; ++k) {
    const std::string nm = Rcpp::as<std::string>(re_names[k]);
    if (!batch_fixef.containsElementNamed(nm.c_str())) {
      Rcpp::stop("batch_fixef must contain element '%s'.", nm.c_str());
    }
    Rcpp::NumericMatrix mat =
      Rcpp::as<Rcpp::NumericMatrix>(batch_fixef[nm.c_str()]);
    if (row_i >= mat.nrow()) {
      Rcpp::stop(
        "chain index %d exceeds nrow(batch_fixef[['%s']]) (%d).",
        chain_i, nm.c_str(), mat.nrow()
      );
    }
    Rcpp::NumericVector row = mat(row_i, Rcpp::_);
    SEXP col_dn = R_NilValue;
    SEXP mat_dn = mat.attr("dimnames");
    if (!Rf_isNull(mat_dn)) {
      Rcpp::List mat_dnl(mat_dn);
      if (mat_dnl.size() >= 2 && !Rf_isNull(mat_dnl[1])) {
        col_dn = mat_dnl[1];
      }
    }
    if (!Rf_isNull(col_dn)) {
      row.attr("names") = col_dn;
    }
    out[k] = row;
  }
  return out;
}

/// All-chains step A: mirror \code{batch$tau2[i, ]} (1-based \code{chain_i};
/// preserve matrix \code{colnames} as vector \code{names}).
// lmerb-unused: candidate for deletion (not on the lmerb path).
Rcpp::NumericVector batch_tau2_chain_row(
    const Rcpp::NumericMatrix& batch_tau2,
    int chain_i
) {
  if (chain_i < 1) {
    Rcpp::stop("'chain_i' must be at least 1.");
  }
  const int row_i = chain_i - 1;
  if (row_i >= batch_tau2.nrow()) {
    Rcpp::stop(
      "chain index %d exceeds nrow(batch_tau2) (%d).",
      chain_i, batch_tau2.nrow()
    );
  }
  Rcpp::NumericVector tau2_i = batch_tau2(row_i, Rcpp::_);
  SEXP dn = batch_tau2.attr("dimnames");
  if (!Rf_isNull(dn)) {
    Rcpp::List dnl(dn);
    if (dnl.size() >= 2 && !Rf_isNull(dnl[1])) {
      tau2_i.attr("names") = dnl[1];
    }
  }
  return tau2_i;
}

/// All-chains step C: mirror \code{batch$b[, , chain_i] <- b_draw}.
/// \code{b_store} must have \code{dim = c(J, p_re, n)} (R column-major).
/// Assignment is positional: \code{b_draw(g, j)} → \code{b[g, j, chain_i]}.
// lmerb-unused: candidate for deletion (not on the lmerb path).
void batch_b_assign_slice(
    Rcpp::NumericVector& b_store,
    int chain_i,
    const Rcpp::NumericMatrix& b_draw
) {
  if (chain_i < 1) {
    Rcpp::stop("'chain_i' must be at least 1.");
  }
  Rcpp::IntegerVector b_dim = b_store.attr("dim");
  if (b_dim.size() != 3) {
    Rcpp::stop("'b_store' must be a 3-dimensional array (J x p_re x n).");
  }
  const int J = b_dim[0];
  const int p_re = b_dim[1];
  const int n = b_dim[2];
  if (chain_i > n) {
    Rcpp::stop("'chain_i' exceeds third dimension of b_store (%d).", n);
  }
  const int chain0 = chain_i - 1;
  if (b_draw.nrow() != J || b_draw.ncol() != p_re) {
    Rcpp::stop(
      "Block 1 slice for chain %d is %d x %d; expected %d x %d.",
      chain_i, b_draw.nrow(), b_draw.ncol(), J, p_re
    );
  }
  for (int j = 0; j < p_re; ++j) {
    for (int g = 0; g < J; ++g) {
      b_store[g + J * (j + p_re * chain0)] = b_draw(g, j);
    }
  }
}

/// Third dimension of \code{batch$b} (\code{n_chains}); requires \code{dim}.
// lmerb-unused: candidate for deletion (not on the lmerb path).
int batch_b_n_chains(const Rcpp::NumericVector& batch_b) {
  Rcpp::IntegerVector b_dim = batch_b.attr("dim");
  if (b_dim.size() != 3) {
    Rcpp::stop("'batch_b' must be a 3-dimensional array (J x p_re x n).");
  }
  return static_cast<int>(b_dim[2]);
}

/// Restore \code{dim = c(J, p_re, n)} and \code{dimnames} on \code{batch$b}
/// (mirror R \code{.two_block_ensure_batch_b_dimnames} / batch init).
// lmerb-unused: candidate for deletion (not on the lmerb path).
void ensure_batch_b_dimnames(
    Rcpp::NumericVector& batch_b,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::CharacterVector& re_names,
    int n_chains
) {
  const int J = group_levels.size();
  const int p_re = re_names.size();
  if (J < 1 || p_re < 1 || n_chains < 1) {
    Rcpp::stop(
      "ensure_batch_b_dimnames: J, p_re, and n_chains must be positive."
    );
  }
  const int expected_len = J * p_re * n_chains;
  if (batch_b.size() != expected_len) {
    Rcpp::stop(
      "batch$b length (%d) != J * p_re * n (%d * %d * %d).",
      batch_b.size(), J, p_re, n_chains
    );
  }
  Rcpp::IntegerVector dim = {J, p_re, n_chains};
  batch_b.attr("dim") = dim;

  Rcpp::List dn(3);
  dn[0] = group_levels;
  dn[1] = re_names;
  dn[2] = R_NilValue;
  batch_b.attr("dimnames") = dn;
}

/// All-chains step D: mirror \code{batch$iters_ranef[chain_i] <-
/// batch$iters_ranef[chain_i] + iters_mean}.
// lmerb-unused: candidate for deletion (not on the lmerb path).
void batch_iters_ranef_add(
    Rcpp::NumericVector& iters_ranef,
    int chain_i,
    double iters_mean
) {
  if (chain_i < 1) {
    Rcpp::stop("'chain_i' must be at least 1.");
  }
  if (chain_i > iters_ranef.size()) {
    Rcpp::stop(
      "'chain_i' (%d) exceeds length(iters_ranef) (%d).",
      chain_i, iters_ranef.size()
    );
  }
  if (!R_finite(iters_mean)) {
    Rcpp::stop("'iters_mean' must be finite.");
  }
  iters_ranef[chain_i - 1] += iters_mean;
}

MuAllBuilder::MuAllBuilder(
    const Rcpp::List& x_hyper,
    const Rcpp::CharacterVector& group_levels
) {
  p_re = x_hyper.size();
  J = group_levels.size();
  X.reserve(p_re);
  row_idx.resize(static_cast<size_t>(p_re));
  for (int i = 0; i < p_re; ++i) {
    Rcpp::NumericMatrix X_k = Rcpp::as<Rcpp::NumericMatrix>(x_hyper[i]);
    X.push_back(X_k);

    Rcpp::CharacterVector rn;
    SEXP dn = X_k.attr("dimnames");
    if (!Rf_isNull(dn)) {
      Rcpp::List dnl(dn);
      if (dnl.size() >= 1 && !Rf_isNull(dnl[0])) {
        rn = Rcpp::CharacterVector(dnl[0]);
      }
    }

    std::vector<int>& idx = row_idx[static_cast<size_t>(i)];
    idx.resize(static_cast<size_t>(J));
    if (rn.size() == 0) {
      if (X_k.nrow() != J) {
        Rcpp::stop(
          "nrow(X_hyper[[%d]]) (%d) must equal length(group_levels) (%d).",
          i + 1, X_k.nrow(), J
        );
      }
      for (int j = 0; j < J; ++j) {
        idx[static_cast<size_t>(j)] = j;
      }
    } else {
      for (int j = 0; j < J; ++j) {
        const std::string lev = Rcpp::as<std::string>(group_levels[j]);
        int found = -1;
        for (int r = 0; r < rn.size(); ++r) {
          if (!Rcpp::CharacterVector::is_na(rn[r]) &&
              Rcpp::as<std::string>(rn[r]) == lev) {
            found = r;
            break;
          }
        }
        if (found < 0) {
          Rcpp::stop(
            "group level(s) not found in rownames(X_hyper[[%d]]): %s",
            i + 1, lev.c_str()
          );
        }
        idx[static_cast<size_t>(j)] = found;
      }
    }
  }
}

Rcpp::NumericMatrix MuAllBuilder::build(
    const std::vector<Rcpp::NumericVector>& fixef
) const {
  Rcpp::NumericMatrix mu_all(p_re, J);
  for (int i = 0; i < p_re; ++i) {
    const Rcpp::NumericMatrix& X_k = X[static_cast<size_t>(i)];
    const Rcpp::NumericVector& gamma_k = fixef[static_cast<size_t>(i)];
    if (gamma_k.size() != X_k.ncol()) {
      Rcpp::stop(
        "length(fixef[[%d]]) (%d) must equal ncol(X_hyper[[%d]]) (%d).",
        i + 1, gamma_k.size(), i + 1, X_k.ncol()
      );
    }
    const std::vector<int>& idx = row_idx[static_cast<size_t>(i)];
    for (int j = 0; j < J; ++j) {
      const int r = idx[static_cast<size_t>(j)];
      double s = 0.0;
      for (int c = 0; c < X_k.ncol(); ++c) {
        s += X_k(r, c) * gamma_k[c];
      }
      mu_all(i, j) = s;
    }
  }
  return mu_all;
}

Rcpp::NumericMatrix two_block_build_mu_all(
    const Rcpp::List& x_hyper,
    const Rcpp::List& fixef,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::CharacterVector& group_levels
) {
  const std::vector<Rcpp::NumericVector> fixef_v =
    fixef_vec_from_list(fixef, re_names);
  MuAllBuilder builder(x_hyper, group_levels);
  Rcpp::NumericMatrix mu_all = builder.build(fixef_v);
  set_matrix_dimnames(mu_all, re_names, group_levels);
  return mu_all;
}

/// Refresh Block~1 prior precision for ING components from chain \code{tau2}
/// (port of \code{.two_block_block1_prior_with_tau2}; v5 \code{any_ing} semantics).
// lmerb-unused: candidate for deletion (not on the lmerb path).
Rcpp::List two_block_groups_prior_with_tau2(
    const Rcpp::List& base_prior,
    const Rcpp::NumericVector& tau2_vec,
    const Rcpp::CharacterVector& ptypes,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::NumericMatrix& mu_all
) {
  const int p_re = re_names.size();
  if (tau2_vec.size() != p_re) {
    Rcpp::stop("length(tau2_vec) must equal length(re_names).");
  }
  if (ptypes.size() != p_re) {
    Rcpp::stop("length(ptypes) must equal length(re_names).");
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("mu") = mu_all,
    Rcpp::Named("dispersion") = optional_list_elt(base_prior, "dispersion"),
    Rcpp::Named("ddef") = optional_list_elt(base_prior, "ddef")
  );

  bool any_ing = false;
  for (int k = 0; k < p_re; ++k) {
    if (Rcpp::as<std::string>(ptypes[k]) == "dIndependent_Normal_Gamma") {
      any_ing = true;
      break;
    }
  }

  if (!any_ing) {
    out["P"] = base_prior["P"];
    return out;
  }

  if (!base_prior.containsElementNamed("P")) {
    Rcpp::stop("base_prior must contain 'P' when refreshing ING precision.");
  }

  Rcpp::NumericMatrix P1 =
    Rcpp::clone(Rcpp::as<Rcpp::NumericMatrix>(base_prior["P"]));
  if (P1.nrow() != p_re || P1.ncol() != p_re) {
    Rcpp::stop("dim(base_prior$P) must be p_re x p_re.");
  }

  for (int k = 0; k < p_re; ++k) {
    if (Rcpp::as<std::string>(ptypes[k]) != "dIndependent_Normal_Gamma") {
      continue;
    }
    for (int c = 0; c < p_re; ++c) {
      P1(k, c) = 0.0;
      P1(c, k) = 0.0;
    }
    P1(k, k) = 1.0 / tau2_vec[k];
  }
  out["P"] = P1;
  return out;
}

/// Mean envelope candidate count across groups from one Block 1 draw
/// (port of .two_block_block1_iters_mean in two_block_batch_gibbs.R).
// lmerb-unused: candidate for deletion (not on the lmerb path).
double two_block_groups_iters_mean(const Rcpp::List& block_out) {
  if (!block_out.containsElementNamed("group_results")) {
    return 1.0;
  }
  Rcpp::List br = Rcpp::as<Rcpp::List>(block_out["group_results"]);
  if (br.size() == 0) {
    return 1.0;
  }
  double sum = 0.0;
  for (int b = 0; b < br.size(); ++b) {
    Rcpp::List out_b = Rcpp::as<Rcpp::List>(br[b]);
    if (!out_b.containsElementNamed("iters") || Rf_isNull(out_b["iters"])) {
      sum += 1.0;
      continue;
    }
    SEXP it = out_b["iters"];
    if (Rf_isMatrix(it)) {
      Rcpp::NumericMatrix im(it);
      sum += im(0, 0);
    } else {
      Rcpp::NumericVector iv(it);
      sum += iv[0];
    }
  }
  return sum / static_cast<double>(br.size());
}

/// Reorder Block~1 coefficient rows to \code{group_levels}
/// (port of R \code{match(group_levels, rownames(b_draw))} in draw_one_chain).
/// When \code{block_ids} is \code{NULL}, returns \code{b_draw} unchanged.
// lmerb-unused: candidate for deletion (not on the lmerb path).
Rcpp::NumericMatrix two_block_reorder_b_to_group_levels(
    Rcpp::NumericMatrix b_draw,
    SEXP block_ids_sexp,
    const Rcpp::CharacterVector& group_levels
) {
  if (Rf_isNull(block_ids_sexp)) {
    return b_draw;
  }
  Rcpp::CharacterVector block_ids(block_ids_sexp);
  if (block_ids.size() == 0) {
    return b_draw;
  }

  const int J = group_levels.size();
  if (b_draw.nrow() != J || block_ids.size() != J) {
    Rcpp::stop("Block 1 group dimension mismatch during reorder.");
  }

  SEXP col_dn = R_NilValue;
  SEXP b_dn = b_draw.attr("dimnames");
  if (!Rf_isNull(b_dn)) {
    Rcpp::List b_dnl(b_dn);
    if (b_dnl.size() >= 2 && !Rf_isNull(b_dnl[1])) {
      col_dn = b_dnl[1];
    }
  }

  bool aligned = true;
  for (int g = 0; g < J; ++g) {
    if (Rcpp::CharacterVector::is_na(block_ids[g]) ||
        Rcpp::CharacterVector::is_na(group_levels[g]) ||
        Rcpp::as<std::string>(block_ids[g]) !=
          Rcpp::as<std::string>(group_levels[g])) {
      aligned = false;
      break;
    }
  }
  if (aligned) {
    Rcpp::NumericMatrix out = Rcpp::clone(b_draw);
    set_matrix_dimnames(out, group_levels, col_dn);
    return out;
  }

  Rcpp::NumericMatrix out(J, b_draw.ncol());
  for (int g = 0; g < J; ++g) {
    const std::string lev = Rcpp::as<std::string>(group_levels[g]);
    int src = -1;
    for (int r = 0; r < J; ++r) {
      if (!Rcpp::CharacterVector::is_na(block_ids[r]) &&
          Rcpp::as<std::string>(block_ids[r]) == lev) {
        src = r;
        break;
      }
    }
    if (src < 0) {
      Rcpp::stop("Block 1 group ids do not match group_levels.");
    }
    out(g, Rcpp::_) = b_draw(src, Rcpp::_);
  }
  set_matrix_dimnames(out, group_levels, col_dn);
  return out;
}

/// Block~1 prep + draw for one replicate chain.
/// Mirrors \code{.two_block_block1_prep_one_chain} then
/// \code{.two_block_block1_draw_one_chain} with all piecewise C++ flags
/// \code{TRUE} (\code{use_cpp_mu_all}, \code{use_cpp_prior_tau2},
/// \code{use_cpp_reorder}, \code{use_cpp_iters}).
// lmerb-unused: candidate for deletion (not on the lmerb path).
Rcpp::List two_block_groups_one_chain_impl(
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
  // ---- prep (.two_block_block1_prep_one_chain; piecewise C++ flags TRUE) ----
  Rcpp::List fixef_i = fixef_list_from_batch_chain(
    batch_fixef, chain_i, re_names
  );
  Rcpp::NumericMatrix mu_all = two_block_build_mu_all(
    x_hyper, fixef_i, re_names, group_levels
  );
  Rcpp::List prior_list = two_block_groups_prior_with_tau2(
    block1_prior, tau2_i, ptypes, re_names, mu_all
  );

  // ---- draw (.two_block_block1_draw_one_chain; piecewise C++ flags TRUE) ----
  Rcpp::List block_out;
  if (is_gaussian) {
    block_out = group::group_rNormalReg_cpp_export(
      1, y, Z, groups, prior_list, R_NilValue, offset, wt,
      f2_gauss, f3_gauss, Gridtype
    );
  } else {
    block_out = group::group_rNormalGLM_cpp_export(
      1, y, Z, groups, prior_list, R_NilValue, offset, wt,
      f2, f3, family, link, Gridtype, n_envopt,
      false, false, false
    );
  }

  Rcpp::List group_info = Rcpp::as<Rcpp::List>(block_out["group_info"]);
  Rcpp::NumericMatrix b_draw =
    Rcpp::as<Rcpp::NumericMatrix>(block_out["coefficients"]);
  const int J = group_levels.size();
  const int p_re = re_names.size();
  if (b_draw.nrow() != J || b_draw.ncol() != p_re) {
    Rcpp::stop(
      "Block 1 returned a %d x %d coefficient matrix; expected %d x %d.",
      b_draw.nrow(), b_draw.ncol(), J, p_re
    );
  }

  set_block_draw_coefficient_dimnames(b_draw, Z, group_info);
  b_draw = two_block_reorder_b_to_group_levels(
    b_draw, matrix_rownames_sexp(b_draw), group_levels
  );
  const double iters_mean = two_block_groups_iters_mean(block_out);

  return Rcpp::List::create(
    Rcpp::Named("mu_all") = mu_all,
    Rcpp::Named("prior_list") = prior_list,
    Rcpp::Named("b") = b_draw,
    Rcpp::Named("iters_mean") = iters_mean
  );
}

/// Mirror R \code{.two_block_block1_one_chain_cpp} offset/weights prep from \code{design}.
void design_offset_wt(
    const Rcpp::List& design,
    Rcpp::NumericVector& offset,
    Rcpp::NumericVector& wt
) {
  const Rcpp::NumericVector y = Rcpp::as<Rcpp::NumericVector>(design["y"]);
  const int l2 = y.size();
  if (l2 < 1) {
    Rcpp::stop("'design$y' must be non-empty.");
  }

  SEXP offset_sexp = R_NilValue;
  if (design.containsElementNamed("offset")) {
    offset_sexp = design["offset"];
  }
  if (Rf_isNull(offset_sexp)) {
    offset = Rcpp::NumericVector(l2);
    std::fill(offset.begin(), offset.end(), 0.0);
  } else {
    offset = Rcpp::as<Rcpp::NumericVector>(offset_sexp);
    if (offset.size() == 1) {
      const double v = offset[0];
      offset = Rcpp::NumericVector(l2);
      std::fill(offset.begin(), offset.end(), v);
    }
  }

  SEXP wt_sexp = R_NilValue;
  if (design.containsElementNamed("weights")) {
    wt_sexp = design["weights"];
  }
  if (Rf_isNull(wt_sexp)) {
    wt = Rcpp::NumericVector(l2);
    std::fill(wt.begin(), wt.end(), 1.0);
  } else {
    wt = Rcpp::as<Rcpp::NumericVector>(wt_sexp);
    if (wt.size() == 1) {
      const double v = wt[0];
      wt = Rcpp::NumericVector(l2);
      std::fill(wt.begin(), wt.end(), v);
    }
  }
}

/// Coerce \code{design$W} like \code{lapply(..., as.matrix)} in R.
Rcpp::List x_hyper_as_matrix(const Rcpp::List& x_hyper) {
  const int p_re = x_hyper.size();
  Rcpp::List out(p_re);
  SEXP nms = x_hyper.names();
  if (!Rf_isNull(nms)) {
    out.names() = nms;
  }
  for (int k = 0; k < p_re; ++k) {
    out[k] = Rcpp::as<Rcpp::NumericMatrix>(x_hyper[k]);
  }
  return out;
}

/// Mirror \code{identical(family$family, "gaussian")} in R.
bool family_is_gaussian(SEXP family) {
  if (Rf_isNull(family)) {
    return false;
  }
  Rcpp::List fam(family);
  if (!fam.containsElementNamed("family")) {
    return false;
  }
  SEXP fam_name = fam["family"];
  if (Rf_isNull(fam_name)) {
    return false;
  }
  return Rcpp::as<std::string>(fam_name) == "gaussian";
}

std::string family_name_string(SEXP family) {
  Rcpp::List fam(family);
  return Rcpp::as<std::string>(fam["family"]);
}

std::string family_link_string(SEXP family) {
  Rcpp::List fam(family);
  return Rcpp::as<std::string>(fam["link"]);
}

/// Block~1 draw for one chain from R-built \code{mu_all} and prior precision \code{P}.
/// Prep (\code{fixef_i} -> \code{mu_all}, \code{tau2_i} -> \code{P}) stays in R.
// lmerb-unused: candidate for deletion (not on the lmerb path).
Rcpp::List two_block_groups_one_chain_from_mu_P_impl(
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
  Rcpp::List prior_list = Rcpp::List::create(
    Rcpp::Named("mu") = mu_all,
    Rcpp::Named("P") = P,
    Rcpp::Named("dispersion") = dispersion,
    Rcpp::Named("ddef") = ddef
  );

  Rcpp::NumericVector offset;
  Rcpp::NumericVector wt;
  design_offset_wt(design, offset, wt);

  const bool is_gaussian = family_is_gaussian(family);
  const std::string fam_str = family_name_string(family);
  const std::string link_str = family_link_string(family);

  const Rcpp::NumericVector y =
    Rcpp::as<Rcpp::NumericVector>(design["y"]);
  const Rcpp::NumericMatrix Z =
    Rcpp::as<Rcpp::NumericMatrix>(design["D"]);
  SEXP groups = design["group"];

  Rcpp::List block_out;
  if (is_gaussian) {
    block_out = group::group_rNormalReg_cpp_export(
      1, y, Z, groups, prior_list, R_NilValue, offset, wt,
      f2_gauss, f3_gauss, 2
    );
  } else {
    block_out = group::group_rNormalGLM_cpp_export(
      1, y, Z, groups, prior_list, R_NilValue, offset, wt,
      f2, f3, fam_str, link_str, 2, 1,
      false, false, false
    );
  }

  Rcpp::List group_info = Rcpp::as<Rcpp::List>(block_out["group_info"]);
  Rcpp::NumericMatrix b_draw =
    Rcpp::as<Rcpp::NumericMatrix>(block_out["coefficients"]);
  const int J = group_levels.size();
  const int p_re = re_names.size();
  if (b_draw.nrow() != J || b_draw.ncol() != p_re) {
    Rcpp::stop(
      "Block 1 returned a %d x %d coefficient matrix; expected %d x %d.",
      b_draw.nrow(), b_draw.ncol(), J, p_re
    );
  }

  set_block_draw_coefficient_dimnames(b_draw, Z, group_info);

  SEXP block_ids = R_NilValue;
  if (group_info.containsElementNamed("ids")) {
    block_ids = group_info["ids"];
  }
  if (Rf_isNull(block_ids)) {
    block_ids = matrix_rownames_sexp(b_draw);
  }
  if (Rf_isNull(block_ids)) {
    Rcpp::stop(
      "Block 1 draw returned no group ids; cannot align rows to group_levels."
    );
  }
  b_draw = two_block_reorder_b_to_group_levels(
    b_draw, block_ids, group_levels
  );
  set_matrix_dimnames(b_draw, group_levels, re_names);

  const double iters_mean = two_block_groups_iters_mean(block_out);

  return Rcpp::List::create(
    Rcpp::Named("b") = b_draw,
    Rcpp::Named("iters_mean") = iters_mean
  );
}


/// Full per-chain Block~1 orchestration (steps A→D).
/// Port of exported R \code{.two_block_block1_one_chain_cpp} before the
/// rcpp_wrappers thin \code{.Call} layer.
// lmerb-unused: candidate for deletion (not on the lmerb path).
Rcpp::List two_block_groups_one_chain_orchestrate_impl(
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
  if (chain_i < 1) {
    Rcpp::stop("'chain_i' must be at least 1.");
  }

  Rcpp::NumericVector offset;
  Rcpp::NumericVector wt;
  design_offset_wt(design, offset, wt);

  const bool is_gaussian = family_is_gaussian(family);
  const std::string fam_str = family_name_string(family);
  const std::string link_str = family_link_string(family);

  Rcpp::NumericVector tau2_i;
  if (use_cpp_tau2_row) {
    tau2_i = batch_tau2_chain_row(batch_tau2, chain_i);
  } else {
    const int row_i = chain_i - 1;
    if (row_i >= batch_tau2.nrow()) {
      Rcpp::stop(
        "chain index %d exceeds nrow(batch_tau2) (%d).",
        chain_i, batch_tau2.nrow()
      );
    }
    tau2_i = batch_tau2(row_i, Rcpp::_);
  }

  const Rcpp::NumericVector y =
    Rcpp::as<Rcpp::NumericVector>(design["y"]);
  const Rcpp::NumericMatrix Z =
    Rcpp::as<Rcpp::NumericMatrix>(design["D"]);
  SEXP groups = design["group"];
  Rcpp::List x_hyper =
    x_hyper_as_matrix(Rcpp::as<Rcpp::List>(design["W"]));

  Rcpp::List draw_out = two_block_groups_one_chain_impl(
    chain_i, batch_fixef, tau2_i, y, Z, groups, offset, wt,
    x_hyper, re_names, group_levels, ptypes, block1_prior,
    is_gaussian, f2, f3, f2_gauss, f3_gauss,
    fam_str, link_str, 2, 1
  );

  Rcpp::NumericMatrix b_draw_group_ordered =
    Rcpp::as<Rcpp::NumericMatrix>(draw_out["b"]);
  const double iters_mean = Rcpp::as<double>(draw_out["iters_mean"]);

  Rcpp::NumericVector batch_b =
    use_cpp_b_slice ? Rcpp::clone(b_store) : b_store;
  batch_b_assign_slice(batch_b, chain_i, b_draw_group_ordered);
  ensure_batch_b_dimnames(
    batch_b, group_levels, re_names, batch_b_n_chains(batch_b)
  );

  Rcpp::NumericVector batch_iters_ranef =
    use_cpp_iters_ranef_add ? Rcpp::clone(iters_ranef) : iters_ranef;
  batch_iters_ranef_add(batch_iters_ranef, chain_i, iters_mean);

  return Rcpp::List::create(
    Rcpp::Named("b") = batch_b,
    Rcpp::Named("iters_ranef") = batch_iters_ranef
  );
}


}  // namespace twoblock
} // namespace sim
} // namespace glmbayes
