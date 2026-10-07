// two_block_population.cpp
// Sweep population update: one chain of gamma and tau2 given b.
// The Gibbs driver in twoBlockGibbs.cpp does not call these functions.
//
// .Call names: inst/TWO_BLOCK_CPP_CALLS.md.

#include "RcppArmadillo.h"
#include "simfuncs_twoblock.h"
#include "package_ns.h"

#include <algorithm>
#include <string>

namespace glmbayes {
namespace sim {
namespace twoblock {

namespace {

NumericVector two_block_coef_draw_from_rglmb(const Rcpp::List& fit_k) {
  SEXP coef_sexp = fit_k["coefficients"];
  if (Rf_isNull(coef_sexp)) {
    Rcpp::stop("rglmb result missing coefficients.");
  }
  if (!Rcpp::is<Rcpp::NumericMatrix>(coef_sexp)) {
    Rcpp::stop("rglmb coefficients must be a matrix.");
  }
  Rcpp::NumericMatrix coef(coef_sexp);
  if (coef.nrow() < 1) {
    Rcpp::stop("rglmb coefficients must have at least one row.");
  }
  Rcpp::NumericVector draw(coef.ncol());
  for (int j = 0; j < coef.ncol(); ++j) {
    draw[j] = coef(0, j);
  }
  if (coef.hasAttribute("dimnames")) {
    Rcpp::List dn(coef.attr("dimnames"));
    if (!Rf_isNull(dn[1])) {
      draw.attr("names") = dn[1];
    }
  }
  return draw;
}


bool charvec_contains_name(
    const Rcpp::CharacterVector& names,
    const std::string& target
) {
  for (int i = 0; i < names.size(); ++i) {
    if (Rcpp::CharacterVector::is_na(names[i])) {
      continue;
    }
    if (Rcpp::as<std::string>(names[i]) == target) {
      return true;
    }
  }
  return false;
}

double charvec_named_value(
    const Rcpp::NumericVector& values,
    const Rcpp::CharacterVector& names,
    const std::string& target
) {
  for (int i = 0; i < names.size(); ++i) {
    if (Rcpp::CharacterVector::is_na(names[i])) {
      continue;
    }
    if (Rcpp::as<std::string>(names[i]) == target) {
      return values[i];
    }
  }
  Rcpp::stop("Group level \"%s\" missing from b.", target.c_str());
  return NA_REAL;
}

Rcpp::CharacterVector matrix_rownames(const Rcpp::NumericMatrix& X_k) {
  if (!X_k.hasAttribute("dimnames")) {
    return Rcpp::CharacterVector();
  }
  Rcpp::List dn(X_k.attr("dimnames"));
  if (Rf_isNull(dn[0])) {
    return Rcpp::CharacterVector();
  }
  return Rcpp::as<Rcpp::CharacterVector>(dn[0]);
}

void two_block_update_fixef_row(
    Rcpp::NumericVector& row,
    const Rcpp::NumericVector& coef_k
) {
  if (coef_k.hasAttribute("names") && !Rf_isNull(coef_k.attr("names"))) {
    Rcpp::CharacterVector row_names = row.attr("names");
    Rcpp::CharacterVector coef_names = coef_k.attr("names");
    for (int i = 0; i < coef_names.size(); ++i) {
      const std::string nm = Rcpp::as<std::string>(coef_names[i]);
      for (int c = 0; c < row_names.size(); ++c) {
        if (Rcpp::CharacterVector::is_na(row_names[c])) {
          continue;
        }
        if (Rcpp::as<std::string>(row_names[c]) == nm) {
          row[c] = coef_k[i];
          break;
        }
      }
    }
    return;
  }
  const int ncopy = std::min(row.size(), coef_k.size());
  for (int c = 0; c < ncopy; ++c) {
    row[c] = coef_k[c];
  }
}


}  // anonymous namespace

Rcpp::NumericVector two_block_align_b_to_xhyper_cpp(
    Rcpp::NumericVector b_vec,
    Rcpp::NumericMatrix X_k,
    Rcpp::CharacterVector group_levels
) {
  const Rcpp::CharacterVector rn = matrix_rownames(X_k);
  if (rn.size() == 0 || Rf_isNull(rn)) {
    if (b_vec.size() != X_k.nrow()) {
      Rcpp::stop(
        "length(b) (%d) must equal nrow(X_hyper) (%d) when X_hyper has no rownames.",
        b_vec.size(), X_k.nrow()
      );
    }
    return b_vec;
  }

  if (b_vec.hasAttribute("names") && !Rf_isNull(b_vec.attr("names"))) {
    Rcpp::CharacterVector bnames = b_vec.attr("names");
    std::string miss;
    for (int r = 0; r < rn.size(); ++r) {
      if (Rcpp::CharacterVector::is_na(rn[r])) {
        continue;
      }
      const std::string lev = Rcpp::as<std::string>(rn[r]);
      if (!charvec_contains_name(bnames, lev)) {
        if (!miss.empty()) {
          miss += ", ";
        }
        miss += lev;
      }
    }
    if (!miss.empty()) {
      Rcpp::stop("Group level(s) missing from b: %s", miss.c_str());
    }
    Rcpp::NumericVector out(rn.size());
    for (int r = 0; r < rn.size(); ++r) {
      if (Rcpp::CharacterVector::is_na(rn[r])) {
        out[r] = NA_REAL;
        continue;
      }
      out[r] = charvec_named_value(
        b_vec, bnames, Rcpp::as<std::string>(rn[r])
      );
    }
    return out;
  }

  if (b_vec.size() != group_levels.size() || b_vec.size() != X_k.nrow()) {
    Rcpp::stop(
      "b and X_hyper row counts disagree (b: %d, X_hyper: %d, group_levels: %d).",
      b_vec.size(), X_k.nrow(), group_levels.size()
    );
  }

  std::string miss;
  for (int r = 0; r < rn.size(); ++r) {
    if (Rcpp::CharacterVector::is_na(rn[r])) {
      continue;
    }
    const std::string lev = Rcpp::as<std::string>(rn[r]);
    if (!charvec_contains_name(group_levels, lev)) {
      if (!miss.empty()) {
        miss += ", ";
      }
      miss += lev;
    }
  }
  if (!miss.empty()) {
    Rcpp::stop(
      "X_hyper rownames do not match group_levels; missing in groups: %s",
      miss.c_str()
    );
  }

  Rcpp::NumericVector out(rn.size());
  for (int r = 0; r < rn.size(); ++r) {
    if (Rcpp::CharacterVector::is_na(rn[r])) {
      out[r] = NA_REAL;
      continue;
    }
    const std::string lev = Rcpp::as<std::string>(rn[r]);
    out[r] = charvec_named_value(b_vec, group_levels, lev);
  }
  return out;
}

// Lazy stats::gaussian() for rglmb Block~2 calls (matches R batch driver).
Rcpp::List two_block_gaussian_family_r() {
  Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats");
  Rcpp::Function gaussian = stats["gaussian"];
  return Rcpp::as<Rcpp::List>(gaussian());
}

// .two_block_rglmb_iter_count() — envelope candidates from one rglmb draw.
int two_block_rglmb_iter_count_r(const Rcpp::List& fit_k) {
  if (!fit_k.containsElementNamed("iters") || Rf_isNull(fit_k["iters"])) {
    return 1;
  }
  SEXP it = fit_k["iters"];
  if (Rcpp::is<Rcpp::NumericMatrix>(it)) {
    Rcpp::NumericMatrix M(it);
    return static_cast<int>(M(0, 0));
  }
  Rcpp::NumericVector v(it);
  return static_cast<int>(v[0]);
}

// Block~2 via rglmb(..., pfamily = pf) — same path as .two_block_block2_one_chain.
Rcpp::NumericVector two_block_population_rglmb_gamma(
    const Rcpp::NumericMatrix& b_i,
    int col_j,
    const Rcpp::NumericMatrix& X_j,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::List& pfamily_j,
    bool is_ing,
    double& tau2_j,
    double& iters_add
) {
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env(GLMBAYES_R_NS);
  Rcpp::Function rglmb = pkg["rglmb"];

  const int J = group_levels.size();
  Rcpp::NumericVector b_col(J);
  for (int g = 0; g < J; ++g) {
    b_col[g] = b_i(g, col_j);
  }
  b_col.attr("names") = group_levels;

  Rcpp::NumericVector y_j = two_block_align_b_to_xhyper_cpp(
    b_col, X_j, group_levels
  );

  Rcpp::List fit_k = rglmb(
    Rcpp::Named("n") = 1,
    Rcpp::Named("y") = y_j,
    Rcpp::Named("x") = X_j,
    Rcpp::Named("family") = two_block_gaussian_family_r(),
    Rcpp::Named("pfamily") = pfamily_j,
    Rcpp::Named("verbose") = false,
    Rcpp::Named("use_parallel") = false
  );

  if (is_ing) {
    Rcpp::NumericVector disp_j =
      Rcpp::as<Rcpp::NumericVector>(fit_k["dispersion"]);
    tau2_j = disp_j[0];
    iters_add = static_cast<double>(two_block_rglmb_iter_count_r(fit_k));
  } else {
    iters_add = 1.0;
  }

  return two_block_coef_draw_from_rglmb(fit_k);
}

// Port of R .two_block_block2_one_chain() for one replicate chain.
Rcpp::List two_block_population_one_chain_cpp_export(
    const Rcpp::NumericMatrix& b_i,
    const Rcpp::List& fixef_rows,
    const Rcpp::NumericVector& tau2_i,
    const Rcpp::NumericVector& iters_i,
    const Rcpp::List& x_hyper,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::List& pfamily_list,
    const Rcpp::CharacterVector& ptypes,
    const Rcpp::CharacterVector& re_names
) {
  const int p_re = re_names.size();
  Rcpp::List fixef_out = Rcpp::clone(fixef_rows);
  Rcpp::NumericVector tau2_out = Rcpp::clone(tau2_i);
  Rcpp::NumericVector iters_out = Rcpp::clone(iters_i);

  for (int k = 0; k < p_re; ++k) {
    const std::string re_k = Rcpp::as<std::string>(re_names[k]);
    Rcpp::NumericMatrix X_k = Rcpp::as<Rcpp::NumericMatrix>(x_hyper[re_k]);
    Rcpp::List pfamily_k = Rcpp::as<Rcpp::List>(pfamily_list[re_k]);
    const std::string ptype_k = Rcpp::as<std::string>(ptypes[k]);
    const bool is_ing = (ptype_k == "dIndependent_Normal_Gamma");

    double tau2_k = tau2_out[k];
    double it_add = 0.0;
    Rcpp::NumericVector gamma_k = two_block_population_rglmb_gamma(
      b_i,
      k,
      X_k,
      group_levels,
      pfamily_k,
      is_ing,
      tau2_k,
      it_add
    );

    Rcpp::NumericVector row_k = Rcpp::as<Rcpp::NumericVector>(fixef_out[k]);
    two_block_update_fixef_row(row_k, gamma_k);
    fixef_out[k] = row_k;
    if (is_ing) {
      tau2_out[k] = tau2_k;
      iters_out[k] = iters_out[k] + it_add;
    } else {
      iters_out[k] = iters_out[k] + 1.0;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("fixef") = fixef_out,
    Rcpp::Named("tau2") = tau2_out,
    Rcpp::Named("iters") = iters_out
  );
}


}  // namespace twoblock
}  // namespace sim
}  // namespace glmbayes
