// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// EnvelopeBuild split around EnvelopeEval: EnvelopeGridSetup (before the
// evaluation) and EnvelopeFinalize (after it). Both copy the matching lines
// of EnvelopeBuild (EnvelopeBuild.cpp) so the result is identical.

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "famfuncs.h"
#include "Envelopefuncs.h"
#include "Envelopefuncs_split.h"
#include "progress_utils.h"
#include "R_interface.h"
#include "openclPort.h"

#include <chrono>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

using namespace Rcpp;

namespace glmbayes {

namespace env {

bool envelope_use_opencl(bool use_opencl, bool verbose) {
#ifdef USE_OPENCL
  (void) verbose;
  return use_opencl;
#else
  if (use_opencl && verbose) {
    Rcpp::Rcout << "[NOTE] OpenCL support was not detected during configuration.\n";
    Rcpp::Rcout << "       Disabling use_opencl and falling back to CPU implementation.\n";
  }
  return false;
#endif
}

namespace {

// EnvelopeOpt (R/simulationpipeline.R) step for step: rank() with average
// ties, prod() accumulated in long double as R does, last index at the
// minimum. Returns dimcount[minindex, ] (1 or 3 per dimension).
NumericVector envelope_opt_cpp(const arma::vec& a1, int n, int core_cnt) {
  if (core_cnt < 1) core_cnt = 1;
  const int l1 = a1.n_elem;

  std::vector<double> v(l1);
  for (int j = 0; j < l1; ++j) v[j] = 1 / (1 + a1[j]);
  std::vector<double> a1rank(l1);
  for (int j = 0; j < l1; ++j) {
    int less = 0, eq = 0;
    for (int k = 0; k < l1; ++k) {
      if (v[k] < v[j]) ++less;
      else if (v[k] == v[j]) ++eq;
    }
    a1rank[j] = less + (eq + 1) / 2.0;
  }

  std::vector<double> dimcount(l1, 1.0);
  std::vector<double> scaleest(l1);
  for (int j = 0; j < l1; ++j) scaleest[j] = std::sqrt(1 + a1[j]);
  std::vector<std::vector<double>> rows;
  rows.reserve(l1 + 1);
  std::vector<double> evalest(l1 + 1);

  for (int i = 1; i <= l1 + 1; ++i) {
    if (i > 1) {
      for (int j = 0; j < l1; ++j) {
        if (a1rank[j] == i - 1) {
          dimcount[j] = 3;
          scaleest[j] = 2 / std::sqrt(M_PI);
        }
      }
    }
    rows.push_back(dimcount);
    const double intest = (i == 1) ? 1.0 : std::pow(3.0, i - 1);
    long double s = 1.0;
    for (int j = 0; j < l1; ++j) s *= scaleest[j];
    const double slopeest = static_cast<double>(s);
    evalest[i - 1] = (intest / core_cnt) + n * slopeest;
  }

  double mn = evalest[0];
  for (int j = 1; j <= l1; ++j) {
    if (evalest[j] < mn) mn = evalest[j];
  }
  int minindex = 0;
  for (int j = 0; j <= l1; ++j) {
    if (evalest[j] == mn) minindex = j;
  }
  return NumericVector(rows[minindex].begin(), rows[minindex].end());
}

struct CppGrid {
  NumericVector gridindex;
  double E_draws;
  NumericMatrix GIndex;
  NumericMatrix G3;
  NumericMatrix G4;
};

// EnvelopeSize's per-dimension pattern and E_draws, then GIndex and G3 in
// expand.grid() row order (first dimension fastest) with the
// as.matrix(expand.grid(.)) dimnames, and G4 = t(G3).
CppGrid grid_setup_cpp(const arma::vec& a, const NumericMatrix& G1,
                       int Gridtype, int n_envopt, bool use_ocl,
                       GridSetupTiming* timing) {
  auto t0 = std::chrono::steady_clock::now();
  const int l1 = a.n_elem;
  CppGrid out;
  out.gridindex = NumericVector(l1);
  if (Gridtype == 2) {
    const int core_cnt =
        use_ocl ? openclPort::opencl_core_count_for_scaling() : 1;
    out.gridindex = envelope_opt_cpp(a, n_envopt, core_cnt);
  }

  std::vector<int> npts(l1);
  double E_draws = 1.0;
  for (int i = 0; i < l1; ++i) {
    bool three;
    if (Gridtype == 1) {
      three = !(std::sqrt(1 + a[i]) <= (2 / std::sqrt(M_PI)));
    } else if (Gridtype == 2) {
      three = !(out.gridindex[i] == 1);
    } else if (Gridtype == 3) {
      three = true;
    } else if (Gridtype == 4) {
      three = false;
    } else {
      Rcpp::stop("EnvelopeGridSetup: unsupported Gridtype %d", Gridtype);
    }
    npts[i] = three ? 3 : 1;
    if (three) {
      E_draws *= (2 / std::sqrt(M_PI));
    } else {
      E_draws *= std::sqrt(1 + a[i]);
    }
  }
  out.E_draws = E_draws;
  if (timing) {
    const auto t1 = std::chrono::steady_clock::now();
    timing->cpp_opt += std::chrono::duration<double>(t1 - t0).count();
    t0 = t1;
  }

  int l2 = 1;
  for (int i = 0; i < l1; ++i) l2 *= npts[i];
  out.GIndex = NumericMatrix(l2, l1);
  out.G3 = NumericMatrix(l2, l1);
  out.G4 = NumericMatrix(l1, l2);
  int stride = 1;
  for (int i = 0; i < l1; ++i) {
    for (int r = 0; r < l2; ++r) {
      const int k = (r / stride) % npts[i];
      const double code = (npts[i] == 1) ? 4.0 : k + 1.0;
      const double val = (npts[i] == 1) ? G1(1, i) : G1(k, i);
      out.GIndex(r, i) = code;
      out.G3(r, i) = val;
      out.G4(i, r) = val;
    }
    stride *= npts[i];
  }
  CharacterVector cn(l1);
  for (int i = 0; i < l1; ++i) cn[i] = "Var" + std::to_string(i + 1);
  out.GIndex.attr("dimnames") = List::create(R_NilValue, cn);
  out.G3.attr("dimnames") = List::create(R_NilValue, cn);
  if (timing) {
    timing->cpp_grid += std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - t0).count();
  }
  return out;
}

std::string type_name(SEXP x) { return Rf_type2char(TYPEOF(x)); }

// Stop unless `c` (C++) equals `r` (R) in type, dim, dimnames, every value
// bit for bit and, with `all_attrib`, every other attribute.
void compare_setup_field(const char* field, SEXP r, SEXP c, bool all_attrib) {
  if (TYPEOF(r) != TYPEOF(c)) {
    Rcpp::stop("setup compare: %s type differs (R %s, C++ %s)", field,
               type_name(r), type_name(c));
  }
  if (!R_compute_identical(Rf_getAttrib(r, R_DimSymbol),
                           Rf_getAttrib(c, R_DimSymbol), 0)) {
    Rcpp::stop("setup compare: %s dim differs", field);
  }
  if (all_attrib &&
      !R_compute_identical(Rf_getAttrib(r, R_DimNamesSymbol),
                           Rf_getAttrib(c, R_DimNamesSymbol), 0)) {
    Rcpp::stop("setup compare: %s dimnames differ", field);
  }
  if (Rf_xlength(r) != Rf_xlength(c)) {
    Rcpp::stop("setup compare: %s length differs (R %d, C++ %d)", field,
               static_cast<long long>(Rf_xlength(r)),
               static_cast<long long>(Rf_xlength(c)));
  }
  if (TYPEOF(r) == REALSXP) {
    const double* pr = REAL(r);
    const double* pc = REAL(c);
    for (R_xlen_t i = 0; i < Rf_xlength(r); ++i) {
      if (std::memcmp(pr + i, pc + i, sizeof(double)) != 0) {
        Rcpp::stop("setup compare: %s value differs at element %d "
                   "(R %.17g, C++ %.17g)",
                   field, static_cast<long long>(i + 1), pr[i], pc[i]);
      }
    }
  }
  if (all_attrib && !R_compute_identical(r, c, 1)) {
    Rcpp::stop("setup compare: %s attributes differ", field);
  }
}

}  // namespace

EnvGrid EnvelopeGridSetup(NumericVector bStar,
                          NumericMatrix A,
                          int Gridtype,
                          int n,
                          int n_envopt,
                          bool use_opencl,
                          bool verbose,
                          GridSetupTiming* timing,
                          int setup_method) {
  if (n_envopt < 0) {
    n_envopt = n;
  }
  auto t0 = std::chrono::steady_clock::now();
  auto lap = [&](double GridSetupTiming::*field) {
    const auto t1 = std::chrono::steady_clock::now();
    if (timing) {
      timing->*field += std::chrono::duration<double>(t1 - t0).count();
    }
    t0 = t1;
  };

  int l1 = A.nrow(), k = A.ncol();
  arma::mat A2(A.begin(), l1, k, false);
  arma::colvec bStar_2(bStar.begin(), bStar.size(), false);

  NumericVector a_1(l1);
  arma::vec a_2(a_1.begin(), a_1.size(), false);

  NumericVector xx_1(3, 1.0);
  NumericVector xx_2 = NumericVector::create(-1.0, 0.0, 1.0);
  NumericVector yy_1(2, 1.0);
  NumericVector yy_2 = NumericVector::create(-0.5, 0.5);
  NumericMatrix G1(3, l1);
  NumericMatrix Lint1(2, l1);
  arma::mat G1b(G1.begin(), 3, l1, false);
  arma::mat Lint(Lint1.begin(), 2, l1, false);

  arma::colvec xx_1b(xx_1.begin(), xx_1.size(), false);
  arma::colvec xx_2b(xx_2.begin(), xx_2.size(), false);
  arma::colvec yy_1b(yy_1.begin(), yy_1.size(), false);
  arma::colvec yy_2b(yy_2.begin(), yy_2.size(), false);
  Rcpp::Function expGrid = glmbayes_R::r_expand_grid();
  Rcpp::Function asMat = glmbayes_R::r_as_matrix();

  a_2 = arma::diagvec(A2);
  arma::vec omega = (sqrt(2) - arma::exp(-1.20491 - 0.7321 * sqrt(0.5 + a_2))) / arma::sqrt(1 + a_2);
  G1b = xx_1b * arma::trans(bStar_2) + xx_2b * arma::trans(omega);
  Lint = yy_1b * arma::trans(bStar_2) + yy_2b * arma::trans(omega);
  lap(&GridSetupTiming::arith);

  EnvGrid grid;
  NumericVector gridindex_r;
  if (setup_method != SETUP_CPP) {
    Rcpp::List size_info = EnvelopeSize(a_2, G1, Gridtype, n, n_envopt,
                                        use_opencl, verbose, timing);
    Rcpp::List G2 = size_info["G2"];
    Rcpp::List GIndex1 = size_info["GIndex1"];
    gridindex_r = size_info["gridindex"];

    grid.E_draws = Rcpp::as<double>(size_info["E_draws"]);
    t0 = std::chrono::steady_clock::now();
    const RObject gi_df = expGrid(GIndex1);
    lap(&GridSetupTiming::gindex_expand);
    grid.GIndex = Rcpp::as<NumericMatrix>(asMat(gi_df));
    lap(&GridSetupTiming::gindex_asmat);
    grid.l2 = grid.GIndex.nrow();
    const RObject g3_df = expGrid(G2);
    lap(&GridSetupTiming::g3_expand);
    grid.G3 = Rcpp::as<NumericMatrix>(asMat(g3_df));
    lap(&GridSetupTiming::g3_asmat);
    grid.G4 = NumericMatrix(grid.G3.ncol(), grid.G3.nrow());
    arma::mat G3b(grid.G3.begin(), grid.G3.nrow(), grid.G3.ncol(), false);
    arma::mat G4b(grid.G4.begin(), grid.G4.nrow(), grid.G4.ncol(), false);
    G4b = trans(G3b);
    lap(&GridSetupTiming::g4);
  }

  if (setup_method != SETUP_R) {
    const CppGrid cg =
        grid_setup_cpp(a_2, G1, Gridtype, n_envopt, use_opencl, timing);
    if (setup_method == SETUP_CPP) {
      grid.E_draws = cg.E_draws;
      grid.GIndex = cg.GIndex;
      grid.G3 = cg.G3;
      grid.G4 = cg.G4;
      grid.l2 = cg.G3.nrow();
    } else {
      if (Gridtype == 2) {
        compare_setup_field("gridindex", gridindex_r, cg.gridindex, false);
      }
      if (std::memcmp(&grid.E_draws, &cg.E_draws, sizeof(double)) != 0) {
        Rcpp::stop("setup compare: E_draws differs (R %.17g, C++ %.17g)",
                   grid.E_draws, cg.E_draws);
      }
      compare_setup_field("GIndex", grid.GIndex, cg.GIndex, true);
      compare_setup_field("G3", grid.G3, cg.G3, true);
      compare_setup_field("G4", grid.G4, cg.G4, true);
      if (timing) timing->compared += 1.0;
    }
  }

  if (verbose) {
    Rcpp::Rcout << "[EnvelopeGridSetup] Grid size (l2) = "
                << glmbayes::progress::format_int_with_commas(grid.l2)
                << "\n";
  }

  grid.Lint1 = Lint1;
  grid.a_1 = a_1;
  grid.l1 = l1;
  return grid;
}

List EnvelopeFinalize(const EnvGrid& grid,
                      NumericVector NegLL,
                      const arma::mat& cbars_eval,
                      NumericMatrix x,
                      bool sortgrid,
                      bool verbose) {
  const int l1 = grid.l1;
  const int l2 = grid.l2;
  NumericMatrix GIndex = grid.GIndex;
  NumericMatrix G3 = grid.G3;
  NumericMatrix Lint1 = grid.Lint1;
  NumericVector a_1 = grid.a_1;
  const double E_draws = grid.E_draws;

  NumericMatrix cbars(l2, l1);
  NumericMatrix Up(l2, l1);
  NumericMatrix Down(l2, l1);
  NumericMatrix logP(l2, 2);
  NumericMatrix logU(l2, l1);
  NumericMatrix loglt(l2, l1);
  NumericMatrix logrt(l2, l1);
  NumericMatrix logct(l2, l1);
  NumericMatrix LLconst(l2, 1);

  arma::mat cbars2(cbars.begin(), l2, l1, false);
  cbars2 = cbars_eval;

  if (verbose) {
    Rcpp::Rcout << "[EnvelopeFinalize:EnvelopeSet_Grid] Entering: "
                << glmbayes::progress::timestamp_cpp() << "\n";
  }
  EnvelopeSet_Grid_C2_pointwise(GIndex, cbars, Lint1, Down, Up, loglt, logrt,
                                logct, logU, logP);
  setlogP_C2(logP, NegLL, cbars, G3, LLconst);

  NumericMatrix::Column logP2 = logP(_, 1);
  double maxlogP = max(logP2);
  NumericVector PLSD = exp(logP2 - maxlogP);
  double sumP = sum(PLSD);
  PLSD = PLSD / sumP;

  bool bad = false;
  for (int i = 0; i < PLSD.size(); ++i) {
    double v = PLSD[i];
    if (!R_finite(v) || v <= 0.0) {
      Rcout << "[ERROR] Invalid PLSD at index " << i
            << "  value=" << v
            << "  logP2=" << logP2[i]
            << "  maxlogP=" << maxlogP
            << "\n";
      bad = true;
      break;
    }
  }

  if (bad) {
    arma::mat X(x.begin(), x.nrow(), x.ncol(), false);
    arma::vec s = arma::svd(X);
    double kappa_X   = s.max() / s.min();
    double kappa_XtX = (s.max() * s.max()) / (s.min() * s.min());
    double maxlogP2 = max(logP2);
    double minlogP2 = min(logP2);
    double spread  = maxlogP2 - minlogP2;
    Rcpp::stop(
      "[EnvelopeFinalize] PLSD construction failed for this model.\n"
      "  - kappa(X)           ~ " + std::to_string(kappa_X)   + "\n"
      "  - kappa(t(X) %*% X)  ~ " + std::to_string(kappa_XtX) + "\n"
      "  - PLSD log-weights span " + std::to_string(spread) + " on the log scale.\n"
      "  The standardized design is highly collinear, so the fixed PLSD\n"
      "  partition cannot remain tight and the PLSD weights blow up.\n"
    );
  }

  if (sortgrid) {
    Rcpp::Function EnvSort = glmbayes_R::r_envelope_sort();
    Rcpp::List outlist = EnvSort(l1, l2, GIndex, G3, cbars, logU, logrt,
                                 loglt, logP, LLconst, PLSD, a_1, E_draws);
    if (outlist.containsElementNamed("sort_ok") &&
        !Rcpp::as<bool>(outlist["sort_ok"])) {
      if (verbose) {
        Rcpp::Rcout << "[EnvelopeFinalize:EnvelopeSort] Using unsorted envelope (memory fallback).\n";
      }
      outlist = Rcpp::List::create(
        Rcpp::Named("GridIndex") = GIndex,
        Rcpp::Named("thetabars") = G3,
        Rcpp::Named("cbars")    = cbars,
        Rcpp::Named("logU")     = logU,
        Rcpp::Named("logrt")    = logrt,
        Rcpp::Named("loglt")    = loglt,
        Rcpp::Named("LLconst")  = LLconst,
        Rcpp::Named("logP")     = logP(_, 0),
        Rcpp::Named("PLSD")     = PLSD,
        Rcpp::Named("a1")       = a_1,
        Rcpp::Named("E_draws")  = E_draws
      );
    }
    return outlist;
  }

  return Rcpp::List::create(Rcpp::Named("GridIndex") = GIndex,
                            Rcpp::Named("thetabars") = G3,
                            Rcpp::Named("cbars") = cbars,
                            Rcpp::Named("logU") = logU,
                            Rcpp::Named("logrt") = logrt,
                            Rcpp::Named("loglt") = loglt,
                            Rcpp::Named("LLconst") = LLconst,
                            Rcpp::Named("logP") = logP(_, 0),
                            Rcpp::Named("PLSD") = PLSD,
                            Rcpp::Named("a1") = a_1,
                            Rcpp::Named("E_draws") = E_draws);
}

} // namespace env

} // namespace glmbayes
