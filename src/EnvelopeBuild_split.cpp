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

#include <string>

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

EnvGrid EnvelopeGridSetup(NumericVector bStar,
                          NumericMatrix A,
                          int Gridtype,
                          int n,
                          int n_envopt,
                          bool use_opencl,
                          bool verbose) {
  if (n_envopt < 0) {
    n_envopt = n;
  }

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

  Rcpp::List size_info = EnvelopeSize(a_2, G1, Gridtype, n, n_envopt,
                                      use_opencl, verbose);
  Rcpp::List G2 = size_info["G2"];
  Rcpp::List GIndex1 = size_info["GIndex1"];

  EnvGrid grid;
  grid.E_draws = Rcpp::as<double>(size_info["E_draws"]);
  grid.GIndex = Rcpp::as<NumericMatrix>(asMat(expGrid(GIndex1)));
  grid.l2 = grid.GIndex.nrow();
  grid.G3 = Rcpp::as<NumericMatrix>(asMat(expGrid(G2)));
  grid.G4 = NumericMatrix(grid.G3.ncol(), grid.G3.nrow());
  arma::mat G3b(grid.G3.begin(), grid.G3.nrow(), grid.G3.ncol(), false);
  arma::mat G4b(grid.G4.begin(), grid.G4.nrow(), grid.G4.ncol(), false);
  G4b = trans(G3b);

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
