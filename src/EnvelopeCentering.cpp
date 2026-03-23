// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"

#include "Envelopefuncs.h"
#include "simfuncs.h"
#include "progress_utils.h"

using namespace Rcpp;
using namespace glmbayes::sim;
using namespace glmbayes::progress;


namespace glmbayes {

namespace env {

List EnvelopeCentering(
    NumericVector y,
    NumericMatrix x,
    NumericVector mu,
    NumericMatrix P,
    NumericVector offset,
    NumericVector wt,
    double shape,
    double rate,
    int Gridtype,
    bool verbose
) {
  const int n_beta_draws = 10000;
  const int n_rss_iter = 10;
  Rcpp::Function lm_wfit("lm.wfit");
  Rcpp::Function gaussian("gaussian");
  Rcpp::Environment glmbayes_ns = Rcpp::Environment::namespace_env("glmbayes");
  Rcpp::Function glmbfamfunc = glmbayes_ns["glmbfamfunc"];

  int n_obs = y.size();
  NumericVector ystar(n_obs);
  for (int i = 0; i < n_obs; i++) {
    ystar[i] = y[i] - offset[i];
  }

  double n_w = 0.0;
  for (int i = 0; i < wt.size(); ++i) n_w += wt[i];

  Rcpp::List fit = lm_wfit(
    Rcpp::_["x"] = x,
    Rcpp::_["y"] = ystar,
    Rcpp::_["w"] = wt
  );

  NumericVector res = fit["residuals"];
  double RSS = 0.0;
  for (int i = 0; i < res.size(); i++) {
    RSS += res[i] * res[i];
  }
  int p = Rcpp::as<int>(fit["rank"]);
  double dispersion2 = RSS / (n_obs - p);

  Rcpp::List famfunc = glmbfamfunc(gaussian());
  Rcpp::Function f2 = famfunc["f2"];
  Rcpp::Function f3 = famfunc["f3"];

  arma::mat X   = Rcpp::as<arma::mat>(x);
  arma::vec Y   = Rcpp::as<arma::vec>(y);
  arma::rowvec y_row = Y.t();
  arma::rowvec off_row = Rcpp::as<arma::rowvec>(offset);
  arma::rowvec wt_row  = Rcpp::as<arma::rowvec>(wt);

  Rcpp::List cpp_out;
  double RSS_Post2 = NA_REAL;

  if (verbose) {
    Rcpp::Rcout << "[EnvelopeCentering] Entering loop: "
                << glmbayes::progress::timestamp_cpp() << "\n";
  }

  for (int j = 0; j < n_rss_iter; ++j) {
    cpp_out = rNormalReg(
      n_beta_draws,
      y, x, mu, P, offset, wt,
      dispersion2,
      f2, f3,
      mu,
      "gaussian",
      "identity",
      Gridtype
    );

    arma::mat beta_draws = Rcpp::as<arma::mat>(cpp_out["coefficients"]);
    arma::mat lp_mat = beta_draws * X.t();
    arma::mat eta_mat = lp_mat.each_row() + off_row;
    arma::mat mu_mat = eta_mat;
    arma::mat diff = mu_mat.each_row() - y_row;
    arma::mat res_sq = diff % diff;
    arma::mat res_sq_weighted = res_sq;
    res_sq_weighted.each_row() %= wt_row;
    arma::vec RSS_temp = arma::sum(res_sq_weighted, 1);
    RSS_Post2 = arma::mean(RSS_temp);

    double shape2 = shape + n_w / 2.0;
    double rate2  = rate  + RSS_Post2 / 2.0;
    dispersion2 = rate2 / (shape2 - 1.0);
  }

  if (verbose) {
    Rcpp::Rcout << "[EnvelopeCentering] Exiting loop: "
                << glmbayes::progress::timestamp_cpp() << "\n";
  }

  return List::create(
    Named("dispersion") = dispersion2,
    Named("RSS_post")   = RSS_Post2
  );
}

}  // namespace env

}  // namespace glmbayes
