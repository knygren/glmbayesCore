// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include "famfuncs.h"
#include "Envelopefuncs.h"
#include <RcppParallel.h>
#include "openclPort.h"
#include "progress_utils.h"
#include "opencl.h"

using namespace Rcpp;
using namespace openclPort;
using namespace glmbayes::fam;
using namespace glmbayes::env;
using namespace glmbayes::opencl;
using namespace glmbayes::progress;

// Internal helper: run OpenCL pilot timing, print diagnostics, and prompt user.
// Not exported to R.


double f2_f3_opencl_pilot(const Rcpp::NumericMatrix& G4,
                        const Rcpp::NumericVector& y,
                        const Rcpp::NumericMatrix& x,
                        const Rcpp::NumericMatrix& mu,
                        const Rcpp::NumericMatrix& P,
                        const Rcpp::NumericVector& alpha,
                        const Rcpp::NumericVector& wt,
                        const std::string& family,
                        const std::string& link,
                        bool use_opencl,
                        bool verbose,
                        double threshold_sec=300) {
  if (!use_opencl) return NA_REAL;
  
  int m1_total = G4.ncol();
  int y_obs    = y.size();
  
  double obs_factor = std::min(1.0, 100.0 / (double)y_obs);
  double frac_A = 0.01 * obs_factor;
  double frac_B = 0.02 * obs_factor;
  
  int m1_pilot_A = std::max(100, (int)(m1_total * frac_A));
  int m1_pilot_B = std::max(m1_pilot_A + 100, (int)(m1_total * frac_B));
  m1_pilot_A = std::min(std::max(1, m1_pilot_A), m1_total);
  m1_pilot_B = std::min(std::max(1, m1_pilot_B), m1_total);
  
  double refined_est_total_sec = NA_REAL;
  double fixed_cost = 0.0;
  double per_grid_cost = 0.0;
  
  if (m1_total > 50000 && use_opencl) {
    auto slice_grid = [&](int m1_pilot) {
      Rcpp::NumericMatrix G4_pilot(G4.nrow(), m1_pilot);
      for (int j = 0; j < m1_pilot; ++j)
        for (int i = 0; i < G4.nrow(); ++i)
          G4_pilot(i, j) = G4(i, j);
      return G4_pilot;
    };
    
    
    
    // Warm-up
    auto G4_pilot_A = slice_grid(m1_pilot_A);
    
    if (verbose) {
      
      Rcpp::Rcout << "Entering f2_F3_opencl warmup: "
                  << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << "\n";
    }
    
    
    Rcpp::List warmup_A = f2_f3_opencl(family, link, G4_pilot_A, y, x, mu, P, alpha, wt, 0);
    
    if (verbose) {
      
      Rcpp::Rcout << "Exiting f2_F3_opencl warmup: "
                  << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << "\n";
    }
    
    
    
    // Timed pilot A
    auto t0A = std::chrono::high_resolution_clock::now();
    Rcpp::List pilot_A = f2_f3_opencl(family, link, G4_pilot_A, y, x, mu, P, alpha, wt, 0);
    auto t1A = std::chrono::high_resolution_clock::now();
    double time_A = std::chrono::duration<double>(t1A - t0A).count();
    
    // Timed pilot B
    auto G4_pilot_B = slice_grid(m1_pilot_B);
    auto t0B = std::chrono::high_resolution_clock::now();
    Rcpp::List pilot_B = f2_f3_opencl(family, link, G4_pilot_B, y, x, mu, P, alpha, wt, 0);
    auto t1B = std::chrono::high_resolution_clock::now();
    double time_B = std::chrono::duration<double>(t1B - t0B).count();
    
    // Estimate costs
    per_grid_cost = (time_B - time_A) / std::max(1.0, (double)(m1_pilot_B - m1_pilot_A));
    fixed_cost    = time_A - m1_pilot_A * per_grid_cost;
    
    if (per_grid_cost <= 0.0) {
      Rcpp::Rcout << "[WARNING] Negative per-grid cost — falling back to Pilot B average.\n";
      per_grid_cost = time_B / m1_pilot_B;
      fixed_cost = 0.0;
    }
    if (fixed_cost < 0.0) {
      Rcpp::Rcout << "[WARNING] Negative fixed cost — overriding to 0 and using Pilot B average.\n";
      per_grid_cost = time_B / m1_pilot_B;
      fixed_cost = 0.0;
    }
    
    // double est_time_sec = fixed_cost + per_grid_cost * m1_total;
    
    int m1_grid = std::max(1, (int)std::ceil(0.01 * (double)m1_total));
    // double est_time_sec_m1 = fixed_cost + per_grid_cost * (double)m1_grid;
    int m2_grid = std::max(1, (int)std::floor((300.0 - fixed_cost) / std::max(1e-12, per_grid_cost)));
    int m_stage_grid = std::min(m1_grid, m2_grid);
    m_stage_grid = std::min(m_stage_grid, m1_total);
    m_stage_grid = std::max(1, m_stage_grid);
    
    
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeBuild:EnvelopeEval:Pilot] Running timed grid slice of size "
                  << m_stage_grid << "...\n";
    }
    
    auto G4_pilot = slice_grid(m_stage_grid);
    auto t0p = std::chrono::high_resolution_clock::now();
    Rcpp::List pilot = f2_f3_opencl(family, link, G4_pilot, y, x, mu, P, alpha, wt, 0);
    auto t1p = std::chrono::high_resolution_clock::now();
    double time_p = std::chrono::duration<double>(t1p - t0p).count();
    
    double per_grid_sec_parallel = time_p / (double)m_stage_grid;
    refined_est_total_sec = per_grid_sec_parallel * m1_total;
    
    auto fmt_hms = [](double seconds) {
      int s = (int)std::round(seconds);
      int h = s / 3600; s %= 3600;
      int m = s / 60;   s %= 60;
      std::ostringstream oss;
      oss << h << "h " << m << "m " << s << "s";
      return oss.str();
    };
    
    if (verbose) {
      
      Rcpp::Rcout << "[EnvelopeBuild:EnvelopeEval:Pilot] Calibration elapsed = " << time_p
                  << " s for " << m_stage_grid
                  << " grid points (" << per_grid_sec_parallel << " s/grid).\n";
      
      Rcpp::Rcout << "[EnvelopeBuild:EnvelopeEval:Pilot] Refined grid build time estimate = "
                  << refined_est_total_sec << " seconds (" << fmt_hms(refined_est_total_sec) << ")\n";
      
    }
    long total = (long)std::round(refined_est_total_sec);
    long h = total / 3600;
    long m = (total % 3600) / 60;
    long s = total % 60;
    
    if (verbose) {
      
      Rcpp::Rcout << "[EnvelopeBuild:EnvelopeEval:Pilot] Estimated full f2_f3 evaluation time = "
                  << refined_est_total_sec << " seconds (" << h << "h " << m << "m " << s << "s)\n";
      
      Rcpp::Rcout << "[EnvelopeBuild:EnvelopeEval:Pilot] Components: fixed=" << fixed_cost
                  << ", per-grid=" << per_grid_cost << "\n";
      
      Rcpp::Rcout << "[EnvelopeBuild:EnvelopeEval:Pilot] Note: estimate is approximate and may vary with system load.\n"; 
      
      
    }
    
  }
  
  if (refined_est_total_sec > threshold_sec) {
    Rcpp::Rcout << "\nEstimated run time exceeds 5 minutes ("
                << refined_est_total_sec << " seconds).\n";
    
    Rcpp::Function r_interactive("interactive");
    bool is_interactive = Rcpp::as<bool>(r_interactive());
    
    if (is_interactive) {
      Rcpp::Function readline("readline");
      std::string response = Rcpp::as<std::string>(
        readline("Do you want to continue? [y/N]: ")
      );
      
      if (response != "y" && response != "Y") {
        Rcpp::Rcout << "Aborting full run per user choice.\n";
        Rcpp::stop("User aborted run due to estimated time exceeding threshold.");
      } else {
        Rcpp::Rcout << "Proceeding with full run...\n";
      }
    } else {
      // Non-interactive (e.g. CI/CRAN): auto-approve
      Rcpp::Rcout << "[NOTE] Non-interactive session: proceeding automatically.\n";
      Rcpp::Rcout << "Proceeding with full run...\n";
    }
  }
  
  return refined_est_total_sec;
}




Rcpp::List f2_f3_non_opencl(
    std::string family,
    std::string link,
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar
) {
  int l2 = b.nrow(); // grid points
  int l1 = b.ncol(); // parameters
  
  Rcpp::NumericVector NegLL(l2);
  Rcpp::NumericMatrix cbars(l2, l1);
  arma::mat cbars2(cbars.begin(), l2, l1, false);
  
  // --- binomial family ---
  if (family == "binomial" && link == "logit") {
    NegLL  = f2_binomial_logit(b, y, x, mu, P, alpha, wt, progbar);
    cbars2 = f3_binomial_logit(b, y, x, mu, P, alpha, wt, progbar);
  }
  else if (family == "binomial" && link == "probit") {
    NegLL  = f2_binomial_probit(b, y, x, mu, P, alpha, wt, progbar);
    cbars2 = f3_binomial_probit(b, y, x, mu, P, alpha, wt, progbar);
  }
  else if (family == "binomial" && link == "cloglog") {
    NegLL  = f2_binomial_cloglog(b, y, x, mu, P, alpha, wt, progbar);
    cbars2 = f3_binomial_cloglog(b, y, x, mu, P, alpha, wt, progbar);
  }
  
  // --- quasibinomial family (reuse binomial kernels) ---
  else if (family == "quasibinomial" && link == "logit") {
    NegLL  = f2_binomial_logit(b, y, x, mu, P, alpha, wt, progbar);
    cbars2 = f3_binomial_logit(b, y, x, mu, P, alpha, wt, progbar);
  }
  else if (family == "quasibinomial" && link == "probit") {
    NegLL  = f2_binomial_probit(b, y, x, mu, P, alpha, wt, progbar);
    cbars2 = f3_binomial_probit(b, y, x, mu, P, alpha, wt, progbar);
  }
  else if (family == "quasibinomial" && link == "cloglog") {
    NegLL  = f2_binomial_cloglog(b, y, x, mu, P, alpha, wt, progbar);
    cbars2 = f3_binomial_cloglog(b, y, x, mu, P, alpha, wt, progbar);
  }
  
  // --- poisson family ---
  else if (family == "poisson") {
    NegLL  = f2_poisson(b, y, x, mu, P, alpha, wt, progbar);
    cbars2 = f3_poisson(b, y, x, mu, P, alpha, wt, progbar);
  }
  
  // --- quasipoisson family (reuse poisson kernels) ---
  else if (family == "quasipoisson") {
    NegLL  = f2_poisson(b, y, x, mu, P, alpha, wt, progbar);
    cbars2 = f3_poisson(b, y, x, mu, P, alpha, wt, progbar);
  }
  
  // --- gamma family ---
  else if (family == "Gamma") {
    NegLL  = f2_gamma(b, y, x, mu, P, alpha, wt, progbar);
    cbars2 = f3_gamma(b, y, x, mu, P, alpha, wt, progbar);
  }
  
  // --- gaussian family ---
  else if (family == "gaussian") {
    NegLL  = f2_gaussian(b, y, x, mu, P, alpha, wt);
    cbars2 = f3_gaussian(b, y, x, mu, P, alpha, wt);
  }
  
  else {
    Rcpp::stop("Unsupported family/link combination in f2_f3_non_opencl: " +
      family + "/" + link);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("qf")   = NegLL,
    Rcpp::Named("grad") = cbars2
  );
}


namespace glmbayes {

namespace env {

Rcpp::List EnvelopeEval(const Rcpp::NumericMatrix& G4,   // grid (parameters × grid points)
                        const Rcpp::NumericVector& y,
                        const Rcpp::NumericMatrix& x,
                        const Rcpp::NumericMatrix& mu,
                        const Rcpp::NumericMatrix& P,
                        const Rcpp::NumericVector& alpha,
                        const Rcpp::NumericVector& wt,
                        const std::string& family,
                        const std::string& link,
                        bool use_opencl ,
                        bool verbose ) {
  int progbar = 0;
  
  // Optional pilot timing for large parameter dimension
  // (originally: if (l1 >= 14) ...; here we use number of columns in G4)
  // --- Pilot timing ---
  if (G4.ncol() >= 14) {
    Timer t_pilot; if (verbose) t_pilot.begin();
    
    if (verbose) {
      
      Rcpp::Rcout << "Entering Run_opencl_pilot: "
                  << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << "\n";
    }
    
    
    double est_time = f2_f3_opencl_pilot(G4, y, x, mu, P, alpha, wt,family, link, use_opencl, verbose);
    if (verbose) {
      print_completed("[EnvelopeBuild::EnvelopeEval] Pilot", t_pilot);
      Rcpp::Rcout << "[EnvelopeBuild::EnvelopeEval] Pilot estimated time = "
                  << est_time << " seconds\n";
    }
  }
  
  
  // --- Dispatch timing ---
  Timer t_dispatch; if (verbose) t_dispatch.begin();
  Rcpp::List prepGrad;
  if (use_opencl) {
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeBuild::EnvelopeEval] Initiating f2_f3_opencl at "
                  << now_hms() << "\n";
    }
    prepGrad = f2_f3_opencl(family, link, G4, y, x, mu, P, alpha, wt, progbar);
  } else {
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeEval] Initiating f2_f3_non_opencl at "
                  << now_hms() << "\n";
    }
    prepGrad = f2_f3_non_opencl(family, link, G4, y, x, mu, P, alpha, wt, progbar);
  }
  if (verbose) {
    print_completed("[EnvelopeBuild::EnvelopeEval] Dispatch", t_dispatch);
  }
  
  // Unpack results
  Rcpp::NumericVector NegLL = prepGrad["qf"];          // negative log likelihood values
  arma::mat cbars = Rcpp::as<arma::mat>(prepGrad["grad"]); // gradient matrix
  
  return Rcpp::List::create(
    Rcpp::Named("NegLL") = NegLL,
    Rcpp::Named("cbars") = cbars
  );
}

}
}