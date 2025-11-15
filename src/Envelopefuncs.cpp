// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include "famfuncs.h"
#include "Set_Grid.h"
#include "Envelopefuncs.h"
#include "kernel_wrappers.h"

using namespace Rcpp;



// Internal helper: run OpenCL pilot timing, print diagnostics, and prompt user.
// Not exported to R.
// [[Rcpp::export]]

double run_opencl_pilot(const Rcpp::NumericMatrix& G4,
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
                        double threshold_sec = 300.0) {
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
    Rcpp::List warmup_A = f2_f3_opencl(family, link, G4_pilot_A, y, x, mu, P, alpha, wt, 0);
    
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
    
    double est_time_sec = fixed_cost + per_grid_cost * m1_total;
    
    int m1_grid = std::max(1, (int)std::ceil(0.01 * (double)m1_total));
    double est_time_sec_m1 = fixed_cost + per_grid_cost * (double)m1_grid;
    int m2_grid = std::max(1, (int)std::floor((300.0 - fixed_cost) / std::max(1e-12, per_grid_cost)));
    int m_stage_grid = std::min(m1_grid, m2_grid);
    m_stage_grid = std::min(m_stage_grid, m1_total);
    m_stage_grid = std::max(1, m_stage_grid);
    
    if (verbose) {
      // Optional debug prints (commented in your original)
      // Rcpp::Rcout << "[GRID EST] m1_grid = " << m1_grid << "\n";
      // Rcpp::Rcout << "[GRID EST] est_time_sec_m1 = " << est_time_sec_m1 << "\n";
      // Rcpp::Rcout << "[GRID EST] m2_grid = " << m2_grid << "\n";
      // Rcpp::Rcout << "[GRID EST] Selected m_stage_grid = " << m_stage_grid << "\n";
      // Rcpp::Rcout << "[GRID EST] m1_pilot_A = " << m1_pilot_A << "\n";
      // Rcpp::Rcout << "[GRID EST] m1_pilot_B = " << m1_pilot_B << "\n";
    }
    
    Rcpp::Rcout << "Running timed grid slice of size " << m_stage_grid << "...\n";
    
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
    
    Rcpp::Rcout << "[GRID CALIB] Calibration elapsed = " << time_p
                << " s for " << m_stage_grid
                << " grid points (" << per_grid_sec_parallel << " s/grid).\n";
    
    Rcpp::Rcout << "Refined grid build time estimate: " << refined_est_total_sec
                << " seconds (" << fmt_hms(refined_est_total_sec) << ").\n";
    
    long total = (long)std::round(refined_est_total_sec);
    long h = total / 3600;
    long m = (total % 3600) / 60;
    long s = total % 60;
    
    Rcpp::Rcout << "Estimated full f2_f3 evaluation time: "
                << refined_est_total_sec << " seconds ("
                << h << "h " << m << "m " << s << "s)"
                << " (fixed=" << fixed_cost
                << ", per-grid=" << per_grid_cost << ").\n"
                << "Note: estimate is approximate and may vary with system load.\n";
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
/*
 EnvelopeBuild_c — envelope grid construction for models with Gaussian priors
 and log-concave likelihoods
 
 Purpose
 Construct an axis-aligned envelope grid and the auxiliary objects required
 by Set_Grid and the envelope sampler. This routine assumes a multivariate
 normal prior (precision-like matrix A) and a log-concave likelihood so that
 the (transformed) posterior is unimodal and well-behaved near the mode.
 The routine is family-agnostic otherwise: envelope placement uses only the
 posterior center (bStar) and curvature-like information (diag(A)).
 
 Lint construction (per-dimension cutpoints)
 - Lint is a 2 × l1 matrix; column i defines the central interval cutpoints
 for dimension i:
 Lint[0,i] := ℓ_{i,1} = θ*_i − 0.5 · ω_i
 Lint[1,i] := ℓ_{i,2} = θ*_i + 0.5 · ω_i
 where θ* = bStar (the posterior mode on the parameter scale used for the
 envelope) and ω is a spread parameter derived from the i-th diagonal of A.
 - Implementation details:
 a_2   = diag(A)            // curvature/precision per dimension
 omega = (sqrt(2) - exp(-1.20491 - 0.7321 * sqrt(0.5 + a_2))) / sqrt(1 + a_2)
 yy_1  = [1, 1]; yy_2 = [-0.5, 0.5]
 Lint  = yy_1 %*% transpose(bStar) + yy_2 %*% transpose(omega)
 The constants in the omega formula are empirical calibrations chosen to
 produce robust interval widths across typical precision regimes.
 
 Why Gaussian-prior + log-concave-likelihood matters
 - Log-concavity of the likelihood combined with a Gaussian prior implies the
 (unnormalized) log-posterior is concave around the mode; a single-mode
 assumption justifies axis-aligned intervals centered at the posterior mode.
 - Using diag(A) to set ω leverages local curvature information: larger
 precision (larger a_i) yields narrower ω_i, smaller a_i yields wider ω_i.
 - These choices ensure the central interval ℓ_{i,1}..ℓ_{i,2} captures the
 high-density mass along each axis while keeping tails handled by the
 left/right intervals.
 
 Gridpoint generation and Gridtype
 - G1 contains candidate per-dimension points {θ* - ω, θ*, θ* + ω}.
 - Gridtype controls whether a dimension uses the three-point set or only the
 mode:
 Gridtype 1: static threshold test using (1 + a_i) ≤ 2/√π → single-point
 Gridtype 2: dynamic selection via EnvelopeOpt(a_2, n) (cost-based)
 Gridtype 3: always three-point
 Gridtype 4: always single-point
 - Rationale: trade-off between build cost and sampling cost; with larger n
 or available parallelism, richer grids can be worthwhile.
 
 Full-grid expansion and cbars
 - Per-dimension G2 lists are combined with expand.grid to form the full set
 of grid locations (G3) and corresponding region codes (GIndex).
 - G4 = transpose(G3); cbars[j,i] = G4[j,i] − bStar[i] is the j-th component's
 offset from the mode for dimension i.
 - In Set_Grid these cbars shift Lint per-row to produce Down[j,i] and Up[j,i],
 the final bounds for truncated-normal evaluations.
 
 Numerical and modelling notes
 - The envelope construction uses only bStar and diag(A); any family-specific
 transforms (e.g., link functions, canonical transforms) should be applied
 upstream so bStar and A are expressed on the scale used by the envelope.
 - Log-concavity guarantees the mode is meaningful for centering intervals;
 if the likelihood is not log-concave the assumptions behind axis-aligned
 modal envelopes break down and the envelope builder may need alternative
 strategies.
 - Keep the Lint construction here as the single source of truth; if ω
 calibration constants are refined, update them only in this function.
 
 OpenCL hinting and verbosity
 - If compiled with OpenCL and use_opencl == true, the function may scale n
 by the detected core count to exploit parallel envelope optimization.
 - If OpenCL is unavailable at compile-time, the function logs a diagnostic
 (when verbose) and disables use_opencl, falling back to CPU.
 
 Outputs used downstream
 - G2: per-dimension gridpoint lists
 - G3/G4: full grid and transpose
 - GIndex: integer region codes per grid component and dimension
 - Lint: two unshifted cutpoints per dimension (lower, upper)
 - cbars: per-component offsets from bStar used to shift Lint in Set_Grid
 - These objects enable Set_Grid to compute per-dimension truncated-normal
 probabilities U_{j,i} and their log-sums logP[j], suitable for envelope
 weighting and acceptance tests under the Gaussian-prior + log-concave
 likelihood assumption.
 */


// [[Rcpp::export]]

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



// [[Rcpp::export]]
Rcpp::List EnvelopeSize(const arma::vec& a,
                        const Rcpp::NumericMatrix& G1,
                        int Gridtype   = 2,
                        int n          = 1000,
                        int n_envopt   = -1,
                        bool use_opencl = false,
                        bool verbose    = false) 
  {
  

  int l1 = a.size();
  Rcpp::List G2(l1);
  Rcpp::List GIndex1(l1);
  double E_draws = 1.0;


  // core count for scaling
  int core_CNT = get_opencl_core_count();
  if (verbose) {
    Rcpp::Rcout << "[INFO] OpenCL core count = " << core_CNT << "\n";
  }

  
  // 3.4.1 - EnvelopOpt
  
  /// If GridType=2, then the Size of the Grid is optimized for performance
  /// while factoring in the tradeoff between a large grid 
  /// (more time consuming /expensive) to build and the acceptance rate 
  /// (which gets better with a larger grid). The number of desired draws
  /// are also factored in as the importance of a high acceptance rate
  /// is more important when the number of draws is greater.
  /// In addition, the call also fators in the number of cores 
  /// since EnvelopeConstruction can occur in parallel. This is treated
  /// as the equivalent of a greater number of draws
  /// so a larger grid is generally constructed when OpenCL is enabled
  /// 
  /// If GridType is not equal to 2 then the size of the Grid is determined 
  /// uniquely by that setting
  

    
  // EnvelopeOpt is an R function
  Rcpp::Function EnvelopeOpt("EnvelopeOpt");
  Rcpp::NumericVector gridindex(l1);
  
  if (Gridtype == 2) {
    if (use_opencl) {
      gridindex = EnvelopeOpt(a, n_envopt, core_CNT);
    } else {
      gridindex = EnvelopeOpt(a, n_envopt, 1);
    }
  }
  

  
  // Loop over dimensions
  for (int i = 0; i < l1; i++) {
    Rcpp::NumericVector Temp1 = G1(_, i);
    double Temp2 = G1(1, i);
    
    if (Gridtype == 1) {
      if (std::sqrt(1 + a[i]) <= (2 / std::sqrt(M_PI))) {
        G2[i] = Rcpp::NumericVector::create(Temp2);
        GIndex1[i] = Rcpp::NumericVector::create(4.0);
        E_draws *= std::sqrt(1 + a[i]);
      } else {
        G2[i] = Rcpp::NumericVector::create(Temp1(0), Temp1(1), Temp1(2));
        GIndex1[i] = Rcpp::NumericVector::create(1.0, 2.0, 3.0);
        E_draws *= (2 / std::sqrt(M_PI));
      }
    }
    else if (Gridtype == 2) {
      if (gridindex[i] == 1) {
        G2[i] = Rcpp::NumericVector::create(Temp2);
        GIndex1[i] = Rcpp::NumericVector::create(4.0);
        E_draws *= std::sqrt(1 + a[i]);
      } else {
        G2[i] = Rcpp::NumericVector::create(Temp1(0), Temp1(1), Temp1(2));
        GIndex1[i] = Rcpp::NumericVector::create(1.0, 2.0, 3.0);
        E_draws *= (2 / std::sqrt(M_PI));
      }
    }
    else if (Gridtype == 3) {
      G2[i] = Rcpp::NumericVector::create(Temp1(0), Temp1(1), Temp1(2));
      GIndex1[i] = Rcpp::NumericVector::create(1.0, 2.0, 3.0);
      E_draws *= (2 / std::sqrt(M_PI));
    }
    else if (Gridtype == 4) {
      G2[i] = Rcpp::NumericVector::create(Temp2);
      GIndex1[i] = Rcpp::NumericVector::create(4.0);
      E_draws *= std::sqrt(1 + a[i]);
    }
  }
  

  return Rcpp::List::create(
    Rcpp::Named("G2")       = G2,
    Rcpp::Named("GIndex1")  = GIndex1,
    Rcpp::Named("E_draws")  = E_draws,
    Rcpp::Named("gridindex")= gridindex
  );
}



// [[Rcpp::export]]
Rcpp::List EnvelopeEval(const Rcpp::NumericMatrix& G4,   // grid (parameters × grid points)
                        const Rcpp::NumericVector& y,
                        const Rcpp::NumericMatrix& x,
                        const Rcpp::NumericMatrix& mu,
                        const Rcpp::NumericMatrix& P,
                        const Rcpp::NumericVector& alpha,
                        const Rcpp::NumericVector& wt,
                        const std::string& family,
                        const std::string& link,
                        bool use_opencl = false,
                        bool verbose = false) {
  int progbar = 0;
  
  // Optional pilot timing for large parameter dimension
  // (originally: if (l1 >= 14) ...; here we use number of columns in G4)
  if (G4.ncol() >= 14) {
    double est_time = run_opencl_pilot(G4, y, x, mu, P, alpha, wt,
                                       family, link, use_opencl, verbose);
    if (verbose) {
      Rcpp::Rcout << "[INFO] OpenCL pilot estimated time = "
                  << est_time << " seconds\n";
    }
  }
  
  
  // Dispatch to OpenCL or CPU evaluation
  Rcpp::List prepGrad;
//  if (use_opencl && family != "gaussian") {
if (use_opencl) {
  if (verbose) {
    Rcpp::Rcout << "Initiating f2_f3_opencl: "
                << Rcpp::as<std::string>(
    Rcpp::Function("format")(Rcpp::Function("Sys.time")()))
    << "\n";
  }
  
  prepGrad = f2_f3_opencl(family, link, G4, y, x, mu, P, alpha, wt, progbar);
  
  // temporary guard: bail out immediately for gaussian
  // if (family == "gaussian") {
  //   Rcpp::stop("OpenCL gaussian path not yet implemented — stopping after prepGrad");
  // }
}
      
         else {
    if (verbose) {
      
      Rcpp::Rcout << "Initiating f2_f3_non_opencl: "                  << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << "\n";
    }
    
    prepGrad = f2_f3_non_opencl(family, link, G4, y, x, mu, P, alpha, wt, progbar);
  }
  
  // Unpack results
  Rcpp::NumericVector NegLL = prepGrad["qf"];          // negative log likelihood values
  arma::mat cbars = Rcpp::as<arma::mat>(prepGrad["grad"]); // gradient matrix
  
  return Rcpp::List::create(
    Rcpp::Named("NegLL") = NegLL,
    Rcpp::Named("cbars") = cbars
  );
}

// [[Rcpp::export("Inv_f3_with_disp")]]
arma::mat Inv_f3_with_disp(Rcpp::List cache,
                           double dispersion,
                           Rcpp::NumericMatrix cbars_small) {
  arma::mat Pmat    = cache["Pmat"];
  arma::mat Pmu     = cache["Pmu"];
  arma::vec base_B0 = cache["base_B0"];
  arma::mat base_A  = cache["base_A"];
  
  // Scale the base terms
  arma::vec B0 = base_B0 / dispersion + Pmu;
  arma::mat A  = Pmat + base_A / dispersion;
  A = 0.5 * (A + A.t());
  
  arma::mat R = arma::chol(A);
  
  // Wrap cbars_small into an Armadillo view
  arma::mat Csmall(cbars_small.begin(), Pmat.n_rows, cbars_small.ncol(), false);
  
  // Use Armadillo's n_cols
  arma::mat Out(Pmat.n_rows, Csmall.n_cols);
  
  for (arma::uword i = 0; i < Csmall.n_cols; i++) {
    arma::vec cbars_i(Csmall.colptr(i), Pmat.n_rows, false);
    arma::vec b = -cbars_i + B0;
    
    arma::vec ytmp = arma::solve(arma::trimatl(R.t()), b);
    arma::vec sol  = arma::solve(arma::trimatu(R), ytmp);
    
    Out.col(i) = -sol;
  }
  
  return Out.t(); // m × p
}





// [[Rcpp::export(".EnvelopeBuild_cpp")]]

List EnvelopeBuild_c(NumericVector bStar,
                     NumericMatrix A, /// Diagonal Precision Matrix for Adjusted Likelihood Function
                     NumericVector y, 
                     NumericMatrix x,
                     NumericMatrix mu,
                     NumericMatrix P, /// Part of the prior precision matrix that is shifted to the likelihood
                     NumericVector alpha,
                     NumericVector wt,
                     std::string family,
                     std::string link,
                     int Gridtype, 
                     int n,
                     int n_envopt,
                     bool sortgrid,
                     bool use_opencl ,        // Enables OpenCL acceleration during envelope construction
                     bool verbose             // Enables diagnostic output
                     
){
  
  if (n_envopt < 0) {
    n_envopt = n;
  }  

  // Handle OpenCL availability at compile/runtime
  // (fall back to CPU if requested but not supported)
  
  

#ifdef USE_OPENCL
  // OpenCL support detected at compile time — proceed as requested
#else
  if (use_opencl) {
    if (verbose) {
      Rcpp::Rcout << "[NOTE] OpenCL support was not detected during configuration.\n";
      Rcpp::Rcout << "       Disabling use_opencl and falling back to CPU implementation.\n";
      Rcpp::Rcout << "       To enable OpenCL, install an OpenCL SDK and ensure CL/cl.h is discoverable.\n";
      Rcpp::Rcout << "       You may need to set OPENCL_HOME or add the SDK to your system PATH.\n";
    }
    use_opencl = false;
  }
#endif
  
  // Print call parameters if verbose

  if (verbose) {
    Rcpp::Rcout << ">>> EnvelopeBuild_c called with:\n"
                << "    Gridtype   = " << Gridtype << "\n"
                << "    n          = " << n << "\n"
                << "    use_opencl = " << use_opencl << "\n"
                << "    sortgrid   = " << sortgrid << "\n";
  }
  
  
  // Basic setup: dimensions, Armadillo views, and working vectors
  // l1 = number of parameters, k = number of predictors
  

  int progbar=0;
  
  int l1 = A.nrow(), k = A.ncol();
  arma::mat A2(A.begin(), l1, k, false);
  arma::colvec bStar_2(bStar.begin(), bStar.size(), false);
  
  
  NumericVector a_1(l1);
  arma::vec a_2(a_1.begin(), a_1.size(), false);
  
  NumericVector xx_1(3, 1.0);
  NumericVector xx_2=NumericVector::create(-1.0,0.0,1.0);
  NumericVector yy_1(2, 1.0);
  NumericVector yy_2=NumericVector::create(-0.5,0.5);
  NumericMatrix G1(3,l1);
  NumericMatrix Lint1(2,l1);
  arma::mat G1b(G1.begin(), 3, l1, false);
  arma::mat Lint(Lint1.begin(), 2, l1, false);
  
  arma::colvec xx_1b(xx_1.begin(), xx_1.size(), false);
  arma::colvec xx_2b(xx_2.begin(), xx_2.size(), false);
  arma::colvec yy_1b(yy_1.begin(), yy_1.size(), false);
  arma::colvec yy_2b(yy_2.begin(), yy_2.size(), false);
  // List G2(a_1.size());
  // List GIndex1(a_1.size());
  Rcpp::Function EnvelopeOpt("EnvelopeOpt");
  Rcpp::Function expGrid("expand.grid");
  Rcpp::Function asMat("as.matrix");
  Rcpp::Function EnvSort("EnvelopeSort");
  
  int i;  
  
  
  // Construct per-dimension tangent points (G1) and linear intercepts (Lint)
  // using diagonal precisions and offsets (omega)
  
  a_2=arma::diagvec(A2);
  arma::vec omega=(sqrt(2)-arma::exp(-1.20491-0.7321*sqrt(0.5+a_2)))/arma::sqrt(1+a_2);
  G1b=xx_1b*arma::trans(bStar_2)+xx_2b*arma::trans(omega);
  Lint=yy_1b*arma::trans(bStar_2)+yy_2b*arma::trans(omega);
  
  
  Rcpp::Rcout << "[DEBUG] 3.0" << std::endl;
  
  // Call EnvelopeSize to determine grid structure and expected draws

    Rcpp::List size_info = EnvelopeSize(a_2, G1,
                                      Gridtype,
                                      n,
                                      n_envopt,
                                      use_opencl,
                                      verbose);
  
  
  // Unpack results
  Rcpp::List G2        = size_info["G2"];
  Rcpp::List GIndex1   = size_info["GIndex1"];
  double E_draws       = Rcpp::as<double>(size_info["E_draws"]);
  Rcpp::NumericVector gridindex = size_info["gridindex"];

// 
// Expand grid indices and candidate points (GIndex, G3, G4)
// l2 = total number of grid combinations


   NumericMatrix GIndex=asMat(expGrid(GIndex1));
   int l2=GIndex.nrow();
   NumericMatrix G3=asMat(expGrid(G2));
   NumericMatrix G4(G3.ncol(),G3.nrow());
   
  
  // print grid size if requested
  if (verbose) {
    Rcpp::Rcout 
    << ">>> EnvelopeBuild_c: grid size (l2) = " 
    << l2 
    << "\n";
  }
  

  
  arma::mat G3b(G3.begin(), G3.nrow(), G3.ncol(), false);
  arma::mat G4b(G4.begin(), G4.nrow(), G4.ncol(), false);
  
  G4b=trans(G3b);
  
  // Allocate containers for evaluation results (cbars, NegLL, logP, etc.)
  
  NumericMatrix cbars(l2,l1);
  NumericMatrix Up(l2,l1);
  NumericMatrix Down(l2,l1);
  NumericMatrix logP(l2,2);
  NumericMatrix logU(l2,l1);
  NumericMatrix loglt(l2,l1);
  NumericMatrix logrt(l2,l1);
  NumericMatrix logct(l2,l1);
  
  NumericMatrix LLconst(l2,1);
  NumericVector NegLL(l2);    
  
  NumericVector NegLL_Alt(l2);    /// Temporary
  
  
  arma::mat cbars2(cbars.begin(), l2, l1, false); 
  arma::mat cbars3(cbars.begin(), l2, l1, false); 
  
  // Note: NegLL_2 only added to allow for QC printing of results 
  
  arma::colvec NegLL_2(NegLL.begin(), NegLL.size(), false);
  
  // Call EnvelopeEval to compute negative log-likelihood and gradients
  // at each grid point
  


  Rcpp::List eval_info = EnvelopeEval(G4, y, x, mu, P, alpha, wt,
                                      family, link, use_opencl, verbose);
  

  // Copy results into cbars/NegLL structures used downstream
  
  NegLL = eval_info["NegLL"];
  cbars2 = Rcpp::as<arma::mat>(eval_info["cbars"]);
  
  // Do a temporary correction here cbars3 should point to correct memory
  // See if this sets cbars
  
  cbars3=cbars2;
  
  // July 2025 - Parallelization Implementation in steps
  
  // 1) Set_Grid_C2_pointwise changes loop to enable parallel processing (suggested by Copilot)
  
  
  //  Rcpp::Rcout << "Entering Set grid C2 pointwise: "
  //              << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
  //              << "\n";
  
  if (verbose) {
    
    Rcpp::Rcout << "Setting Grid: "
                << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                << "\n";
  }
  
  // Set Grid
  
  //  Set_Grid_C2(GIndex, cbars, Lint1,Down,Up,loglt,logrt,logct,logU,logP);
  Set_Grid_C2_pointwise(GIndex, cbars, Lint1,Down,Up,loglt,logrt,logct,logU,logP);
  
  
  
  //    Rcpp::Rcout << "Entering setlogP_C2: "
  //                << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
  //                << "\n";
  
  
// Set LOG P

  if (verbose) {
    
    Rcpp::Rcout << "Setting logP: "
                << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                << "\n";
  }
  
  
  setlogP_C2(logP,NegLL,cbars,G3,LLconst);
  
  // Normalize probabilities (PLSD) for likelihood subgradient densities
  
  NumericMatrix::Column logP2 = logP( _, 1);
  double  maxlogP=max(logP2);

  NumericVector PLSD=exp(logP2-maxlogP);
  
  double sumP=sum(PLSD);
  
  PLSD=PLSD/sumP;
  
  // Optionally sort grid for efficiency if sortgrid = TRUE
  
  if(sortgrid==true){
    
    if (verbose) {
      
      Rcpp::Rcout << "Sorting Grid: "
                  << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << "\n";
    }
    
    
    Rcpp::List outlist=EnvSort(l1,l2,GIndex,G3,cbars,logU,logrt,loglt,logP,LLconst,PLSD,a_1,E_draws);
    
    return(outlist);
    
  }
  
  // Return assembled envelope components as a list
  
  
  return Rcpp::List::create(Rcpp::Named("GridIndex")=GIndex,
                            Rcpp::Named("thetabars")=G3,
                            Rcpp::Named("cbars")=cbars,
                            Rcpp::Named("logU")=logU,
                            Rcpp::Named("logrt")=logrt,
                            Rcpp::Named("loglt")=loglt,
                            Rcpp::Named("LLconst")=LLconst,
                            Rcpp::Named("logP")=logP(_,0),
                            Rcpp::Named("PLSD")=PLSD,
                            Rcpp::Named("a1")=a_1,
                            Rcpp::Named("E_draws")=E_draws
                              );
  
}





// [[Rcpp::export(".EnvelopeBuild_Ind_Normal_Gamma")]]

List EnvelopeBuild_Ind_Normal_Gamma(NumericVector bStar,NumericMatrix A,
                                    NumericVector y, 
                                    NumericMatrix x,
                                    NumericMatrix mu,
                                    NumericMatrix P,
                                    NumericVector alpha,
                                    NumericVector wt,
                                    std::string family="binomial",
                                    std::string link="logit",
                                    int Gridtype=2, 
                                    int n=1,
                                    int n_envopt=-1,
                                    bool sortgrid=false,
                                    bool use_opencl    = false,
                                    bool verbose       = false
){
  
  
  //  int progbar=0;
  
  int l1 = A.nrow(), k = A.ncol();
  arma::mat A2(A.begin(), l1, k, false);
  arma::colvec bStar_2(bStar.begin(), bStar.size(), false);
  
  
  NumericVector a_1(l1);
  arma::vec a_2(a_1.begin(), a_1.size(), false);
  
  NumericVector xx_1(3, 1.0);
  NumericVector xx_2=NumericVector::create(-1.0,0.0,1.0);
  NumericVector yy_1(2, 1.0);
  NumericVector yy_2=NumericVector::create(-0.5,0.5);
  NumericMatrix G1(3,l1);
  NumericMatrix Lint1(2,l1);
  arma::mat G1b(G1.begin(), 3, l1, false);
  arma::mat Lint(Lint1.begin(), 2, l1, false);
  
  arma::colvec xx_1b(xx_1.begin(), xx_1.size(), false);
  arma::colvec xx_2b(xx_2.begin(), xx_2.size(), false);
  arma::colvec yy_1b(yy_1.begin(), yy_1.size(), false);
  arma::colvec yy_2b(yy_2.begin(), yy_2.size(), false);
  List G2(a_1.size());
  List GIndex1(a_1.size());
  Rcpp::Function EnvelopeOpt("EnvelopeOpt");
  Rcpp::Function expGrid("expand.grid");
  Rcpp::Function asMat("as.matrix");
  Rcpp::Function EnvSort("EnvelopeSort");
  
  int i;  
  
  a_2=arma::diagvec(A2);
  arma::vec omega=(sqrt(2)-arma::exp(-1.20491-0.7321*sqrt(0.5+a_2)))/arma::sqrt(1+a_2);
  G1b=xx_1b*arma::trans(bStar_2)+xx_2b*arma::trans(omega);
  Lint=yy_1b*arma::trans(bStar_2)+yy_2b*arma::trans(omega);
  
  // Second row in G1b here is the posterior mode
  
  NumericVector gridindex(l1);
  
  if(Gridtype==2){
    gridindex=EnvelopeOpt(a_2,n);
  }
  
  NumericVector Temp1=G1( _, 0);
  double Temp2;
  
  // Should write a small note with logic behind types 1 and 2
  
  for(i=0;i<l1;i++){
    
    if(Gridtype==1){
      
      // For Gridtype==1, small 1+a[i]<=(2/sqrt(M_PI) yields grid over full line
      // Can check speed for simulation when Gridtype=1 vs. Gridtyp=2 or 3     
      
      if((1+a_2[i])<=(2/sqrt(M_PI))){ 
        Temp2=G1(1,i);
        G2[i]=NumericVector::create(Temp2);
        GIndex1[i]=NumericVector::create(4.0);
      }
      if((1+a_2[i])>(2/sqrt(M_PI))){
        Temp1=G1(_,i);
        G2[i]=NumericVector::create(Temp1(0),Temp1(1),Temp1(2));
        GIndex1[i]=NumericVector::create(1.0,2.0,3.0);
      }    
    }  
    if(Gridtype==2){
      if(gridindex[i]==1){
        Temp2=G1(1,i);
        G2[i]=NumericVector::create(Temp2);
        GIndex1[i]=NumericVector::create(4.0);
      }
      if(gridindex[i]==3){
        Temp1=G1(_,i);
        G2[i]=NumericVector::create(Temp1(0),Temp1(1),Temp1(2));
        GIndex1[i]=NumericVector::create(1.0,2.0,3.0);
      }
    }
    
    if(Gridtype==3){
      Temp1=G1(_,i);
      G2[i]=NumericVector::create(Temp1(0),Temp1(1),Temp1(2));
      GIndex1[i]=NumericVector::create(1.0,2.0,3.0);
    }
    
    if(Gridtype==4){
      Temp2=G1(1,i);
      G2[i]=NumericVector::create(Temp2);
      GIndex1[i]=NumericVector::create(4.0);
    }
    
    
    
  }
  
  NumericMatrix G3=asMat(expGrid(G2));
  NumericMatrix GIndex=asMat(expGrid(GIndex1));
  NumericMatrix G4(G3.ncol(),G3.nrow());
  int l2=GIndex.nrow();
  
  arma::mat G3b(G3.begin(), G3.nrow(), G3.ncol(), false);
  arma::mat G4b(G4.begin(), G4.nrow(), G4.ncol(), false);
  
  G4b=trans(G3b);
  
  NumericMatrix cbars(l2,l1);
  NumericMatrix cbars_slope(l2,l1);
  NumericMatrix Up(l2,l1);
  NumericMatrix Down(l2,l1);
  NumericMatrix logP(l2,2);
  NumericMatrix logU(l2,l1);
  NumericMatrix loglt(l2,l1);
  NumericMatrix logrt(l2,l1);
  NumericMatrix logct(l2,l1);
  
  NumericMatrix LLconst(l2,1);
  NumericVector NegLL(l2);    
  NumericVector NegLL_slope(l2);    
  NumericVector RSS_Out(l2);
  arma::mat cbars2(cbars.begin(), l2, l1, false); 
  arma::mat cbars3(cbars.begin(), l2, l1, false); 
  
  arma::mat cbars_slope2(cbars_slope.begin(), l2, l1, false); 
  arma::mat cbars_slope3(cbars_slope.begin(), l2, l1, false); 
  
  
  // Note: NegLL_2 only added to allow for QC printing of results 
  
  arma::colvec NegLL_2(NegLL.begin(), NegLL.size(), false);
  
  //    G4b.print("tangent points");
  
  //  Rcpp::Rcout << "Gridtype is :"  << Gridtype << std::endl;
  //  Rcpp::Rcout << "Number of Variables in model are :"  << l1 << std::endl;
  //  Rcpp::Rcout << "Number of points in Grid are :"  << l2 << std::endl;
  
  
  
  
  if(family=="gaussian" ){
    //Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    
    // Adjust the slope calculations to split into several terms:
    // (i) Terms from shifted "prior" that does not depend on the dispersion
    // (ii) Constant terms from the actual LL that do not depend on dispersion or beta
    // (iii) Term from the LL that depends on the dispersion but not beta
    // (iv) Term from the LL that depends on beta and the dispersion (scaled RSS)
    
    Rcpp::List eval_info = EnvelopeEval(G4, y, x, mu, P, alpha, wt,
                                        family, link, use_opencl, verbose);
    
    
    // Copy results into cbars/NegLL structures used downstream
    
    NegLL = eval_info["NegLL"];
    cbars2 = Rcpp::as<arma::mat>(eval_info["cbars"]);
    
    
//    NegLL=f2_gaussian(G4,y, x, mu, P, alpha, wt);  
    //Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
//    cbars2=f3_gaussian(G4,y, x,mu,P,alpha,wt);

Rcpp::List eval_info2 = EnvelopeEval(G4, y, x, mu,0*P, alpha, wt,
                                    family, link, use_opencl, verbose);


NegLL_slope   = eval_info2["NegLL"];
cbars_slope2  = Rcpp::as<arma::mat>(eval_info2["cbars"]);


if (verbose) {
  Rcpp::Rcout << "Finished assigning NegLL_slope and cbars_slope2" << std::endl;
}

//    NegLL_slope=f2_gaussian(G4,y, x, mu, 0*P, alpha, wt);  
//    cbars_slope2=f3_gaussian(G4,y, x,mu,0*P,alpha,wt);
    RSS_Out=RSS(y, x,G4,alpha,wt); // Note currenly includes the dispersion in the weight
    
  }
  
  
  //  Rcpp::Rcout << "Finished Log-posterior evaluations:" << std::endl;
  
  // Do a temporary correction here cbars3 should point to correct memory
  // See if this sets cbars
  
  cbars3=cbars2;
  cbars_slope3=cbars_slope2;

  if (verbose) {
    Rcpp::Rcout << "Entering Set_Grid_C2" << std::endl;
  }
  

  
  Set_Grid_C2(GIndex, cbars, Lint1,Down,Up,loglt,logrt,logct,logU,logP);

  if (verbose) {
    Rcpp::Rcout << "Entering Set_logP_C2" << std::endl;
  }
  
  
  setlogP_C2(logP,NegLL,cbars,G3,LLconst);
  
  
  if (verbose) {
    Rcpp::Rcout << "Computing PLSD" << std::endl;
  }
  
    
  NumericMatrix::Column logP2 = logP( _, 1);

    
  
  double  maxlogP=max(logP2);
  
  NumericVector PLSD=exp(logP2-maxlogP);
  
  double sumP=sum(PLSD);
  
  PLSD=PLSD/sumP;
  
  // Add sorting step back later after modifying EnvSort function
  // Should accomodate ready List
  
  //  if(sortgrid==true){
  //    Rcpp::List outlist=EnvSort(l1,l2,GIndex,G3,cbars,logU,logrt,loglt,logP,LLconst,PLSD,a_1);
  //    return(outlist);
  //  }
  
  
  if (verbose) {
    Rcpp::Rcout << "Finished Computing PLSD" << std::endl;
  }
  
  
  return Rcpp::List::create(Rcpp::Named("GridIndex")=GIndex,
                            Rcpp::Named("thetabars")=G3,
                            Rcpp::Named("cbars")=cbars,
                            Rcpp::Named("cbars_slope")=cbars_slope,
                            Rcpp::Named("NegLL")=NegLL,
                            Rcpp::Named("NegLL_slope")=NegLL_slope,
                            Rcpp::Named("Lint1")=Lint1,
                            Rcpp::Named("RSS_Out")=RSS_Out,
                            Rcpp::Named("logU")=logU,
                            Rcpp::Named("logrt")=logrt,
                            Rcpp::Named("loglt")=loglt,
                            Rcpp::Named("LLconst")=LLconst,
                            Rcpp::Named("logP")=logP(_,0),
                            Rcpp::Named("PLSD")=PLSD,
                            Rcpp::Named("a1")=a_1
  );
  
  
}





// [[Rcpp::export]]
Rcpp::List Inv_f3_precompute_disp(NumericMatrix cbars,
                                  NumericVector y,
                                  NumericMatrix x,
                                  NumericMatrix mu,
                                  NumericMatrix P,
                                  NumericVector alpha,
                                  NumericVector wt) {
  int n = x.nrow();
  int p = x.ncol();
  int m = cbars.ncol();
  
  arma::mat X(x.begin(), n, p, false);
  arma::mat Xt = X.t();
  arma::vec yv(y.begin(), n, false);
  arma::vec alphav(alpha.begin(), n, false);
  arma::vec xb = alphav - yv;
  
  arma::mat Pmat(P.begin(), p, p, false);
  Pmat = 0.5 * (Pmat + Pmat.t());
  
  arma::mat Mu(mu.begin(), p, 1, false);
  arma::mat Pmu = Pmat * Mu;
  
  arma::vec wv(wt.begin(), n, false);
  
  arma::vec base_B0 = Xt * (wv % xb);
  arma::mat base_A  = Xt * (X.each_col() % wv);
  
  arma::mat C(cbars.begin(), p, m, false);
  
  return Rcpp::List::create(
    Rcpp::Named("Pmat")    = Pmat,
    Rcpp::Named("Pmu")     = Pmu,
    Rcpp::Named("base_B0") = base_B0,
    Rcpp::Named("base_A")  = base_A,
    Rcpp::Named("C")       = C
  );
}



// // [[Rcpp::export]]
// 
// double rss_face_at_disp(const Rcpp::List cache,
//                                       double dispersion,
//                                       const Rcpp::NumericVector cbars_j) {
//   // Extract cached matrices
//   arma::mat Pmat    = cache["Pmat"];
//   arma::mat Pmu     = cache["Pmu"];
//   arma::vec base_B0 = cache["base_B0"];
//   arma::mat base_A  = cache["base_A"];
//   
//   // Scale terms by dispersion
//   arma::vec B0 = base_B0 / dispersion + Pmu;
//   arma::mat A  = Pmat + base_A / dispersion;
//   A = 0.5 * (A + A.t());                // ensure symmetry
//   
//   arma::mat R = arma::chol(A);          // Cholesky
//   
//   // Wrap cbars_j as Armadillo vector
// //  arma::vec c_j(cbars_j.begin(), cbars_j.size(), false);
//   arma::vec c_j = Rcpp::as<arma::vec>(cbars_j);
//   
//   
//   // Solve A^{-1}(-c_j + B0)
//   arma::vec b    = -c_j + B0;
//   arma::vec ytmp = arma::solve(arma::trimatl(R.t()), b);
//   arma::vec sol  = arma::solve(arma::trimatu(R), ytmp);
//   
//   // RSS is squared norm of the solution
//   return arma::dot(sol, sol);
// }

// [[Rcpp::export("rss_face_at_disp")]]

double rss_face_at_disp(double dispersion,
                               Rcpp::List cache,
                               Rcpp::NumericVector cbars_j,
                               Rcpp::NumericVector y,
                               Rcpp::NumericMatrix x,
                               Rcpp::NumericVector alpha,
                               Rcpp::NumericVector wt) {
  // Build 1×l1 matrix, then transpose to l1×1 for Inv_f3_with_disp
  int l1 = cbars_j.size();
  Rcpp::NumericMatrix cbars_small(1, l1);
  for (int k = 0; k < l1; ++k) cbars_small(0, k) = cbars_j[k];
  
  arma::mat theta_row = Inv_f3_with_disp(cache, dispersion, Rcpp::transpose(cbars_small));
  arma::vec beta = theta_row.t(); // 1×l1 -> l1×1
  
  arma::vec y2(y.begin(), y.size(), false);
  arma::vec a2(alpha.begin(), alpha.size(), false);
  arma::mat X(x.begin(), x.nrow(), x.ncol(), false);
  arma::vec w(wt.begin(), wt.size(), false);
  
  arma::vec resid = (y2 - a2 - X * beta) % arma::sqrt(w);
  return arma::as_scalar(resid.t() * resid);
}


// [[Rcpp::export]]
double UB2(double dispersion,
           Rcpp::List cache,
           Rcpp::NumericVector cbars_j,
           Rcpp::NumericVector y,
           Rcpp::NumericMatrix x,
           Rcpp::NumericVector alpha,
           Rcpp::NumericVector wt,
           double rss_min_global) {
  
  // Call the existing RSS function
  double rss_val = rss_face_at_disp(dispersion, cache, cbars_j, y, x, alpha, wt);
  
  // Compute UB2
  double UB2_val = (1.0 / dispersion) * (rss_val - rss_min_global);
  
  return UB2_val;
}



// Utility: safe max for NumericVector
static inline double max_vec(const NumericVector& v) {
  double m = R_NegInf;
  for (int i = 0; i < v.size(); ++i) if (v[i] > m) m = v[i];
  return m;
}


NumericVector EnvBuildLinBound_cpp(NumericMatrix thetabars,
                                   NumericMatrix cbars,
                                   NumericVector y,
                                   NumericMatrix x,
                                   NumericMatrix P,
                                   NumericVector alpha,
                                   double dispstar) {
  // Convert to Armadillo
  arma::mat thetabarsA = as<arma::mat>(thetabars);
  arma::mat cbarsA     = as<arma::mat>(cbars);
  arma::vec yA         = as<arma::vec>(y);
  arma::mat xA         = as<arma::mat>(x);
  arma::mat PA         = as<arma::mat>(P);
  arma::vec alphaA     = as<arma::vec>(alpha);
  
  int gs = cbarsA.n_rows;
  
  arma::mat XtX   = xA.t() * xA;
  arma::vec rhs   = xA.t() * (yA - alphaA);
  arma::mat M     = XtX + dispstar * PA;
  arma::mat Minv  = arma::inv(M);           // match R's solve(M)
  arma::mat H1    = -Minv * PA * Minv;
  
  arma::mat V = -thetabarsA * PA + cbarsA;                 // gs x p
  arma::mat Minv_cbars = cbarsA * Minv.t();                // gs x p
  arma::vec term1 = arma::sum(V % Minv_cbars, 1);
  
  arma::mat rhs_mat = arma::repmat(rhs, 1, gs);            // p x gs
  arma::mat H1_rhs  = (H1 * (rhs_mat + dispstar * cbarsA.t())).t(); // gs x p
  arma::vec term2 = arma::sum(V % H1_rhs, 1);
  
  arma::vec result = term1 + term2;
  
  // Return explicitly as NumericVector
  NumericVector out(gs);
  std::copy(result.begin(), result.end(), out.begin());
  return out;
}


NumericVector thetabar_const_cpp(NumericMatrix P,
                                 NumericMatrix cbars,
                                 NumericMatrix thetabars) {
  arma::mat PA         = as<arma::mat>(P);
  arma::mat cbarsA     = as<arma::mat>(cbars);
  arma::mat thetabarsA = as<arma::mat>(thetabars);
  
  int gs = cbarsA.n_rows;
  arma::vec thetaconst(gs);
  
  for (int j = 0; j < gs; ++j) {
    arma::vec theta_temp = thetabarsA.row(j).t();
    arma::vec cbars_temp = cbarsA.row(j).t();
    thetaconst[j] = -0.5 * arma::as_scalar(theta_temp.t() * PA * theta_temp)
      + arma::as_scalar(cbars_temp.t() * theta_temp);
  }
  
  NumericVector out(gs);
  std::copy(thetaconst.begin(), thetaconst.end(), out.begin());
  return out;
}



// [[Rcpp::export]]
List EnvelopeDispersionBuild_cpp(
    List Env,
    double Shape,
    double Rate,
    NumericMatrix P,
    NumericVector y,
    NumericMatrix x,
    NumericVector alpha,
    int n_obs,
    double RSS_post,
    double RSS_ML,
    NumericMatrix mu,         // ← new
    NumericVector wt,         // ← new
    double max_disp_perc = 0.99,
    Nullable<double> disp_lower = R_NilValue,
    Nullable<double> disp_upper = R_NilValue,
    bool verbose = false
)
  {
  // Step 1: Posterior Gamma parameters (precision prior)
  double shape2 = Shape + static_cast<double>(n_obs) / 2.0;
  double rate3  = Rate  + RSS_post / 2.0;
  
  // Step 2: Dispersion bounds (on sigma^2)
  double low, upp;
  if (disp_lower.isNull() || disp_upper.isNull()) {
    // Call R's qgamma for tail quantiles, then invert to get sigma^2 bounds
    Function qgamma("qgamma");
    NumericVector q_low = qgamma(
      Named("p")     = max_disp_perc,
      Named("shape") = shape2,
      Named("rate")  = rate3
    );
    NumericVector q_upp = qgamma(
      Named("p")     = 1.0 - max_disp_perc,
      Named("shape") = shape2,
      Named("rate")  = rate3
    );
    low = 1.0 / q_low[0];
    upp = 1.0 / q_upp[0];
  } else {
    low = as<double>(disp_lower);
    upp = as<double>(disp_upper);
    if (!R_finite(low) || !R_finite(upp))
      stop("disp_lower/disp_upper must be finite.");
    if (low <= 0.0 || upp <= 0.0)
      stop("disp_lower/disp_upper must be positive.");
    if (upp <= low)
      stop("disp_upper must be strictly greater than disp_lower.");
  }
  
  // Step 3: Extract envelope faces
  NumericMatrix cbars     = Env["cbars"];      // gs x l1
  NumericMatrix thetabars = Env["thetabars"];  // gs x l1 (grid of tangencies)
  NumericVector logP1     = Env["logP"];       // length gs
  int gs = cbars.nrow();
  int l1 = cbars.ncol();
  
  /// Step 3B: Precompute elements for finding inverse function for cbars
  
  Rcpp::List cache = Inv_f3_precompute_disp(cbars, y, x, mu, P, alpha, wt);
  
  // Step 3C: Minimize RSS over dispersion for each face (optional diagnostics / UB2 prep)
  // Strategy A (pure C++): call a Brent/golden-section minimizer using rss_face_at_disp()
  // Strategy B (R-side): call optim("Brent") on [low, upp] — easier to prototype
  
  // Step 3C: Minimize RSS over dispersion for each face
  Rcpp::Function optim("optim");
  Rcpp::Function rss_fn("rss_face_at_disp");
  

  NumericVector rss_min(gs), disp_min(gs);
  for (int j = 0; j < gs; ++j) {
    // Extract row j of cbars into a NumericVector
    NumericVector cbars_j(l1);
    for (int k = 0; k < l1; ++k) cbars_j[k] = cbars(j, k);
    
    double par0 = 0.5 * (low + upp);
    
    // Call R's optim, passing the exported C++ function object
    Rcpp::List res = optim(
      Rcpp::Named("par")     = par0,
      Rcpp::Named("fn")      = rss_fn,
      Rcpp::Named("method")  = "Brent",
      Rcpp::Named("lower")   = low,
      Rcpp::Named("upper")   = upp,
      Rcpp::Named("cache")   = cache,
      Rcpp::Named("cbars_j") = cbars_j,
      Rcpp::Named("y")       = y,
      Rcpp::Named("x")       = x,
      Rcpp::Named("alpha")   = alpha,
      Rcpp::Named("wt")      = wt
    );    
    disp_min[j] = res["par"];
    rss_min[j]  = res["value"];
    
  }
  
  
  // Optionally: keep the best across faces
  double rss_min_global = R_PosInf;
  double disp_min_global = NA_REAL;
  int j_best = -1;
  for (int j = 0; j < gs; ++j) {
    if (rss_min[j] < rss_min_global) {
      rss_min_global = rss_min[j];
      disp_min_global = disp_min[j];
      j_best = j;
    }
  }  
  

  //////////////////////////////////////
  
  // Assume UB2 has been exported as shown earlier
  Rcpp::Function ub2_fn("UB2");
  

  // --- UB2 minimization loop ---
  NumericVector ub2_min(gs), disp_min_ub2(gs);
  
  for (int j = 0; j < gs; ++j) {
    // Extract row j of cbars into a NumericVector
    NumericVector cbars_j(l1);
    for (int k = 0; k < l1; ++k) cbars_j[k] = cbars(j, k);
    
    double par0 = 0.5 * (low + upp);
    
    // Call R's optim, passing the UB2 function object
    Rcpp::List res = optim(
      Rcpp::Named("par")     = par0,
      Rcpp::Named("fn")      = ub2_fn,
      Rcpp::Named("method")  = "Brent",
      Rcpp::Named("lower")   = low,
      Rcpp::Named("upper")   = upp,
      Rcpp::Named("cache")   = cache,
      Rcpp::Named("cbars_j") = cbars_j,
      Rcpp::Named("y")       = y,
      Rcpp::Named("x")       = x,
      Rcpp::Named("alpha")   = alpha,
      Rcpp::Named("wt")      = wt,
      Rcpp::Named("rss_min_global") = rss_min_global
    );
    
    disp_min_ub2[j] = res["par"];
    ub2_min[j]      = res["value"];
    
  }
  
  // Find global UB2 minimum
  double ub2_min_global = R_PosInf;
  double disp_min_global_ub2 = NA_REAL;
  int j_best_ub2 = -1;
  for (int j = 0; j < gs; ++j) {
    if (ub2_min[j] < ub2_min_global) {
      ub2_min_global = ub2_min[j];
      disp_min_global_ub2 = disp_min_ub2[j];
      j_best_ub2 = j;
    }
  }
  

  
  
  
  
  // Step 4: Base face constants via R helper (keep R version for now)
//  Function thetabar_const_R("thetabar_const");
//  NumericVector thetabar_const_base =
//    thetabar_const_R(P, cbars, thetabars);   // length gs
  
  NumericVector thetabar_const_base =
    thetabar_const_cpp(P, cbars, thetabars);
  
    
  // Step 5: initial anchor (posterior mean; optional)
  // Note: consistency with external rate2: here rate3=Rate+RSS_post/2 as in R function
  double dispstar = rate3 / (shape2 - 1.0);
  

  // Step 6: Face slopes at dispstar via R helper
//  Function EnvBuildLinBound_R("EnvBuildLinBound");
//  NumericVector New_LL_Slope =
//    EnvBuildLinBound_R(thetabars, cbars, y, x, P, alpha, dispstar); // length gs
  

  NumericVector New_LL_Slope =
    EnvBuildLinBound_cpp(thetabars, cbars, y, x, P, alpha, dispstar);
  

  // Step 7: Linear extrapolation of face constants to bounds
  NumericVector thetabar_const_upp_apprx(gs), thetabar_const_low_apprx(gs);
  for (int j = 0; j < gs; ++j) {
    thetabar_const_upp_apprx[j] = thetabar_const_base[j] + (upp - dispstar) * New_LL_Slope[j];
    thetabar_const_low_apprx[j] = thetabar_const_base[j] + (low - dispstar) * New_LL_Slope[j];
  }
  
  // Step 8: Global upper line geometry (match original mean-slope correction)
  double max_low = max_vec(thetabar_const_low_apprx);
  double max_upp = max_vec(thetabar_const_upp_apprx);
  
  // No-op in original; keep for parity via mean slope correction
  double m_New_LL_Slope = Rcpp::mean(New_LL_Slope);
  double max_low_mean   = max_upp - m_New_LL_Slope * (upp - low);
  max_low = max_low_mean;
  
  double new_slope = (max_upp - max_low) / (upp - low);
  double new_int   = max_low - new_slope * low;
  
  // Step 9a: Dispersion anchor (exactly as in original: b1/(-c1))
  double b1 = (upp - low);
  double c1 = -std::log(upp / low);
  dispstar  = b1 / (-c1);  // equivalently (upp - low)/log(upp/low)
  

  // Step 9: Mixture weights per face (match original)
  NumericVector New_logP2(gs);
  NumericVector prob_factor(gs);
  NumericVector prob_factor2(gs);
  for (int j = 0; j < gs; ++j) {
    // cbars_temp is row j (length l1)
    double norm2 = 0.0;
    for (int k = 0; k < l1; ++k) {
      double cjk = cbars(j, k);
      norm2 += cjk * cjk;
    }
    New_logP2[j] = logP1[j] + 0.5 * norm2;
    
    double pf_upp = thetabar_const_upp_apprx[j] - max_upp;
    double pf_low = thetabar_const_low_apprx[j] - max_low;
    prob_factor[j] = (pf_upp > pf_low ? pf_upp : pf_low);
    prob_factor2[j] =prob_factor[j]-ub2_min[j];

    
  }
  

  // Log-space prob factors (kept separate for UB_list, as in R)
  NumericVector lg_prob_factor = clone(prob_factor);
  NumericVector lg_prob_factor2 = clone(prob_factor2);
  
  
  
  
  
  // Normalize weights (PLSD)
  NumericVector prob_factor_exp(gs);
  NumericVector prob_factor_exp2(gs);
  for (int j = 0; j < gs; ++j){
    prob_factor_exp[j] = std::exp(New_logP2[j] + prob_factor[j]);
    prob_factor_exp2[j] = std::exp(New_logP2[j] + prob_factor2[j]);
    
  }
  double sumP = std::accumulate(prob_factor_exp.begin(), prob_factor_exp.end(), 0.0);
  double sumP2 = std::accumulate(prob_factor_exp2.begin(), prob_factor_exp2.end(), 0.0);
  for (int j = 0; j < gs; ++j){
    prob_factor_exp[j] /= sumP;
    prob_factor_exp2[j] /= sumP2;
    
  }   
  // Step 10: Envelope constants for dispersion and gamma tilt
  double lm_log2 = new_slope * dispstar;
  double lm_log1 = new_int + new_slope * dispstar - new_slope * std::log(dispstar);
  double shape3  = shape2 - lm_log2;
  
  // Step 11: Package outputs
//  Env["PLSD"] = prob_factor_exp;
  Env["PLSD"] = prob_factor_exp2;
  
  List gamma_list = List::create(
    Named("shape3")     = shape3,
//    Named("rate2")      = Rate + RSS_ML / 2.0,  // matches original definition
    Named("rate2")      = Rate + rss_min_global / 2.0,  // matches original definition
    Named("disp_upper") = upp,
    Named("disp_lower") = low
  );
  
  List UB_list = List::create(
    Named("RSS_ML")         = RSS_ML,               // not RSS_post
    Named("RSS_Min")        = rss_min_global,       // Minimum across faces
    Named("max_New_LL_UB")  = max_upp,
    Named("max_LL_log_disp")= lm_log1 + lm_log2 * std::log(upp),
    Named("lm_log1")        = lm_log1,
    Named("lm_log2")        = lm_log2,
    Named("lg_prob_factor") = lg_prob_factor,
    Named("lmc1")           = new_int,
    Named("lmc2")           = new_slope,
    Named("UB2min")           = ub2_min
  
  );
  
  List diagnostics = List::create(
    Named("dispstar")     = dispstar,
    Named("New_LL_Slope") = New_LL_Slope,
    Named("shape2")       = shape2,
    Named("rate3")        = rate3,
    Named("shape3")       = shape3,
    Named("max_low")      = max_low,
    Named("max_upp")      = max_upp,
    Named("new_slope")    = new_slope,
    Named("new_int")      = new_int,
    Named("prob_factor")  = prob_factor_exp,
    Named("UB2min")           = ub2_min
//  Named("prob_factor2")  = prob_factor_exp2
  );
  
  if (verbose) {
    Rcout << "EnvelopeDispersionBuild diagnostics:\n";
    Rcout << "  dispstar      = " << dispstar << "\n";
    Rcout << "  new_slope     = " << new_slope << "\n";
    Rcout << "  new_int       = " << new_int << "\n";
    Rcout << "  lm_log1       = " << lm_log1 << "\n";
    Rcout << "  lm_log2       = " << lm_log2 << "\n";
    Rcout << "  shape3        = " << shape3 << "\n";
    Rcout << "  max_low       = " << max_low << "\n";
    Rcout << "  max_upp       = " << max_upp << "\n";
    Rcout << "  RSS_ML       = " << RSS_ML << "\n";
    Rcout << "  RSS_Min       = " << rss_min_global << "\n";
    Rcout << "  disp_lower       = " << low << "\n";
    Rcout << "  disp_upper       = " << upp << "\n";
  }
  
  return List::create(
    Named("Env_out")    = Env,
    Named("gamma_list") = gamma_list,
    Named("UB_list")    = UB_list,
    Named("diagnostics")= diagnostics
  );
}





