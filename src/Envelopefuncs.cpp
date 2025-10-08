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
  
  
  if (verbose) {
    Rcpp::Rcout << ">>> EnvelopeBuild_c called with:\n"
                << "    Gridtype   = " << Gridtype << "\n"
                << "    n          = " << n << "\n"
                << "    use_opencl = " << use_opencl << "\n"
                << "    sortgrid   = " << sortgrid << "\n";
  }
  
  
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
  
  int core_CNT=get_opencl_core_count();
  
  if (verbose) {
    Rcpp::Rcout << "[INFO] OpenCL core count = " << core_CNT << "\n";
  }
  
  // Second row in G1b here is the posterior mode
  
  
  
  NumericVector gridindex(l1);
  
  if(Gridtype==2){
    // Temporarily scale n to encourage richer envelopes when GPU parallelism is available
    int scaled_n = std::max(1, n * core_CNT);
    
//    if (verbose) {
//      Rcpp::Rcout << "[INFO] Scaling n from " << n << " to " << scaled_n
//                  << " for envelope optimization.\n";
//    }
    
//    gridindex = EnvelopeOpt(a_2, scaled_n,core_CNT);

  if(use_opencl==0 ){
     gridindex=EnvelopeOpt(a_2,n_envopt,1);}
    else{
    gridindex = EnvelopeOpt(a_2, n_envopt,core_CNT);
    }
    
  }
  
  NumericVector Temp1=G1( _, 0);
  double Temp2;
  
  
  if (verbose) {
    
    Rcpp::Rcout << "Entering Grid Loop: "
                << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                << "\n";
  }
  
  // Should write a small note with logic behind types 1 and 2
  
  /*
   GRIDTYPE LOGIC (Nygren & Nygren 2006)
   
   Let a_i = posterior precision diagonal for dimension i (i.e. A2(i,i)).
   Let ω_i = width parameter computed below.
   
   1) Gridtype 1: static threshold test
   • If sqrt(1 + a_i) ≤ 2/√π  ≈ 1.128379 ⇒
   the theoretical upper bound on expected candidates per draw
   for a full normal envelope is ≤ 1, so building more points
   doesn’t reduce rejections.
   ⇒ use a single-point envelope at the mode:
   G2[i] = {θ⋆_i},  GIndex1[i] = {4}
   
   • Else
   ⇒ build a three-point envelope at
   {θ⋆_i − ω_i, θ⋆_i, θ⋆_i + ω_i}, GIndex1[i] = {1,2,3}
   
   2) Gridtype 2: dynamic envelope via `EnvelopeOpt(a, n)`
   • Rather than a fixed-bound test, we solve (per dimension) the
   cost-tradeoff between:
   – T_build(g_i) ∝ g_i  (linear grid-construction cost)
   – T_sample(n, acc_i(g_i)) ∝ n / acc_i(g_i)
   where acc_i(g_i) ≈ acceptance rate given g_i grid points
   • `EnvelopeOpt` returns gridindex[i] ∈ {1, 3} by minimizing
   T_build + T_sample under approximations from subgradient theory.
   • This lets us invest in more envelope points up-front when
   n is large, because reduced rejections during sampling
   more than pay for the extra build cost.
   
   3) Gridtype 3: always three-point grid (regardless of a_i or n)
   4) Gridtype 4: always single-point grid (mode only)
   */
  
  
  double E_draws=1.0L;
  
//  Rcpp::Rcout << "[DEBUG] half=" << i 
//              << " E_draws=" << E_draws << "\n";
  
  
    
  for(i=0;i<l1;i++){
    
    if(Gridtype==1){
      
      // For Gridtype==1, small 1+a[i]<=(2/sqrt(M_PI) yields grid over full line
      // Can check speed for simulation when Gridtype=1 vs. Gridtyp=2 or 3     
      
      if(sqrt(1+a_2[i])<=(2/sqrt(M_PI))){ 
        Temp2=G1(1,i);
        G2[i]=NumericVector::create(Temp2);
        GIndex1[i]=NumericVector::create(4.0);
        E_draws=E_draws*sqrt(1+a_2[i]);
        
      }
      if(sqrt(1+a_2[i])>(2/sqrt(M_PI))){
        Temp1=G1(_,i);
        G2[i]=NumericVector::create(Temp1(0),Temp1(1),Temp1(2));
        GIndex1[i]=NumericVector::create(1.0,2.0,3.0);
        E_draws=E_draws*(2/sqrt(M_PI));
                               
      }    
    }  
    if(Gridtype==2){
      if(gridindex[i]==1){
        Temp2=G1(1,i);
        G2[i]=NumericVector::create(Temp2);
        GIndex1[i]=NumericVector::create(4.0);
        E_draws=E_draws*sqrt(1+a_2[i]);
      }
      if(gridindex[i]==3){
        Temp1=G1(_,i);
        G2[i]=NumericVector::create(Temp1(0),Temp1(1),Temp1(2));
        GIndex1[i]=NumericVector::create(1.0,2.0,3.0);
        E_draws=E_draws*(2/sqrt(M_PI));
      }
    }
    
    if(Gridtype==3){
      Temp1=G1(_,i);
      G2[i]=NumericVector::create(Temp1(0),Temp1(1),Temp1(2));
      GIndex1[i]=NumericVector::create(1.0,2.0,3.0);
      E_draws=E_draws*(2/sqrt(M_PI));
    }
    
    if(Gridtype==4){
      Temp2=G1(1,i);
      G2[i]=NumericVector::create(Temp2);
      GIndex1[i]=NumericVector::create(4.0);
      E_draws=E_draws*sqrt(1+a_2[i]);
      
    }
  
  
  // after updating E_draws, print it
//  Rcpp::Rcout << "[DEBUG] i=" << i 
//              << " E_draws=" << E_draws << "\n";
  
    
  }
  
  
  
  
  NumericMatrix G3=asMat(expGrid(G2));
  NumericMatrix GIndex=asMat(expGrid(GIndex1));
  NumericMatrix G4(G3.ncol(),G3.nrow());
  int l2=GIndex.nrow();
  
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
  
  
  //////////////////  Pilot Phase

  
  int m1_total = G4.ncol();
  int y_obs = y.size();  // or however you extract observation count
  
  double obs_factor = std::min(1.0, 100.0 / (double)y_obs);  // dampens pilot size for large n
//  double frac_A = 0.01 * obs_factor;  // base fraction scaled down
//  double frac_B = 0.02 * obs_factor;
  
  double frac_A = 0.01 * obs_factor;  // base fraction scaled down
  double frac_B = 0.02 * obs_factor;
  
  
    int m1_pilot_A = std::max(100, (int)(m1_total * frac_A));
  int m1_pilot_B = std::max(m1_pilot_A + 100, (int)(m1_total * frac_B));
  

  double est_time_sec = NA_REAL;
  
  if (m1_total > 50000 && use_opencl) {
    auto slice_grid = [&](int m1_pilot) {
      Rcpp::NumericMatrix G4_pilot(G4.nrow(), m1_pilot);
      for (int j = 0; j < m1_pilot; ++j)
        for (int i = 0; i < G4.nrow(); ++i)
          G4_pilot(i, j) = G4(i, j);
      return G4_pilot;
    };
    
    // Run Pilot A (warm-up, discard timing)
    auto G4_pilot_A = slice_grid(m1_pilot_A);
    Rcpp::List warmup_A = f2_f3_opencl(family, link, G4_pilot_A, y, x, mu, P, alpha, wt, 0);


    // Run Pilot A again (timed)
    auto t0A = std::chrono::high_resolution_clock::now();
    Rcpp::List pilot_A = f2_f3_opencl(family, link, G4_pilot_A, y, x, mu, P, alpha, wt, 0);
    auto t1A = std::chrono::high_resolution_clock::now();
    double time_A = std::chrono::duration<double>(t1A - t0A).count();
    
    // Run Pilot B (timed)
    auto G4_pilot_B = slice_grid(m1_pilot_B);
    

    auto t0B = std::chrono::high_resolution_clock::now();
    Rcpp::List pilot_B = f2_f3_opencl(family, link, G4_pilot_B, y, x, mu, P, alpha, wt, 0);
    auto t1B = std::chrono::high_resolution_clock::now();
    double time_B = std::chrono::duration<double>(t1B - t0B).count();
    
    // Estimate per-grid cost and fixed cost
    double per_grid_cost = (time_B - time_A) / std::max(1.0, (double)(m1_pilot_B - m1_pilot_A));
    double fixed_cost = time_A - m1_pilot_A * per_grid_cost;
    
    if (per_grid_cost <= 0.0) {
      Rcpp::Rcout << "[WARNING] Negative per-grid cost — falling back to Pilot B average.\n";
      per_grid_cost = time_B / m1_pilot_B;
      fixed_cost = 0.0;
    }
    
    
    est_time_sec = fixed_cost + per_grid_cost * m1_total;
    
    
    long total = static_cast<long>(std::round(est_time_sec));
    long h = total / 3600;
    long m = (total % 3600) / 60;
    long s = total % 60;
    
    Rcpp::Rcout << "Estimated full f2_f3 evaluation time: "
                << est_time_sec << " seconds ("
                << h << "h " << m << "m " << s << "s)"
                << " (fixed=" << fixed_cost
                << ", per-grid=" << per_grid_cost << ").\n"
                << "Note: estimate is approximate and may vary with system load.\n";
    
  }
  
  
  
  if (est_time_sec > 300.0) {
    Rcpp::Rcout << "\nEstimated run time exceeds 5 minutes (" << est_time_sec << " seconds).\n";
    Rcpp::Rcout << "Do you want to continue? [y/N]: ";
    Rcpp::Function readline("readline");
    std::string response = Rcpp::as<std::string>(readline("Do you want to continue? [y/N]: "));
    
    if (response != "y" && response != "Y") {
      Rcpp::Rcout << "Aborting full run per user choice.\n";
      Rcpp::stop("User aborted run due to estimated time exceeding threshold.");
    } else {
      Rcpp::Rcout << "Proceeding with full run...\n";
    }
  }
  
  
  
  //    G4b.print("tangent points");
  
  //  Rcpp::Rcout << "Gridtype is :"  << Gridtype << std::endl;
  //  Rcpp::Rcout << "Number of Variables in model are :"  << l1 << std::endl;
  //  Rcpp::Rcout << "Number of points in Grid are :"  << l2 << std::endl;
  
  if( family=="binomial" && link=="logit"){
    
    if (verbose) {
      
      Rcpp::Rcout << "Initiating NegLL Calculations: "
                  << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << "\n";
    }
    
    if(use_opencl==0 ){
      NegLL=f2_binomial_logit(G4,y, x, mu, P, alpha, wt,progbar);  
      
      
      
      if (verbose) {
        Rcpp::Rcout << "Initiating Gradient Evaluations: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      
      cbars2=f3_binomial_logit(G4,y, x,mu,P,alpha,wt,progbar);
      
      
    }
    
    else{
      
      
      
      if (verbose) {
        Rcpp::Rcout << "Initiating f2_f3_opencl: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      
      
      Rcpp::List prepGrad_v3= f2_f3_opencl(
        family,
        link,
        G4,          // NumericMatrix b
        y,           // NumericVector y
        x,           // NumericMatrix x
        mu,          // NumericMatrix mu 
        P,           // NumericMatrix P
        alpha,       // NumericVector alpha
        wt,          // NumericVector wt
        progbar     // int progbar
      );
      
      
      NegLL = prepGrad_v3["qf"];
      cbars2 = Rcpp::as<arma::mat>(prepGrad_v3["grad"]);
      
      
      
    }
    
    
  }
  if(family=="binomial"  && link=="probit"){
    
    if(use_opencl==0 ){
      
      if (verbose) {
        
        Rcpp::Rcout << "Initiating NegLL Calculations: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      
      NegLL=f2_binomial_probit(G4,y, x, mu, P, alpha, wt,progbar);  
      
      
      
      if (verbose) {
        Rcpp::Rcout << "Initiating Gradient Evaluations: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      cbars2=f3_binomial_probit(G4,y, x,mu,P,alpha,wt,progbar);
      
      
    }
    
    else{
      
      
      if (verbose) {
        Rcpp::Rcout << "Initiating f2_f3_opencl: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      
      
      Rcpp::List prepGrad_v3= f2_f3_opencl(
        family,
        link,
        G4,          // NumericMatrix b
        y,           // NumericVector y
        x,           // NumericMatrix x
        mu,          // NumericMatrix mu 
        P,           // NumericMatrix P
        alpha,       // NumericVector alpha
        wt,          // NumericVector wt
        progbar     // int progbar
      );
      
      
      NegLL = prepGrad_v3["qf"];
      cbars2 = Rcpp::as<arma::mat>(prepGrad_v3["grad"]);
      
      
      
    }    
    
    
  }
  if(family=="binomial"   && link=="cloglog"){
    
    if(use_opencl==0 ){
      
      
      if (verbose) {
        
        Rcpp::Rcout << "Initiating NegLL Calculations: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      NegLL=f2_binomial_cloglog(G4,y, x, mu, P, alpha, wt,progbar);  
      
      
      if (verbose) {
        Rcpp::Rcout << "Initiating Gradient Evaluations: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      cbars2=f3_binomial_cloglog(G4,y, x,mu,P,alpha,wt,progbar);
      
      
      
    }
    
    else{
      
      
      
      if (verbose) {
        Rcpp::Rcout << "Initiating f2_f3_opencl: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      
      
      Rcpp::List prepGrad_v3= f2_f3_opencl(
        family,
        link,
        G4,          // NumericMatrix b
        y,           // NumericVector y
        x,           // NumericMatrix x
        mu,          // NumericMatrix mu 
        P,           // NumericMatrix P
        alpha,       // NumericVector alpha
        wt,          // NumericVector wt
        progbar     // int progbar
      );
      
      
      
      NegLL = prepGrad_v3["qf"];
      cbars2 = Rcpp::as<arma::mat>(prepGrad_v3["grad"]);
      
      
      
    }
    
  }
  
  if(family=="quasibinomial"  && link=="logit"){
    
    
    if (verbose) {
      
      Rcpp::Rcout << "Initiating NegLL Calculations: "
                  << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << "\n";
    }
    
    if(use_opencl==0 ){
      NegLL=f2_binomial_logit(G4,y, x, mu, P, alpha, wt,progbar);  
      
      
      
      if (verbose) {
        Rcpp::Rcout << "Initiating Gradient Evaluations: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      
      cbars2=f3_binomial_logit(G4,y, x,mu,P,alpha,wt,progbar);
      
      
    }
    
    else{
      
      
      
      if (verbose) {
        Rcpp::Rcout << "Initiating f2_f3_opencl: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      
      
      Rcpp::List prepGrad_v3= f2_f3_opencl(
        family,
        link,
        G4,          // NumericMatrix b
        y,           // NumericVector y
        x,           // NumericMatrix x
        mu,          // NumericMatrix mu 
        P,           // NumericMatrix P
        alpha,       // NumericVector alpha
        wt,          // NumericVector wt
        progbar     // int progbar
      );
      
      
      NegLL = prepGrad_v3["qf"];
      cbars2 = Rcpp::as<arma::mat>(prepGrad_v3["grad"]);
      
      
      
    }
    
    
    
  }
  if(family=="quasibinomial" && link=="probit"){
    
    if(use_opencl==0 ){
      
      if (verbose) {
        
        Rcpp::Rcout << "Initiating NegLL Calculations: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      
      NegLL=f2_binomial_probit(G4,y, x, mu, P, alpha, wt,progbar);  
      
      
      
      if (verbose) {
        Rcpp::Rcout << "Initiating Gradient Evaluations: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      cbars2=f3_binomial_probit(G4,y, x,mu,P,alpha,wt,progbar);
      
      
    }
    
    else{
      
      
      if (verbose) {
        Rcpp::Rcout << "Initiating f2_f3_opencl: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      
      
      Rcpp::List prepGrad_v3= f2_f3_opencl(
        family,
        link,
        G4,          // NumericMatrix b
        y,           // NumericVector y
        x,           // NumericMatrix x
        mu,          // NumericMatrix mu 
        P,           // NumericMatrix P
        alpha,       // NumericVector alpha
        wt,          // NumericVector wt
        progbar     // int progbar
      );
      
      
      NegLL = prepGrad_v3["qf"];
      cbars2 = Rcpp::as<arma::mat>(prepGrad_v3["grad"]);
      
      
      
    }    
    
    
  }
  
  if(family=="poisson" ){
    
    if(use_opencl==0 ){
      
      if (verbose) {
        
        Rcpp::Rcout << "Initiating NegLL Calculations: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      NegLL=f2_poisson(G4,y, x, mu, P, alpha, wt,progbar);  
      
      
      if (verbose) {
        Rcpp::Rcout << "Initiating Gradient Evaluations: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      cbars2=f3_poisson(G4,y, x,mu,P,alpha,wt,progbar);
      
      
      
    }
    
    else{
      
      
      if (verbose) {
        Rcpp::Rcout << "Initiating f2_f3_opencl: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      
      
      Rcpp::List prepGrad_v3= f2_f3_opencl(
        family,
        link,
        G4,          // NumericMatrix b
        y,           // NumericVector y
        x,           // NumericMatrix x
        mu,          // NumericMatrix mu 
        P,           // NumericMatrix P
        alpha,       // NumericVector alpha
        wt,          // NumericVector wt
        progbar     // int progbar
      );
      
      if (verbose) {
        
        Rcpp::Rcout << "Exiting f2_f3_opencl: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      
      
      NegLL = prepGrad_v3["qf"];
      cbars2 = Rcpp::as<arma::mat>(prepGrad_v3["grad"]);
      
    }
    
    
  }
  
  if(family=="quasipoisson" ){
    
    if(use_opencl==0 ){
      
      if (verbose) {
        
        Rcpp::Rcout << "Initiating NegLL Calculations: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      NegLL=f2_poisson(G4,y, x, mu, P, alpha, wt,progbar);  
      
      
      if (verbose) {
        Rcpp::Rcout << "Initiating Gradient Evaluations: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      cbars2=f3_poisson(G4,y, x,mu,P,alpha,wt,progbar);
      
      
      
    }
    
    else{
      
      
      if (verbose) {
        Rcpp::Rcout << "Initiating f2_f3_opencl: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      
      
      Rcpp::List prepGrad_v3= f2_f3_opencl(
        family,
        link,
        G4,          // NumericMatrix b
        y,           // NumericVector y
        x,           // NumericMatrix x
        mu,          // NumericMatrix mu 
        P,           // NumericMatrix P
        alpha,       // NumericVector alpha
        wt,          // NumericVector wt
        progbar     // int progbar
      );
      
      
      NegLL = prepGrad_v3["qf"];
      cbars2 = Rcpp::as<arma::mat>(prepGrad_v3["grad"]);
      
    }
    
  }
  
  if(family=="Gamma" ){
    
    if(use_opencl==0 ){
      
      if (verbose) {
        
        Rcpp::Rcout << "Initiating NegLL Calculations: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      NegLL=f2_gamma(G4,y, x, mu, P, alpha, wt,progbar);  
      if (verbose) {
        Rcpp::Rcout << "Initiating Gradient Evaluations: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      cbars2=f3_gamma(G4,y, x,mu,P,alpha,wt,progbar);
      
      
      
    } // End use_opencl
    
    
    else{
      
      
      
      if (verbose) {
        Rcpp::Rcout << "Initiating f2_f3_opencl: "
                    << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                    << "\n";
      }
      
      
      Rcpp::List prepGrad_v3= f2_f3_opencl(
        family,
        link,
        G4,          // NumericMatrix b
        y,           // NumericVector y
        x,           // NumericMatrix x
        mu,          // NumericMatrix mu 
        P,           // NumericMatrix P
        alpha,       // NumericVector alpha
        wt,          // NumericVector wt
        progbar     // int progbar
      );
      
      NegLL = prepGrad_v3["qf"];
      cbars2 = Rcpp::as<arma::mat>(prepGrad_v3["grad"]);
      
      
      
    } 
    
    
  }
  
  if(family=="gaussian" ){
    if (verbose) {
      
      Rcpp::Rcout << "Initiating NegLL Calculations: "
                  << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << "\n";
    }
    NegLL=f2_gaussian(G4,y, x, mu, P, alpha, wt);  
    if (verbose) {
      Rcpp::Rcout << "Initiating Gradient Evaluations: "
                  << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << "\n";
    }
    cbars2=f3_gaussian(G4,y, x,mu,P,alpha,wt);
  }
  
  
  //  Rcpp::Rcout << "Finished cbars Calculations: "
  //              << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
  //              << "\n";
  
  
  //  Rcpp::Rcout << "Finished Log-posterior evaluations:" << std::endl;
  
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
  
  
  
  //  Set_Grid_C2(GIndex, cbars, Lint1,Down,Up,loglt,logrt,logct,logU,logP);
  Set_Grid_C2_pointwise(GIndex, cbars, Lint1,Down,Up,loglt,logrt,logct,logU,logP);
  
  
  
  //    Rcpp::Rcout << "Entering setlogP_C2: "
  //                << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
  //                << "\n";
  
  
  if (verbose) {
    
    Rcpp::Rcout << "Setting logP: "
                << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                << "\n";
  }
  
  
  setlogP_C2(logP,NegLL,cbars,G3,LLconst);
  
  
  //      Rcpp::Rcout << "Exiting setlogP_C2: "
  //                  << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
  //                  << "\n";
  
  
  
  NumericMatrix::Column logP2 = logP( _, 1);
  
  
  double  maxlogP=max(logP2);
  
  NumericVector PLSD=exp(logP2-maxlogP);
  
  double sumP=sum(PLSD);
  
  PLSD=PLSD/sumP;
  
  //  Rcout << "Entering Enveloped sort: " << Rcpp::as<std::string>(Rcpp::Function("Sys.time")()) << "\n";
  
  
  if(sortgrid==true){
    
    if (verbose) {
      
      Rcpp::Rcout << "Sorting Grid: "
                  << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")())) 
                  << "\n";
    }
    
    
    Rcpp::List outlist=EnvSort(l1,l2,GIndex,G3,cbars,logU,logrt,loglt,logP,LLconst,PLSD,a_1,E_draws);
    
    return(outlist);
    
  }
  
  
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
                                    bool sortgrid=false
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
    
    NegLL=f2_gaussian(G4,y, x, mu, P, alpha, wt);  
    NegLL_slope=f2_gaussian(G4,y, x, mu, 0*P, alpha, wt);  
    //Rcpp::Rcout << "Finding Value of Gradients at Log-posteriors:" << std::endl;
    cbars2=f3_gaussian(G4,y, x,mu,P,alpha,wt);
    cbars_slope2=f3_gaussian(G4,y, x,mu,0*P,alpha,wt);
    RSS_Out=RSS(y, x,G4,alpha,wt); // Note currenly includes the dispersion in the weight
    
  }
  
  
  //  Rcpp::Rcout << "Finished Log-posterior evaluations:" << std::endl;
  
  // Do a temporary correction here cbars3 should point to correct memory
  // See if this sets cbars
  
  cbars3=cbars2;
  cbars_slope3=cbars_slope2;
  
  Set_Grid_C2(GIndex, cbars, Lint1,Down,Up,loglt,logrt,logct,logU,logP);
  
  setlogP_C2(logP,NegLL,cbars,G3,LLconst);
  
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
