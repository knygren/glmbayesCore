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
#include "R_interface.h"

#include <chrono>


using namespace Rcpp;
using namespace openclPort;
using namespace glmbayes::env;

// no export tag

namespace {
double seconds_since_(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
      .count();
}
}  // namespace

namespace glmbayes {

namespace env {
Rcpp::List EnvelopeSize(const arma::vec& a,
                        const Rcpp::NumericMatrix& G1,
                        int Gridtype   ,
                        int n          ,
                        int n_envopt   ,
                        bool use_opencl ,
                        bool verbose    ,
                        GridSetupTiming* timing) 
{
  
  
  int l1 = a.size();
  Rcpp::List G2(l1);
  Rcpp::List GIndex1(l1);
  double E_draws = 1.0;
  
  
  // core count for scaling
  auto t0 = std::chrono::steady_clock::now();
  int core_CNT = opencl_core_count_for_scaling();
  if (timing) timing->core_count += seconds_since_(t0);
  if (verbose) {
    Rcpp::Rcout << "[EnvelopeBuild:EnvelopeSize] OpenCL core count = "
                << core_CNT << "\n";
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
  
  
  
  Rcpp::Function EnvelopeOpt = glmbayes_R::r_envelope_opt();
  Rcpp::NumericVector gridindex(l1);
  
  t0 = std::chrono::steady_clock::now();
  if (Gridtype == 2) {
    if (use_opencl) {
      gridindex = EnvelopeOpt(a, n_envopt, core_CNT);
    } else {
      gridindex = EnvelopeOpt(a, n_envopt, 1);
    }
  }
  if (timing) timing->envelope_opt += seconds_since_(t0);
  
  t0 = std::chrono::steady_clock::now();
  
  
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
  if (timing) timing->size_lists += seconds_since_(t0);
  
  
  return Rcpp::List::create(
    Rcpp::Named("G2")       = G2,
    Rcpp::Named("GIndex1")  = GIndex1,
    Rcpp::Named("E_draws")  = E_draws,
    Rcpp::Named("gridindex")= gridindex
  );
}


} //envelopefuncs
} //glmbayes