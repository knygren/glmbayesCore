// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "nmath_local.h"
#include "dpq_local.h"


#include "rng_utils.h"
#include <random>



// Thread-local RNG and distribution
thread_local std::mt19937 safe_rng_engine(std::random_device{}());
thread_local std::uniform_real_distribution<> safe_rng_dist(0.0, 1.0);

using namespace glmbayes::rng;


// Safe inverse-gamma CDF using nmath/rmath pgamma
double p_inv_gamma_safe(double dispersion,
                        double shape,
                        double rate) {
  // For X ~ InvGamma(shape, rate), Y = 1/X ~ Gamma(shape, rate)
  // So P(X <= d) = P(Y >= 1/d) = 1 - F_Y(1/d)
  double y = 1.0 / dispersion;
  
  // Call the ported pgamma (not R::pgamma)
  // Arguments: x, shape, scale, lower_tail, log_p
  double Fy = pgamma_local(y, shape, 1.0 / rate, /*lower_tail=*/1, /*log_p=*/0);
  
  return 1.0 - Fy;
}


double q_inv_gamma_safe(double p,
                        double shape,
                        double rate,
                        double disp_upper,
                        double disp_lower) {
  // Compute probabilities at the bounds using safe pgamma
  double p_upp = p_inv_gamma_safe(disp_upper, shape, rate);
  double p_low = p_inv_gamma_safe(disp_lower, shape, rate);
  
  // Map uniform p into [p_low, p_upp]
  double p1 = p_low + p * (p_upp - p_low);
  double p2 = 1.0 - p1;
  
  // Invert via safe qgamma (ported from nmath/rmath)
  return 1.0 / qgamma_local(p2, shape, 1.0 / rate, /*lower_tail=*/1, /*log_p=*/0);
}






namespace glmbayes {

namespace rng {

// Core sampling function
double safe_runif() {
  return safe_rng_dist(safe_rng_engine);
}


// 
// // Declaration (e.g. in a header if needed)
// // double r_invgamma_safe(double shape, double rate,
// //                        double disp_upper, double disp_lower);
// 
// // Definition (in your .cpp file)
double r_invgamma_safe(double shape,
                       double rate,
                       double disp_upper,
                       double disp_lower) {
  // draw uniform(0,1) from thread‑local RNG
  double p = safe_runif();
  
  // invert CDF at p to get inverse‑gamma draw
  // q_inv_gamma must be pure C++ math, no R calls
  return q_inv_gamma_safe(p, shape, rate, disp_upper, disp_lower);
}

}
}
