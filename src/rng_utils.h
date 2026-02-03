// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#ifndef RNG_UTILS_H
#define RNG_UTILS_H

namespace glmbayes{

namespace rng {

// Thread-safe uniform RNG [0, 1)
double safe_runif();

double r_invgamma_safe(double shape,
                       double rate,
                       double disp_upper,
                       double disp_lower);



double ctrnorm_cpp(double lgrt,double lglt,double mu,double sigma);

double r_invgamma(double shape,double rate,double disp_upper,double disp_lower);


}
}
#endif
