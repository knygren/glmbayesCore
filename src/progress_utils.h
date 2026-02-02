#ifndef PROGRESS_UTILS_H
#define PROGRESS_UTILS_H

#include <string>
#include <tuple>
#include <chrono>
#include <ctime>
#include <RcppArmadillo.h>


// Dependencies:

// 1) progress_utils.cpp

// 2) EnvelopeBuild_Ind_Normal_Gamma.cpp
// 3) EnvelopeEval.cpp
// 4) famfuncs_Gamma.cpp
// 5) famfuncs_binomial.cpp
// 6) famfuncs_poisson.cpp 
// 7) rnnorm_reg_cpp.cpp

namespace glmbayes {

namespace progress {

inline std::string now_hms() {
  std::time_t t = std::time(nullptr);
  char buf[16];
  std::strftime(buf, sizeof(buf), "%H:%M:%S", std::localtime(&t));
  return std::string(buf);
}

struct Timer {
  std::chrono::steady_clock::time_point start;
  void begin() { start = std::chrono::steady_clock::now(); }
  std::tuple<int,int,int> hms() const {
    auto dur = std::chrono::duration_cast<std::chrono::seconds>(
      std::chrono::steady_clock::now() - start
    ).count();
    int h = static_cast<int>(dur / 3600);
    int m = static_cast<int>((dur - h*3600) / 60);
    int s = static_cast<int>(dur - h*3600 - m*60);
    return {h,m,s};
  }
};

inline void print_completed(const char* prefix, const Timer& tm) {
  auto [h,m,s] = tm.hms();
  Rcpp::Rcout << prefix << " completed in: " << h << "h " << m << "m " << s << "s.\n";
}


void progress_bar2(double x, double N);

}
}

#endif
