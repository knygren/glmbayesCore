#include "progress_utils.h"


using namespace Rcpp;
using namespace glmbayes::progress;

namespace glmbayes {

namespace progress {

void progress_bar2(double x, double N)
{
  // how wide you want the progress meter to be
  int totaldotz=40;
  double fraction = x / N;
  // part of the progressmeter that's already "full"
  int dotz = round(fraction * totaldotz);
  
  Rcpp::Rcout.precision(3);
  Rcout << "\r                                                                 " << std::flush ;
  Rcout << "\r" << std::flush ;
  Rcout << std::fixed << fraction*100 << std::flush ;
  Rcout << "% [" << std::flush ;
  int ii=0;
  for ( ; ii < dotz;ii++) {
    Rcout << "=" << std::flush ;
  }
  // remaining part (spaces)
  for ( ; ii < totaldotz;ii++) {
    Rcout << " " << std::flush ;
  }
  // and back to line begin 
  
  Rcout << "]" << std::flush ;
  
  // and back to line begin 
  
  Rcout << "\r" << std::flush ;
  
}

}
}