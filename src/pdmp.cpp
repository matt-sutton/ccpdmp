#include <RcppArmadillo.h>
using namespace Rcpp;
#include "types.h"
#include "piecewise.h"
#include "poly.h"
#include "zigzag_cpp.h"
#include "bps.h"


// [[Rcpp::export]]
List sim_rate_poly(arma::vec eval_times, arma::vec eval_rates, int poly_order ){

  piecewise f_pw;
  poly poly_approx;

  poly_approx.make_poly(eval_times, eval_rates, poly_order);

  // // Add points for f_pw
  for(int i = 0; i < eval_times.size(); ++i) {
    f_pw.push_tf( eval_times(i), (poly_approx.evalpoly(eval_times(i))).head(3));
  }

  arma::vec tausim(2);
  tausim(0) = 0.0; tausim(1) = R::rexp(1);

  double upper = 0;
  arma::vec poly_rate(3);

  while(tausim(1) > 0.0){
    tausim = f_pw.simt(tausim);
    if(tausim(0) > eval_times.max()-1e9){
      break;
    }
    upper = f_pw.evalpw(tausim(0));
    poly_rate = poly_approx.evalpoly(tausim(0));
    bool flip = (R::runif(0,1) < (poly_rate(0) + poly_rate(1))/upper & tausim(1) < 1e-9);
    f_pw.push_tf( tausim(0), poly_rate);
    if(!flip){
      tausim(1) = R::rexp(1);
    }
  }

  if(tausim(1)< 1e-9){
    upper = f_pw.evalpw(tausim(0));
  }

  List ret ;
  ret["t"] = tausim(0) ;
  ret["u"] = tausim(1) ;
  ret["f_evall"] = upper;
  return(ret) ;
}


