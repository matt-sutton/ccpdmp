#include <RcppArmadillo.h>
using namespace Rcpp;
#include "types.h"
#include "piecewise.h"
#include "poly.h"
#include "zigzag_cpp.h"
#include "bps_cpp.h"
#include "hmc.h"
#include "mala.h"


// [[Rcpp::export]]
List sim_rates(arma::vec eval_times, arma::mat eval_rates, int poly_order = 0, int n_points = -1, double u = -1 ){

  piecewise f_pw;
  poly poly_approx;

  bool no_poly = all(arma::abs(eval_rates.row(3)) < 1e-10);

  // If No poly
  if(no_poly){
    for(int i = 0; i < eval_times.size(); ++i) {
      f_pw.push_tf( eval_times(i), eval_rates.col(i).head(3) );
    }
  } else {
    poly_approx.make_poly(eval_times, eval_rates.row(3).t(), poly_order);
    // // Add points for f_pw
    for(int i = 0; i < eval_times.size(); ++i) {
      // fu(0), fn(1), gn(2), poly(3), exact(4)
      f_pw.push_tf( eval_times(i), eval_rates.col(i).head(3) + (poly_approx.evalpoly(eval_times(i))).head(3));
    }
  }

  arma::vec tausim(2);
  tausim(0) = 0.0;
  if(u < 0) {
    tausim(1) = R::rexp(1);
    } else {
    tausim(1) = u;
  };
  tausim = f_pw.simt(tausim);

  double upper = 0;
  if(tausim(1)< 1e-9){
    upper = f_pw.evalpw(tausim(0));
  }

  List ret ;
  ret["t"] = tausim(0) ;
  ret["u"] = tausim(1) ;
  ret["f_evall"] = upper;

  if(n_points > 0){
    arma::vec range = arma::linspace<arma::vec>(1e-10, eval_times.max()-1e-10, n_points);
    arma::vec upper_range(n_points);
    for(int i = 0; i < range.size(); ++i) {
      upper_range(i) = f_pw.evalpw(range(i));
    }
    ret["f_range"] = range;
    ret["f_upper_range"] = upper_range;
  }

  return(ret) ;
}

// [[Rcpp::export]]
List cc_sim(arma::vec eval_times, arma::mat eval_rates, int n_points = -1, double u = -1 ){

  piecewise f_pw;

  for(int i = 0; i < eval_times.size(); ++i) {
    f_pw.push_tf( eval_times(i), eval_rates.col(i).head(3) );
  }

  arma::vec tausim(2);
  tausim(0) = 0.0;
  if(u < 0) {
    tausim(1) = R::rexp(1);
  } else {
    tausim(1) = u;
  };
  tausim = f_pw.simt(tausim);

  double upper = 0;

  if(abs(tausim(0) - eval_times.max())> 1e-13){
    upper = f_pw.evalpw(tausim(0));
  }

  List ret ;
  ret["t"] = tausim(0) ;
  ret["u"] = tausim(1) ;
  ret["f_evall"] = upper;

  // Additional (optional) plotting returns
  if(n_points > 0){
    arma::vec range = arma::linspace<arma::vec>(1e-10, eval_times.max()-1e-10, n_points);
    arma::vec upper_range(n_points);
    for(int i = 0; i < range.size(); ++i) {
      upper_range(i) = f_pw.evalpw(range(i));
    }
    ret["range"] = range;
    ret["upper_range"] = upper_range;
  }

  return(ret) ;
}

// [[Rcpp::export]]
List sim_rate_poly(arma::vec eval_times, arma::vec eval_rates, int poly_order, int n_points = -1 ){

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
    if(tausim(0) > (eval_times.max()-1e-9)){
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

  // Additional (optional) plotting returns
  if(n_points > 0){
    arma::vec range = arma::linspace<arma::vec>(1e-10, eval_times.max()-1e-10, n_points);
    arma::vec upper_range(n_points);
    for(int i = 0; i < range.size(); ++i) {
      upper_range(i) = f_pw.evalpw(range(i));
    }
    ret["range"] = range;
    ret["upper_range"] = upper_range;
  }

  return(ret) ;
}


