#ifndef INVERSION_H
#define INVERSION_H
#include <RcppArmadillo.h>

// Simulate from a + bs
// [[Rcpp::export]]
arma::vec linear_inv_t(double a, double b, double u, double tmax) {
  arma::vec res(2);
  res.zeros();
  double tsim, t_zero, eps = 1e-14;

  if(std::abs(u) < eps){
    return(res);
  }
  // Checked
  if( a < eps & b > eps){
    t_zero = -a/b;
    if( t_zero > tmax){
      res(0) = tmax;
    } else {
      tsim = -a/b+sqrt(2.0*u/b);
      if(tsim < tmax){
        res(0) = tsim; res(1) = -u;
      } else {
        res(0) = tmax; res(1) = -b*pow(tmax + a/b,2.0)/2.0;
      }
    }
    return(res);
  } else if ( a > eps & b < eps ){
    t_zero = -a/b;
    if(t_zero > tmax){
      // always positive
      tsim = -a/b -sqrt(pow(a/b,2.0) + 2.0*u/b);
      if(tsim < tmax){
        res(0) = tsim; res(1) = -u;
      } else {
        res(0) = tmax; res(1) = -(a*tmax + b*pow(tmax,2.0)/2.0);
      }
    } else {
      // switch on interval
      if(a*t_zero + b*pow(t_zero,2.0)/2.0 >= u ){
        tsim = -a/b -sqrt(pow(a/b,2.0) + 2.0*u/b);
        if(tsim < tmax){
          res(0) = tsim; res(1) = -u;
        } else {
          res(0) = tmax; res(1) = -(a*tmax + b*pow(tmax,2.0)/2.0);
        }
      } else {
        res(0) = tmax;
        res(1) = -(a*t_zero + b*pow(t_zero,2.0)/2.0);
      }
    }
    return(res);
  } else if ( a >= eps & b > eps ){
    tsim = -a/b + sqrt(pow(a/b,2) + 2.0*u/b);
    if(tsim < tmax){
      res(0) = tsim; res(1) = -u;
    } else {
      res(0) = tmax; res(1) = -(a*tmax + b*pow(tmax,2.0)/2.0);
    }
    return(res);

  } else if ( a > eps & abs(b) <= eps ){
    tsim = u/a;
    if(tsim < tmax){
      res(0) = tsim; res(1) = -u;
    } else {
      res(0) = tmax; res(1) = -a*tmax;
    }
    return(res);

  } else if ( a <= eps & b <= eps ){
    res(0) = tmax;
    return(res);
  }
}

// [[Rcpp::export]]
double exp_inv_t(double a, double b, double u) {
  double tsim = 0.0, eps = 1e-14;
  // Inverse of integrated max(0, bexp(a + bs))
  if(std::abs(u) < eps){
    return(tsim);
  }
  if( b > 0.0 & exp(a) >= log(u)){
    return((log(exp(a) - log(u)) - a)/b);
  } else {
    return(10e4);
  }
}

#endif
