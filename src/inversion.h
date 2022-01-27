#ifndef INVERSION_H
#define INVERSION_H
#include <RcppArmadillo.h>


// Simulate from a + bs
// [[Rcpp::export]]
arma::vec linear_inv_t(double a, double b, double u, double tmax) {
  arma::vec res(2);
  res.zeros();
  double tsim, t_zero, eps = 1e-13;

  if(std::abs(u) < eps){
    return(res);

  } else if ( a > eps & abs(b) < eps ){
    tsim = u/a;
    if(tsim < tmax){
      res(0) = tsim; res(1) = -u;
    } else {
      res(0) = tmax; res(1) = -a*tmax;
    }
    return(res);

  } else if( a < eps & b > eps){
    t_zero = -a/b;

    if( t_zero > tmax){
      res(0) = tmax;
    } else {
      tsim = t_zero + sqrt(2*u/b);
      if(tsim < tmax){
        res(0) = tsim; res(1) = -u;
      } else {
        res(0) = tmax; res(1) = -(tmax - t_zero)*(a+b*tmax)/2.0;//-b*pow(tmax + a/b,2)/2;
      }
    }
    return(res);
  } else if ( a > eps & b < eps ){
    if(-pow(a,2)/b + pow(-a,2)/(2*b) >= u ){
      tsim = -a/b -sqrt(pow(a/b,2) + 2*u/b); // tsim lands in triangle
      if(tsim < tmax){
        res(0) = tsim; res(1) = -u;
      } else {
        res(0) = tmax; res(1) = -(a*tmax + b*pow(tmax,2)/2);
      }
    } else {
      res(0) = tmax;
      if(tmax > -a/b){
        res(1) = pow(a,2)/b - pow(-a,2)/(2*b);
      } else {
        res(1) = -(a*tmax + b*pow(tmax,2)/2);
      }
    }
    return(res);

  } else if ( a >= eps & b > 0 ){
    tsim = -a/b + sqrt(pow(a/b,2) + 2*u/b);
    if(tsim < tmax){
      res(0) = tsim; res(1) = -u;
    } else {
      res(0) = tmax; res(1) = -(a*tmax + b*pow(tmax,2)/2);
    }
    return(res);

  } else if ( a <= 0 & b <= 0 ){
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
