#ifndef INVERSION_H
#define INVERSION_H
#include <RcppArmadillo.h>


// Simulate from a + bs
// [[Rcpp::export]]
arma::vec linear_inv_t(double a, double b, double u, double tmax) {
  arma::vec res(2);
  res.zeros();
  double tsim, t_zero;

  if(std::abs(u) < 1e-14){
    return(res);
  }
  if( a < 0 & b > 0){//ddd
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
  } else if ( a > 0 & b < 0 ){
    if(-pow(a,2)/b + pow(-a,2)/(2*b) >= u ){
      tsim = -a/b -sqrt(pow(a/b,2) + 2*u/b);
      if(tsim < tmax){
        res(0) = tsim; res(1) = -u;
      } else {
        res(0) = tmax; res(1) = -(a*tmax + b*pow(tmax,2)/2);
      }
    } else {
      res(0) = tmax;
      res(1) = pow(a,2)/b - pow(-a,2)/(2*b);
    }
    return(res);

  } else if ( a >= 0 & b > 0 ){
    tsim = -a/b + sqrt(pow(a/b,2) + 2*u/b);
    if(tsim < tmax){
      res(0) = tsim; res(1) = -u;
    } else {
      res(0) = tmax; res(1) = -(a*tmax + b*pow(tmax,2)/2);
    }
    return(res);

  } else if ( a > 0 & b == 0 ){
    tsim = u/a;
    if(tsim < tmax){
      res(0) = tsim; res(1) = -u;
    } else {
      res(0) = tmax; res(1) = -a*tmax;
    }
    return(res);

  } else if ( a <= 0 & b <= 0 ){
    res(0) = tmax;
    return(res);
  }
}

// // Simulate from a + bs
// // [[Rcpp::export]]
// arma::vec linear_inv_t(double a, double b, double u, double tmax) {
//
//   arma::vec res(2);
//   res.zeros();
//   double tsim, t_zero, eps = 1e-15;
//
//   // if B is (effectivley) zero
//   if ( (a > 0.0) & (abs(b) < eps) ){
//     // If constant rate i.e. b=0
//     tsim = exp(log(u) - log(a));
//     if(tsim < tmax){
//       res(0) = tsim; res(1) = -u;
//     } else {
//       res(0) = tmax; res(1) = -exp(log(a)+log(tmax));
//     }
//     return(res);
//   }
//
//   if( (a < 0.0) & (b > 0.0) ){
//     t_zero = - exp(log(a)-log(b));
//
//     if( t_zero > tmax){
//       res(0) = tmax;
//     } else {
//       tsim = t_zero + sqrt(exp(log(2)+log(u)-log(b)));
//       if(tsim < tmax){
//         res(0) = tsim; res(1) = -u;
//       } else {
//         res(0) = tmax; res(1) = -b*pow(tmax + a/b,2)/2; //-(tmax - t_zero)*(a+b*tmax)/2.0;//-exp( log(b) + 2.0*log(tmax-t_zero) - log(2.0) );//  -b*pow(tmax + a/b,2.0)/2.0;
//       }
//     }
//     return(res);
//   } else if ( (a > 0.0) & (b < 0.0) ){
//     t_zero = - exp(log(a)-log(b));
//     if(t_zero > tmax){
//       // always positive
//       if(a*t_zero + b*pow(t_zero,2.0)/2.0 >= u){
//         // check if it doesn't go to infty
//         tsim = t_zero - exp(0.5*log(pow(t_zero,2.0) + exp(log(2.0)+log(u)-log(b))));
//       } else {
//         tsim = t_zero; // Technically goes to Inf
//       }// Maybe don't need???....
//       if(tsim < tmax){
//         res(0) = tsim; res(1) = -u;
//       } else {
//         res(0) = tmax; res(1) = -(exp(log(a)+log(tmax)) + exp(log(b)+2.0*log(tmax)-log(2.0)));
//       }
//     } else {
//       // switch on interval
//       if(exp(log(a)+log(t_zero)) + exp(log(b)+2.0*log(t_zero)-log(2.0)) >= u ){
//         // Doesn't run to end...
//         tsim = t_zero -exp( 0.5*log(pow(t_zero,2.0) + 2.0*exp(log(u)-log(b))));//sqrt(pow(a/b,2.0) + 2.0*u/b);
//         if(tsim < tmax){
//           res(0) = tsim; res(1) = -u;
//         } else {
//           res(0) = tmax; res(1) = -(exp(log(a)+log(tmax)) + exp(log(b)+2.0*log(tmax)-log(2.0)));
//         }
//       } else {
//         res(0) = tmax;
//         res(1) = -(exp(log(a)+log(t_zero)) + exp(log(b)+2.0*log(t_zero)-log(2.0)));
//       }
//     }
//     return(res);
//   } else if ( (a >= 0.0) & (b > 0.0) ){
//     t_zero = - exp(log(a)-log(b));
//     tsim = t_zero + sqrt(pow(a/b,2) + 2.0*u/b);
//     if(tsim < tmax){
//       res(0) = tsim; res(1) = -u;
//     } else {
//       res(0) = tmax; res(1) = -(exp(log(a)+log(tmax)) + exp(log(b)+2.0*log(tmax)-log(2.0)));
//     }
//     return(res);
//   } else if ( (a <= 0.0) & (b <= 0.0) ){
//     // Rcout << "1";
//     res(0) = tmax;
//     return(res);
//   }
// }

// // Simulate from a + bs
// // [[Rcpp::export]]
// arma::vec linear_inv_t(double a, double b, double u, double tmax) {
//   // Rcout<< "a:" << a << "b:" <<b;
//   arma::vec res(2);
//   res.zeros();
//   double tsim, t_zero, eps = 1e-15;
//
//   if(std::abs(u) < eps){
//     return(res);
//   }
//
//   // if always negative
//   if ( a <= eps & b <= eps ){
//     // Rcout << "1";
//     res(0) = tmax;
//     return(res);
//
//   } else if ( a > eps & abs(b) <= eps ){
//     // Rcout << "2";
//     // If constant rate
//     tsim = u/a;
//     if(tsim < tmax){
//       res(0) = tsim; res(1) = -u;
//     } else {
//       res(0) = tmax; res(1) = -a*tmax;
//     }
//     return(res);
//   }
//
//   // Checked
//   if( a < eps & b > eps){
//     // Rcout << "3";
//     t_zero = -a/b;
//     if( t_zero > tmax){
//       res(0) = tmax;
//     } else {
//       tsim = -a/b+sqrt(2.0*u/b);
//       if(tsim < tmax){
//         res(0) = tsim; res(1) = -u;
//       } else {
//         res(0) = tmax; res(1) = -b*pow(tmax + a/b,2.0)/2.0;
//       }
//     }
//     return(res);
//   } else if ( a > eps & b < eps ){
//     // Rcout << "4";
//     t_zero = -a/b;
//     if(t_zero > tmax){
//       // always positive
//       tsim = -a/b -sqrt(pow(a/b,2.0) + 2.0*u/b);
//       if(tsim < tmax){
//         res(0) = tsim; res(1) = -u;
//       } else {
//         res(0) = tmax; res(1) = -(a*tmax + b*pow(tmax,2.0)/2.0);
//       }
//     } else {
//       // switch on interval
//       if(a*t_zero + b*pow(t_zero,2.0)/2.0 >= u ){
//         tsim = -a/b -sqrt(pow(a/b,2.0) + 2.0*u/b);
//         if(tsim < tmax){
//           res(0) = tsim; res(1) = -u;
//         } else {
//           res(0) = tmax; res(1) = -(a*tmax + b*pow(tmax,2.0)/2.0);
//         }
//       } else {
//         res(0) = tmax;
//         res(1) = -(a*t_zero + b*pow(t_zero,2.0)/2.0);
//       }
//     }
//     return(res);
//   } else if ( a >= eps & b > eps ){
//     // Rcout << "5";
//     tsim = -a/b + sqrt(pow(a/b,2) + 2.0*u/b);
//     if(tsim < tmax){
//       res(0) = tsim; res(1) = -u;
//     } else {
//       res(0) = tmax; res(1) = -(a*tmax + b*pow(tmax,2.0)/2.0);
//     }
//     return(res);
//   }
// }

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
