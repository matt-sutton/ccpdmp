#ifndef POLY_H
#define POLY_H

#include <RcppArmadillo.h>
using namespace Rcpp;

class poly {
public:

  void make_poly(arma::vec tau, arma::vec f_eval, int poly_order){
    // assumes that tau is on correct scale.
    arma::vec coef = arma::polyfit(tau, f_eval, poly_order);

    coef(arma::find(arma::abs(coef) < 1e-10)).zeros();

    arma::vec zero_coef = arma::zeros(coef.size());
    coefu = arma::max(zero_coef, coef);
    coefn = arma::min(zero_coef, coef);
    coefgn = coefn.head(coefn.size() - 1)%arma::linspace<arma::vec>(coefn.size()-1,1,coefn.size()-1);
  }

  arma::vec evalpoly(double t_ev) {
    arma::vec rate_eval(3);
    arma::vec t_evs(1);
    t_evs(0) = t_ev;
    rate_eval(0) = arma::as_scalar(arma::polyval(coefu, t_evs));
    rate_eval(1) = arma::as_scalar(arma::polyval(coefn, t_evs));
    rate_eval(2) = arma::as_scalar(arma::polyval(coefgn, t_evs));
    return(rate_eval);
  }

  double evalpoly_full(double t_ev) {
    arma::vec t_evs(1);
    t_evs(0) = t_ev;
    return(arma::as_scalar(arma::polyval(coefu+coefn, t_evs)));
  }

protected:
  arma::vec f_eval, tau, coefu, coefn, coefgn;
  int poly_order;
  double eps = 0.1, tmax;
};

// class poly {
// public:
//
//   arma::vec get_polyCoef(){
//     int N = cheb_coef.size();
//     double sv;
//     arma::vec poly_coef(N), poly_coeff(N);
//     poly_coef.zeros(); poly_coeff.zeros();
//     poly_coef(0)=cheb_coef(N-1);
//     for (int j=N-2;j>0;j--) {
//       for (int k=N-j;k>0;k--) {
//         sv=poly_coef(k);
//         poly_coef(k)=2.0*poly_coef(k-1)-poly_coeff(k);
//         poly_coeff(k)=sv;
//       }
//       sv=poly_coef(0);
//       poly_coef(0) = -poly_coeff(0)+cheb_coef(j);
//       poly_coeff(0) = sv;
//     }
//     for (int j=N-1;j>0;j--) {
//       poly_coef(j)=poly_coef(j-1)-poly_coeff(j);
//     }
//     poly_coef[0] = -poly_coeff(0)+0.5*cheb_coef(0);
//     return(poly_coef);
//   }
//
//   arma::vec polyshift(arma::vec polyc){
//     double cnst = 2.0/tmax;
//     double fac = cnst;
//     int n = polyc.size();
//     for (int j=1;j<n;j++) {
//       polyc(j) *= fac;
//       fac *= cnst;
//     }
//     cnst=0.5*tmax;
//     for (int j=0;j<=n-2;j++){
//       for (int k = n-2;k>=j;k--){
//         polyc(k) -= cnst*polyc(k+1);
//       }
//     }
//     return(polyc);
//   }
//
//   void get_cheb_coef(){
//     int n = tau_shifted.size();
//     double c0 = 2.0/n;
//     cheb_coef.zeros(n);
//     arma::vec ts = (arma::linspace<arma::vec>(0.0, n-1, n) + 0.5)/n;
//     for(int j = 0; j < n; j++){
//       cheb_coef(j) = c0*arma::dot(f_eval,arma::cos(arma::datum::pi*j*ts));
//     }
//     eps = 2*abs(cheb_coef(poly_order));
//     arma::uvec zero_coef = arma::find(arma::abs(cheb_coef) < 1e-5);
//     cheb_coef(zero_coef).zeros();
//     cheb_coef = cheb_coef(arma::linspace<arma::uvec>(0, poly_order-1, poly_order));
//     poly_order -= 1;
//   }
//
//   void make_poly(arma::vec tau0, arma::vec f_eval0, double tmax0){
//     // assumes that tau is on correct scale.
//     tmax = tmax0;
//     tau_shifted = tau0;
//     tau = adjust_range(tau_shifted, tmax);
//     f_eval = f_eval0;
//     poly_order = tau0.size() - 1;
//
//     get_cheb_coef();
//     arma::vec coef = get_polyCoef();
//     coef = polyshift(coef);
//     coef = arma::reverse(coef);
//
//     arma::vec zero_coef = arma::zeros(coef.size());
//     coefu = arma::max(zero_coef, coef);
//     coefn = arma::min(zero_coef, coef);
//     coefgn = coefn.head(coefn.size() - 1)%arma::linspace<arma::vec>(coefn.size()-1,1,coefn.size()-1);
//
//     coefu(poly_order) += eps;
//   }
//
//   // add new (t, f) --> f = (fu, fn, gn)
//   void push_tf(arma::vec tau0, arma::vec f_eval0) {
//     poly_order = tau0.size() - 1;
//     f_eval = f_eval0;
//     tau = tau0;
//   }
//   void push_coef(double t_new, double f_new){
//     arma::vec tn(1), fn(1);
//     tn(0) = t_new; fn(0) = f_new;
//     tau.insert_rows(tau.n_rows, tn);
//     f_eval.insert_rows(f_eval.n_rows, fn);
//
//     // f_eval.shed_rows(arma::find(tau < t_new-1e-10));
//     // tau.shed_rows(arma::find(tau < t_new-1e-10));
//     // poly_order = tau.size() -1;
//
//     //poly_order += 1;
//     arma::vec coef = arma::polyfit(tau, f_eval, poly_order);
//
//     arma::vec zero_coef = arma::zeros(coef.size());
//     coefu = arma::max(zero_coef, coef);
//     coefn = arma::min(zero_coef, coef);
//     coefgn = coefn.head(coefn.size() - 1)%arma::linspace<arma::vec>(coefn.size()-1,1,coefn.size()-1);
//
//     //coefu(poly_order) += eps; // abs(cheb_coef(poly_order))
//   }
//
//   arma::vec evalpoly(double t_ev) {
//     arma::vec rate_eval(3);
//     arma::vec t_evs(1);
//     t_evs(0) = t_ev;
//     rate_eval(0) = arma::as_scalar(arma::polyval(coefu, t_evs));
//     rate_eval(1) = arma::as_scalar(arma::polyval(coefn, t_evs));
//     rate_eval(2) = arma::as_scalar(arma::polyval(coefgn, t_evs));
//     return(rate_eval);
//   }
//   double evalpoly_full(double t_ev) {
//     arma::vec t_evs(1);
//     t_evs(0) = t_ev;
//     return(arma::as_scalar(arma::polyval(coefu+coefn, t_evs)));
//   }
//
//   arma::vec get_poly(){
//     return(coefu + coefn);
//   }
//
// protected:
//   arma::vec f_eval, tau, tau_shifted, cheb_coef, coefu, coefn, coefgn;
//   int poly_order;
//   double eps = 0.1, tmax;
// };

#endif
