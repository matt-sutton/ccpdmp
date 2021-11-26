#ifndef POLY_H
#define POLY_H

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::vec get_chebspace_u(int n_pnts){
  arma::vec ts = (arma::linspace<arma::vec>(0.0, n_pnts-1, n_pnts) + 0.5)/n_pnts;
  arma::vec tau_shifted = arma::cos(ts*arma::datum::pi);
  return(tau_shifted);
}

arma::vec adjust_range(arma::vec tau_shifted, double tmax){
  arma::vec tau = (tau_shifted*tmax + tmax)/2.0;
  return(tau);
}

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

#endif
