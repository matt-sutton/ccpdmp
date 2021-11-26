#ifndef CAVE_VEX_RATE_H
#define CAVE_VEX_RATE_H
#include "poly.h"
#include "rate.h"

class cave_vex_rate: public rate_ref {
public:
  ratePtr get_rate_fun;
  // Return fu(0), fn(1), gn(2), poly(3), exact(4)

  cave_vex_rate(arma::vec& x0, arma::vec& theta0, arma::uvec local_update0,
                    arma::vec& grad_vals0, SEXP fun_, const arma::vec& y0,
                    const List& Data0, arma::vec& t_old0, arma::uvec rate_index0)
    : rate_ref{ x0, theta0, local_update0 }, grad_vals{grad_vals0}, Data{Data0}, y{y0}, t_old{t_old0}
    {
      rate_index = rate_index0;
      XPtr<ratePtr> xpfun(fun_);
      get_rate_fun = *xpfun;
    }

  void init(arma::vec s_eval, double t_curr){
    int n_pnt = s_eval.size();

    arma::mat rate_mat(5,n_pnt);
    for(int i = 0; i < n_pnt; ++i) {
      rate_mat.col(i) = get_rate_fun(t_curr + s_eval(i), t_old, x, theta, y, Data, rate_index, false);
    }
    no_poly = all(arma::abs(rate_mat.row(3)) < 1e-10);

    if(no_poly){
      // If No poly
      for(int i = 0; i < s_eval.size(); ++i) {
        f_pw.push_tf( s_eval(i), rate_mat.col(i).head(3) );
      }

    } else {
      // Poly ~~~ Builds poly on Cheb points...
      poly_approx.make_poly(s_eval, rate_mat.row(3).t(), s_eval.size()-1);

      // // Add points for f_pw
      for(int i = 0; i < s_eval.size(); ++i) {
        f_pw.push_tf( s_eval(i), rate_mat.col(i).head(3) + (poly_approx.evalpoly(s_eval(i))).head(3));
      }
    }
    tausim(0) = 0.0; tausim(1) = R::rexp(1);
  }

  void refresh(arma::vec s_eval,double t_curr){
    f_pw.refresh();
    init(s_eval, t_curr);
  }

  arma::uvec get_state_update(){
    return(rate_index);
  }

  void bounce_bps(){
    theta(rate_index) -= 2*arma::dot(theta(rate_index), grad_vals(rate_index))/arma::dot(grad_vals(rate_index),grad_vals(rate_index))*grad_vals(rate_index);
  }
  void bounce_zigzag(){
    theta(rate_index) = -theta(rate_index);
  }

  bool thinning(double tau_val, double tau_rel = 0.0, double t_curr = 0.0){
    arma::vec fun_eval = get_rate_fun(t_curr + tau_val, t_old, x, theta, y, Data, rate_index, true);
    double upper = f_pw.evalpw(tau_val + tau_rel);
    if(fun_eval(4) > upper + 2*eps){
      Rcout << "\n Invalid thining for Rate : " << rate_index(0) << " at time t= "<< tau_val + tau_rel<< " rate:" << fun_eval(4) << " upper-bound:" << upper;
    }
    bool flip = (R::runif(0,1) < fun_eval(4)/upper & tausim(1) < 2*eps);

    if(!flip){
      if(no_poly){
        f_pw.push_tf( tau_val + tau_rel, fun_eval.head(3));
      } else {
        arma::vec poly_eval = poly_approx.evalpoly(tau_val + tau_rel);
        arma::vec rate_eval = poly_eval.head(3) + fun_eval.head(3);
        f_pw.push_tf( tau_val + tau_rel, rate_eval);
      }
    } else{
      // update gradient for flip..
      grad_vals(rate_index) = fun_eval.tail(fun_eval.size() - 5);
    }
    tausim(1) = R::rexp(1);
    return( flip );
  }

protected:
  arma::vec& t_old, grad_vals;
  const arma::vec& y;
  const List& Data;
  arma::uvec rate_index;
  double eps = 1e-10;
  bool no_poly;
  poly poly_approx;
};

#endif
