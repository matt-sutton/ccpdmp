#ifndef MALA_H
#define MALA_H
#include "types.h"

double log_q_mala(arma::vec& x, arma::vec& y, arma::vec& x_grad, double eps){
  arma::vec res = y - x + (eps/2.0)*x_grad;
  return(-1.0/(2.0*eps) * arma::dot(res,res));
}

// [[Rcpp::export]]
List mala_pnt(double maxTime, SEXP post_f, SEXP grad_f, const List& Data, const arma::vec& y,
              arma::vec x0, double epsilon = 0.1, int nmax = 10^6,
              int burn = -1, int thin = 1){

  XPtr<mcmcGradPtr> xpgrad(grad_f);
  mcmcGradPtr get_grad_fun = *xpgrad;

  XPtr<mcmcPostPtr> xppost(post_f);
  mcmcPostPtr get_post_fun = *xppost;

  int p = x0.size(), nEvent = 1;
  double eps = 1e-10, acc_prob, Ux, Ux_prop, Uz, Uz_prop;

  arma::mat sk_points(p,nmax);
  arma::vec acc_count(nmax), nlogpi(nmax), x = x0, x_prop(p), z_vals(p), z_prop(p);

  if( burn < 0){
    burn = 0;   sk_points.col(0) = x0;
  }

  Ux = get_post_fun(x, Data, y);
  nlogpi(0) = Ux;
  arma::vec x_grad(p), x_prop_grad(p);
  x_grad = get_grad_fun(x, Data, y);

  arma::wall_clock timer;
  double current_time = 0.0;
  timer.tic();

  while( nEvent < nmax + burn){
    z_vals = rnorm(p);// p_curr
    x_prop = x - (epsilon/2.0)  * x_grad + sqrt(epsilon) * z_vals;

    Ux_prop = get_post_fun(x_prop, Data, y);
    x_prop_grad = get_grad_fun(x_prop, Data, y);

    acc_prob = std::min(1.0, exp(Ux - Ux_prop +
      log_q_mala(x_prop, x, x_prop_grad, epsilon) -
      log_q_mala(x, x_prop, x_grad, epsilon)));

    if( R::runif(0,1) < acc_prob ){
      x = x_prop;
      Ux = Ux_prop;
      x_grad = x_prop_grad;
    }
    if( nEvent >= burn){
      sk_points.col(nEvent-burn) = x;
      nlogpi(nEvent-burn) = Ux;
      acc_count(nEvent-burn) = acc_prob;
    }
    nEvent++;

    current_time = timer.toc();
    if(current_time > maxTime){
      if(nEvent < burn){
        Rcout << "Sampler still in burnin phase - set a longer runtime";
      } else {
        sk_points.shed_cols(nEvent-burn, nmax-1);
        acc_count.shed_rows(nEvent-burn, nmax-1);
        nlogpi.shed_rows(nEvent-burn, nmax-1);
      }
      break;
    }
  }

  List ret ;
  ret["samples"] = sk_points ;
  ret["acc_probs"] = acc_count ;
  ret["nlogpi"] = nlogpi ;
  ret["Etime"] = current_time;
  return(ret) ;
}

#endif
