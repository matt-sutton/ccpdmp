#ifndef HMC_H
#define HMC_H
#include "types.h"

// [[Rcpp::export]]
List hmc(double maxTime, SEXP post_f, SEXP grad_f, const List& Data, const arma::vec& y,
              arma::vec x0, double epsilon = 0.1, double stoch_time = 1, int nmax = 10^6,
              int burn = -1, int thin = 1){

  XPtr<mcmcGradPtr> xpgrad(grad_f);
  mcmcGradPtr get_grad_fun = *xpgrad;

  XPtr<mcmcPostPtr> xppost(post_f);
  mcmcPostPtr get_post_fun = *xppost;

  int p = x0.size(), nEvent = 1;
  int L = std::ceil(stoch_time/epsilon);
  double eps = 1e-10, acc_prob, Ux, Ux_prop, Uz, Uz_prop;

  arma::mat sk_points(p,nmax);
  arma::vec acc_count(nmax), nlogpi(nmax), x = x0, x_prop(p), z_vals(p), z_prop(p);

  if( burn < 0){
    burn = 0;   sk_points.col(0) = x0;
  }
  Ux = get_post_fun(x, Data, y);
  nlogpi(0) = Ux;
  arma::vec x_grad(p);

  arma::wall_clock timer;
  double current_time = 0.0;
  timer.tic();

  while( nEvent < nmax + burn){
    z_vals = rnorm(p);// p_curr
    x_prop = x; z_prop = z_vals; // q
    x_grad = get_grad_fun(x_prop, Data, y);

    // half step momentum first
    z_prop -= epsilon/2.0 * x_grad;

    for( int l_val = 0; l_val < L; l_val++){
      // full position move
      x_prop += epsilon * z_prop;
      // full momentum move
      x_grad = get_grad_fun(x_prop, Data, y);
      if( l_val < L-1){
        z_prop -= epsilon * x_grad;
      }
    }
    // final half step
    z_prop -= epsilon/2.0 * x_grad;
    // flip velocity
    z_prop = - z_prop;

    Uz_prop = arma::dot(z_prop,z_prop)/2.0;
    Uz = arma::dot(z_vals,z_vals)/2.0;
    Ux_prop = get_post_fun(x_prop, Data, y);

    acc_prob = std::min(1.0, exp(Ux - Ux_prop + Uz - Uz_prop));
    if( R::runif(0,1) < acc_prob ){
      x = x_prop;
      Ux = Ux_prop;
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
