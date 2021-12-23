#ifndef BPS_CPP_H
#define BPS_CPP_H
#include <queue>

#include "types.h"
#include "poly.h"
#include "cave_vex_rate.h"
#include "piecewise.h"
#include "inversion.h"
#include "rate.h"

// [[Rcpp::export]]
List bps_cpp(double maxTime, SEXP rate_f, const List& factors, const List& local_updates,
         const List& Data, const arma::vec y, arma::vec x0, arma::vec theta0,
         arma::uvec trac_coords, double tmax = 1.0, int poly_order = 0,
         int nmax = 10^6, int burn = -1, double ref_rate = 0.1, bool sphere = false){

  int p = x0.size(), nEvent= 1, numRates = factors.size();
  double eps = 1e-10, t = 0.0, tau_val, tau_rel;
  bool check_event = true, autoTmax = false, event = false;

  // Information on states and velocities
  arma::vec sk_times(nmax), theta = theta0, x = x0, grad_vals(p), theta_times(p);
  arma::vec clock_taus(numRates), clock_old(numRates);
  arma::ivec sk_flip(nmax);

  if(sphere){
    theta /= sqrt(arma::dot(theta, theta));
  }

  arma::mat sk_points(trac_coords.size(), nmax);
  sk_points.col(0) = x(trac_coords);

  theta_times.zeros(); clock_taus.zeros(); clock_old.zeros();
  sk_times.zeros(); sk_flip.zeros();
  if( burn < 0){
    burn = 0;
  }

  int nShadow = 0, nRef = 0, nIterall = 0;

  std::vector<cave_vex_rate> local_rates;
  local_rates.reserve(numRates);

  typedef std::pair<double, int> Elt;
  std::priority_queue<Elt, std::vector<Elt>, std::greater<Elt> > tau_queue;

  arma::vec s_eval(poly_order+2);
  arma::vec u_vals(poly_order);
  if(poly_order > 1e-4){
    u_vals = get_chebspace_u(poly_order);
    s_eval(arma::span(1, poly_order)) = adjust_range(u_vals, tmax);
  } else {
    u_vals.zeros();
  }
  s_eval(0) = 0.0; s_eval(poly_order+1) = tmax;

  arma::uvec factor_updates, state_updates;
  Elt elt;
  for (int i = 0; i < numRates; i++){
    arma::uvec state_ind = factors[i];
    arma::uvec local_update = local_updates[i];
    local_rates.push_back( cave_vex_rate(x, theta, local_update, grad_vals,
                                             rate_f, y, Data, theta_times, state_ind) );
    local_rates[i].init(s_eval, t);
    local_rates[i].sim_time();
    clock_taus(i) = local_rates[i].get_tau();

    elt = Elt(clock_taus(i), i);
    tau_queue.push(elt);
  }
  //Add refreshment rate
  elt = Elt(-log(R::runif(0,1))/ref_rate, numRates);
  tau_queue.push(elt);

  arma::wall_clock timer;
  double current_time = 0.0;
  timer.tic();

  int mini, nev = 0, rt = 0;
  while( nEvent < nmax + burn){
    nIterall++;
    elt = tau_queue.top();
    mini = elt.second;
    tau_queue.pop();

    if(elt.second == factors.size()){
      // Refreshment
      t = elt.first;
      x = x + (t - theta_times)%theta;
      theta.randn();
      if(sphere){
        theta /= sqrt(arma::dot(theta, theta));
      }
      theta_times = t*arma::ones(p);

      // Propose new times //
      for( int i = 0; i < factors.size(); i++ ){
        clock_old(i) = t;
        local_rates[i].refresh(s_eval, t);
        local_rates[i].sim_time();
        clock_taus(i) = t + local_rates[i].get_tau();
        elt = Elt(clock_taus(i), i);
        tau_queue.push(elt);
      }
      elt = Elt(t -log(R::runif(0,1))/ref_rate, numRates);
      tau_queue.push(elt);

    } else if(std::abs(clock_taus(elt.second) - elt.first) < eps){

      tau_val = clock_taus(mini) - t;
      tau_rel = t - clock_old(mini);

      // Not found on increment //
      if(tau_val + tau_rel >= tmax - 2*eps){
        clock_old(mini) = t + tau_val;

        // Propose new time //
        local_rates[mini].refresh(s_eval, t + tau_val);
        local_rates[mini].sim_time();
        clock_taus(mini) = clock_old(mini) + local_rates[mini].get_tau();

        elt = Elt(clock_taus(mini), mini);
        tau_queue.push(elt);
        nRef++;

      } else {
        // Bounce //
        // Check if time was an event //
        check_event = local_rates[mini].thinning(tau_val, tau_rel, t);
        if( check_event ){
          event = true;
          clock_old(mini) = t;

          // Update states
          t += tau_val;
          factor_updates = local_rates[mini].get_local_update();
          for (arma::uvec::iterator f = factor_updates.begin(); f != factor_updates.end(); ++f){
            int fac = *f;

            // Update states associated to factor
            state_updates = local_rates[fac].get_state_update();
            for(arma::uvec::iterator st = state_updates.begin(); st != state_updates.end(); ++st){
              int state = *st;
              x(state) += (t - theta_times(state))*theta(state);
              theta_times(state) = t;
            }
          }

          // Update Theta
          local_rates[mini].bounce_bps();

          // Update Times
          for (arma::uvec::iterator f = factor_updates.begin(); f != factor_updates.end(); ++f){
            int fac = *f;
            // Update clock //
            local_rates[fac].refresh(s_eval, t);
            local_rates[fac].sim_time();
            clock_taus(fac) = t + local_rates[fac].get_tau();
            clock_old(fac) = t;

            elt = Elt(clock_taus(fac), fac);
            tau_queue.push(elt);
          }

          if( nEvent >= burn ){
            sk_times(nEvent-burn) = t;
            sk_flip(nEvent-burn) = mini;
            sk_points.col(nEvent-burn) = x(trac_coords) + (t - theta_times(trac_coords))%theta(trac_coords);;
            if(nEvent == burn){
              x0 = x;    theta0 = theta;
            }
          }
          nEvent++;

        } else {
          //
          clock_old(mini) = t+tau_val;
          local_rates[mini].refresh(s_eval, clock_old(mini));
          //
          local_rates[mini].sim_time();
          clock_taus(mini) = clock_old(mini) + local_rates[mini].get_tau();

          elt = Elt(clock_taus(mini), mini);
          tau_queue.push(elt);
          nShadow++;
        }
      }
    }

    current_time = timer.toc();
    if(current_time > maxTime){
      if(nEvent < burn){
        Rcout << "Sampler still in burnin phase - set a longer runtime";
      } else {
        sk_flip.shed_rows(nEvent-burn, nmax-1);
        sk_times.shed_rows(nEvent-burn, nmax-1);
        sk_points.shed_cols(nEvent-burn, nmax-1);
      }
      break;
    }
  }
  sk_times -= sk_times(0);

  List ret ;
  ret["times"] = sk_times ;
  ret["flips"] = sk_flip ;
  ret["x0"] = x0 ;
  ret["theta0"] = theta0 ;
  ret["positions"] = sk_points;
  ret["nRef"] = nRef;
  ret["nShadow"] = nShadow;
  ret["nIter"] = nIterall;
  ret["Etime"] = current_time;
  return(ret) ;
}


#endif
