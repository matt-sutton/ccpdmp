#ifndef PIECEWISE_H
#define PIECEWISE_H

#include <RcppArmadillo.h>
using namespace Rcpp;
#include "inversion.h"

class piecewise {
public:

  // Store the function evaluations
  std::map< double, arma::vec > fEvaluations;
  double eps = 1e-14;

  // add new (t, f) pair --> f = (fu, fn, gn)
  void push_tf(double t_in, arma::vec f_in) {
    fEvaluations.insert({t_in, f_in});
  }

  void refresh() {
    fEvaluations.clear();
  }

  void update(){
    // Pulls everything back and saves the final element to reuse
    arma::vec f_in = (*std::prev(fEvaluations.end(),1)).second;
    fEvaluations.clear();
    fEvaluations.insert({0.0, f_in});
  }

  // Get intercept (a) and gradients (b)
  arma::vec get_ab(std::pair<double, arma::vec> f0, std::pair<double, arma::vec> f1) {
    // ab_values = (a0,b0,t0,a1,b1,t1)
    arma::vec ab_values(6);
    ab_values.zeros();
    // Convex part
    double b_u = (f1.second[0] - f0.second[0])/(f1.first - f0.first);
    double a_u = f1.second[0] - b_u*f1.first;

    // Concave part 1
    double a_n1 = f0.second[1] - f0.second[2]*f0.first;
    double b_n1 = f0.second[2];

    // Upper bound
    ab_values(0) = a_u + a_n1;
    ab_values(1) = b_u + b_n1;

    double delta_grad = f1.second[2] - f0.second[2];
    double a_n2 = f1.second[1] - f1.second[2]*f1.first;
    double b_n2 = f1.second[2];

    if( std::abs(delta_grad) < eps ){
      ab_values(2) = f1.first;
      if( delta_grad >= 0 ){
        ab_values(0) = a_u + a_n2;
        ab_values(1) = b_u + b_n2;
      } else {
        ab_values(0) = a_u + a_n1;
        ab_values(1) = b_u + b_n1;
      }
    } else {
      // Concave part 2
      ab_values(3) = a_u + a_n2;
      ab_values(4) = b_u + b_n2;

      // time concave gradients intersect
      double t_intersect = (a_n1 - a_n2)/(b_n2 - b_n1);
      t_intersect = std::min(t_intersect,f1.first); // uperbound ?
      t_intersect = std::max(t_intersect,f0.first); // lowerbound ?

      ab_values(2) = t_intersect;
      ab_values(5) = f1.first;
    }
    return(ab_values);
  }

  // evaluate f_pw at t
  double evalpw(double t_ev) {
    // find interval on t_ev
    std::map< double, arma::vec >::iterator it_t0, it_t1;
    // it_t0 = fEvaluations.lower_bound(t_ev-1e-10);
    // it_t1 = std::next(it_t0,1);
    it_t1 = fEvaluations.upper_bound(t_ev);
    it_t0 = std::prev(it_t1,1);

    // calc intercepts and gradients (a0,b0,t0,a1,b1,t1)
    arma::vec ab_values = get_ab(*it_t0, *it_t1);

    // Check which interval t_ev is on
    if( t_ev <= ab_values[2] | ab_values[5] <= eps){ // minus tau ???
      return(ab_values[0] + ab_values[1]*t_ev);
    } else {
      return(ab_values[3] + ab_values[4]*t_ev);
    }
  }

  // Simulate a time on the current range
  arma::vec simt(arma::vec tu) {

    std::map< double, arma::vec >::iterator it_t0, it_t1;
    it_t0 = fEvaluations.lower_bound(tu(0)-1e-14); // Need to subtract for rounding precision.
    double tmax = (*std::prev(fEvaluations.end(),1)).first;

    arma::vec ab_values(6);
    while( tu(0) < tmax ) {
      //Calculate the (intercept,  gradient)
      it_t1 = std::next(it_t0,1);
      if(it_t1 == fEvaluations.end()){
        tu(0) = tmax;
        return(tu);
      }
      ab_values = get_ab(*it_t0, *it_t1);

      // tu += linear_inv_t(ab_values[0]+exp(log(ab_values[1])+log(tu(0))), ab_values[1], tu(1), ab_values[2] - tu(0));
      tu += linear_inv_t(ab_values[0]+ab_values[1]*tu(0), ab_values[1], tu(1), ab_values[2] - tu(0));
      // Check if time simulated
      if(tu(1) < eps | tu(0) > tmax-eps){
        return(tu);
      }
      if( ab_values[5] -tu(0) > 1e-10 ){
        // tu += linear_inv_t(ab_values[3]+exp(log(ab_values[4])+log(tu(0))), ab_values[4], tu(1), ab_values[5] - tu(0));
        tu += linear_inv_t(ab_values[3]+ab_values[4]*tu(0), ab_values[4], tu(1), ab_values[5] - tu(0));
        // Check if time simulated
        if(tu(1) < eps | tu(0) > tmax-eps){
          return(tu);
        }
      }
      // Increment the iterator
      it_t0 = it_t1;
    }
    return(tu);
  }

};

#endif
