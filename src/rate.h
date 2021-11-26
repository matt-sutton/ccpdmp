#ifndef rate_H
#define rate_H

#include "types.h"
#include "piecewise.h"
#include "inversion.h"

class rate_ref {
public:
  rate_ref(arma::vec& x0, arma::vec& theta0, arma::uvec local_update0)
    : x{x0}, theta{theta0}, local_update{local_update0}
    {
      if(local_update(0) >= x.size()){
        update_all = true;
      } else {
        update_all = false;
      }
    }

  void linear_move(double tau_val, int localu){
    x(localu) += tau_val*theta(localu);
  }
  void sim_time(){
    tausim = f_pw.simt(tausim);
  }
  double get_tau(){
    return(tausim(0));
  }
  void refresh_v(){
    theta.randn();
    theta = theta/sqrt(arma::dot(theta,theta));
  }
  void bounce_zigzag(int flip){
    theta(flip) = -theta(flip);
  }
  arma::uvec get_local_update(){
    if(update_all){
      return(arma::regspace<arma::uvec>(0,x.size()-1));
    } else {
      return(local_update);
    }
  }

protected:
  // store position and velocity
  piecewise f_pw;
  arma::vec& x;
  arma::vec& theta;
  arma::uvec local_update;
  bool update_all = false;
  arma::vec tausim = arma::zeros(2);
};


#endif
