#ifndef types_H
#define types_H

// Define types so Cpp arguments may be passed to functions

// Pointer for function returning rate and possibly gradients
// Arguments t, t_old, x, theta, y, Data, rate_index, grad
typedef arma::vec (*ratePtr)(double, arma::vec&, arma::vec&,
                   arma::vec&, const arma::vec&, const List&,
                   arma::uvec, bool);

// Pointer for function returning rate and possibly gradients
// Arguments t, t_old, x, theta, y, Data, Factors, Factor_updates, grad
typedef arma::mat (*rateMPtr)(double, arma::vec&, arma::vec&,
                   arma::vec&, const arma::vec&, const List&,
                   const List&, arma::uvec, bool);

// Pointer for function returning -log pi and grad -log pi
// Arguments x, Data, y
typedef arma::vec (*mcmcGradPtr)(arma::vec&, const List&, const arma::vec&);
typedef double (*mcmcPostPtr)(arma::vec&, const List&, const arma::vec&);

#endif
