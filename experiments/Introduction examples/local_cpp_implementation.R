library(RcppXPtrUtils)
local_rate <- cppXPtr("arma::vec get_rate(double t, arma::vec& t_old, arma::vec& x, arma::vec& theta, const arma::vec& y, const List& Data, arma::uvec rate_inds, bool grad) {
                     const double rho = Data[0];
                     int num_terms = rate_inds.size();

                     int rsize = 5;
                     if(grad){
                     rsize += num_terms;
                     }
                     arma::vec res(rsize);
                     res.zeros();

                     int d = x.size();
                     double rho2 = std::pow(rho,2.0);
                     double x_t, grd;

                     for(int i = 0; i < num_terms; i++){
                          int term = rate_inds(i);
                          grd = 0.0;
                          x_t = x(term) + theta(term)*(t-t_old(term));
                          res(3) -= theta(term)*y(term);
                          // Likelihood terms
                          if( theta(term) < 0 ) {
                              res(1) += theta(term)*exp(x_t);
                              res(2) += theta(term)*theta(term)*exp(x_t);
                          } else {
                              res(0) += theta(term)*(exp(x_t));
                          }
                          grd += exp(x_t) - y(term);
                          // Latent Gaussian terms
                          for(int rt = std::max(0,term-1); rt < std::min(term+2,d); rt ++){
                              x_t = x(rt) + theta(rt)*(t-t_old(rt));
                              if(rt == term){
                                 res(3) += theta(term)*(1+rho2)*x_t;
                                 grd += (1+rho2)*x_t;
                              } else {
                                 res(3) += theta(term)*rho*x_t;
                                 grd += rho*x_t;
                              }
                          }
                          if(term == d-1){
                             x_t = x(term) + theta(term)*(t-t_old(term));
                             res(3) += theta(term)*x_t-theta(term)*(1+rho2)*x_t;
                             grd += x_t - (1+rho2)*x_t;
                          }
                          if(grad){
                          res(5 + i) = grd;
                          }
                     }
                     res(4) = res(0) + res(1) + res(3);
                     return(res);
                     }", depends = c("RcppArmadillo"))

precision_banded = function(N, rho, band_length = 1, sparse=T){
  if(length(rho) < band_length+1) rho <- rep(-rho, band_length+1)
  diags <- lapply(0:band_length, function(i) rep(-rho[1+i], N))
  cholV <- Matrix::bandSparse(N, k = c(0:band_length),
                              diagonals = diags, giveCsparse = T)
  diag(cholV) <- 1
  Q <- Matrix::crossprod(cholV)
  Q[N,] <- Q[,N] <- Q[1,N:1]
  Q[1,1] <- Q[N,N] <- 1
  return(Q)
}
library(Matrix)
library(coda)
set.seed(0)
d <- 28
rho <- -0.95
V <- precision_banded(d, rho = rho, band_length = 1)
V[1,1] <- V[2,2]
Vs <- as(V, "sparseMatrix")

xs <- as.matrix(solve(chol(Vs), rnorm(d)))
y <- rpois(d, lambda = as.numeric(exp(xs)))
nmax <- 10^5
x0 <- xs
ref_rate <- 0.05
tmax <- 15
Datann <- list(rho=rho)

## Standard BPS -- A single factor consisting of all variables
Neighbourhoods <- list(c(0))
Factors <- list(c(0:(d-1)))
Datann <- list(rho=Vs[1,2])

theta <- sample(c(1,-1),d, replace = T)
theta <- rnorm(d)
library(ccpdmp)
timing_bfl <- system.time({set.seed(1);bps_global <- ccpdmp::bps(maxTime = tmax,
                                                        trac_coords = c(0,1,d-1),
                                                        rate_f = local_rate,
                                                        factors = Factors,
                                                        local_updates = Neighbourhoods,
                                                        Data = Datann, y = y[1:d],
                                                        nmax = d*3*10^2, tmax = 1,
                                                        x0 = x0[1:d],
                                                        theta0 = theta,
                                                        ref_rate = 1e-19)})
length(bps_global$times)
plot_pdmp(bps_global, pch = '.',coords = c(1,2,3), nsamples = 1e3, inds = 1:1e3)




