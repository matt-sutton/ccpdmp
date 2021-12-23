library(RcppXPtrUtils)
local_rate <- cppXPtr("arma::vec get_rate(double t, arma::vec& t_old, arma::vec& x, arma::vec& theta,
                      const arma::vec& y, const List& Data, arma::uvec rate_inds, bool grad) {

                     // Standard header //
                     int num_terms = rate_inds.size(), rsize = 5;
                     double x_t, grd;
                     if(grad){
                        rsize += num_terms;
                     }
                     arma::vec res = arma::zeros(rsize);

                     // Evaluate the rate: (f_u, f_n, f_n', p, f) for the factor
                     // Return additional gradient terms if grad = TRUE

                     for(int i = 0; i < num_terms; i++){
                          int partial = rate_inds(i);
                          x_t = x(partial) + theta(partial)*(t-t_old(partial));
                          grd = x_t;
                          res(0) += theta(partial)*grd;

                          if(grad){
                          res(5 + i) = grd;
                          }
                     }
                     res(4) = res(0) + res(1) + res(3);
                     return(res);
                     }", depends = c("RcppArmadillo"))

library(Matrix)
library(coda)
set.seed(0)
d <- 1000
x0<-xs <- rnorm(d)
y <- c(1)
nmax <- 10^5
ref_rate <- 0.05
tmax <- 5
Datann <- list(rho=0)

## Standard BPS -- A single factor consisting of all variables
Neighbourhoods <- list(c(0))
Factors <- list(c(0:(d-1)))

theta <- rnorm(d)
library(ccpdmp)
system.time({set.seed(1);bps_global <- bps(maxTime = tmax,
                                                        trac_coords = c(0,1,d-1),
                                                        rate_f = local_rate,
                                                        factors = Factors,
                                                        local_updates = Neighbourhoods,
                                                        Data = Datann, y = y,
                                                        nmax = d*3*10^2, tmax = 1,
                                                        x0 = x0[1:d], poly_order = 0,
                                                        theta0 = theta,
                                                        ref_rate = 1)})
length(bps_global$times)
plot_pdmp(bps_global, pch = '.',coords = c(1,2,3), nsamples = 1e3, inds = 1:1e3)

## local BPS -- 10 factors consisting of d/10 variables each
Factors <- lapply(0:(d/10-1), function(i) seq(from = 10*i, to = 10*i+9))
Neighbourhoods <- lapply(0:(d/10-1), function(i) c(i))

system.time({set.seed(1);bps_local <- bps(maxTime = tmax,
                                                         trac_coords = c(0,1,d-1),
                                                         rate_f = local_rate,
                                                         factors = Factors,
                                                         local_updates = Neighbourhoods,
                                                         Data = Datann, y = y[1:d],
                                                         nmax = d*3*10^2, tmax = 1,
                                                         x0 = x0[1:d], poly_order = 0,
                                                         theta0 = theta,
                                                         ref_rate = 1)})
length(bps_local$times)
plot_pdmp(bps_local, pch = '.',coords = c(1,2,3), nsamples = 1e3, inds = 1:1e3)


## Standard ZigZag -- A factor for each variable (no assumption on sparsity)
Neighbourhoods <- lapply(0:(d-1), function(i) c(0:(d-1)))
Factors <- lapply(0:(d-1), function(i) c(i))
theta <- sample(c(1,-1),d, replace = T)
system.time({set.seed(1);zigzag_global <- zigzag_cpp(maxTime = tmax,
                                                         trac_coords = c(0,1,d-1),
                                                         rate_f = local_rate,
                                                         factors = Factors,
                                                         local_updates = Neighbourhoods,
                                                         Data = Datann, y = y[1:d],
                                                         nmax = d*3*10^2, tmax = 1,
                                                         x0 = x0[1:d], poly_order = 0,
                                                         theta0 = theta)})
length(zigzag_global$times)
plot_pdmp(zigzag_global, pch = '.',coords = c(1,2,3), nsamples = 1e3, inds = 1:1e3)

## Using full sparsity knowledge:
Neighbourhoods <- lapply(0:(d-1), function(i) c(i))
system.time({set.seed(1);zigzag_local <- zigzag_cpp(maxTime = tmax,
                                                     trac_coords = c(0,1,d-1),
                                                     rate_f = local_rate,
                                                     factors = Factors,
                                                     local_updates = Neighbourhoods,
                                                     Data = Datann, y = y[1:d],
                                                     nmax = d*3*10^2, tmax = 1,
                                                     x0 = x0[1:d], poly_order = 0,
                                                     theta0 = theta)})
length(zigzag_local$times)
plot_pdmp(zigzag_local, pch = '.',coords = c(1,2,3), nsamples = 1e3, inds = 1:1e3)




