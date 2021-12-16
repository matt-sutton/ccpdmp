library(RcppXPtrUtils)

local_rate <- cppXPtr("arma::vec get_rate(double t, arma::vec& t_old, arma::vec& x, arma::vec& theta, const arma::vec& y, const List& Data, arma::uvec rate_inds, bool grad) {

                     // Standard header //
                     int num_terms = rate_inds.size(), rsize = 5;
                     double x_t, grd;
                     if(grad){
                        rsize += num_terms;
                     }
                     arma::vec res = arma::zeros(rsize);

                     // Evaluate the rate: (f_u, f_n, f_n', p, f) for the factor
                     // Return additional gradient terms if grad = TRUE

                     int d = x.size();
                     const double rho = Data[0];
                     double rho2 = std::pow(rho,2.0);

                     for(int i = 0; i < num_terms; i++){
                          int term = rate_inds(i);
                          grd = 0.0;
                          x_t = x(term) + theta(term)*(t-t_old(term));

                          // Likelihood terms theta*(exp(x) - y)
                          res(0) -= theta(term)*y(term);
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
                                 res(0) += theta(term)*(1+rho2)*x_t;
                                 grd += (1+rho2)*x_t;
                              } else {
                                 res(0) -= theta(term)*rho*x_t;
                                 grd -= rho*x_t;
                              }
                          }
                          if(term == d-1){
                             x_t = x(term) + theta(term)*(t-t_old(term));
                             res(0) += theta(term)*x_t - theta(term)*(1+rho2)*x_t;
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
d_vals <- c(2^c(15:11))
rho <- -0.5
d <- max(d_vals)
V <- precision_banded(d, rho = rho, band_length = 1)
V[1,1] <- V[2,2]
Vs <- as(V, "sparseMatrix")

xs <- as.matrix(solve(chol(Vs), rnorm(d)))
y <- rpois(d, lambda = as.numeric(exp(xs)))

x0 <- xs
ref_rate <- 0.05
Datann <- list(rho=rho)

ref_rate <- 1
lhmc <- 0.97
tmax <- .5
d <- 3
nmax <- 1e3
library(ccpdmp)

# timinghml <- system.time({set.seed(1);hml <- hmc(maxTime = tmax, post_f = post_cpp,grad_f = grad_cpp,
#                                                  Data = Datann,y = y[1:d],
#                                                  nmax = nmax, epsilon = lhmc*d^(-1/4),
#                                                  x0 = x0[1:d], burn = 1)})

# rate_updates <- lapply(0:(d-1), function(i) c(max(0, i-1):min((d-1),i+1)))
d <- 500
nm <- 1e5
rate_updates <- lapply(0:(d-1), function(i) c(max(0, i-1):min((d-1),i+1)))
# rate_updates <- lapply(0:(d-1), function(i) c(i))
#rate_updates <- lapply(0:(d-1), function(i) c(0:d-1))
timing_z <- system.time({set.seed(1);zz <- zigzag_cpp(maxTime = 20,
                                                      trac_coords = c(0,1,d-1),
                                                      local_updates = rate_updates,
                                                      rate_f = local_rate, Data = Datann, y = y[1:d],
                                                      nmax = nm, tmax = 10,
                                                      x0 = x0[1:d], theta0 = sample(c(1,-1),d, replace = T))})

length(zz$times)
timing_z <- system.time({set.seed(1);zz_c <- zigzag_cpp(maxTime = 20,
                                                      trac_coords = c(0,1,d-1),
                                                      local_updates = rate_updates,
                                                      rate_f = local_rate, Data = Datann, y = y[1:d],
                                                      nmax = nm, tmax = .1,
                                                      x0 = x0[1:d], theta0 = sample(c(1,-1),d, replace = T))})
length(zz_c$times)

plot_pdmp_multiple(list(zz, zz_c), coords = c(1,2,3), inds = 1:500, pch = '.', nsamples = 1e4)
zz_c$times
plot(zz$times, zz$positions, type = 'l')
lines(zz_c$times, zz_c$positions, col = 'red')


return_rates <- function(x, theta, tau_grid, y){
  rates_eval <- matrix(0, nrow = 5, ncol = length(tau_grid))
  for( i in 1:length(tau_grid) ){
    f_u <- f_n <- f_n_d <- 0
    x_g <- x+tau_grid[i]*theta

    for(d in 1:length(x)){
      f_u <- f_u + theta[d]*(x_g[d] - y[d])
      if(theta[d] < 0){
        f_n <- f_n + theta[d]*exp(x_g[d])
        f_n_d <- f_n_d + theta[d]^2*exp(x_g[d])
      } else {
        f_u <- f_u + theta[d]*exp(x_g[d])
      }
    }
    rates_eval[,i] <- c(f_u,f_n,f_n_d,0, f_u + f_n)
  }
  return(rates_eval)
}
bps_cc <- function(nEvents,lambda_max = 10, x0, v0, y, tau_max, tau_length, poly_order=2, verb = FALSE, adapt_tau =T){

  t = 0; eps = 1e-10;
  x = x0;
  nvel <- 1#length(x)

  # Velo set
  taus = rep(Inf, nvel)
  u_s = rep(Inf, nvel)
  f_s = rep(Inf, nvel)
  theta = v0;

  tau_grid = seq(from = 0, to = tau_max, length.out = tau_length)
  rates_eval <- return_rates(x, theta, tau_grid, y)
  tus = ccpdmp::sim_rates(eval_times = tau_grid, eval_rates = rates_eval, poly_order)
  taus = tus$t;  u_s = tus$u;  f_s = tus$f_evall

  samples = x;   ts = 0;   v_samples = theta
  v = 0;  j = 0;  num_evts = 0; nE = 1
  n_p <- 0; n_total <- 0
  n_props <- rep(0, nEvents)

  while(nE < nEvents){
    n_total <- n_total + 1
    mini = which.min(taus)[1]
    x = x + taus[mini]*theta
    t = t + taus[mini]

    if(u_s[mini] < 1e-9){
      n_p <- n_p + 1
      grad <- (x + exp(x) - y)

      acc_prb <-  sum(theta*grad)/f_s[mini]
      if(acc_prb > 1){
        print("oops? ");        print(c(sum(theta*grad), f_s[mini]))
        print(nE)
      }
      if(runif(1) <= acc_prb){
        theta = theta - 2*sum(theta*grad)/sum(grad*grad)*grad
        samples = cbind(samples, x)
        v_samples = cbind(v_samples, theta)
        ts = c(ts, t)
        nE = nE + 1
        if(nE %% 1e2 ==0){
          if(verb){print(nE) }
          if( adapt_tau ){
            tau_max <- quantile(diff(ts), 0.9)[[1]]
            tau_grid = seq(from = 0, to = tau_max, length.out = tau_length)
          }
        }
      }
    }else{
      # print("extend the timelength")
    }
    # # Simulate new time
    rates_eval <- return_rates(x, theta, tau_grid, y)
    tus = ccpdmp::sim_rates(eval_times = tau_grid, eval_rates = rates_eval, poly_order)
    taus = tus$t;  u_s = tus$u;  f_s = tus$f_evall
    n_props[nE] <- n_p
  }
  return (list(xx=samples,v_samples=v_samples,tt=ts, n_prop = n_props, n_total = n_total))
}

b <- bps_cc(nEvents = 1e3, y = y[1:d],x0 = x0[1:d], v0=1, tau_max = 10, tau_length = 2, adapt_tau = F)
lines(b$tt, b$xx, col = 'green')

# nplot <- 1e3
# plot_pdmp(zz, coords = c(1,2,3), inds = 1:2,
#           nsamples = 2*nplot, pch = '.', mcmc_samples = t(hml$samples[c(1,2,d),1:nplot]))
