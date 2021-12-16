library(RcppXPtrUtils)

grad_cpp <- cppXPtr("arma::vec get_grad(arma::vec x, const List& Data, const arma::vec& y) {
                     const double rho = Data[0];
                     int d = x.size();
                     arma::vec res(d);
                     arma::vec xmy(d);
                     xmy = x - y;
                     double rho2 = std::pow(rho,2.0);

                     res(0) = (1+rho2)*xmy(0) - rho*xmy(1);
                     res(d-1) = xmy(d-1) - rho*xmy(d-2);
                     for(int i =1; i < d-1; i++){
                        res(i) = (1+rho2)*xmy(i) - rho*xmy(i-1) - rho*xmy(i+1);
                     }
                     return(res);
                     }", depends = c("RcppArmadillo"))

post_cpp <- cppXPtr("double get_post(arma::vec x, const List& Data, const arma::vec& y) {
                     const double rho = Data[0];
                     int d = x.size();
                     arma::vec res(d);
                     arma::vec xmy(d);
                     xmy = x - y;
                     double rho2 = std::pow(rho,2.0);

                     res(0) = (1+rho2)*xmy(0) - rho*xmy(1);
                     res(d-1) = xmy(d-1) - rho*xmy(d-2);
                     for(int i =1; i < d-1; i++){
                        res(i) = (1+rho2)*xmy(i) - rho*xmy(i-1) - rho*xmy(i+1);
                     }
                     return(arma::dot(xmy, res)/2.0);
                     }", depends = c("RcppArmadillo"))

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
                     double x_t, x_t_m1, x_t_p1, grd;

                     for(int i = 0; i < num_terms; i++){

                          grd = 0.0;
                          int term = rate_inds(i);
                          x_t = x(term) + theta(term)*(t-t_old(term)) - y(term);

                          // Latent Gaussian terms
                          if(term == d-1){
                             x_t_m1 = x(term-1) + theta(term-1)*(t-t_old(term-1)) - y(term-1);
                             res(0) += theta(term)*(x_t - rho*x_t_m1);
                             grd += x_t - rho*x_t_m1;
                          }
                          if(term == 0){
                             x_t_p1 = x(term+1) + theta(term+1)*(t-t_old(term+1)) - y(term+1);
                             res(0) += theta(term)*((1+rho2)*x_t - rho*x_t_p1);
                             grd += (1+rho2)*x_t - rho*x_t_p1;
                          }
                          if((0 < term) & (term < d-1)){
                            x_t_m1 = x(term-1) + theta(term-1)*(t-t_old(term-1)) -y(term-1);
                            x_t_p1 = x(term+1) + theta(term+1)*(t-t_old(term+1)) -y(term+1);

                            res(0) += theta(term)*((1+rho2)*x_t - rho*x_t_p1- rho*x_t_m1);
                            grd += (1+rho2)*x_t - rho*x_t_p1- rho*x_t_m1;
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
nmax <- 10^5
x0 <- xs
ref_rate <- 0.05
tmax <- 15
Datann <- list(rho=rho)

tmax <- 5
d <- 100
library(ccpdmp)
timinghml <- system.time({set.seed(1);hml <- hmc(maxTime = tmax, post_f = post_cpp,grad_f = grad_cpp,
                                                 Data = Datann,y = y[1:d],
                                                 nmax = 1e4, epsilon = 1*d^(-1/4),
                                                 x0 = x0[1:d], burn = 1)})
factors <- lapply(0:(d-1), function(i) i)
rate_updates <- lapply(0:(d-1), function(i) c(max(0, i-1):min((d-1),i+1)))
timing_bfl <- system.time({set.seed(1);zz <- bps(maxTime = tmax,
                                                 trac_coords = c(0,1,d-2),
                                                 factors = factors,
                                                 local_updates = rate_updates,
                                                 rate_f = local_rate, Data = Datann, y = y[1:d],
                                                 nmax = d*10^3, ref_rate = 1e-16,
                                                 x0 = x0[1:d], theta0 = sample(c(1,-1),d, replace = T))})
plot_pdmp(zz, coords = c(1,2,3), inds = 1:2,
          nsamples = 5000, pch = '.', mcmc_samples = t(hml$samples[c(1,2,d),]))

rate_updates <- list(c(0))
factors <- list(c(0:(d-1)))
timing_b <- system.time({set.seed(r);bps_1 <- bps(maxTime = tmax,
                                                  trac_coords = c(0,1,d-1),
                                                  factors = factors,
                                                  local_updates = rate_updates,
                                                  rate_f = local_rate, Data = Datann, y = y[1:d],
                                                  nmax = d*1e4, ref_rate = 1,
                                                  x0 = x0[1:d], theta0 = rnorm(d))})

nrep <- 40
nit_hml <- time_hml <- ESS_hml <- matrix(0, nrow = nrep, ncol = length(d_vals))
nit_ml <- time_ml <- ESS_ml <- matrix(0, nrow = nrep, ncol = length(d_vals))
nit_zz <- time_zz <- ESS_zz <- matrix(0, nrow = nrep, ncol = length(d_vals))
nit_bps <- time_bps <- ESS_bps <- matrix(0, nrow = nrep, ncol = length(d_vals))
nit_bps_4 <- time_bps_4 <- ESS_bps_4 <- matrix(0, nrow = nrep, ncol = length(d_vals))
nit_bps_8 <- time_bps_8 <- ESS_bps_8 <- matrix(0, nrow = nrep, ncol = length(d_vals))
nit_bps_16 <- time_bps_16 <- ESS_bps_16 <- matrix(0, nrow = nrep, ncol = length(d_vals))
h_vals <- rep(0,length(d_vals))


tune_mcmc <- function(n_eval=10, args, args_tune, opt_acc, method = 'mala', warmstart = FALSE){
  acc_vals <- c()
  for(i in 1:n_eval){
    method_eval <- do.call(method, args)
    acc <- mean(method_eval$acc_probs)
    acc_vals <- c(acc_vals, acc)
    args[[args_tune]] = args[[args_tune]]*(1 + (acc - opt_acc))
    if(warmstart){
      args[["x0"]] = method_eval$samples[,length(method_eval$acc_probs)]
    }
  }
  return(list(acc = acc_vals, args = args))
}

d <- 4000
Datann <- list(rho=rho);
args <- list(maxTime = 7, y = y[1:d],nmax = 10^4, post_f = post_cpp,Data = Datann,
             grad_f = grad_cpp, epsilon = 0.4*d^(-1/3), x0 = x0[1:d], burn = 10^3)
r1 <- tune_mcmc(30, args = args, args_tune = 'epsilon', opt_acc = 0.574, method = 'mala_pnt')
print(r1$acc)
print(r1$args$epsilon*d^(1/3)) ## ell

args <- list(maxTime = 7, y = y[1:d],nmax = 10^4, post_f = post_cpp,Data = Datann,
             grad_f = grad_cpp, epsilon = 1*d^(-1/4), x0 = x0[1:d], burn = 10^2)
r2 <- tune_mcmc(30, args = args, args_tune = 'epsilon', opt_acc = 0.651,
                method = 'hmc', warmstart = T)
print(r2$acc)
print(r2$args$epsilon*d^(1/4))

lmala <- 0.4 #(maybe round to .45)
lhmc <- 0.97
tmax <- 40
d_index <- 4
library(ccpdmp)

for(d_index in 1:length(d_vals)){
  set.seed(1)
  d <- d_vals[d_index]
  d <- 100
  Datann <- list(rho=rho)
  print(paste0(" dim = ", d))
  r <- 1

  for( r in 1:nrep ){

    ## MALA
    timingmala <- system.time({set.seed(r);mala <- mala_pnt(maxTime = tmax,
                                                            post_f = post_cpp,grad_f = grad_cpp,
                                                            Data = Datann, y = y[1:d],
                                                            nmax = floor(6e3*nmax/d), epsilon = lmala*d^(-1/3),
                                                            x0 = x0[1:d], burn = 1)})
    time_ml[r,d_index] <- mala$Etime
    print("mala")
    print(time_ml[r,d_index])
    print(mean(mala$acc_probs))
    print(length(mala$acc_probs))
    ESS_ml[r,d_index] <- effectiveSize(mala$samples[1,])
    nit_ml[r,d_index] <- length(mala$acc_probs)
    print(ESS_ml[r,d_index])
    rm(mala)
    gc()

    ## HMC
    timinghml <- system.time({set.seed(r);hml <- hmc(maxTime = tmax, post_f = post_cpp,grad_f = grad_cpp,
                                                          Data = Datann,y = y[1:d],
                                                          nmax = floor(1e2*nmax/d), epsilon = lhmc*d^(-1/4),
                                                          x0 = x0[1:d], burn = 1)})
    time_hml[r,d_index] <- hml$Etime
    print("hmc")
    print(time_hml[r,d_index])
    print(mean(hml$acc_probs))
    print(length(hml$acc_probs))
    ESS_hml[r,d_index] <- effectiveSize(hml$samples[1,])
    nit_hml[r,d_index] <- length(hml$acc_probs)
    print(ESS_hml[r,d_index])
    rm(hml)
    gc()

    ## Zig-Zag (local)
    factors <- lapply(0:(d-1), function(i) i)
    rate_updates <- lapply(0:(d-1), function(i) c(max(0, i-1):min((d-1),i+1)))
    timing_bfl <- system.time({set.seed(1);zz <- bps(maxTime = tmax,
                                                         trac_coords = c(0,1,d-1),
                                                         factors = factors,
                                                         local_updates = rate_updates,
                                                         rate_f = local_rate, Data = Datann, y = y[1:d],
                                                         nmax = d*10^3, ref_rate = 1e-16,
                                                         x0 = x0[1:d], theta0 = sample(c(1,-1),d, replace = T))})
    plot_pdmp(zz, coords = c(1,2,3), inds = 1:2,
              nsamples = 5000, pch = '.', mcmc_samples = t(hml$samples))

    # rate_updates <- lapply(0:(d-1), function(i) c(max(0, i-1):min((d-1),i+1)))
    # timing_z <- system.time({set.seed(r);zz <- zigzag_cpp(maxTime = tmax,
    #                                                       trac_coords = c(0,1,d-1),
    #                                                       local_updates = rate_updates,
    #                                                       rate_f = local_rate, Data = Datann, y = y[1:d],
    #                                                       nmax = d*1e4, tmax = 1,
    #                                                       x0 = x0[1:d], theta0 = sample(c(1,-1),d, replace = T))})
    time_zz[r,d_index] <- zz$Etime
    print("ZZ")
    print(time_zz[r,d_index])
    nit_zz[r,d_index] <- length(zz$times)
    print(nit_zz[r,d_index])
    samp_zz <- gen_samples(zz$positions,zz$times,nsample = max(10^3,10*nit_hml[r,d_index]))
    ESS_zz[r,d_index] <- effectiveSize(samp_zz$xx[1,])
    print(ESS_zz[r,d_index])
    rm(zz); rm(samp_zsk)
    gc()

    plot_pdmp(zz, coords = c(1,2,3), inds = 1:2,
                       nsamples = 5000, pch = '.', mcmc_samples = t(hml$samples))



    ## BPS (global) factors = 1
    rate_updates <- list(c(0))
    factors <- list(c(0:(d-1)))
    timing_b <- system.time({set.seed(r);bps_1 <- bps(maxTime = tmax,
                                                      trac_coords = c(0,1,d-1),
                                                      factors = factors,
                                                      local_updates = rate_updates,
                                                      rate_f = local_rate, Data = Datann, y = y[1:d],
                                                      nmax = d*1e4, ref_rate = 1,
                                                      x0 = x0[1:d], theta0 = rnorm(d))})
    time_bps[r,d_index] <- bps_1$Etime
    nit_bps[r,d_index] <- length(bps_1$times)
    print("BPS Global")
    print(time_bps[r,d_index])
    print(nit_bps[r,d_index])
    samp_bps <- gen_samples(bps_1$positions,bps_1$times,nsample = max(10^3,10*nit_hml[r,d_index]))
    ESS_bps[r,d_index] <- effectiveSize(samp_bps$xx[1,])
    print(ESS_bps[r,d_index])
    rm(bps_1); rm(samp_bps)
    gc()

    plot_pdmp(bps_1, coords = c(1,2,3), inds = 1:2,
              nsamples = 5000, pch = '.', mcmc_samples = t(hml$samples))


    ## BPS (local) Factors = d/4
    factors <- lapply(0:(d/4-1), function(i) c(4*i, 4*i+1, 4*i+2, 4*i+3))
    rate_updates <- lapply(0:(d/4-1), function(i) c(max(0, i-1):min((d/4-1),i+1)))
    timing_bfl <- system.time({set.seed(1);bps_4 <- bps(maxTime = tmax,
                                                        trac_coords = c(0,1,d-1),
                                                        factors = factors,
                                                        local_updates = rate_updates,
                                                        rate_f = local_rate, Data = Datann, y = y[1:d],
                                                        nmax = d*10^3, ref_rate = .1,
                                                        x0 = x0[1:d], theta0 = rnorm(d))})
    samp_bps_4 <- gen_samples(bps_4$positions,bps_4$times,nsample = max(10^3,10*nit_hml[r,d_index]))
    time_bps_4[r,d_index] <- bps_4$Etime
    nit_bps_4[r,d_index] <- length(bps_4$times)
    print("BPS Loc 4")
    print(time_bps_4[r,d_index])
    print(nit_bps_4[r,d_index])

    samp_bps <- gen_samples(bps_4$positions,bps_4$times,nsample = max(1e3,10*nit_hml[r,d_index]))
    ESS_bps_4[r,d_index] <- effectiveSize(samp_bps$xx[1,])
    print(ESS_bps_4[r,d_index])
    rm(bps_4); rm(samp_bps)
    gc()

    ## BPS (local) Factors = d/8
    factors <- lapply(0:(d/8-1), function(i) 8*i +0:7)
    rate_updates <- lapply(0:(d/8-1), function(i) c(max(0, i-1):min((d/8-1),i+1)))
    timing_bfl <- system.time({set.seed(1);bps_8 <- bps(maxTime = tmax,
                                                        trac_coords = c(0,1,d-1),
                                                        factors = factors,
                                                        local_updates = rate_updates,
                                                        rate_f = local_rate, Data = Datann, y = y[1:d],
                                                        nmax = d*10^3, tmax = 1,ref_rate = .1,
                                                        x0 = x0[1:d], theta0 = rnorm(d))})
    samp_bps_8 <- gen_samples(bps_8$positions,bps_8$times,nsample = max(10^3,10*nit_hml[r,d_index]))
    time_bps_8[r,d_index] <- bps_8$Etime
    nit_bps_8[r,d_index] <- length(bps_8$times)
    print("BPS Loc 8")
    print(time_bps_8[r,d_index])
    print(nit_bps_8[r,d_index])

    samp_bps <- gen_samples(bps_8$positions,bps_8$times,nsample = max(10^3,10*nit_hml[r,d_index]))
    ESS_bps_8[r,d_index] <- effectiveSize(samp_bps$xx[1,])
    print(ESS_bps_8[r,d_index])
    rm(bps_8); rm(samp_bps)
    gc()

    plot_pdmp_multiple(list(zz,bps_8), coords = c(1,2,3), inds = 1:2,
                       nsamples = 5000, pch = '.', mcmc_samples = t(hml$samples))

    ## BPS (local) Factors = d/16
    factors <- lapply(0:(d/16-1), function(i) 16*i +0:15)
    rate_updates <- lapply(0:(d/16-1), function(i) c(max(0, i-1):min((d/16-1),i+1)))
    timing_bfl <- system.time({set.seed(1);bps_16 <- bps(maxTime = tmax,
                                                        trac_coords = c(0,1,d-1),
                                                        factors = factors,
                                                        local_updates = rate_updates,
                                                        rate_f = local_rate, Data = Datann, y = y[1:d],
                                                        nmax = d*10^3, tmax = 1,ref_rate = .1,
                                                        x0 = x0[1:d], theta0 = rnorm(d))})
    samp_bps_16 <- gen_samples(bps_16$positions,bps_16$times,nsample = max(10^3,10*nit_hml[r,d_index]))
    time_bps_16[r,d_index] <- bps_16$Etime
    nit_bps_16[r,d_index] <- length(bps_16$times)
    print("BPS Loc 16")
    print(time_bps_16[r,d_index])
    print(nit_bps_16[r,d_index])

    samp_bps <- gen_samples(bps_16$positions,bps_16$times,nsample = max(10^3,10*nit_hml[r,d_index]))
    ESS_bps_16[r,d_index] <- effectiveSize(samp_bps$xx[1,])
    print(ESS_bps_16[r,d_index])
    rm(bps_16); rm(samp_bps)
    gc()


    plot_pdmp_multiple(list(zz,bps_16), coords = c(1,2,3), inds = 1:2,
                       nsamples = 3000, pch = '.',
                       mcmc_samples = t(hml$samples[c(1,2,d),]) )



    # #rjpdmp::plot_pdmp(zsk, coords = c(1,2,3), inds = 1:2, nsamples = 2*nit_hml[r,d_index],mcmc_samples = t(hml$samples)[,c(1,2,d)] )


    rjpdmp::plot_pdmp(zsk, coords = c(1,2,3), inds = 1:2, nsamples = 2*length(hml$acc_probs),
                      mcmc_samples = sampSTAN$x[,c(1,2,d)], pch = '.' )

    # mls <- t(mala$samples[c(1,2,d),seq(from = 1,to = floor(length(mala$acc_probs)/length(hml$acc_probs))*length(hml$acc_probs), length.out = length(hml$acc_probs))])
    # rjpdmp::plot_pdmp(zsk, coords = c(2,3), inds = 1:2, nsamples = length(hml$acc_probs),
    #                   mcmc_samples = mls, pch = '.' )

    # timing_b <- system.time({set.seed(2);bps <- bps_cave_vex_cpp(maxTime = tmax,rate_f = ratez_cpp,
    #                                                               grad_f = gradn_cpp, Data = Datan,
    #                                                               s_eval = c(0,1),
    #                                                               standardise = T, theta0 = rnorm(d),
    #                                                               nmax = nmax, x0 = x0[1:d], ref = 2*pi/d)})
    # time_bps[r,d_index] <- timing_b[3]
    # nit_bps[r,d_index] <- length(bps$times)
    # samp_bps <- rjpdmp::gen_sample(bps$positions,bps$times,nsample = nit_hml[r,d_index], theta = bps$theta)
    # ESS_bps[r,d_index] <- effectiveSize(samp_bps$x[1,])
    # print("bps")
    # print(length(bps$times))
    # print(ESS_bps[r,d_index])
    # rm(bps); rm(samp_bps)
    # gc()
    #
    # rjpdmp::plot_pdmp(bps, coords = c(1,2,d), inds = 1:2, nsamples = nit_hml[r,d_index],mcmc_samples = t(hml$samples) )
  }
}


