library(RcppXPtrUtils)
grad_cpp <- cppXPtr("arma::vec get_grad(arma::vec x, const List& Data, const arma::vec& y) {
                     const double rho = Data[0];
                     int d = x.size();
                     arma::vec res(d);
                     double rho2 = std::pow(rho,2.0);

                     res(0) = (1+rho2)*x(0) - rho*x(1) + exp(x(0)) - y(0);
                     res(d-1) = x(d-1) - rho*x(d-2) + exp(x(d-1)) - y(d-1);
                     for(int i =1; i < d-1; i++){
                        res(i) = (1+rho2)*x(i) - rho*x(i-1) - rho*x(i+1) + exp(x(i)) - y(i);
                     }
                     return(res);
                     }", depends = c("RcppArmadillo"))

post_cpp <- cppXPtr("double get_post(arma::vec x, const List& Data, const arma::vec& y) {
                     const double rho = Data[0];
                     int d = x.size();
                     arma::vec res(d);
                     double rho2 = std::pow(rho,2.0);

                     res(0) = (1+rho2)*x(0) - rho*x(1);
                     res(d-1) = x(d-1) - rho*x(d-2);
                     for(int i =1; i < d-1; i++){
                        res(i) = (1+rho2)*x(i) - rho*x(i-1) - rho*x(i+1);
                     }
                     return(arma::dot(x, res)/2.0 + arma::sum(exp(x) - y%x));
                     }", depends = c("RcppArmadillo"))

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
  diags <- lapply(0:band_length, function(i) rep(rho[1+i], N))
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
rho <- 0.5
d <- max(d_vals)
V <- precision_banded(d, rho = rho, band_length = 1)
V[1,1] <- V[2,2]
Vs <- as(V, "sparseMatrix")

xs <- as.matrix(solve(chol(Vs), rnorm(d)))
y <- rpois(d, lambda = as.numeric(exp(xs)))
nmax <- 10^5
x0 <- xs
Datann <- list(rho=rho)


nrep <- 50
nit_hml <- time_hml <- ESS_hml <- matrix(0, nrow = nrep, ncol = length(d_vals))
nit_ml <- time_ml <- ESS_ml <- matrix(0, nrow = nrep, ncol = length(d_vals))
nit_zz <- time_zz <- ESS_zz <- matrix(0, nrow = nrep, ncol = length(d_vals))
nit_bps <- time_bps <- ESS_bps <- matrix(0, nrow = nrep, ncol = length(d_vals))
nit_bps_4 <- time_bps_4 <- ESS_bps_4 <- matrix(0, nrow = nrep, ncol = length(d_vals))
nit_bps_8 <- time_bps_8 <- ESS_bps_8 <- matrix(0, nrow = nrep, ncol = length(d_vals))
nit_bps_16 <- time_bps_16 <- ESS_bps_16 <- matrix(0, nrow = nrep, ncol = length(d_vals))
h_vals <- rep(0,length(d_vals))


# tune_mcmc <- function(n_eval=10, args, args_tune, opt_acc, method = 'mala', warmstart = FALSE){
#   acc_vals <- c()
#   for(i in 1:n_eval){
#     method_eval <- do.call(method, args)
#     acc <- mean(method_eval$acc_probs)
#     acc_vals <- c(acc_vals, acc)
#     args[[args_tune]] = args[[args_tune]]*(1 + (acc - opt_acc))
#     if(warmstart){
#       args[["x0"]] = method_eval$samples[,length(method_eval$acc_probs)]
#     }
#   }
#   return(list(acc = acc_vals, args = args))
# }
#
# d <- 4000
# Datann <- list(rho=rho);
# args <- list(maxTime = 7, y = y[1:d],nmax = 10^4, post_f = post_cpp,Data = Datann,
#              grad_f = grad_cpp, epsilon = 0.4*d^(-1/3), x0 = x0[1:d], burn = 10^3)
# r1 <- tune_mcmc(30, args = args, args_tune = 'epsilon', opt_acc = 0.574, method = 'mala_pnt')
# print(r1$acc)
# print(r1$args$epsilon*d^(1/3)) ## ell
#
# args <- list(maxTime = 7, y = y[1:d],nmax = 10^4, post_f = post_cpp,Data = Datann,
#              grad_f = grad_cpp, epsilon = 1*d^(-1/4), x0 = x0[1:d], burn = 10^2)
# r2 <- tune_mcmc(30, args = args, args_tune = 'epsilon', opt_acc = 0.651,
#                 method = 'hmc', warmstart = T)
# print(r2$acc)
# print(r2$args$epsilon*d^(1/4))


# Parameters found from commented section above
ref_rate <- 0.5
lmala <- 0.4
lhmc <- 0.98
tmax <- 90
library(ccpdmp)

for(d_index in 1:length(d_vals)){
  set.seed(1)
  d <- d_vals[d_index]
  Datann <- list(rho=rho)
  print(paste0(" dim = ", d))
  r <- 1
  plot_it <- T

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
                                                     nmax = floor(1e3*nmax/d), epsilon = lhmc*d^(-1/4),
                                                     x0 = x0[1:d], burn = 1)})
    time_hml[r,d_index] <- hml$Etime
    print("hmc")
    print(time_hml[r,d_index])
    print(mean(hml$acc_probs))
    print(length(hml$acc_probs))
    ESS_hml[r,d_index] <- effectiveSize(hml$samples[1,])
    nit_hml[r,d_index] <- length(hml$acc_probs)
    print(ESS_hml[r,d_index])
    nplot <- length(hml$acc_probs)
    if(!plot_it){
      rm(hml)
      gc()
    }

    ## Zig-Zag (local)
    rate_updates <- lapply(0:(d-1), function(i) c(max(0, i-1):min((d-1),i+1)))
    timing_z <- system.time({set.seed(r);zz <- zigzag_cpp(maxTime = tmax,
                                                          trac_coords = c(0,1,d-1),
                                                          local_updates = rate_updates,
                                                          rate_f = local_rate, Data = Datann, y = y[1:d],
                                                          nmax = d*1e4, tmax = 1,
                                                          x0 = x0[1:d], theta0 = sample(c(1,-1),d, replace = T))})
    time_zz[r,d_index] <- zz$Etime
    print("ZZ")
    print(time_zz[r,d_index])
    nit_zz[r,d_index] <- length(zz$times)
    print(nit_zz[r,d_index])
    samp_zz <- gen_samples(zz$positions,zz$times,nsample = max(10^3,nit_hml[r,d_index]))
    ESS_zz[r,d_index] <- effectiveSize(samp_zz$xx[1,])
    print(ESS_zz[r,d_index])

    if(plot_it){
      plot_pdmp(zz, coords = c(1,2,3), inds = 1:2,
                nsamples = nplot, pch = '.', mcmc_samples = t(hml$samples[c(1,2,d),1:nplot]))
      title(paste("d=",d, "ZZ"))
    }
    rm(zz); rm(samp_zz)
    gc()

    ## BPS (global) factors = 1
    rate_updates <- list(c(0))
    factors <- list(c(0:(d-1)))
    timing_b <- system.time({set.seed(r);bps_1 <- bps(maxTime = tmax,
                                                      trac_coords = c(0,1,d-1),
                                                      factors = factors,
                                                      local_updates = rate_updates,
                                                      rate_f = local_rate, Data = Datann, y = y[1:d],
                                                      nmax = d*10^3, tmax = .1, ref_rate = ref_rate,
                                                      x0 = x0[1:d], theta0 = rnorm(d))})
    time_bps[r,d_index] <- bps_1$Etime
    nit_bps[r,d_index] <- length(bps_1$times)
    print("BPS Global")
    print(time_bps[r,d_index])
    print(nit_bps[r,d_index])
    samp_bps <- gen_samples(bps_1$positions,bps_1$times,nsample = max(10^3,10*nit_hml[r,d_index]))
    ESS_bps[r,d_index] <- effectiveSize(samp_bps$xx[1,])
    print(ESS_bps[r,d_index])
    if(plot_it){
      plot_pdmp(bps_1, coords = c(1,2,3), inds = 1:2,
                nsamples = nplot, pch = '.', mcmc_samples = t(hml$samples[c(1,2,d),1:nplot]))
      title(paste("d=",d, "BPS"))
    }
    rm(bps_1); rm(samp_bps)
    gc()

    ## BPS (local) Factors = d/4
    factors <- lapply(0:(d/4-1), function(i) 4*i + 0:3)
    rate_updates <- lapply(0:(d/4-1), function(i) c(max(0, i-1):min((d/4-1),i+1)))
    ref_rate_local <- ref_rate*(4-1)/(d-1)
    timing_bfl <- system.time({set.seed(1);bps_4 <- bps(maxTime = tmax,
                                                        trac_coords = c(0,1,d-1),
                                                        factors = factors, tmax = 1,
                                                        local_updates = rate_updates,
                                                        rate_f = local_rate, Data = Datann, y = y[1:d],
                                                        nmax = d*10^3, ref_rate = ref_rate_local,
                                                        x0 = x0[1:d], theta0 = rnorm(d))})

    samp_bps_4 <- gen_samples(bps_4$positions,bps_4$times,nsample = max(10^3,10*nit_hml[r,d_index]))
    time_bps_4[r,d_index] <- bps_4$Etime
    nit_bps_4[r,d_index] <- length(bps_4$times)
    print("BPS Loc 4")
    print(time_bps_4[r,d_index])
    print(nit_bps_4[r,d_index])

    ESS_bps_4[r,d_index] <- effectiveSize(samp_bps_4$xx[1,])
    print(ESS_bps_4[r,d_index])
    if(plot_it){
      plot_pdmp(bps_4, coords = c(1,2,3), inds = 1:2,
                nsamples = nplot, pch = '.', mcmc_samples = t(hml$samples[c(1,2,d),1:nplot]))
      title(paste("d=",d, "BPS (4)"))
    }
    rm(bps_4); rm(samp_bps_4)
    gc()

    ## BPS (local) Factors = d/8
    factors <- lapply(0:(d/8-1), function(i) 8*i +0:7)
    rate_updates <- lapply(0:(d/8-1), function(i) c(max(0, i-1):min((d/8-1),i+1)))
    ref_rate_local <- ref_rate*(8-1)/(d-1)
    timing_bfl <- system.time({set.seed(1);bps_8 <- bps(maxTime = tmax,
                                                        trac_coords = c(0,1,d-1),
                                                        factors = factors,
                                                        local_updates = rate_updates,
                                                        rate_f = local_rate, Data = Datann, y = y[1:d],
                                                        nmax = d*10^3, tmax = 1, ref_rate = ref_rate_local,
                                                        x0 = x0[1:d], theta0 = rnorm(d))})
    samp_bps_8 <- gen_samples(bps_8$positions,bps_8$times,nsample = max(10^3,10*nit_hml[r,d_index]))
    time_bps_8[r,d_index] <- bps_8$Etime
    nit_bps_8[r,d_index] <- length(bps_8$times)
    print("BPS Loc 8")
    print(time_bps_8[r,d_index])
    print(nit_bps_8[r,d_index])
    ESS_bps_8[r,d_index] <- effectiveSize(samp_bps_8$xx[1,])
    print(ESS_bps_8[r,d_index])
    if(plot_it){
      plot_pdmp(bps_8, coords = c(1,2,3), inds = 1:2,
                nsamples = nplot, pch = '.', mcmc_samples = t(hml$samples[c(1,2,d),1:nplot]))
      title(paste("d=",d, "BPS (8)"))
    }
    rm(bps_8); rm(samp_bps_8)
    gc()

    ## BPS (local) Factors = d/16
    factors <- lapply(0:(d/16-1), function(i) 16*i +0:15)
    rate_updates <- lapply(0:(d/16-1), function(i) c(max(0, i-1):min((d/16-1),i+1)))
    ref_rate_local <- ref_rate*(16-1)/(d-1)
    timing_bfl <- system.time({set.seed(1);bps_16 <- bps(maxTime = tmax,
                                                         trac_coords = c(0,1,d-1),
                                                         factors = factors,
                                                         local_updates = rate_updates,
                                                         rate_f = local_rate, Data = Datann, y = y[1:d],
                                                         nmax = d*10^3, tmax = 1,ref_rate = ref_rate_local,
                                                         x0 = x0[1:d], theta0 = rnorm(d))})
    samp_bps_16 <- gen_samples(bps_16$positions,bps_16$times,nsample = max(10^3,10*nit_hml[r,d_index]))
    time_bps_16[r,d_index] <- bps_16$Etime
    nit_bps_16[r,d_index] <- length(bps_16$times)
    print("BPS Loc 16")
    print(time_bps_16[r,d_index])
    print(nit_bps_16[r,d_index])

    ESS_bps_16[r,d_index] <- effectiveSize(samp_bps_16$xx[1,])
    print(ESS_bps_16[r,d_index])
    if(plot_it){
      plot_pdmp(bps_16, coords = c(1,2,3), inds = 1:2,
                nsamples = nplot, pch = '.', mcmc_samples = t(hml$samples[c(1,2,d),1:nplot]))
      title(paste("d=",d, "BPS (16)"))
    }
    rm(bps_16); rm(samp_bps_16)
    gc()

    save(time_ml,time_hml,time_zz,time_bps,time_bps_4,time_bps_8,time_bps_16,
         ESS_ml,ESS_hml,ESS_zz,ESS_bps,ESS_bps_4,ESS_bps_8,ESS_bps_16,
         nit_ml,nit_hml,nit_zz,nit_bps,nit_bps_4,nit_bps_8,nit_bps_16, file = "MALA_HMC_Comp.Rdata")
    plot_it <- F
  }
}

inds <- 1:(length(d_vals))
types <- c("Computational Efficiency", "Statistical Efficiency", "Overall Efficiency")
text_width <- 1.2
par(mfrow = c(1,3),cex.main=text_width, cex.lab=text_width, cex.axis=text_width, lwd = text_width, mar = c(4,4,2,1.2))
tv <- 1
for ( tv in 1:3 ){
  type <- types[tv]
  if(type == "Overall Efficiency"){
    # ESS per Sec
    calc_ml <- ESS_ml/time_ml
    calc_hml <- ESS_hml/time_hml
    calc_zz <- ESS_zz/time_zz
    calc_zl <- ESS_bps_4/time_bps_4
    calc_zml <- ESS_bps_8/time_bps_8
    calc_zeml <- ESS_bps_16/time_bps_16
    calc_zb <- ESS_bps/time_bps

  } else if(type == "Computational Efficiency"){
    # Iteration per Sec
    calc_hml <- nit_hml/time_hml
    calc_ml <- nit_ml/time_ml
    calc_zz <- nit_zz/time_zz
    calc_zl <- nit_bps_4/time_bps_4
    calc_zml <- nit_bps_8/time_bps_8
    calc_zeml <- nit_bps_16/time_bps_16
    calc_zb <- nit_bps/time_bps
  } else if(type == "Statistical Efficiency"){
    #  ESS per iteration
    calc_ml <- ESS_ml/nit_ml
    calc_hml <- ESS_hml/nit_hml
    calc_zz <- ESS_zz/nit_zz
    calc_zl <- ESS_bps_4/nit_bps_4
    calc_zml <- ESS_bps_8/nit_bps_8
    calc_zeml <- ESS_bps_16/nit_bps_16
    calc_zb <- ESS_bps/nit_bps
  }

  calcs <- c(calc_hml,calc_ml,calc_zb,calc_zz,calc_zl,calc_zml,calc_zeml )
  ylim <- range(log2(calcs)[log2(calcs) <10^5], na.rm = T)

  calc_mat <- rbind(colMeans(calc_hml),
                    colMeans(calc_ml),
                    colMeans(calc_zb),
                    colMeans(calc_zz),
                    colMeans(calc_zl),
                    colMeans(calc_zml),
                    colMeans(calc_zeml))
  sds_mat <- rbind(apply(calc_hml, 2, sd),
                   apply(calc_ml, 2, sd),
                   apply(calc_zb, 2, sd),
                   apply(calc_zz, 2, sd),
                   apply(calc_zl, 2, sd),
                   apply(calc_zml, 2, sd),
                   apply(calc_zeml, 2, sd))
  methods <- c("hmc","mala","bps","zigzag","local (4)","local (8)","local (16)")

  require("Hmisc")
  #par(mfrow = c(1,1))
  leg_label <- NULL
  print(paste("type: ", type))
  for( m in 1:length(methods)){
    errbar(log2(d_vals[inds]), log(calc_mat[m,inds],2), log(calc_mat[m,inds]+sds_mat[m,inds],2),
           log(calc_mat[m,inds]-sds_mat[m,inds],2), add = (m > 1),ylim = round(ylim,2), col = m, errbar.col = m,
           pch=20, cap=.015, ann = FALSE, ylab = type, sub = "???")
    X <- matrix(NaN,nrow=length(inds),ncol=2)
    X[,1] <- rep(1,length(inds))
    X[,2] <- log2(d_vals[inds])
    y <- log(calc_mat[m,],2)

    beta_method <- solve(t(X)%*%X)%*%t(X)%*%y
    abline(a = beta_method[1], b = beta_method[2], col = m, lwd = 1.5)
    leg_label <- c(leg_label, paste0(methods[m],": ", round(beta_method[1],2), ' ', round(beta_method[2],2), 'd'))
    print(paste0(methods[m],": ", round(beta_method[1],2), ' ', round(beta_method[2],2), 't'))
  }

  df <- data.frame(calc_mat)

  title(type, xlab = 'Dimension log 2', ylab = type)
  legend('bottomleft', leg_label, col = 1:m, lwd = 1)
}
