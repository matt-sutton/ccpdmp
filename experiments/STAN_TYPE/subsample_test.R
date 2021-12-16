## Stan simple Normal example: ---
library(rstan)
library(ccpdmp)
stanmodelcode <- "
data {
  int<lower=0> N;
  real y[N];
}
parameters {
  real mu;
}
model {
  target += normal_lpdf(mu | 0, 10);
  target += normal_lpdf(y  | mu, 1);
}
"
G <- 5
N <- 1000
y <- rnorm(N*G)
y_split <- matrix(y, ncol = G)
model <- stan_model(model_code = stanmodelcode, model_name = "example")

## A control variate is a guess at where the distribution is most concentrated
control_variate <- mean(y)
dat <- list(N = N*G, y = y);
fit_full <- sampling(model, warmup = 10,
        data = dat, iter = 1e5, chains = 1, verbose = FALSE)
stan_samples <- extract(fit_full, 'mu')
grad_cv <- -1 * grad_log_prob(fit_full, control_variate)

dimension <- 1
fits <- list()
grad_cv_subsample <- matrix(0, dimension, G)
for(g in 1:G){
  dat <- list(N = N, y = y_split[,g]);
  fits[[g]] <- sampling(model, warmup = 0,
                        data = dat, iter = 10, chains = 1, verbose = FALSE)
  grad_cv_subsample[,g] <- -1 * grad_log_prob(fits[[g]], control_variate)
}

## Return the rates
return_rates <- function(x, theta, tau_grid, dnlogpi, rate_updates){

  tau_length <- length(tau_grid);
  n_rates <- length(rate_updates)
  subsample <<- sample(1:G, size = 1) ## Assign global so can be used in grad function

  rates_eval <- matrix(0, n_rates, tau_length)
  C <- 5500

  for(i in 1:n_rates){
    ## Calculate f(t)
    partial_i <- rate_updates[i]

    rates_eval[i,] <- max(0,theta[partial_i]*grad_cv[partial_i]) +
      C*(norm(x - control_variate, "2") + tau_grid*norm(theta,"2"))
  }
  return(rates_eval)
}

dnlogpi <- function(x, partial){
  grad_ss <- -1 * grad_log_prob(fits[[subsample]], x) ## return the gradient of the (subsampled) negative log posterior
  res <- G*(grad_ss[partial] - grad_cv_subsample[partial, subsample]) + grad_cv[partial] # Control variate idea
  return(res)
}
library(ccpdmp)
system.time(zigzag_fit <- zigzag(1e3, dnlogpi, x0 = c(0), poly_order = 1, return_rates = return_rates)) ## Using 1/5 the number of iterations.
samples <- gen_samples(nsample = 1e4, positions = zigzag_fit$positions, times = zigzag_fit$times)
plot(density(samples$xx), col = 1, main = "Density plot for Zig Zag (black) and STAN (red)")
lines(density(stan_samples$mu), col = 2)


