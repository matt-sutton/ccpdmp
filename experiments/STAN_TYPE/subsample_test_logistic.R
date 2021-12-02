## Stan simple Normal example: ---
library(rstan)
library(ccpdmp)
stanmodelcode <- "
data {
  int<lower = 0> N;
  vector[N] x;
  int y[N];
}
parameters {
  real alpha;
  real beta;
}
model {
  alpha ~ normal(0, 2);
  beta ~ normal(0, 2);
  y ~ bernoulli_logit(alpha + x * beta);
}
"
G <- 5
N <- 5000
inv_logit <- function(u) 1 / (1 + exp(-u))
x <- matrix(rnorm(N * G), N, G)
alpha <- 1
beta <- -1
y_split <- matrix(rbinom(G * N, 1, inv_logit(alpha + beta * as.vector(x))), N, G)
y <- as.vector(y_split)
model <- stan_model(model_code = stanmodelcode, model_name = "example")

control_variate <- c(alpha, beta)
dat <- list(N = N*G, x=as.vector(x), y = y);
fit_full <- sampling(model, warmup = 10,
        data = dat, iter = 5e3, chains = 1, verbose = FALSE)
stan_samples <- extract(fit_full)
grad_cv <- -1 * grad_log_prob(fit_full, control_variate)

dimension <- 2
fits <- list()
grad_cv_subsample <- matrix(0, dimension, G)
for(g in 1:G){
  dat <- list(N = N, x = x[,g],y = y_split[,g]);
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
  C <- 5000

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
system.time(zigzag_fit <- zigzag(5e3, dnlogpi, x0 = control_variate, poly_order = 1, return_rates = return_rates)) ## Using 1/5 the number of iterations.
samples <- gen_samples(nsample = 1e4, positions = zigzag_fit$positions, times = zigzag_fit$times)

plot_pdmp(zigzag_fit, coords = c(1,2), burn = .1,inds = 200:2e3,
          mcmc_samples = cbind(stan_samples$alpha, stan_samples$beta), nsamples = 1e4)

fit_full
coda::effectiveSize(samples$xx[1,])
coda::effectiveSize(samples$xx[2,])

lp <- apply(samples$xx, 2, function(s) log_prob(fit_full, s))
coda::effectiveSize(lp)
plot(lp)
plot(stan_samples$lp__)
