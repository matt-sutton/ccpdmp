## Example Code for logistic regression
## More exotic GLMs can be implemented with different phi and f
##
library(ccpdmp)
generate.logistic.data <- function(beta, n.obs, siginv) {
  p <- length(beta)
  dataX <- mvtnorm::rmvnorm(n=n.obs,sigma = solve(siginv))
  vals <- dataX%*%beta
  generateY <- function(p) { rbinom(1, 1, p)}
  dataY <- sapply(1/(1 + exp(-vals)), generateY)
  return(list(dataX = dataX, dataY = dataY))
}

set.seed(1)
n_ev <- 1e3
p <- 5
n <- 1000
sigma2_prior <- 1e10
beta <- c(c(-1.25, 0.5), rep(-0.4,p-2))
siginv <- diag(1, p,p)
siginv[1,2] <- siginv[2,1] <- 0.9
data <- generate.logistic.data(beta, n, siginv)
X <- data$dataX
y <- data$dataY

get_grad <- function(x, index, subsample = TRUE){
  if(subsample){
    ss_index <- sample(1:n, size = 1)
    expa <- exp(sum(X[ss_index,]*x))
    phi_1 <- expa/(1+expa) - y[ss_index]

    expa_star <- exp(sum(X[ss_index,]*control_variate))
    phi_1_star <- expa_star/(1+expa_star) - y[ss_index]

    grad <- rep(0, length(index))
    grad_cv_ss <- rep(0, length(index))
    for(i in 1:length(index)){
      grad[i] <- phi_1*X[ss_index,index[i]]
      grad_cv_ss[i] <- phi_1_star*X[ss_index,index[i]]
    }
    res <- n*(grad - grad_cv_ss) + grad_cv[index]# Control variate idea
    res <- res + x[index]/sigma2_prior # Add in the prior term
    return(res)

  } else {
    # Grad likelihood with prior removed for subsampling
    expa <- exp(X%*%x)
    phi_1 <- expa/(1+expa) - y

    grad <- rep(0, length(index))
    for(i in 1:length(index)){
      grad[i] <- sum(phi_1*X[,index[i]])
    }
    return(grad)
  }
}

control_variate <- beta
grad_cv <- get_grad(beta, 1:length(beta), subsample = FALSE)

C_i_vals <- rep(0, p)
for( i in 1:p){
  C_i_vals[i] <- n*0.25*max(apply(X, 1, function(s) norm(s,"2"))*abs(X[,i]))
}
C_i_vals

## Return the rates for an upper-bound on all subsampled rates
return_rates_subsample <- function(x, theta, tau_grid, dnlogpi, rate_updates){

  tau_length <- length(tau_grid);
  n_rates <- length(rate_updates)

  rates_eval <- matrix(0, n_rates, tau_length)

  for(i in 1:n_rates){
    partial_i <- rate_updates[i]

    rates_eval[i,] <- max(0,theta[partial_i]*grad_cv[partial_i]) +
      C_i_vals[partial_i]*(norm(x - control_variate, "2") + tau_grid*norm(theta,"2"))

    rates_eval[i,] <-  rates_eval[i,] + (x[partial_i] + tau_grid*theta[partial_i])/sigma2_prior

  }
  return(rates_eval)
}

n_ev <- 1e3
set.seed(1);z_ss <- zigzag(max_events = n_ev,  return_rates = return_rates_subsample, dnlogpi = get_grad,
                          x0 = beta, tau_max = 1, poly_order = 1, adapt_tau_max = T)

plot_pdmp(z_ss, pch = '.',coords = c(1,2,p), nsamples = 1e3, inds = 1:1e2)
