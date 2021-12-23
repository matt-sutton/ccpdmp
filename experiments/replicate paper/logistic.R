## Example Code for logistic regression
## More exotic GLMs can be implemented with different phi and f
##

generate.logistic.data <- function(beta, n.obs, siginv) {
  p <- length(beta)
  dataX <- mvtnorm::rmvnorm(n=n.obs*p,sigma = solve(siginv))
  vals <- dataX%*%beta
  generateY <- function(p) { rbinom(1, 1, p)}
  dataY <- sapply(1/(1 + exp(-vals)), generateY)
  return(list(dataX = dataX, dataY = dataY))
}

set.seed(1)
n_ev <- 1e3
p <- 5
n <- 500
sigma2_prior <- 1000
beta <- c(c(-1.25, 0.5), rep(-0.4,p-2))
siginv <- diag(1, p,p)
siginv[1,2] <- siginv[2,1] <- 0.9
data <- generate.logistic.data(beta, n, siginv)
X <- data$dataX
y <- data$dataY

get_grad <- function(x, index){
  expa <- exp(X%*%x)
  phi_1 <- expa/(1+expa) - y

  grad <- rep(0, length(index))
  for(i in 1:length(index)){
    grad[i] <- sum(phi_1*X[,index[i]]) + x[index[i]]
  }
  return(grad)
}

## Return the rates for an upper-bounding 3rd order polynomial:
return_rates_3 <- function(x, theta, tau_grid, dnlogpi, rate_updates){

  tau_length <- length(tau_grid);
  n_rates <- length(rate_updates)

  rates_eval <- matrix(0, n_rates, tau_length)

  a <- X%*%x
  expa <- exp(a)
  da_dt <- X%*%theta

  # First 3 derivatives of phi(a) = log(1+exp(a)) - ya
  phi_1 <- expa/(1+expa) - y
  phi_2 <- expa/(1+expa)^2
  phi_3 <- -expa*(expa -1)/(expa+1)^3
  phi_4_bound <- 1/8

  for(i in 1:n_rates){
    ## Calculate f(t) and derivatives of f(t) for Taylor expansion
    partial_i <- rate_updates[i]
    f <- theta[partial_i]*sum(phi_1*X[,partial_i]) + theta[partial_i]*x[partial_i]/1
    f_1 <- theta[partial_i]*sum(phi_2*X[,partial_i]*da_dt) + theta[partial_i]^2/1
    f_2 <- theta[partial_i]*sum(phi_3*X[,partial_i]*da_dt^2)
    f_3 <- sum(phi_4_bound*abs(X[,partial_i]*da_dt^3))

    poly <- c(f_3/factorial(3), f_2/factorial(2), f_1, f)

    rates_eval[i,] <- pracma::polyval(poly, tau_grid)
  }
  return(rates_eval)
}

n_ev <- 1e4
set.seed(1);z_3 <- zigzag(max_events = n_ev,  return_rates = return_rates_3, dnlogpi = get_grad,
                          x0 = beta, tau_max = 1, poly_order = 3, adapt_tau_max = T)
plot_pdmp(o1 = z_3,pch = '.',coords = c(1,2,p), nsamples = 1e3, inds = 1:1e3)
