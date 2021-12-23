## Return evaluations of the rates
## Will evaluate the rate on the tau_grid points
## dnlogpi returns the partials of the negative log target

return_rates_zigzag <- function(x, theta, tau_grid, dnlogpi, rate_updates){

  tau_length <- length(tau_grid);
  n_rates <- length(rate_updates)

  rates_eval <- matrix(0, n_rates, tau_length)

  for( i in 1:tau_length ){
    x_g <- x + tau_grid[i]*theta
    rates_eval[,i] <- dnlogpi(x_g, rate_updates)*theta[rate_updates]
  }
  return(rates_eval)
}

## Return evaluations of the global BPS rate
## Will evaluate the rate on the tau_grid points
## dnlogpi returns all partials of the negative log target
return_rates_bps <- function(x, theta, tau_grid, dnlogpi){

  tau_length <- length(tau_grid);

  rates_eval <- matrix(0, 1, tau_length)

  for( i in 1:tau_length ){
    x_g <- x + tau_grid[i]*theta
    rates_eval[i] <- sum(dnlogpi(x_g)*theta)
  }
  return(rates_eval)
}

## Vectorised code for simulating from lambda = exp(x)
vexp<- Vectorize(exp_inv_t,SIMPLIFY = T)
