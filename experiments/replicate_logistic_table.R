## Replicate Table 1 Increasing poly order for logistic regression
## See logistic.R for a simple example
##
generate.logistic.data <- function(beta, n.obs, siginv) {
  p <- length(beta)
  dataX <- mvtnorm::rmvnorm(n=n.obs*p,sigma = solve(siginv))
  vals <- dataX%*%beta
  generateY <- function(p) { rbinom(1, 1, p)}
  dataY <- sapply(1/(1 + exp(-vals)), generateY)
  return(list(dataX = dataX, dataY = dataY))
}

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
return_rates_poly <- function(x, theta, tau_grid, dnlogpi, rate_updates, poly){

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

  for(i in 1:n_rates){
    ## Calculate f(t) and derivatives of f(t) for Taylor expansion
    partial_i <- rate_updates[i]
    f <- theta[partial_i]*sum(phi_1*X[,partial_i]) + theta[partial_i]*x[partial_i]/1

    if(poly == 1){
      f_1 <- sum(1/4 * abs(X[,partial_i]*da_dt))+ theta[partial_i]^2/1
      poly_coef <- c(f_1, f)
    }
    if(poly == 2){
      f_1 <- theta[partial_i]*sum(phi_2*X[,partial_i]*da_dt) + theta[partial_i]^2/1
      f_2 <- sum(1/(6*sqrt(3))*abs(X[,partial_i]*da_dt^2))
      poly_coef <- c(f_2/factorial(2), f_1, f)
    }
    if(poly == 3){
      f_1 <- theta[partial_i]*sum(phi_2*X[,partial_i]*da_dt) + theta[partial_i]^2/1
      f_2 <- theta[partial_i]*sum(phi_3*X[,partial_i]*da_dt^2)
      f_3 <- sum(1/8*abs(X[,partial_i]*da_dt^3))

      poly_coef <- c(f_3/factorial(3), f_2/factorial(2), f_1, f)
    }

    rates_eval[i,] <- pracma::polyval(poly_coef, tau_grid)
  }
  return(rates_eval)
}

return_rates <- function(x, theta, tau_grid, dnlogpi, rate_updates){
  return(return_rates_poly(x, theta, tau_grid, dnlogpi, rate_updates, 1))
}

corrs <- round(log(seq(exp(0), exp(0.95), length.out = 7)),2)#seq(0, 0.95, length.out = 7)
corrs <- c(0.00, 0.25, 0.5, 0.65, 0.75, 0.85, 0.95)

leff <- list()
lefft <- list()
for( reps in 1:20){
  efft <- eff <- matrix(0, nrow = 3, ncol = length(corrs))
  for( cor_ind in 1:length(corrs)){

    set.seed(reps)
    n_ev <- 5*1e3
    p <- 5
    n <- 200
    sigma2_prior <- 100
    beta <- c(c(-1.25, 0.5), rep(-0.4,p-2))
    siginv <- diag(1, p,p)
    siginv[1,2] <- siginv[2,1] <- corrs[cor_ind]
    data <- generate.logistic.data(beta, n, siginv)
    X <- data$dataX
    y <- data$dataY

    return_rates <- function(x, theta, tau_grid, dnlogpi, rate_updates){
      return(return_rates_poly(x, theta, tau_grid, dnlogpi, rate_updates, 1))
    }
    set.seed(1);z_1_order <- zigzag(max_events = n_ev,  return_rates = return_rates, dnlogpi = get_grad,
                                    x0 = beta, tau_max = 1, poly_order = 1)
    eff[1, cor_ind] <- length(z_1_order$times)/max(z_1_order$n_prop)
    efft[1, cor_ind] <- length(z_1_order$times)/max(z_1_order$n_iterations)

    return_rates <- function(x, theta, tau_grid, dnlogpi, rate_updates){
      return(return_rates_poly(x, theta, tau_grid, dnlogpi, rate_updates, 2))
    }
    set.seed(1);z_2_order <- zigzag(max_events = n_ev,  return_rates = return_rates, dnlogpi = get_grad,
                                    x0 = beta, tau_max = 1, poly_order = 2)
    eff[2, cor_ind] <- length(z_2_order$times)/max(z_2_order$n_prop)
    efft[2, cor_ind] <- length(z_2_order$times)/max(z_2_order$n_iterations)

    return_rates <- function(x, theta, tau_grid, dnlogpi, rate_updates){
      return(return_rates_poly(x, theta, tau_grid, dnlogpi, rate_updates, 3))
    }
    set.seed(1);z_3_order <- zigzag(max_events = n_ev,  return_rates = return_rates, dnlogpi = get_grad,
                                    x0 = beta, tau_max = 1, poly_order = 3)
    eff[3, cor_ind] <- length(z_3_order$times)/max(z_3_order$n_prop)
    efft[3, cor_ind] <- length(z_3_order$times)/max(z_3_order$n_iterations)

    print("\n")
    print(efft)
    print("\n")
    plot_pdmp_multiple(list(o1 = z_1_order,o2 = z_2_order,o3 = z_3_order),pch = '.',
                       coords = c(1,2,p), nsamples = 1e3, inds = 1:1e3)

  }
  leff[[reps]] <- eff
  lefft[[reps]] <- efft
}

lmean_eff <- eff*0
for( reps in 1:length(leff)){
  lmean_eff <- lmean_eff + leff[[reps]]
}
lmean_eff <- lmean_eff/length(leff)
colnames(lmean_eff) <- corrs

lmean_efft <- efft*0
for( reps in 1:length(lefft)){
  lmean_efft <- lmean_efft + lefft[[reps]]
}
lmean_efft <- lmean_efft/length(lefft)
colnames(lmean_efft) <- corrs

# save(lmean_eff, lmean_efft, file = "Logit_poly.RData")
