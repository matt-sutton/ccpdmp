## Base Zig-Zag Sampler
example_dnlogpi <- function(x, index){
  return(x[index])
}

## Return the rates for each dimension evaluated at times tau_grid ahead
return_rates <- function(x, theta, tau_grid, dnlogpi, index = NULL){
  l_tau <- length(tau_grid);
  nx <- length(x)
  if(is.null(index)){
    index <- 1:nx
  }
  rates_eval <- matrix(0, nx, l_tau)

  for( i in 1:l_tau ){
    x_g <- x + tau_grid[i]*theta
    grads <- dnlogpi(x_g, index)

    for(j in 1:length(index)){
      partial_index <- index[j]
      rates_eval[j,i] <- grads[j]*theta[partial_index]
    }
  }
  return(rates_eval)
}

zigzag <- function(max_events, dnlogpi, x0, tau_max = 0.5, poly_order = 2, echo = FALSE){
  t = 0; eps = 1e-10; nits <- 0
  nvel <- length(x0)
  x = x0; theta = rep(1, nvel)

  taus = rep(Inf, nvel)
  u_s = rexp(nvel)
  f_s = rep(Inf, nvel)

  # Simulate times
  tau_grid <- seq(0, to = tau_max,  length.out = poly_order + 1)
  rates <- return_rates(x, theta, tau_grid, dnlogpi)
  for( i in 1:nvel){
    tus = sim_rate_poly(eval_times = tau_grid, eval_rates = rates[i,], poly_order)
    taus[i] = tus$t
    u_s[i] = tus$u
    f_s[i] = tus$f_evall
  }

  v_samples <- positions <- matrix(0, nrow = nvel, ncol = max_events);
  ts = rep(0,max_events);
  v_samples[,1] <- theta; positions[,1] <- x;

  v = 0;  j = 0;  num_evts = 1
  event = FALSE

  while(num_evts < max_events){

    mini_x <- which.min(taus)
    tau <- taus[mini_x]

    x = x + tau*theta
    t = t + tau
    if(u_s[mini_x] < 1e-10){
      grad <- dnlogpi(x, mini_x)
      rate <- grad*theta[mini_x]
      acc_prb <- rate/f_s[mini_x]
      if(acc_prb > 1.001){
        print(paste("oops x:",acc_prb))
      }
      if(runif(1) <= acc_prb){
        theta[mini_x] = -theta[mini_x]
        event = TRUE
      }
    }

    if(event){
      # If there was an event store this information
      # Re-simulate rates for all dimensions

      num_evts = num_evts + 1
      v_samples[,num_evts] <- theta; positions[,num_evts] <- x
      ts[num_evts] = t
      event = FALSE
      if(echo){
        print(num_evts)
      }
      rates <- return_rates(x, theta, tau_grid, dnlogpi)
      for( i in 1:nvel){
        tus = sim_rate_poly(eval_times = tau_grid, eval_rates = rates[i,], poly_order)
        taus[i] = tus$t
        u_s[i] = tus$u
        f_s[i] = tus$f_evall
      }
    } else {
      # If there was no event
      # Adjust simulated times
      taus <- taus-tau

      # Re-simulate times for all taus less than zero:
      update_rates <- which(taus <= 0)
      rates <- return_rates(x, theta, tau_grid, dnlogpi, update_rates)
      for( j in update_rates){
        tus = sim_rate_poly(eval_times = tau_grid, eval_rates = rates[j,], poly_order)
        taus[j] = tus$t
        u_s[j] = tus$u
        f_s[j] = tus$f_evall
      }
    }
    nits = nits +1
  }
  return (list(positions=positions,theta=v_samples,times=ts, nits = nits))
}
