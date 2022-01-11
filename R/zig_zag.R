## Base Zig-Zag Sampler
zigzag <- function(max_events, dnlogpi, x0, tau_max = 0.5, poly_order = 2, adapt_tau_max = T,
                   return_rates = return_rates_zigzag, add_interp = 0){

  t = 0; n_prop <- n_iterations <- 0; n_rates <- length(x0)
  x = x0; theta = rep(1, n_rates)

  tau_s = rep(Inf, n_rates)
  u_s = rexp(n_rates)
  f_s = rep(Inf, n_rates)

  # Simulate times
  tau_grid <- seq(from = 0, to = tau_max,  length.out = poly_order + 1)
  rates <- return_rates(x, theta, tau_grid, dnlogpi, 1:n_rates) + add_interp

  for( i in 1:n_rates){
    event_sim = sim_rate_poly(eval_times = tau_grid, eval_rates = rates[i,], poly_order)

    tau_s[i] = event_sim$t
    u_s[i] = event_sim$u
    f_s[i] = event_sim$f_evall
  }

  velocities <- positions <- matrix(0, nrow = n_rates, ncol = max_events);
  times = rep(0,max_events);

  velocities[,1] <- theta; positions[,1] <- x;

  num_evts = 1
  event = FALSE

  while(num_evts < max_events){
    n_iterations = n_iterations +1

    mini_x <- which.min(tau_s)
    tau <- tau_s[mini_x]

    x = x + tau*theta
    t = t + tau
    if(u_s[mini_x] < 1e-10){
      n_prop <- n_prop + 1
      grad <- dnlogpi(x, mini_x)
      rate <- grad*theta[mini_x]
      acc_prb <- rate/f_s[mini_x]

      if(acc_prb > 1 +1e-10){
        print(paste("Invalid thinning probability:",acc_prb))
      }
      if(runif(1) <= acc_prb){
        theta[mini_x] = - theta[mini_x]
        event = TRUE
      }
    }

    if(event){
      # If there was an event store this information
      # Re-simulate rates for all dimensions (non-local implementation)

      num_evts = num_evts + 1
      velocities[,num_evts] <- theta; positions[,num_evts] <- x
      times[num_evts] = t

      rates <- return_rates(x, theta, tau_grid, dnlogpi, 1:n_rates) + add_interp

      for( i in 1:n_rates){
        event_sim = sim_rate_poly(eval_times = tau_grid, eval_rates = rates[i,], poly_order)

        tau_s[i] = event_sim$t
        u_s[i] = event_sim$u
        f_s[i] = event_sim$f_evall
      }
      event = FALSE

      ## Adapt tau_max
      if((num_evts %% 1e2 == 0) & adapt_tau_max) {
        tau_max <- quantile(diff(times[1:num_evts]), 0.8)[[1]]
        tau_grid <- seq(from = 0, to = tau_max,  length.out = poly_order + 1)
      }


    } else {
      # If there was no event
      # Adjust simulated times
      tau_s <- tau_s - tau

      # Re-simulate times for all tau_s less than zero:
      update_rates <- which(tau_s <= 0)
      rates <- return_rates(x, theta, tau_grid, dnlogpi, update_rates) + add_interp
      for( i in 1:length(update_rates)){
        event_sim = sim_rate_poly(eval_times = tau_grid, eval_rates = rates[i,], poly_order)
        j <- update_rates[i]
        tau_s[j] = event_sim$t
        u_s[j] = event_sim$u
        f_s[j] = event_sim$f_evall
      }
    }
  }
  return (list(positions=positions,theta=velocities, times=times, n_iterations = n_iterations, n_prop = n_prop))
}

