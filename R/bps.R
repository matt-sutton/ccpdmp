## Base bps Sampler
bps <- function(max_events, dnlogpi, x0, tau_max = 0.5, poly_order = 2, adapt_tau_max = T,
                   return_rates = return_rates_bps, ref_rate = .5){

  t = 0; n_prop <- n_iterations <- 0; p <- length(x0)
  x = x0; theta = rnorm(p)

  n_rates <- 2
  tau_s <- rep(Inf, n_rates)

  # Simulate times
  tau_grid <- seq(from = 0, to = tau_max,  length.out = max(2, poly_order + 1))
  rates <- return_rates(x, theta, tau_grid, dnlogpi)

  type_sim = if(nrow(rates) == 5) "CC" else "poly"

  # Simulate event time
  if(type_sim == "poly"){
    event_sim = sim_rate_poly(eval_times = tau_grid, eval_rates = rates, poly_order = poly_order)
  } else{
    event_sim = sim_rates(eval_times = tau_grid, eval_rates = rates, poly_order = poly_order)
  }

  tau_s[1] = event_sim$t
  u_s = event_sim$u
  f_s = event_sim$f_evall

  # Simulate refreshment event
  tau_s[2] <- rexp(1, rate = ref_rate)

  velocities <- positions <- matrix(0, nrow = p, ncol = max_events);
  times = rep(0,max_events);

  velocities[,1] <- theta; positions[,1] <- x;

  num_evts = 1; n_ref = 0;
  event = FALSE

  while(num_evts < max_events){
    n_iterations = n_iterations +1

    mini_x <- which.min(tau_s)
    tau <- tau_s[mini_x]

    x = x + tau*theta
    t = t + tau

    if(mini_x == 2){
      # Refreshment
      theta <- rnorm(p)
      n_ref = n_ref + 1
      event = TRUE
    } else {
      # Bounce
      if(u_s < 1e-10){
        n_prop <- n_prop + 1
        grad <- dnlogpi(x)
        rate <- sum(grad*theta)
        acc_prb <- rate/f_s

        if(acc_prb > 1 +1e-10){
          print(paste("Invalid thinning probability:",acc_prb))
        }
        if(runif(1) <= acc_prb){
          theta = theta - 2*sum(theta*grad)/sum(grad*grad)*grad
          event = TRUE
        }
      }
    }

    if(event){
      # If there was an event store this information
      # Re-simulate rates for all dimensions (non-local implementation)

      num_evts = num_evts + 1
      velocities[,num_evts] <- theta; positions[,num_evts] <- x
      times[num_evts] = t

      rates <- return_rates(x, theta, tau_grid, dnlogpi)

      if(type_sim == "poly"){
        event_sim = sim_rate_poly(eval_times = tau_grid, eval_rates = rates, poly_order = poly_order)
      } else{
        event_sim = sim_rates(eval_times = tau_grid, eval_rates = rates, poly_order = poly_order)
      }
      tau_s[1] = event_sim$t
      u_s = event_sim$u
      f_s = event_sim$f_evall

      # Simulate refreshment event
      tau_s[2] <- rexp(1, rate = ref_rate)

      event = FALSE

      ## Adapt tau_max
      if((num_evts %% 1e2 == 0) & adapt_tau_max) {
        tau_max <- quantile(diff(times[1:num_evts]), 0.8)[[1]]
        tau_grid <- seq(from = 0, to = tau_max,  length.out = max(2, poly_order + 1))
      }


    } else {
      # If there was no event
      # Adjust simulated times
      tau_s <- tau_s - tau

      # Re-simulate times for all tau_s less than zero
      if(tau_s[1] <= 0){
        rates <- return_rates(x, theta, tau_grid, dnlogpi)
        if(type_sim == "poly"){
          event_sim = sim_rate_poly(eval_times = tau_grid, eval_rates = rates, poly_order = poly_order)
        } else{
          event_sim = sim_rates(eval_times = tau_grid, eval_rates = rates, poly_order = poly_order)
        }
        tau_s[1] = event_sim$t
        u_s = event_sim$u
        f_s = event_sim$f_evall
      }
      if(tau_s[2] <= 0){
        tau_s[2] <- rexp(1, rate = ref_rate)
      }
    }
  }
  return (list(positions=positions,theta=velocities, times=times,
               n_iterations = n_iterations, n_prop = n_prop, n_ref = n_ref))
}

