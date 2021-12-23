library(ccpdmp)

## Superposition method to simulate from lambda = max(0, x - y + exp(x))
superpos_rates <- function(x, theta, y){

  # Simulate time for rate lambda1 = max(0, x - y )
  tau_linear <- ccpdmp::linear_inv_t(a = sum(theta*(x-y)), b = sum(theta^2), u = rexp(1), tmax = Inf)[1]

  # Simulate times for rates lambda2 = max(0, exp(x))
  tau_exp <- vexp(x, theta, runif(length(x)))

  # take the min to simulate from lambda1 + lambda2
  tau <- min(c(tau_linear, tau_exp))

  # Evaluate the upper bound at the simulated time
  x_t <- x+theta*tau
  ub_eval <- max(0, sum(theta*(x_t-y))) + sum(pmax(0, theta*exp(x_t)))

  res <- list(t = tau, f_evall = ub_eval)
  return(res)
}

bps_superposition <- function(max_events, y, x0, theta0, verb = FALSE, ref_rate = .5){

  t = 0; eps = 1e-10; x = x0;
  theta = theta0;

  tus = superpos_rates(x, theta, y)
  taus = c(tus$t, rexp(1, ref_rate));
  f_s = tus$f_evall

  positions = matrix(0, nrow = d, ncol = max_events);   times = rep(0,max_events);
  velocities = matrix(0, nrow = d, ncol = max_events)

  positions[,1] <- x; velocities[,1] <- theta

  num_evts = 0; nE = 1
  n_ref <- n_iterations <- 0

  while(nE < max_events){
    n_iterations <- n_iterations + 1
    mini = which.min(taus)[1]
    x = x + taus[mini]*theta
    t = t + taus[mini]

    if(mini == 1){
      ## Check event
      grad <- (x + exp(x) - y)
      acc_prb <-  sum(theta*grad)/f_s[mini]
      if(acc_prb > 1 + 1e-12){
        print(paste("Invalid thinning probability:",acc_prb))
      }
      if(runif(1) <= acc_prb){
        nE = nE + 1
        theta = theta - 2*sum(theta*grad)/sum(grad*grad)*grad
        positions[,nE] <- x;
        velocities[,nE] <- theta;
        times[nE] <- t

        if(nE %% 1e2 ==0 & verb){
          print(nE)
        }
      }
    } else {
      ## Refreshment
      n_ref = n_ref + 1
      theta <- rnorm(d)
    }

    # # Simulate new time
    tus = superpos_rates(x, theta, y)
    taus = c(tus$t, rexp(1, ref_rate));  f_s = tus$f_evall
  }
  return (list(positions=positions,theta=theta,times=times, n_iterations = n_iterations, n_ref = n_ref))
}


return_rates <- function(x, theta, tau_grid, dnlogpi){

  rates_eval <- matrix(0, nrow = 5, ncol = length(tau_grid))
  for( i in 1:length(tau_grid) ){
    f_u <- f_n <- f_n_d <- 0
    x_g <- x+tau_grid[i]*theta

    for(d in 1:length(x)){
      f_u <- f_u + theta[d]*(x_g[d] - y[d])
      if(theta[d] < 0){
        f_n <- f_n + theta[d]*exp(x_g[d])
        f_n_d <- f_n_d + theta[d]^2*exp(x_g[d])
      } else {
        f_u <- f_u + theta[d]*exp(x_g[d])
      }
    }
    rates_eval[,i] <- c(f_u,f_n,f_n_d,0, f_u + f_n)
  }
  return(rates_eval)
}
get_grad <- function(x, index){
  grad <- exp(x) + x - y
  return(grad)
}

(d_length <- 2^c(2:6))
nrep <- 20
n_iter <- array(0, dim = c(2,nrep,length(d_length)))
secs <- array(0, dim = c(2,nrep,length(d_length)))

d <- max(d_length)
set.seed(1)
xs <- rnorm(d)
y_full <- NULL
for(r in 1:nrep){
  y_full <- rbind(y_full,rpois(d, lambda = as.numeric(exp(xs))))
}

n_ev <- 1e3
for(i in 1:length(d_length)){
  d <- d_length[i]
  print(paste("###### Now at d = ",d))

  for( r in 1:nrep){

    ## Update y for each iteration
    y <- y_full[r,1:d]
    t1 <- system.time({set.seed(r);b <- bps(max_events = n_ev,  return_rates = return_rates, dnlogpi = get_grad,
                         x0 = xs[1:d], tau_max = 1, poly_order = 1, adapt_tau_max = T, ref_rate = 1e-10)})

    t2 <- system.time({set.seed(r);b_s <- bps_superposition(max_events = n_ev, y = y, x0 = xs[1:d],
                                                            ref_rate = 1e-10, theta0 = b$theta[,1])})

    # Store compute time and numer of algorithm iterations
    secs[1,r,i] <- t1[3]; n_iter[1,r,i] <- b$n_iterations
    secs[2,r,i] <- t2[3]; n_iter[2,r,i] <- b_s$n_iterations
  }
  save(secs, n_iter, file = "compare_super.Rdata")
}
n_iter_a <- apply(n_iter, c(1,3), mean) ## row 1 = CC, row 2 = super position
secs_a <- apply(secs, c(1,3), mean)

text_width <- 1.2
par(cex.main=text_width, cex.lab=text_width, cex.axis=text_width, lwd = text_width, mar = c(4,4,1.2,1.2))
par(mfrow = c(1,2))
plot(x = 1:length(d_length), y = 1e3/n_iter_a[2,], ylim = c(0,max(1e3/n_iter_a[1,])),
     xlim = c(1,5), xaxt = "n", type = 'l', lwd = 1.7,
     ylab = "Events / Iterations", xlab = "Dimension (n)")
axis(1, at = 1:length(d_length),
     labels = c(expression(2^2),expression(2^3),expression(2^4),expression(2^5), expression(2^6)))
lines(1e3/n_iter_a[1,],lwd = 1.7,col = 2)
grid()

plot(x = 1:5, y = secs_a[2,], ylim = c(min(secs_a[1,]),max(secs_a[2,])), xlim = c(1,5), xaxt = "n", type = 'l',
     ylab = "Seconds", xlab = "Dimension (n)", lwd = 1.7)
axis(1, at = 1:length(d_length),
     labels = c(expression(2^2),expression(2^3),expression(2^4),expression(2^5), expression(2^6)))
lines(secs_a[1,],col = 2, lwd = 1.7)
grid()
legend("topleft", legend = c("Superposition  ","CCPDMP  "), col = c(1,2), lwd = 1.7, cex = 1.05)


