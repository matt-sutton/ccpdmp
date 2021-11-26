gen_samples <- function(positions, times, theta=NULL,
                        nsample = 10^3, burn = 1,
                        dynamics = 'lin', xstar = NULL){

  if(is.null(dim(positions))) positions <- matrix(positions, nrow = 1)

  positions <- positions[,burn:length(times), drop = F]
  if(dynamics != 'lin') {theta <- theta[,burn:length(times), drop = F]}
  times <- times[burn:length(times)] - times[burn]
  nsteps <- length(times)

  Tmax <- times[nsteps]
  dt <- Tmax/(nsample+2)
  t = dt
  t0 = times[1]
  x0 = positions[,1]
  samples <- matrix(0, nrow = length(x0), ncol = nsample)
  n <- 0

  for(i in 2:nsteps){
    x1 = positions[,i]
    t1 = times[i]
    theta0 = if(dynamics == 'lin') (x1 - x0)/(t1 - t0) else theta[,i-1]#(x1 - x0*cos(t1 - t0))/sin(t1 - t0)
    while(t + dt < t1 && n < nsample){
      n <- n+1
      t <- t + dt
      if(dynamics == 'lin'){
        x_s <- x0 + (t-t0)*theta0
      } else{
        x_s <- xstar + (x0-xstar)*cos(t-t0)+theta0*sin(t-t0)
      }
      samples[,n] <- x_s
    }
    x0 = x1; t0 = t1
  }

  return(list(xx = samples))
}
