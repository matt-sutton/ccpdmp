## Effect of Tau_max

kappa <- 1
get_grad <- function(x, index){
  x1 <- x[1]; x2 <- x[2]
  grad <- c(2*(x1-1) + 4*kappa*(x1^2-x2)*x1,
            2*kappa*(x2-x1^2))
  return(grad[index])
}

tau_max_vals <- seq(0.15, to = 4, length.out = 30)
set.seed(1)
n_ev <- 1e4
eff_p <-eff <- rep(0,length(tau_max_vals))
i<-1
for( i in 1:length(tau_max_vals)){
  set.seed(1);zz <- zigzag(max_events = n_ev, get_grad, x0 = c(1,1), tau_max = tau_max_vals[i], poly_order = 3, adapt_tau_max = F)
  eff[i] <- n_ev/zz$n_iterations
  eff_p[i] <- n_ev/zz$n_prop
  plot(t(zz$positions), type = 'l')
}

## The adaptive approach is initalised at a poor starting tau_max value
set.seed(1);zz <- zigzag(max_events = n_ev, get_grad, x0 = c(1,1), poly_order = 3,tau_max = 4, adapt_tau_max = T)
n_ev/zz$n_iterations
n_ev/max(zz$n_prop)
eff
text_width <- 1
par(cex.main=text_width, cex.lab=text_width, cex.axis=text_width, lwd = text_width, mar = c(4,4,1.2,1.2))
par(mfrow = c(1,1))
library(latex2exp)
text_width <- 1
par(cex.main=text_width, cex.lab=1.3, cex.axis=text_width, lwd = text_width, mar = c(4,5,1.2,1.2))
plot(tau_max_vals, eff, type = 'l', xlab = TeX('$\\tau_\\max$'), ylab = 'Efficiency', ylim = c(0,1))
abline(h = n_ev/zz$n_iterations, col = 'red')
grid()

