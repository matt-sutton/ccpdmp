kappa <- 1
get_grad <- function(x, index){
  x1 <- x[1]; x2 <- x[2]
  grad <- c(2*(x1-1) + 4*kappa*(x1^2-x2)*x1,
            2*kappa*(x2-x1^2))
  return(grad[index])
}

tau_max_vals <- seq(0.15, to = 4, length.out = 30)
set.seed(1)
n_ev <- 5*1e3
eff_p <-eff <- rep(0,length(tau_max_vals))
i<-1
for( i in 1:length(tau_max_vals)){
  set.seed(1);z_j <- zigzag(max_events = n_ev, get_grad, x0 = c(1,1), tau_max = tau_max_vals[i], poly_order = 3, adapt_tau_max = F)
  eff[i] <- n_ev/z_j$n_iterations
  eff_p[i] <- n_ev/z_j$n_prop
  plot(t(z_j$positions), type = 'l')
}
set.seed(1);z_j <- zigzag(max_events = n_ev, get_grad, x0 = c(1,1), tau_max = 4, poly_order = 3, adapt_tau_max = T)
n_ev/z_j$n_iterations
n_ev/max(z_j$n_prop)
eff
eff_p
text_width <- 1
par(cex.main=text_width, cex.lab=text_width, cex.axis=text_width, lwd = text_width, mar = c(4,4,1.2,1.2))
library(latex2exp)
plot(tau_max_vals, eff, type = 'l', xlab = TeX('$\\tau_\\max$'), ylab = 'Efficiency', ylim = c(0,1))
lines(tau_max_vals, eff_p, col = 'green')
abline(h = n_ev/z_j$n_total, col = 'red')

text_width <- 1
par(cex.main=text_width, cex.lab=1.3, cex.axis=text_width, lwd = text_width, mar = c(4,5,1.2,1.2))
plot(tau_max_vals, eff, type = 'l', xlab = TeX('$\\tau_\\max$'), ylab = 'Efficiency', ylim = c(0,1))
abline(h = n_ev/z_j$n_iterations, col = 'red')
grid()

