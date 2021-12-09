theta <- 0
v <- -1
b <- -1
library(ccpdmp)
library(latex2exp)
## function to thin
fun <- function(t){
  return( -t^3 + 3*t^2 -3*t + 0.95 )
}
x_seq <- seq(0, 1, length.out = 100)
ply <- pracma::polyfit(x_seq, fun(x_seq), n = 3)

## Decomposition 1 (default)
f_u <- function(t) 3*t^2 + 0.95
f_n <- function(t) -t^3 - 3*t
f_n_d <- function(t) -3*t^2 - 3

x_eval <- c(0,1) ## Abscessisa for creating the upper bound

rates_eval_f1 <- rbind(f_u(x_eval),f_n(x_eval),f_n_d(x_eval), c(0,0), f_u(x_eval) + f_n(x_eval))

sim_decomp_1 = sim_rates(eval_times = x_eval, eval_rates = rates_eval_f1,
                poly_order = 3, n_points = 100)

t_range <- sim_decomp_1$f_range
par(mfrow=c(2,2), mar = c(2,2,2,1))
decomp_1_col <- 2
plot(t_range, f_u(t_range), type = 'l', xlab = 't',
     main = 'f_u',ylab='');grid()
lines(t_range, f_u(0) + (f_u(1) - f_u(0))*t_range, col = decomp_1_col)

plot(t_range, f_n(t_range), type = 'l',  xlab = 't',
     main = 'f_n',ylab='');grid()

lines(t_range, f_n(0) + f_n_d(0)*t_range,
      col = decomp_1_col)
lines(t_range, f_n(0) + f_n_d(1)*t_range,
      col = decomp_1_col)

a + f_n_d(1)*1 =



## Decomposition 2 (optimal)
f_u <- function(t) -t^3 + 3*t^2 -3*t + 1
f_n <- function(t) rep(-0.05, length(t))
f_n_d <- function(t) rep(0, length(t))
rates_eval_f2 <- rbind(f_u(x_eval),f_n(x_eval),f_n_d(x_eval),c(0,0),f_u(x_eval) + f_n(x_eval))

sim_decomp_2 = sim_rates(eval_times = x_eval, eval_rates = rates_eval_f2,
                         poly_order = 3, n_points = 100)
decomp_2_col <- 3
plot(sim_decomp_1$f_range, f_u(sim_decomp_1$f_range), type = 'l', xlab = 't',
     main = 'f_u',ylab='', col = decomp_2_col);grid()
plot(sim_decomp_1$f_range, f_n(sim_decomp_1$f_range), type = 'l', xlab = 't',
     main = 'f_u',ylab='', col = decomp_2_col);grid()
plot(sim_decomp_1$f_range, sim_decomp_1$f_upper_range, type = 'l', col=3)

plot(tus$f_range, tus$f_upper_range, type = 'l', ylim = c(-0.1,1), col = 'red', xlab = "t", ylab = "f(t)")
lines(x_seq, fun(x_seq), col = 'black')


plot(tus$f_range, tus$f_upper_range, type = 'l', ylim = c(-0.1,1), col = 'red', xlab = "t", ylab = "f(t)")
lines(x_seq, fun(x_seq), col = 'black')
tus = sim_rates(eval_times = x_eval, eval_rates = rates_eval_f2,
                           poly_order = 3, n_points = 100)
lines(tus$f_range, tus$f_upper_range, col = 'green')
grid()
legend('bottomleft', legend = c("Decomposition 1", "Decomposition 2", "Exact"), col = c("red", "green", "black"), lwd = 1)
