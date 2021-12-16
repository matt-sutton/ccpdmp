library(ccpdmp)
library(latex2exp)

## function to thin
Const <- 0.5
fun <- function(t){
  return( -t^3 + 3*t^2 -3*t + Const )
}
x_seq <- seq(0, 1, length.out = 100)
ply <- pracma::polyfit(x_seq, fun(x_seq), n = 3)

## Decomposition 1 (default)
f_u <- function(t) 3*t^2 + Const
f_n <- function(t) -t^3 - 3*t
f_n_d <- function(t) -3*t^2 - 3

x_eval <- c(0,1) ## Abscessisa for creating the upper bound

rates_eval_f1 <- rbind(f_u(x_eval),f_n(x_eval),f_n_d(x_eval), c(0,0), f_u(x_eval) + f_n(x_eval))

set.seed(1);sim_decomp_1 = sim_rates(eval_times = x_eval, eval_rates = rates_eval_f1,
                                     poly_order = 0,n_points = 500, u = .1)
print(sim_decomp_1$t)
t_range <- sim_decomp_1$f_range

decomp_1_col <- 2
par(mfrow=c(2,3), mar = c(2,2,2,1))

# plot(t_range, f_u(t_range), type = 'l', xlab = 't',
#      main = 'f_u',ylab='');grid()


# lines(t_range, f_u(0) + (f_u(1) - f_u(0))*t_range, col = decomp_1_col)
# plot(t_range, f_n(t_range), type = 'l',  xlab = 't',
#      main = 'f_n',ylab='');grid()
#
# intersect_fn <- (f_n(1)-f_n_d(1) - f_n(0))/(f_n_d(0) - f_n_d(1))
# lines(t_range[t_range<=intersect_fn], f_n(0) + f_n_d(0)*t_range[t_range<intersect_fn],
#       col = decomp_1_col)
# lines(t_range[t_range>intersect_fn], (f_n(1) - f_n_d(1)) + f_n_d(1)*t_range[t_range>intersect_fn],
#       col = decomp_1_col)
#
# plot(t_range, fun(t_range), type = 'l', xlab = 't',
#      main = 'f = f_u + f_n',ylab='');grid()
# lines(t_range, sim_decomp_1$f_upper_range, col = decomp_1_col)
#
# abline(v = sim_decomp_1$t)


############ Actually simulate
gen_a_rate <- function(tau_max, int_length = 1, uval = 0.5){

  tau <- 0
  tau_m <- tau_max

  while(uval > 1e-14 & tau < max(int_length)){
    x_eval <- c(tau, tau_m)
    rates_eval_f1 <- rbind(f_u(x_eval),f_n(x_eval),f_n_d(x_eval), c(0,0), f_u(x_eval) + f_n(x_eval))

    tus = sim_rates(eval_times = x_eval-tau, eval_rates = rates_eval_f1,
                               poly_order = 0, u = uval, n_points = 100)
    print(x_eval)
    uval <- tus$u
    # tau_m <- tau + tau_max
    plot(tus$f_range+tau, fun(tus$f_range+tau), type = 'l', col = 1,ylim = range(c(fun(tus$f_range+tau), tus$f_upper_range)))
    lines(tus$f_range+tau, tus$f_upper_range, col = 2)

    tau <- tau + tus$t
    abline(v= tau)
    ## Additional acc/rej
    if(uval < 1e-14){
      prb <- max(0,fun(tau)/tus$f_evall)
      print(prb)
      if( runif(1) < prb){
        break
      } else{
        uval <- rexp(1)
      }
    }
  }

  # if( uval > 0 ){
  #   tau <- int_length
  # }
  return(list(tau = tau, u_val = uval))
}

Const <- 1
dev.off();
par(mfrow=c(2,3), mar = c(2,2,2,1));rt <- gen_a_rate(2,int_length = 2, uval = 0.5)
rt$tau


## Decomposition 2 (optimal)
f_u <- function(t) -t^3 + 3*t^2 -3*t + 1
f_n <- function(t) rep(-0.05, length(t))
f_n_d <- function(t) rep(0, length(t))
rates_eval_f2 <- rbind(f_u(x_eval),f_n(x_eval),f_n_d(x_eval),c(0,0),f_u(x_eval) + f_n(x_eval))

sim_decomp_2 = sim_rates(eval_times = x_eval, eval_rates = rates_eval_f2,
                         poly_order = 3, n_points = 500)
decomp_2_col <- 3
plot(t_range, f_u(t_range), type = 'l', xlab = 't',
     main = 'f_u',ylab='');grid()
lines(t_range, f_u(0) + (f_u(1) - f_u(0))*t_range, col = decomp_2_col)
plot(t_range, f_n(t_range), type = 'l',  xlab = 't',
     main = 'f_n',ylab='');grid()

intersect_fn <- 1
lines(t_range[t_range<=intersect_fn], f_n(0) + f_n_d(0)*t_range[t_range<intersect_fn],
      col = decomp_2_col)

plot(t_range,fun(t_range),  type = 'l', ylim = c(-0.1,1), xlab = "t", main = "f = f_n + f_u")
lines(t_range, sim_decomp_2$f_upper_range, col = decomp_2_col);grid()
legend('topright', legend = c("Upper-bound 1", "Upper-bound 2", "Exact"), col = c("red", "green", "black"), lwd = 1)

