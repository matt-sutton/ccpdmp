# Concave-convex PDMP

## Description

This is a trimmed version of the code I have which focuses on simulating the rate using the CC-PDMP approach and the algorithm for Zig-Zag in R. The c++ code is in the src directory and consists of classes to simulate a rate based on a concave-convex decomposition, evaluate and construct a polynomial and simulate a rate based on a polynomial bound. 

The code in R contains the functions for implementing the Zig-Zag. 

## Install

First clone the repo locally. Install R and Rstudio. If you click the ccpdmp.Rproj it should open up Rstudio in a project. On the right you'll see all the files and information of the package. Before installing this package you may need to run the following to get some required packages:

```{r}
install.packages("Rcpp")
install.packages("RcppArmadillo")
```
To compile and run this package press control shift B. You should be able to run the code below.

## A Quick example of the method

We are interested in sampling a density:
$$
\pi(x) \propto \exp(-U(x))
$$
A simple example is the normal distribution with mean zero and variance 1: $U(x) = \frac{1}{2}x^2$. We need a function for the derivative of this:
$$
\frac{dU(x)}{dx_i} = x_i
$$
The rate for the i-th component of the Zig-Zag is $\lambda_i(t) = \max(0, \frac{dU(x+vt)}{dx_i}) = \max(0, x_i + tv_i)$. This event rate (ignoring the max part) is linear. The code below shows how this can be simulated using the Zig-Zag process in this package. 

```{r}
example_dnlogpi # Example potential calculation

z <- zigzag(1e3, example_dnlogpi, x0 = c(0,0), tau_max = 1, poly_order = 1)
plot_pdmp(z, mcmc_samples = matrix(rnorm(2*1e3), ncol = 2))
```

The example_dnlogpi is a function taking x, and index as arguments and returning a vector with elements $\frac{dU(x+vt)}{dx_i} = x_i$ for all $i$ in the index argument. The polynomial order is 1 indicating a linear function is used to simulate the rate (which is exact here). The tau_max parameter indicates the maximum length of time considered ahead of the process. 


