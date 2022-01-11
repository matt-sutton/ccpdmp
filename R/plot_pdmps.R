plot_pdmp <- function(pdmp_res, coords = 1:2, inds = 1:10^3, nsamples = 10^3,
                      burn = 0.1, mcmc_samples=NULL, pch = 20){
  ndim <- length(coords)
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  par(mfrow = c(ndim,ndim),
      mar = rep(2,4))
  samples <- gen_samples(pdmp_res$positions, pdmp_res$times, pdmp_res$theta,
                         nsample = nsamples, burn = burn*length(pdmp_res$times))$samples
  for( i in 1:ndim ){
    for(j in 1:ndim ){
      if(i == j & nsamples > 0){
        ds <- density(samples[coords[i],])
        plot(ds, main='',xlab='', ylab='', col = 'red')
        if(!is.null(mcmc_samples)){
          lines(density(mcmc_samples[,coords[i]], bw = ds$bw), col = 'blue')
        }
      }
      if( i != j ){
        xrange <- range(c(samples[coords[i],], pdmp_res$positions[coords[i],inds]))
        yrange <- range(c(samples[coords[j],], pdmp_res$positions[coords[j],inds]))
        plot(pdmp_res$positions[coords[i],inds], pdmp_res$positions[coords[j],inds],
             xlim = xrange, ylim = yrange, type = 'l', main = '')
        points(samples[coords[i],], samples[coords[j],], col = 'red', pch = pch)
        if(!is.null(mcmc_samples)){
          points(mcmc_samples[,coords[i]], mcmc_samples[,coords[j]], col = 'blue', pch = pch)
        }
      }
    }
  }
  legend('topleft',legend = c("pdmp","mcmc"), col = c("red","blue"), lwd = 1 )
}
plot_pdmp_multiple <- function (list_pdmp, coords = 1:2, inds = 1:10^3, nsamples = 10^3,
                                burn = 0.1, mcmc_samples = NULL, pch = 19, cols = NULL) {
  ndim <- length(coords)
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  par(mfrow = c(ndim, ndim), mar = rep(2, 4))
  nres <- length(list_pdmp)
  if(is.null(list_pdmp)){
    names(list_pdmp) <- 1:nres
  }
  samples <- lapply(list_pdmp, function(res) gen_samples(res$positions,
                                                         res$times, res$theta, nsample = nsamples, burn = burn *
                                                           length(res$times))$samples)
  colres <- if (is.null(cols)) {
    1:(nres + !is.null(mcmc_samples))
  }
  else {
    cols
  }
  for (i in 1:ndim) {
    for (j in 1:ndim) {
      if (i == j & nsamples > 0) {
        max_l <- max(sapply(samples, function(vecs) {
          ds <- density(vecs[coords[i], ])
          return(max(ds$y))
        }))
        ds <- density(samples[[1]][coords[i], ])
        max_l <- if (!is.null(mcmc_samples)) {
          max(c(max_l, max(density(mcmc_samples[, coords[i]])$y)))
        }
        else {
          max_l
        }
        plot(ds, main = "", xlab = "", ylab = "",
             col = colres[1], ylim = c(0, max_l))
        for (r in 2:nres) {
          lines(density(samples[[r]][coords[i], ],
                        bw = ds$bw), col = colres[r])
        }
        if (!is.null(mcmc_samples)) {
          lines(density(mcmc_samples[, coords[i]]), col = colres[r +
                                                                   1])
        }
      }
      if (i != j) {
        xrange <- range(sapply(1:nres, function(rs) range(c(samples[[rs]][coords[i],
        ], list_pdmp[[rs]]$positions[coords[i], inds]))))
        yrange <- range(sapply(1:nres, function(rs) range(c(samples[[rs]][coords[j],
        ], list_pdmp[[rs]]$positions[coords[j], inds]))))
        plot(list_pdmp[[1]]$positions[coords[i], inds],
             list_pdmp[[1]]$positions[coords[j], inds],
             xlim = xrange, ylim = yrange, type = "l",
             main = "", col = colres[1])
        points(samples[[1]][coords[i], ], samples[[1]][coords[j],
        ], col = colres[1], pch = pch)
        for (r in 2:nres) {
          lines(list_pdmp[[r]]$positions[coords[i], inds],
                list_pdmp[[r]]$positions[coords[j], inds],
                col = colres[r])
          points(samples[[r]][coords[i], ], samples[[r]][coords[j],
          ], col = colres[r], pch = pch)
        }
        if (!is.null(mcmc_samples)) {
          points(mcmc_samples[, coords[i]], mcmc_samples[,
                                                         coords[j]], col = colres[r + 1], pch = pch)
        }
      }
    }
  }
  legend('topleft',legend = c(names(list_pdmp),"mcmc"), col = colres, lwd = 1 )
}
