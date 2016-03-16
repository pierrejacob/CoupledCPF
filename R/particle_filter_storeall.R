#'@export
## particle filter that takes the whole randomness as an argument
## and stores all the generated particles, weights and ancestor 
particle_filter_storeall <- function(nparticles, model, theta, observations, randomness){
  datalength <- nrow(observations)
  # initialization
  xparticles <- model$rinit(nparticles, theta, randomness$init)
  normweights <- rep(1/nparticles, nparticles)
  ll <- 0
  #
  xhistory <- rep(list(matrix(nrow = nparticles, ncol = ncol(xparticles))), datalength + 1)
  whistory <- rep(list(rep(0, nparticles)), datalength + 1)
  ahistory <- rep(list(rep(0, nparticles)), datalength)
  
  xhistory[[1]] <- xparticles
  whistory[[1]] <- normweights
  # step t > 1
  for (time in 1:datalength){
    ancestors <- systematic_resampling_n(normweights, nparticles, runif(1))
    ahistory[[time]] <- ancestors
    xparticles <- xparticles[ancestors,]
    xparticles <- model$rtransition(xparticles, theta, time, randomness$transition[,time])
    logw <- model$dmeasurement(xparticles, theta, observations[time,])
    maxlw <- max(logw)
    w <- exp(logw - maxlw)
    # update log likelihood estimate
    ll <- ll + maxlw + log(mean(w))
    normweights <- w / sum(w)
    #
    xhistory[[time+1]] <- xparticles
    whistory[[time+1]] <- normweights
  }
  return(list(ll = ll, xhistory = xhistory, whistory = whistory, ahistory = ahistory))
}
