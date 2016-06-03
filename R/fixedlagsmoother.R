#'@rdname fixedlagsmoother
#'@title Fixed-lag smoother
#'@description Estimates the mean of the smoothing distributio at all times, by fixed-lag smoothing
#'@param nparticles a number of particles
#'@param lag a lag value
#'@param model a list representing a model, for instance as given by \code{\link{get_ar}}.
#'@param theta a parameter to give to the model functions
#'@param observations a matrix of observations of size datalength x dimension(observation)
#'@return a list with the estimator of the smoothing means, and the second moments of the smoothing distributions
#'@export
fixedlagsmoother <- function(nparticles, lag, model, theta, observations){
  datalength <- nrow(observations)
  if (lag > datalength) lag <- datalength
  smoothing_means <- matrix(nrow = datalength + 1, ncol = model$dimension)
  smoothing_secondmoments <- matrix(nrow = datalength + 1, ncol = model$dimension)
  # create tree representation of the trajectories
  Tree <- new(TreeClass, nparticles, 10*nparticles*model$dimension, model$dimension)
  # initialization
  model_precomputed <- model$precompute(theta)
  xparticles <- model$rinit(nparticles, theta, model$rinit_rand(nparticles, theta), model_precomputed)
  Tree$init(xparticles)
  #  
  normweights <- rep(1/nparticles, nparticles)
  indexsmoothing <- 1
  if (lag == 0){
    smoothing_means[indexsmoothing,] <- apply(xparticles, 1, function(v) sum(v * normweights))
    smoothing_secondmoments[indexsmoothing,] <- apply(xparticles, 1, function(v) sum(v^2 * normweights))
    indexsmoothing <- indexsmoothing + 1
  }
  # step t > 1
  for (time in 1:datalength){
    ancestors <- multinomial_resampling_n(normweights, nparticles)
    # if no observation or first time, no resampling
    if (time == 1 || (time > 1 && is.na(observations[time-1,1]))){
      ancestors <- 1:nparticles
    }
    xparticles <- xparticles[,ancestors]
    if (is.null(dim(xparticles))) xparticles <- matrix(xparticles, nrow = dimension)
    xparticles <- model$rtransition(xparticles, theta, time, model$rtransition_rand(nparticles, theta), model_precomputed)
    #    
    logw <- model$dmeasurement(xparticles, theta, observations[time,], model_precomputed)
    maxlw <- max(logw)
    w <- exp(logw - maxlw)
    normweights <- w / sum(w)
    #
    # print(apply(xparticles, 1, function(v) sum(v * normweights)))
    Tree$update(xparticles, ancestors - 1)
    # smoothing calculations
    if (time - lag >= 0){
      # cat("current time: ", time, "\n")
      # cat("smoothing mean at time ", time - lag, "\n")
      xpast <- Tree$retrieve_xgeneration(lag)
      smoothing_means[indexsmoothing,] <- apply(xpast, 1, function(v) sum(v * normweights))
      smoothing_secondmoments[indexsmoothing,] <- apply(xpast, 1, function(v) sum(v^2 * normweights))
      indexsmoothing <- indexsmoothing + 1
    }
  }
  if (lag >= 1){
    for (reversetime in lag:1){
      xpast <- Tree$retrieve_xgeneration(reversetime-1)
      smoothing_means[indexsmoothing,] <- apply(xpast, 1, function(v) sum(v * normweights))
      smoothing_secondmoments[indexsmoothing,] <- apply(xpast, 1, function(v) sum(v^2 * normweights))
      indexsmoothing <- indexsmoothing + 1
    }
  }
  return(list(smoothing_means = smoothing_means, smoothing_secondmoments = smoothing_secondmoments))
}