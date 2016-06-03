#'@rdname particlefilter
#'@title Particle Filter
#'@description Particle filter with storage of trajectories
#'@param nparticles a number of particles
#'@param model a list representing a model, for instance as given by \code{\link{get_ar}}.
#'@param theta a parameter to give to the model functions
#'@param observations a matrix of observations of size datalength x dimension(observation)
#'@return A list with all the trajectories, stored in an array of dimension dimx x (T+1) x N,
#'and the associated normalized weights
#'@export
particlefilter <- function(nparticles, model, theta, observations){
  datalength <- nrow(observations)
  # create tree representation of the trajectories
  Tree <- new(TreeClass, nparticles, 10*nparticles*model$dimension, model$dimension)
  # initialization
  model_precomputed <- model$precompute(theta)
  xparticles <- model$rinit(nparticles, theta, model$rinit_rand(nparticles, theta), model_precomputed)
  Tree$init(xparticles)
  #  
  normweights <- rep(1/nparticles, nparticles)
  # step t > 1
  for (time in 1:datalength){
    ancestors <- multinomial_resampling_n(normweights, nparticles)
    # if no observation or first time, no resampling
    if (time == 1 || (time > 1 && is.na(observations[time-1,1]))){
      ancestors <- 1:nparticles
    }
    xparticles <- xparticles[,ancestors]
    if (is.null(dim(xparticles))) xparticles <- matrix(xparticles, nrow = model$dimension)
    xparticles <- model$rtransition(xparticles, theta, time, model$rtransition_rand(nparticles, theta), model_precomputed)
    #    
    logw <- model$dmeasurement(xparticles, theta, observations[time,], model_precomputed)
    maxlw <- max(logw)
    w <- exp(logw - maxlw)
    normweights <- w / sum(w)
    #
    Tree$update(xparticles, ancestors - 1)
  }
  trajectories <- array(dim = c(model$dimension, datalength + 1, nparticles))
  for (k in 0:(nparticles-1)){
    trajectories[,,k+1] <- Tree$get_path(k)
  }
  return(list(trajectories = trajectories, weights = normweights))
}