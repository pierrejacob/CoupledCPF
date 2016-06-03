#'@rdname CPF
#'@title Conditional Particle Filter
#'@description Runs a conditional particle filter, with or without ancestor sampling
#'@param nparticles number of particles
#'@param model a list representing a model, for instance as given by \code{\link{get_ar}}.
#'@param theta a parameter to give to the model functions
#'@param observations a matrix of observations of size datalength x dimension(observation)
#'@param ref_trajectory a reference trajectory, of size dimension(process) x datalength; if missing, runs a standard
#'particle filter.
#'@param with_as whether ancestor sampling should be used (TRUE/FALSE)
#'@return A new trajectory.
#'@export

CPF <- function(nparticles, model, theta, observations, ref_trajectory = NULL, with_as = FALSE){
  datalength <- nrow(observations)
  # create tree representation of the trajectories
  Tree <- new(TreeClass, nparticles, 10*nparticles*model$dimension, model$dimension)
  # initialization
  model_precomputed <- model$precompute(theta)
  xparticles <- model$rinit(nparticles, theta, model$rinit_rand(nparticles, theta), model_precomputed)
  if (!is.null(ref_trajectory)){
    xparticles[,nparticles] <- ref_trajectory[,1]
  }
  Tree$init(xparticles)
  #  
  normweights <- rep(1/nparticles, nparticles)
  # if ancestor sampling, needs to keep the last generation of particles at hand
  if (with_as){
    x_last <- xparticles
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
    if (!is.null(ref_trajectory)){
      xparticles[,nparticles] <- ref_trajectory[,time+1]
      if (with_as){
        # Ancestor sampling
        logm <- model$dtransition(ref_trajectory[,time+1], x_last, theta, time, model_precomputed)
        logm <- log(normweights) + logm
        w_as <- exp(logm - max(logm))
        w_as <- w_as / sum(w_as)
        ancestors[nparticles] <- systematic_resampling_n(w_as, 1, runif(1))
        x_last <- xparticles
      } else {
        ancestors[nparticles] <- nparticles
      }
    }
    #    
    logw <- model$dmeasurement(xparticles, theta, observations[time,], model_precomputed)
    maxlw <- max(logw)
    w <- exp(logw - maxlw)
    normweights <- w / sum(w)
    #
    Tree$update(xparticles, ancestors - 1)    
  }
  new_trajectory <- Tree$get_path(systematic_resampling_n(normweights, 1, runif(1))-1)
  return(new_trajectory)
}