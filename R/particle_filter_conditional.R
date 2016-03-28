#'@rdname conditional_particle_filter
#'@title conditional_particle_filter
#'@description runs a conditional particle filter, with or without ancestor sampling
#'@export

conditional_particle_filter <- function(nparticles, model, theta, observations, randomness, TreeClass, ref_trajectory, with_as = FALSE){
  datalength <- nrow(observations)
  # create tree representation of the trajectories
  Tree <- new(TreeClass, nparticles, 10*nparticles, model$dimension)
  # initialization
  xparticles <- model$rinit(nparticles, theta, randomness$init)
  if (!missing(ref_trajectory)){
    xparticles[nparticles,] <- ref_trajectory[1,]
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
    ancestors <- systematic_resampling_n(normweights, nparticles, randomness$resampling[time,1])
    xparticles <- xparticles[ancestors,]
    xparticles <- model$rtransition(xparticles, theta, time, randomness$transition[,time])
    if (!missing(ref_trajectory)){
      xparticles[nparticles,] <- ref_trajectory[time+1,]
      if (with_as){
        # % Ancestor sampling
        logm <- model$dtransition(ref_trajectory[time+1,], x_last, theta, time)
        logm <- log(normweights) + logm
        w_as <- exp(logm - max(logm))
        w_as <- w_as / sum(w_as)
        ancestors[nparticles] <- systematic_resampling_n(w_as, 1, randomness$resampling[time,2])
        x_last <- xparticles
      } else {
        ancestors[nparticles] <- nparticles
      }
    }
    #    
    logw <- model$dmeasurement(xparticles, theta, observations[time,])
    maxlw <- max(logw)
    w <- exp(logw - maxlw)
    normweights <- w / sum(w)
    #
    Tree$update(xparticles, ancestors - 1)    
  }
  new_trajectory <- Tree$get_path(systematic_resampling_n(normweights, 1, randomness$resampling[datalength+1,1])-1)
  return(new_trajectory)
}