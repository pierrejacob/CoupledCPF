#'@rdname CPF_coupled
#'@title Coupled Conditional Particle Filter
#'@description Runs a coupled conditional particle filter, with or without ancestor sampling
#'@param nparticles number of particles
#'@param model a list representing a model, for instance as given by \code{\link{get_ar}}.
#'@param theta a parameter to give to the model functions
#'@param observations a matrix of observations of size datalength x dimension(observation)
#'@param ref_trajectory1 a first reference trajectory, of size dimension(process) x datalength
#'@param ref_trajectory2 a second reference trajectory, of size dimension(process) x datalength
#'@param coupled_resampling a coupled resampling scheme, such as \code{\link{CR_indexmatching}}.
#'@param with_as whether ancestor sampling should be used (TRUE/FALSE)
#'@return A pair of new trajectories.
#'@export

CPF_coupled <- function(nparticles, model, theta, observations, ref_trajectory1, ref_trajectory2, 
                        coupled_resampling, with_as = FALSE){
  #
  datalength <- nrow(observations)
  # create tree representation of the trajectories
  Tree1 <- new(TreeClass, nparticles, 10*nparticles*model$dimension, model$dimension)
  Tree2 <- new(TreeClass, nparticles, 10*nparticles*model$dimension, model$dimension)
  # initialization
  model_precomputed <- model$precompute(theta)
  init_rand <- model$rinit_rand(nparticles, theta)
  xparticles1 <- model$rinit(nparticles, theta, init_rand, model_precomputed)
  xparticles1[,nparticles] <- ref_trajectory1[,1]
  Tree1$init(xparticles1)
  normweights1 <- rep(1/nparticles, nparticles)
  #  
  xparticles2 <- model$rinit(nparticles, theta, init_rand, model_precomputed)
  xparticles2[,nparticles] <- ref_trajectory2[,1]
  Tree2$init(xparticles2)
  normweights2 <- rep(1/nparticles, nparticles)
  #
  # if ancestor sampling, needs to keep the last generation of particles at hand
  if (with_as){
    x_last1 <- xparticles1
    x_last2 <- xparticles2
  }
  # step t > 1
  for (time in 1:datalength){
    ancestors <- coupled_resampling(xparticles1, xparticles2, normweights1, normweights2)
    ancestors1 <- ancestors[,1]
    ancestors2 <- ancestors[,2]
    # if no observation or first time, no resampling
    if (time == 1 || (time > 1 && is.na(observations[time-1,1]))){
      ancestors1 <- 1:nparticles
      ancestors2 <- 1:nparticles
    }
    #
    xparticles1 <- xparticles1[,ancestors1]
    xparticles2 <- xparticles2[,ancestors2]
    
    if (is.null(dim(xparticles1))) xparticles1 <- matrix(xparticles1, nrow = dimension)
    if (is.null(dim(xparticles2))) xparticles2 <- matrix(xparticles2, nrow = dimension)
    #
    transition_rand <- model$rtransition_rand(nparticles, theta)
    xparticles1 <- model$rtransition(xparticles1, theta, time, transition_rand, model_precomputed)
    xparticles2 <- model$rtransition(xparticles2, theta, time, transition_rand, model_precomputed)
    if (is.null(dim(xparticles1))) xparticles1 <- matrix(xparticles1, nrow = dimension)
    if (is.null(dim(xparticles2))) xparticles2 <- matrix(xparticles2, nrow = dimension)
    #
    xparticles1[,nparticles] <- ref_trajectory1[,time+1]
    xparticles2[,nparticles] <- ref_trajectory2[,time+1]
    if (with_as){
      # % Ancestor sampling
      logm1 <- model$dtransition(ref_trajectory1[,time+1], x_last1, theta, time, model_precomputed)
      logm1 <- log(normweights1) + logm1
      w_as1 <- exp(logm1 - max(logm1))
      w_as1 <- w_as1 / sum(w_as1)
      unif_resampling_as <- runif(1)
      ancestors1[nparticles] = systematic_resampling_n(w_as1, 1, unif_resampling_as)
      x_last1 <- xparticles1
      #
      logm2 <- model$dtransition(ref_trajectory2[,time+1], x_last2, theta, time, model_precomputed)
      logm2 <- log(normweights2) + logm2
      w_as2 <- exp(logm2 - max(logm2))
      w_as2 <- w_as2 / sum(w_as2)
      ancestors2[nparticles] = systematic_resampling_n(w_as2, 1, unif_resampling_as)
      x_last2 <- xparticles2
    } else {
      ancestors1[nparticles] <- nparticles
      ancestors2[nparticles] <- nparticles
    }
    #    
    logw1 <- model$dmeasurement(xparticles1, theta, observations[time,], model_precomputed)
    logw2 <- model$dmeasurement(xparticles2, theta, observations[time,], model_precomputed)
    #
    maxlw1 <- max(logw1)
    w1 <- exp(logw1 - maxlw1)
    normweights1 <- w1 / sum(w1)
    #
    maxlw2 <- max(logw2)
    w2 <- exp(logw2 - maxlw2)
    normweights2 <- w2 / sum(w2)    
    #
    Tree1$update(xparticles1, ancestors1 - 1)    
    Tree2$update(xparticles2, ancestors2 - 1)    
  }
  u <- runif(1)
  k_path1 <- systematic_resampling_n(normweights1, 1, u)
  k_path2 <- systematic_resampling_n(normweights2, 1, u)
  ##
  new_trajectory1 <- Tree1$get_path(k_path1 - 1)
  new_trajectory2 <- Tree2$get_path(k_path2 - 1)
  return(list(new_trajectory1 = new_trajectory1, new_trajectory2 = new_trajectory2))
}
