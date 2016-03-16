#'@export
L2_distance <- function(observations, model, theta, niterations_max, nparticles, 
                        coupled_resampling, final_resampling, resampling_parameters, with_as,
                        TreeClass){
  datalength <- nrow(observations)
  rand <- model$generate_randomness(nparticles, datalength)
  rand$resampling <- matrix(runif(2*(datalength+1)), ncol = 2)
  xref <- conditional_particle_filter(nparticles, model, theta, observations, rand, TreeClass)
  rand <- model$generate_randomness(nparticles, datalength)
  rand$resampling <- matrix(runif(2*(datalength+1)), ncol = 2)
  xref_tilde <- conditional_particle_filter(nparticles, model, theta, observations, rand, TreeClass)
  
  iteration <- 0
  distance <- matrix(0, ncol = niterations_max + 1, nrow = datalength + 1)
  d <- (xref - xref_tilde)^2
  distance[,1] <- sqrt(apply(d, 1, sum))
  
  while (iteration < niterations_max){
    iteration <- iteration + 1
    rand <- model$generate_randomness(nparticles, nrow(observations))
    rand$resampling <- matrix(runif(2*(datalength+1)), ncol = 2)
    res <- coupled_conditional_particle_filter(nparticles, model, theta, observations, rand,
                                               xref, xref_tilde,
                                               coupled_resampling, final_resampling, resampling_parameters,
                                               with_as, TreeClass)
    xref <- res$new_trajectory1
    xref_tilde <- res$new_trajectory2
    d <- (xref - xref_tilde)^2
    distance[,iteration+1] <- sqrt(apply(d, 1, sum))
    if (isTRUE(all.equal(xref, xref_tilde))){
      break
    }
  }
  return(distance)
}
