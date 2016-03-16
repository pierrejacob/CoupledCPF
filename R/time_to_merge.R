#'@export
time_to_merge <- function(observations, model, theta, niterations_max, nparticles, 
                                coupled_resampling, resampling_parameters, with_as,
                                TreeClass){
  datalength <- nrow(observations)
  rand <- model$generate_randomness(nparticles, datalength)
  rand$resampling <- matrix(runif(2*(datalength+1)), ncol = 2)
  xref <- conditional_particle_filter(nparticles, model, theta, observations, rand, TreeClass)
  rand <- model$generate_randomness(nparticles, datalength)
  rand$resampling <- matrix(runif(2*(datalength+1)), ncol = 2)
  xref_tilde <- conditional_particle_filter(nparticles, model, theta, observations, rand, TreeClass)
  
  iteration <- 0
  while (iteration < niterations_max){
    iteration <- iteration + 1
    rand <- model$generate_randomness(nparticles, nrow(observations))
    rand$resampling <- matrix(runif(2*(datalength+1)), ncol = 2)
    
    res <- coupled_conditional_particle_filter(nparticles, model, theta, observations, rand, TreeClass, 
                                               xref, xref_tilde,
                                               coupled_resampling, resampling_parameters,
                                               with_as = with_as)
    xref <- res$new_trajectory1
    xref_tilde <- res$new_trajectory2
    if (isTRUE(all.equal(xref, xref_tilde))){
      break
    }
  }
  return(iteration)
}