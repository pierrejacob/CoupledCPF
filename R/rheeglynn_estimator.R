#'@export
geom_unbiased_estimator <- function(observations, model, theta, nparticles,
                               proba_geom, coupled_resampling, final_resampling, resampling_parameters, with_as, TreeClass){
  datalength <- nrow(observations)
  geom_niterations <- rgeom(1, proba_geom)
    survival_probabilities <- 1 - pgeom(1:geom_niterations, prob = proba_geom)
  rand <- model$generate_randomness(nparticles, datalength)
  rand$resampling <- matrix(runif(2*(datalength+1)), ncol = 2)
  xref <- conditional_particle_filter(nparticles, model, theta, observations, rand, TreeClass)
  rand <- model$generate_randomness(nparticles, datalength)
  rand$resampling <- matrix(runif(2*(datalength+1)), ncol = 2)
  xref_tilde <- conditional_particle_filter(nparticles, model, theta, observations, rand, TreeClass)
  # xref2 <- conditional_particle_filter(nparticles, model, theta, observations, rand, TreeClass, ref_trajectory = xref2,
  # with_as = with_as)
  iteration <- 0
  estimate <- xref
  if (geom_niterations == 0){
    return(list(estimate = estimate, iteration = iteration, geom_niterations = geom_niterations))
  }
  iteration <- iteration + 1
  rand <- model$generate_randomness(nparticles, datalength)
  rand$resampling <- matrix(runif(2*(datalength+1)), ncol = 2)
  xref <- conditional_particle_filter(nparticles, model, theta, observations, rand, TreeClass,
                                      xref, with_as)
  estimate <- estimate + (xref_tilde - xref) / survival_probabilities[1]
  if (geom_niterations == 1){
    return(list(estimate = estimate, iteration = iteration, geom_niterations = geom_niterations))
  }
  while (iteration < geom_niterations){
    iteration <- iteration + 1
    rand <- model$generate_randomness(nparticles, nrow(observations))
    rand$resampling <- matrix(runif(2*(datalength+1)), ncol = 2)
    res <- coupled_conditional_particle_filter(nparticles, model, theta, observations, rand,
                                               xref, xref_tilde,
                                               coupled_resampling, final_resampling, resampling_parameters,
                                               with_as,  TreeClass)
    xref <- res$new_trajectory1
    xref_tilde <- res$new_trajectory2
    estimate <- estimate + (xref - xref_tilde) / survival_probabilities[iteration]
    if (isTRUE(all.equal(xref, xref_tilde))){
      break
    }
  }
  return(list(estimate = estimate, iteration = iteration, geom_niterations = geom_niterations))
}