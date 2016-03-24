# pmmh_parameters should contain nparticles, mcmciterations and proposalsd, rho_perturb
# resampling_scheme and resampling_parameters
# model should contain perturb_randomness, taking randomness and rho_perturb as arguments
#'@export
coupled_pmmh <- function(pmmh_parameters, model, theta_init, observations){
  current_theta <- theta_init
  theta_dim <- length(theta_init)
  nparticles <- pmmh_parameters$nparticles
  mcmciterations <- pmmh_parameters$mcmciterations
  proposal_covariance <- pmmh_parameters$proposal_covariance
  rho_perturb <- pmmh_parameters$rho_perturb
  coupled_resampling <- pmmh_parameters$resampling_scheme
  resampling_parameters <- pmmh_parameters$resampling_parameters
  datalength <- nrow(observations)
  #
  current_randomness <- model$generate_randomness(nparticles = nparticles, datalength = datalength)
  current_system <- particle_filter_storeall(nparticles, model, current_theta, observations, 
                                             current_randomness)
  current_ll <- current_system$ll
  current_posterior <- current_ll + model$dprior(current_theta)
  pmmh_naccepts <- 0
  pmmh_chain <- matrix(nrow = mcmciterations, ncol = theta_dim)
  pmmh_chain[1,] <- current_theta
  loglikelihoods <- rep(0, mcmciterations)
  loglikelihoods[1] <- current_ll
  for (iteration in 2:mcmciterations){
    if (iteration %% 100 == 1){
      cat("iteration: ", iteration, " / ", mcmciterations, "\n")
      cat("acceptance rate: ", pmmh_naccepts / iteration * 100, "%\n")
    }
    proposal <- current_theta + as.numeric(fast_rmvnorm(1, rep(0, theta_dim), proposal_covariance))
    proposal_randomness <- model$perturb_randomness(current_randomness, rho_perturb)
    proposal_system <- try(particle_filter_given_ref(nparticles, model, proposal, observations, proposal_randomness,
                                                 coupled_resampling, resampling_parameters, current_system))
    if (inherits(proposal_system, "try-error")){
      proposal_ll <- -Inf
    } else {
      proposal_ll <- proposal_system$ll
    }
    proposal_posterior <- proposal_ll + model$dprior(proposal)
    if (log(runif(1)) < (proposal_ll - current_ll)){
      current_theta <- proposal
      current_ll <- proposal_ll
      current_posterior <- proposal_posterior
      current_randomness <- proposal_randomness
      current_system <- proposal_system
      pmmh_naccepts <- pmmh_naccepts + 1
    }
    pmmh_chain[iteration,] <- current_theta
    loglikelihoods[iteration] <- current_ll
  }
  return(list(chain = pmmh_chain, acceptance_rate = pmmh_naccepts / mcmciterations,
              loglikelihoods = loglikelihoods))
}
