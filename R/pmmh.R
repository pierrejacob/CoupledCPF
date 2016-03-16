# pmmh_parameters should contain nparticles, mcmciterations and proposalsd
#'@export
pmmh <- function(pmmh_parameters, model, theta_init, observations){
  current_theta <- theta_init
  theta_dim <- length(theta_init)
  nparticles <- pmmh_parameters$nparticles
  mcmciterations <- pmmh_parameters$mcmciterations
  proposal_covariance <- pmmh_parameters$proposal_covariance
  datalength <- nrow(observations)
  #
  randomness <- model$generate_randomness(nparticles = nparticles, datalength = datalength)
  current_ll <- particle_filter_storeall(nparticles, model, current_theta, observations, randomness)$ll
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
    proposal <- current_theta + fast_rmvnorm(1, rep(0, theta_dim), proposal_covariance)
    randomness <- model$generate_randomness(nparticles = nparticles, datalength = datalength)
    proposal_ll <- try(particle_filter_storeall(nparticles, model, proposal, observations, randomness)$ll)
    if (inherits(proposal_ll, "try-error")){
      proposal_ll <- -Inf
    }
    proposal_posterior <- proposal_ll + model$dprior(proposal)
    if (log(runif(1)) < (proposal_posterior - current_posterior)){
      current_theta <- proposal
      current_ll <- proposal_ll
      current_posterior <- proposal_posterior
      pmmh_naccepts <- pmmh_naccepts + 1
    }
    pmmh_chain[iteration,] <- current_theta
    loglikelihoods[iteration] <- current_ll
  }
  return(list(chain = pmmh_chain, acceptance_rate = pmmh_naccepts / mcmciterations,
              loglikelihoods = loglikelihoods))
}