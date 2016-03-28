
particle_filter_given_uresampling <- function(nparticles, model, theta, observations, randomness, u_resampling){
  uniforms <- pnorm(u_resampling)
  datalength <- nrow(observations)
  # initialization
  xparticles <- model$rinit(nparticles, theta, randomness$init)
  normweights <- rep(1/nparticles, nparticles)
  ll <- 0
  # step t > 1
  for (time in 1:datalength){
    horder <- hilbert_order(xparticles)
    nw_sorted <- normweights[horder]
    ancestors <- systematic_resampling_n(nw_sorted, nparticles, uniforms[time])
    ancestors <- horder[ancestors]
    #
    xparticles <- matrix(xparticles[ancestors,], ncol = model$dimension)
    xparticles <- model$rtransition(xparticles, theta, time, randomness$transition[,time])
    logw <- model$dmeasurement(xparticles, theta, observations[time,])
    maxlw <- max(logw)
    w <- exp(logw - maxlw)
    # update log likelihood estimate
    ll <- ll + maxlw + log(mean(w))
    normweights <- w / sum(w)
    #
  }
  return(ll)
}

#'@export
hilbert_pmmh <- function(pmmh_parameters, model, theta_init, observations){
  current_theta <- theta_init
  theta_dim <- length(theta_init)
  nparticles <- pmmh_parameters$nparticles
  mcmciterations <- pmmh_parameters$mcmciterations
  proposal_covariance <- pmmh_parameters$proposal_covariance
  rho_perturb <- pmmh_parameters$rho_perturb
  v_perturb <- sqrt(1 - rho_perturb^2)
  datalength <- nrow(observations)
  current_u_resampling <- rnorm(datalength)
  #
  current_randomness <- model$generate_randomness(nparticles = nparticles, datalength = datalength)
  current_ll <- particle_filter_given_uresampling(nparticles, model, current_theta, observations,
                                                  current_randomness, current_u_resampling)
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
    proposal <- current_theta + fast_rmvnorm(1, rep(0, theta_dim), proposal_covariance)[1,]
    proposal_prior <- model$dprior(proposal)
    if (!is.infinite(proposal_prior)){
      proposal_randomness <- model$perturb_randomness(current_randomness, rho_perturb)
      proposal_u_resampling <- rho_perturb * current_u_resampling + v_perturb * rnorm(datalength)
      proposal_ll <- try(particle_filter_given_uresampling(nparticles, model, proposal, observations, 
                                                           proposal_randomness, proposal_u_resampling))
      if (inherits(proposal_ll, "try-error")){
        proposal_ll <- -Inf
      } else {
        if (is.na(proposal_ll)){
          proposal_ll <- -Inf
        }
      }
      proposal_posterior <- proposal_ll + proposal_prior
      if (log(runif(1)) < (proposal_posterior - current_posterior)){
        current_theta <- proposal
        current_ll <- proposal_ll
        current_posterior <- proposal_posterior
        current_randomness <- proposal_randomness
        current_u_resampling <- proposal_u_resampling
        pmmh_naccepts <- pmmh_naccepts + 1
      }
    }
    pmmh_chain[iteration,] <- current_theta
    loglikelihoods[iteration] <- current_ll
  }
  return(list(chain = pmmh_chain, acceptance_rate = pmmh_naccepts / mcmciterations,
              loglikelihoods = loglikelihoods))
  
}
