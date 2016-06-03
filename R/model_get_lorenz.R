# lorenz model
# theta = (mu_alpha, sd_alpha, c, e, ml, mq)
# transformed so that all parameters are in R
# theta = (mu_alpha, log sd_alpha, logit c, logit e, logit ml, logit mq)
#'@rdname get_lorenz
#'@title Lorenz 96 model
#'@description This function returns a list with objects such as
#'* rinit  to sample from the initial distribution
#'* rtransition to sample from the transition
#'* dtransition to evaluate the transition density
#'* dmeasurement to evaluate the measurement density
#'* generate_randomness to evaluate the measurement density
#'* perturb_randomness to evaluate the measurement density
#'* dimension, which represents the dimension of the latent process
#'@return A list 
#'@export
get_lorenz <- function(){
  rinit <- function(nparticles, theta, rand, ...){
    unif <- -1 + 4 * pnorm(rand)
    return(matrix(unif, ncol = nparticles))
  }
  rinit_rand <- function(nparticles, theta){
    return(rnorm(8*nparticles))
  }
  #  
  rtransition <- function(xparticles, theta, time, rand, ...){
    nparticles <- ncol(xparticles)
    xparticles <- lorenz_transition(xparticles, time_start = 0.05*(time - 1), time_end = 0.05*time, dt = 0.05, theta[1])
    xparticles <- xparticles + matrix(theta[2] * rand, ncol = nparticles)
    return(xparticles)
  }
  #
  rtransition_rand <- function(nparticles, theta){
    return(rnorm(8*nparticles))
  }
  #
  dmeasurement <- function(xparticles, theta, observation, ...) {
    if (is.na(observation[1])){
      return(rep(0, ncol(xparticles)))
    } else {
      return(fast_dmvnorm_transpose_cholesky(xparticles[1:4,], observation[1:4], diag(0.5, 4, 4)))
    }
  }
  #
  precompute <- function(...) NULL
  #
  lorenz_model <- list(rinit = rinit, rinit_rand = rinit_rand,
                       rtransition = rtransition, rtransition_rand = rtransition_rand,
                   dmeasurement = dmeasurement, precompute = precompute,
                   dimension = 8)
  return(lorenz_model)
}
