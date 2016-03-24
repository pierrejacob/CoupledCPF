# PZ model
# theta = (mu_alpha, sd_alpha, c, e, ml, mq)
# transformed so that all parameters are in R
# theta = (mu_alpha, log sd_alpha, logit c, logit e, logit ml, logit mq)

#'@export
get_pz <- function(){
  logit <- function(z) log(z / (1 - z))
  expit <- function(z) 1 / (1 + exp(-z))
  
  pz_generate_randomness <- function(nparticles, datalength){
    init_rand <- rnorm(2*nparticles, 0, 1)
    transition_rand <- matrix(rnorm(nparticles * datalength, 0, 1), nrow = nparticles, ncol = datalength)
    return(list(init = init_rand, transition = transition_rand))
  }
  
  pz_rinit <- function(nparticles, theta, rand){
    return(exp(matrix(log(2) + rand, ncol = 2)))
  }
  
  pz_rtransition <- function(xparticles, theta, time, rand){
    # nparticles <- nrow(xparticles)
    alphas <- theta[1] + exp(theta[2]) * rand
    xparticles <- pz_transition(xparticles, alphas, time-1, expit(theta[3:6]))
    return(xparticles)
  }
  
  pz_dmeasurement <- function(xparticles, theta, observation) {
    return(dnorm(x = observation, mean = log(xparticles[,1]), sd = 0.2, log = TRUE))
  }
  pzmodel <- list(rinit = pz_rinit,
                  rtransition = pz_rtransition,
                  dmeasurement = pz_dmeasurement,
                  generate_randomness = pz_generate_randomness,
                  dimension = 2)
  return(pzmodel)
}
