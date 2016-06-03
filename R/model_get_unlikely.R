#'@export
get_unlikely <- function(){
  #
  rinit <- function(nparticles, theta, rand, precomputed, ...){
    return(matrix(theta[4] * rand, nrow = 1))
  }
  rinit_rand <- function(nparticles, theta){
    return(rnorm(nparticles))
  }
  #
  rtransition <- function(xparticles, theta, time, rand, precomputed, ...){
    return(theta[1] * xparticles + theta[2] * rand)
  }
  #
  rtransition_rand <- function(nparticles, theta){
    return(rnorm(nparticles))
  }
  #
  dtransition <- function(next_x, xparticles, theta, time, precomputed, ...){
    return(dnorm(next_x, theta[1] * xparticles, sd = theta[2], log = TRUE))
  }
  
  dmeasurement <- function(xparticles, theta, observation, precomputed, ...){
    if (is.na(observation[1])){
      return(rep(0, ncol(xparticles)))
    } else {
      return(dnorm(xparticles, observation, sd = theta[3], log = TRUE))
    }
  }
  
  precompute <- function(theta){
    return(list())
  }
  
  ar_model <- list(rinit = rinit, rinit_rand = rinit_rand, rtransition = rtransition,
                   rtransition_rand = rtransition_rand,
                   dtransition = dtransition, 
                   dmeasurement = dmeasurement, precompute = precompute, dimension = 1)
  return(ar_model)
}
