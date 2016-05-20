#'@rdname get_gss
#'@title Nonlinear growth model as in Gordon Salmond Smith 1993
#'@description This function returns a list with objects such as
#'* rinit, rinit_rand to sample from the initial distribution
#'* rtransition, rtransition_rand to sample from the transition
#'* dtransition to evaluate the transition density
#'* dmeasurement to evaluate the measurement density
#'* dimension, which represents the dimension of the latent process
#'@return A list
#'@export
get_gss <- function(){
  ## Gordon Salmond Smith model
  rinit <- function(nparticles, theta, rand, ...){
    return(matrix(sqrt(2) * rand, nrow = 1))
  }
  #
  rinit_rand <- function(nparticles, theta){
    return(rnorm(nparticles))
  }
  #
  rtransition <- function(xparticles, theta, time, rand, ...){
    means <- 0.5 * xparticles + 25 * xparticles / (1 + xparticles^2) + 8 * cos(1.2*(time - 1))
    return(matrix(means + sqrt(10) * rand, nrow = 1))
  }
  #
  rtransition_rand <- function(nparticles, theta){
    return(rnorm(nparticles))
  }
  #
  dtransition <- function(next_x, xparticles, theta, time, ...){
    means <- 0.5 * xparticles + 25 * xparticles / (1 + xparticles^2) + 8 * cos(1.2*(time - 1))
    return(dnorm(next_x, mean = means, sd = sqrt(10), log = TRUE))
  }
  
  dmeasurement <- function(xparticles, theta, observation, ...){
    return(dnorm(observation, mean = (xparticles^2) / 20, sd = 1, log = TRUE))
  }
  #
  precompute <- function(...){
    return(list())
  }
  #
  gss <- list(rinit = rinit, rinit_rand = rinit_rand, rtransition = rtransition,
              rtransition_rand = rtransition_rand,
              dtransition = dtransition, 
              dmeasurement = dmeasurement, precompute = precompute, dimension = 1)
  return(gss)
}
