#'@rdname get_ar
#'@title Hidden auto-regressive model
#'@description This function returns a list with objects such as
#'* rinit, rinit_rand to sample from the initial distribution
#'* rtransition, rtransition_rand to sample from the transition
#'* dtransition to evaluate the transition density
#'* dmeasurement to evaluate the measurement density
#'* dimension, which represents the dimension of the latent process
#'@return A list 
#'@export
get_ar <- function(dimension){
  #
  rinit <- function(nparticles, theta, rand, precomputed, ...){
    return(matrix(rand, nrow = dimension))
  }
  rinit_rand <- function(nparticles, theta){
    return(rnorm(nparticles * dimension))
  }
  #
  rtransition <- function(xparticles, theta, time, rand, precomputed, ...){
    return(precomputed$A %*% xparticles + rand)
  }
  #
  rtransition_rand <- function(nparticles, theta){
    return(rnorm(nparticles * dimension))
  }
  #
  dtransition <- function(next_x, xparticles, theta, time, precomputed, ...){
    return(dmvnorm_transpose_cholesky(precomputed$A %*% xparticles, next_x, precomputed$di))
  }
  
  dmeasurement <- function(xparticles, theta, observation, precomputed, ...){
    return(dmvnorm_transpose_cholesky(xparticles, observation, precomputed$di))
  }
  
  precompute <- function(theta){
    A <- create_A(theta, dimension)
    return(list(A = A, di = diag(1, dimension, dimension)))
  }
  
  ar_model <- list(rinit = rinit, rinit_rand = rinit_rand, rtransition = rtransition,
                   rtransition_rand = rtransition_rand,
                   dtransition = dtransition, 
                   dmeasurement = dmeasurement, precompute = precompute, dimension = dimension)
  return(ar_model)
}
