#'@export
get_ar5 <- function(){
  ###
  rinit <- function(nparticles, theta, rand){
    return(matrix(rand, ncol = 5))
  }
  rtransition <- function(xparticles, theta, time, rand){
    A <- create_A(theta, 5)
    return(xparticles %*% A + matrix(rand, ncol = 5))
  }
  
  dtransition <- function(next_x, xparticles, theta, time){
    A <- create_A(theta, 5)
    return(fast_dmvnorm(xparticles %*% A, next_x, diag(1, 5, 5)))
  }
  
  dmeasurement <- function(xparticles, theta, observation){
    return(fast_dmvnorm(xparticles, observation, diag(1, 5, 5)))
  }
  
  generate_randomness <- function(nparticles, datalength){
    init_rand <- rnorm(nparticles * 5, 0, 1)
    transition_rand <- matrix(rnorm(nparticles * datalength * 5, 0, 1), nrow = nparticles * 5, ncol = datalength)
    return(list(init = init_rand, transition = transition_rand))
  }
  ar5 <- list(rinit = rinit, rtransition = rtransition, dtransition = dtransition, 
                          dmeasurement = dmeasurement, generate_randomness = generate_randomness,
                          dimension = 5)
  return(ar5)
}
