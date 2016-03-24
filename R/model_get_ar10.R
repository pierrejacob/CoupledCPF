#'@export
get_ar10 <- function(){
  ###
  rinit <- function(nparticles, theta, rand){
    return(matrix(rand, ncol = 10))
  }
  rtransition <- function(xparticles, theta, time, rand){
    A <- create_A(theta, 10)
    return(xparticles %*% A + matrix(rand, ncol = 10))
  }
  
  dtransition <- function(next_x, xparticles, theta, time){
    A <- create_A(theta, 10)
    return(fast_dmvnorm(xparticles %*% A, next_x, diag(1, 10, 10)))
  }
  
  dmeasurement <- function(xparticles, theta, observation){
    return(fast_dmvnorm(xparticles, observation, diag(1, 10, 10)))
  }
  
  generate_randomness <- function(nparticles, datalength){
    init_rand <- rnorm(nparticles * 10, 0, 1)
    transition_rand <- matrix(rnorm(nparticles * datalength * 10, 0, 1), nrow = nparticles * 10, ncol = datalength)
    return(list(init = init_rand, transition = transition_rand))
  }
  
  ar10 <- list(rinit = rinit, rtransition = rtransition, dtransition = dtransition, 
                          dmeasurement = dmeasurement, generate_randomness = generate_randomness,
                          dimension = 10)
  return(ar10)
}
