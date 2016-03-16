# hidden AR model
#'@export
get_ar1 <- function(){
  ar1 <- list(rinit = function(nparticles, theta, rand){
    return(matrix(rand, ncol = 1))
  },
  rtransition = function(xparticles, theta, time, rand){
    return(theta * xparticles + matrix(rand, ncol = 1))
  },
  dtransition = function(next_x, xparticles, theta, time){
    return(dnorm(next_x, mean = theta * xparticles, sd = 1, log = TRUE))
  },
  dmeasurement = function(xparticles, theta, observation){
    return(dnorm(observation, mean = xparticles, sd = 1, log = TRUE))
  },
  generate_randomness = function(nparticles, datalength){
    init_rand <- rnorm(nparticles, 0, 1)
    transition_rand <- matrix(rnorm(nparticles * datalength, 0, 1), nrow = nparticles, ncol = datalength)
    return(list(init = init_rand, transition = transition_rand))
  },
  dimension = 1)
  return(ar1)
}
