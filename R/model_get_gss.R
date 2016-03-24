#'@export
get_gss <- function(){
  ## Gordon Salmond Smith model
  rinit <- function(nparticles, theta, rand){
    return(matrix(sqrt(2) * rand, ncol = 1))
  }
  rtransition <- function(xparticles, theta, time, rand){
    means <- 0.5 * xparticles + 25 * xparticles / (1 + xparticles^2) + 8 * cos(1.2*(time - 1))
    return(matrix(means + sqrt(10) * rand, ncol = 1))
  }
  dtransition <- function(next_x, xparticles, theta, time){
    means <- 0.5 * xparticles + 25 * xparticles / (1 + xparticles^2) + 8 * cos(1.2*(time - 1))
    return(dnorm(next_x, mean = means, sd = sqrt(10), log = TRUE))
  }
  dmeasurement <- function(xparticles, theta, observation){
    return(dnorm(observation, mean = (xparticles^2) / 20, sd = 1, log = TRUE))
  }
  generate_randomness <- function(nparticles, datalength){
    init_rand <- rnorm(nparticles, 0, 1)
    transition_rand <- matrix(rnorm(nparticles * datalength, 0, 1), nrow = nparticles, ncol = datalength)
    return(list(init = init_rand, transition = transition_rand))
  }
  gss <- list(rinit = rinit, rtransition = rtransition, dtransition = dtransition, dmeasurement = dmeasurement,
              generate_randomness = generate_randomness, dimension = 1)
  return(gss)
}
