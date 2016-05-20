#'@rdname get_levydriven
#'@title Levy-driven stochastic volatility model as in Chopin Jacob Papaspiliopoulos 2013
#'@description This function returns a list with objects such as
#'* rinit, rinit_rand to sample from the initial distribution
#'* rtransition, rtransition_rand to sample from the transition
#'* dtransition to evaluate the transition density
#'* dmeasurement to evaluate the measurement density
#'* dimension, which represents the dimension of the latent process
#'@return A list
#'@export
## Stochastic volatility : one-factor model
#
# Y_t = mu + beta * v_t + v_t^(0.5) * epsilon_t
# X_t = (v_t, z_t)
# v_t+1 = lambda^(-1) ( z_t - z_{t+1} + sum_j=1^k e_j )
# z_t+1 = e^(-lambda) * z_t + sum_j=1^k exp(-lambda(t + 1 - c_j)) e_j
# k ~ Poisson(lambda * xi^2 / omega^2)
# c_{1:k} ~ Uniform(t, t + 1)
# e_{1:k} ~ Exp(xi / omega^2) (rate parameter)
#
# v_0 does not matter
# z_0 ~ Gamma(xi^2 / omega^2, xi / omega^2) (rate parameter)
#
# parameters[0, :] = mu
# parameters[1, :] = beta
# parameters[2, :] = xi
# parameters[3, :] = omega^2
# parameters[4, :] = lambda
# 0 0 0.5 0.0625 0.01
#
get_levydriven <- function(){
  # theta <- c(0, 0, 0.5, 0.0625, 0.01)
  # mu <- theta[1]
  # beta <- theta[2]
  # xi <- theta[3]
  # omega2 <- theta[4]
  # lambda <- theta[5]
  #
  rinit <- function(nparticles, theta, rand, ...){
    xparticles <- matrix(0, ncol = nparticles, nrow = 2)
    xparticles[2,] <- rand
    return(xparticles)
  }
  #
  rinit_rand <- function(nparticles, theta){
    return(rgamma(nparticles, shape = theta[3] * theta[3] / theta[4], scale = theta[4]/theta[3]))
  }
  #  
  rtransition <- function(xparticles, theta, time, rand, ...){
    new_z <- exp(-theta[5]) * xparticles[2,] + rand$sum_weighted_e 
    new_v <- (1/theta[5]) * (xparticles[2,] - new_z + rand$sum_e)
    xparticles[1,] <- new_v
    xparticles[2,] <- new_z
    return(xparticles)
  }
  #
  rtransition_rand <- function(nparticles, theta){
    return(levydriven_rtransition_rand_cpp(nparticles, theta))
  }
  #
  dmeasurement <- function(xparticles, theta, observation, ...) {
    return(dnorm(observation, mean = theta[1] + theta[2] * xparticles[1,], sd = sqrt(xparticles[1,]), log = TRUE))
  }
  #
  precompute <- function(...){
    return(list())
  }
  #
  levydriven_model <- list(rinit = rinit, rinit_rand = rinit_rand, rtransition = rtransition,
                           rtransition_rand = rtransition_rand,
                           dmeasurement = dmeasurement, precompute = precompute, dimension = 2)
  return(levydriven_model)
}