#' #'@rdname get_pzwithar
#' #'@title Phytoplankton-zooplankton model as in Jones, Parslow, Murray 2010, with a twist
#' #'@description This function returns a list with objects such as
#' #'* rinit, rinit_rand to sample from the initial distribution
#' #'* rtransition, rtransition_rand to sample from the transition
#' #'* dtransition to evaluate the transition density
#' #'* dmeasurement to evaluate the measurement density
#' #'* dimension, which represents the dimension of the latent process
#' #'@return A list
#' #'@export
#' # PZ model with auto-regressive process for alpha
#' # theta = (mu_alpha, sd_alpha, c, e, ml, mq, rho_alpha)
#' # transformed so that all parameters are in R (except rho_alpha at this point)
#' # theta = (mu_alpha, log sd_alpha, logit c, logit e, logit ml, logit mq, rho_alpha)
#' get_pzwithar <- function(){
#'   # logit <- function(z) log(z / (1 - z))
#'   expit <- function(z) 1 / (1 + exp(-z))
#'   
#'   rinit <- function(nparticles, theta, rand, ...){
#'     randx <- rand[1:(2*nparticles)]
#'     randalpha <- rand[(2*nparticles+1):(3*nparticles)]
#'     x <- matrix(nrow = 3, ncol = nparticles)
#'     x[1:2,] <- exp(matrix(log(2) + randx, nrow = 2))
#'     x[3,] <- theta[1] + exp(theta[2]) * randalpha
#'     return()
#'   }
#'   rinit_rand <- function(nparticles, theta){
#'     return(rnorm(3*nparticles))
#'   }
#'   #  
#'   rtransition <- function(xparticles, theta, time, rand, ...){
#'     alphas <- theta[1] + exp(theta[2]) * rand
#'     alphas <- xparticles[1:2,]
#'     xparticles[1:2,] <- pz_transition(xparticles[1:2,], alphas, time-1, expit(theta[3:6]))
#'     return(xparticles)
#'   }
#'   #
#'   rtransition_rand <- function(nparticles, theta){
#'     return(rnorm(nparticles))
#'   }
#'   #
#'   dmeasurement <- function(xparticles, theta, observation, ...) {
#'     return(dnorm(x = observation, mean = log(xparticles[1,]), sd = 0.2, log = TRUE))
#'   }
#'   #
#'   precompute <- function(...){
#'     return(list())
#'   }
#'   #
#'   pz_model <- list(rinit = rinit, rinit_rand = rinit_rand, rtransition = rtransition,
#'                    rtransition_rand = rtransition_rand,
#'                    dmeasurement = dmeasurement, precompute = precompute, dimension = 2)
#'   return(pz_model)
#' }
