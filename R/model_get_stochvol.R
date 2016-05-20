#'@rdname get_stochvol
#'@title Multivariate stochastic volatility model as in Chopin Gerber 2015
#'@description This function returns a list with objects such as
#'* rinit, rinit_rand to sample from the initial distribution
#'* rtransition, rtransition_rand to sample from the transition
#'* dtransition to evaluate the transition density
#'* dmeasurement to evaluate the measurement density
#'* dimension, which represents the dimension of the latent process
#'@return A list
## Multivariate Stochastic volatility 
#
# Y_t = S_t^{1/2} epsilon_t 
# X_t = mu + Phi * (X_{t-1} - mu) + Psi^{1/2} * eta_t
# X_0 follows Normal(mu, C0)
#
# epsilon_t follows Normal(0, Cy)
# eta_t follows Normal(0, Cx)
# Psi^{1/2} C_x Psi^{1/2}  is denoted V_x such that
# X_t = mu + Phi * (X_{t-1} - mu) + V_x^{1/2} * Normal(0, I)

# theta is a list containing
# mu, C0, Cx, Cy
#'@export
get_stochvol <- function(dimension){
  #
  rinit <- function(nparticles, theta, rand, ...){
    return(rand)
  }
  #
  rinit_rand <- function(nparticles, theta){
    return(fast_rmvnorm_transpose_cholesky(nparticles, theta$mu, theta$cholC0))
  }
  #  
  rtransition <- function(xparticles, theta, time, rand, ...){
    return(theta$mu + theta$Phi %*% (xparticles - theta$mu) + rand)
  }
  #
  rtransition_rand <- function(nparticles, theta){
    return(fast_rmvnorm_transpose_cholesky(nparticles, rep(0, dimension), theta$cholVx))
  }
  #
  dmeasurement <- function(xparticles, theta, observation, ...) {
    # return(fast_dmvnorm_transpose_cholesky(exp(-0.5*xparticles) * observation, rep(0, dimension), theta$cholCy) - 0.5*apply(xparticles, 2, sum))
    return(stochvol_dmeasurement(xparticles, theta, observation, dimension))
  }
  #
  dtransition <- function(next_x, xparticles, theta, time, precomputed, ...){
    return(fast_dmvnorm_transpose_cholesky(theta$mu + theta$Phi %*% (xparticles - theta$mu), next_x, theta$cholVx))
  }
  #
  precompute <- function(...){
    return(list())
  }
  #
  stochvol_model <- list(rinit = rinit, rinit_rand = rinit_rand, rtransition = rtransition,
                           rtransition_rand = rtransition_rand, dtransition = dtransition,
                           dmeasurement = dmeasurement, precompute = precompute, dimension = dimension)
  return(stochvol_model)
}