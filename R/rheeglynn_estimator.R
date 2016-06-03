#'@rdname rheeglynn_estimator
#'@title Rhee--Glynn smoothing estimator 
#'@description Estimates the mean of the smoothing distribution at all times.
#'@param observations a matrix of observations of size datalength x dimension(observation)
#'@param model a list representing a model, for instance as given by \code{\link{get_ar}}.
#'@param theta a parameter to give to the model functions
#'@param algoparameters a list containing the algorithmic parameters, that is, nparticles, coupled_resampling,
#' with_as and  lambda. For the first three, see \code{\link{CPF_coupled}}. Lambda refers to the parameter
#' of the Geometric truncation variable. Set to zero to by-pass this choice.
#'@param force_niterations forces the algorithm to run for a set number of iterations (default to zero, which disables this feature)
#'@param max_niterations makes the algorithm stop at the set number of iterations (default to zero, which disables this features)
#'@return a list with the estimator of the smoothing mean, the number of iterations it took, and the truncation variable if any
#'@export
rheeglynn_estimator <- function(observations, model, theta, algoparameters, force_niterations = 0, max_niterations = 0){
  # number of particles
  nparticles <- algoparameters$nparticles
  # coupled resampling scheme
  coupled_resampling <- algoparameters$coupled_resampling
  # probability of the Geometric variable used for the truncation (mean is 1/lambda)
  lambda <- algoparameters$lambda
  # ancestor Sampling (TRUE or FALSE)
  with_as <- algoparameters$with_as
  #
  truncation <- Inf
  survival_probabilities <- 1
  #
  if (lambda > 0){
    truncation <- rgeom(1, lambda)
    survival_probabilities <- (1-lambda)^(1:max(truncation, force_niterations))
      #FALSE: 1 - pgeom(1:max(truncation, force_niterations), prob = lambda)
  }
  xref <- CPF(nparticles, model, theta, observations)
  xref_tilde <- CPF(nparticles, model, theta, observations)
  iteration <- 0
  estimate <- xref
  if (truncation == 0 && force_niterations == 0){
    return(list(estimate = estimate, iteration = iteration, truncation = truncation))
  }
  iteration <- iteration + 1
  xref <- CPF(nparticles, model, theta, observations, xref, with_as)
  estimate <- estimate + (xref - xref_tilde) / survival_probabilities[1]
  if (truncation == 1 && force_niterations == 0){
    return(list(estimate = estimate, iteration = iteration, truncation = truncation))
  }
  
  while ((force_niterations == 0 && iteration < truncation) || (force_niterations != 0 && iteration < force_niterations)){
    iteration <- iteration + 1
    res <- CPF_coupled(nparticles, model, theta, observations, xref, xref_tilde, coupled_resampling, with_as)
    xref <- res$new_trajectory1
    xref_tilde <- res$new_trajectory2
    if (lambda > 0){
      estimate <- estimate + (xref - xref_tilde) / survival_probabilities[iteration]
    } else {
      estimate <- estimate + (xref - xref_tilde)
    }
    if (isTRUE(all.equal(xref, xref_tilde)) && force_niterations == 0){
      break
    }
    if (max_niterations > 0 && iteration >= max_niterations){
      break
    }
  }
  return(list(estimate = estimate, iteration = iteration, truncation = truncation))
}


