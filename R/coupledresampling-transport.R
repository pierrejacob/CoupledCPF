#'@export
CR_transport <- function(xparticles1, xparticles2, normweights1, normweights2, 
                         parameters = list(epsilon = 0.1, desired_alpha = 0.9)){
  nparticles <- nrow(xparticles1)
  ancestors <- transport_cpp(xparticles1, xparticles2, normweights1, normweights2, runif(nparticles + 2), parameters$epsilon, parameters$desired_alpha)
  return(ancestors + 1)
}