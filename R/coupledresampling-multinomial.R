#'@rdname CR_multinomial
#'@title Coupled Resampling: multinomial
#'@description This function performs independent multinomial resampling on each system
#'@return Two vectors of ancestors, column-binded in a matrix.
#'@export
CR_multinomial <- function(xparticles1, xparticles2, normweights1, normweights2, ...){
  nparticles <- ncol(xparticles1)
  ancestors <- coupled_multinomial_resampling_n_(normweights1, normweights2, nparticles)
  return(ancestors + 1)
  # ancestors1 <- sample(x = 1:nparticles, size = nparticles, replace = TRUE, prob = normweights1)
  # ancestors2 <- sample(x = 1:nparticles, size = nparticles, replace = TRUE, prob = normweights2)
  # return(cbind(ancestors1, ancestors2))
}