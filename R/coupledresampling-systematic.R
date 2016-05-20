#'@rdname CR_systematic
#'@title Coupled Resampling: systematic
#'@description This function performs systematic resampling with the same uniform on each system.
#'@return Two vectors of ancestors, column-binded in a matrix.
#'@export
CR_systematic <- function(xparticles1, xparticles2, normweights1, normweights2, ...){
  u <- runif(1)
  nparticles <- ncol(xparticles1)
  ancestors1 <- systematic_resampling_n(normweights1, nparticles, u)
  ancestors2 <- systematic_resampling_n(normweights2, nparticles, u)
  return(cbind(ancestors1, ancestors2))
}
