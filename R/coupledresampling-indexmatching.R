#'@rdname CR_indexmatching
#'@title Coupled Resampling: index-matching
#'@description This function performs coupled resampling based on index-matching.
#'@return Two vectors of ancestors, column-binded in a matrix.
#'@export
CR_indexmatching <- function(xparticles1, xparticles2, normweights1, normweights2, ...){
  nparticles <- ncol(xparticles1)
  uniforms <- runif(nparticles + 2)
  return(indexmatching_cpp(nparticles, normweights1, normweights2, uniforms) + 1)
}