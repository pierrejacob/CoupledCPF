#'@export
CR_systematic_given <- function(xparticles1, xparticles2, normweights1, normweights2, 
                                              parameters = list(),
                                              ancestors_ref, uniforms){
  # u <- runif(1)
  nparticles <- nrow(xparticles1)
  ancestors2 <- systematic_resampling_n(normweights2, nparticles, uniforms[1])
  return(ancestors2)
}
