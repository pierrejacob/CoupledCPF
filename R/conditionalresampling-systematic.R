#'@export
ind_systematic_resampling_given_ref <- function(xparticles1, xparticles2, normweights1, normweights2, 
                                              parameters = list(),
                                              ancestors_ref){
  # u <- runif(1)
  nparticles <- nrow(xparticles1)
  ancestors2 <- systematic_resampling_n(normweights2, nparticles, runif(1))
  return(ancestors2)
}
