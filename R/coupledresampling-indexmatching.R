# index matching resampling function in C++
# the uniform vector must be of size (ndraws + 2)
#'@export
CR_indexmatching <- function(xparticles1, xparticles2, normweights1, normweights2, ...){
  nparticles <- nrow(xparticles1)
  uniforms <- runif(nparticles + 2)
  return(indexmatching_cpp(nparticles, nw1, nw2, uniforms) + 1)
}