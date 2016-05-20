#'@export
stochvol_dmeasurement <- function(xparticles, theta, observation, dim){
  return(stochvol_dmeas_(xparticles, theta, observation, dim))
}
