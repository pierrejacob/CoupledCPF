#'@export
levydriven_rtransition_rand <- function(nparticles, theta){
  return(levydriven_rtransition_rand_cpp(nparticles, theta))
}
#' 
#' #'@export
#' levydriven_rtransition <- function(nparticles, theta){
#'   return(levydriven_rtransition_rand_cpp(nparticles, theta))
#' }
