#'@rdname CR_indexmatching
#'@title Coupled Resampling: index-matching
#'@description This function performs coupled resampling based on index-matching. It takes
#' vectors of normalized weights, of same sizes, and number of desired samples ('ntrials').
#' If 'ntrials' is not given, default to size of normweights1.
#'@return Two vectors of ancestors, column-binded in a matrix. Takes values between 1 and length(normweights1)
#'@export
CR_indexmatching <- function(xparticles1, xparticles2, normweights1, normweights2, ntrials){
  if (missing(ntrials)){
    ntrials <- length(normweights1)
  }
  nu <- pmin(normweights1, normweights2)
  alpha <- sum(nu)
  mu <- nu / alpha
  R1 <- (normweights1-nu)/(1-alpha)
  R2 <- (normweights2-nu)/(1-alpha)  
  uniforms <- runif(ntrials, 0, 1)
  coupled <- (uniforms < alpha)
  ncoupled <- sum(coupled)
  ancestors <- matrix(nrow = ntrials, ncol = 2)
  if (ncoupled > 0){
    ancestors[coupled] <- CoupledCPF:::multinomial_resampling_n_(mu, ncoupled) + 1
  }
  if (ncoupled < ntrials){
    ancestors[!coupled] <- CoupledCPF:::coupled_multinomial_resampling_n_(R1, R2, ntrials - ncoupled) + 1
  }
  return(ancestors)
  # nparticles <- ncol(xparticles1)
  # return(indexmatching_cpp(normweights1, normweights2, nparticles) + 1)
}
