#'@export
CR_indexmatching_given <- function(xparticles1, xparticles2, normweights1, normweights2, 
                                   parameters = list(),
                                   ancestors_ref, uniforms){
  nparticles <- nrow(xparticles1)
  # uniforms <- runif(nparticles + 1)
  return(indexmatching_given_cpp(nparticles, normweights1, normweights2, uniforms, ancestors_ref - 1) + 1)
}

# CR_indexmatching_given <- function(xparticles1, xparticles2, normweights1, normweights2, 
#                                              parameters = list(),
#                                              ancestors_ref){
#   nparticles <- nrow(xparticles1)
#   # common measure  
#   nu <- pmin(normweights1, normweights2)
#   alpha <- sum(nu)
#   # check if the weight vectors are equal, in which case we don't need to sweat too much
#   if (alpha > 1-1e-20){
#     ancestors2 <- ancestors_ref
#     return(ancestors2)
#   }
#   mu <- nu / alpha 
#   # residuals
#   R2 <- normweights2 -  nu
#   R2 <- R2 / (1 - alpha)
#   #
#   coupled <- (runif(nparticles) < alpha)
#   ncoupled <- sum(coupled)
#   ancestors2 <- rep(0, nparticles)
#   if (ncoupled > 0){
#     ancestors2[coupled] <- ancestors_ref[coupled]
#     if (ncoupled < nparticles){
#       u <- runif(1)
#       ancestors2[!coupled] <- systematic_resampling_n(R2, nparticles - ncoupled, u)
#     }
#   } else {
#     u <- runif(1)
#     ancestors2 <- systematic_resampling_n(R2, nparticles, u)
#   }
#   return(ancestors2)
# }
