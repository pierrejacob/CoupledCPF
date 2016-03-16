# matching: each row is associated with one column
# (but then one column can be associated with multiple rows)

# algorithm to construct the decomposition
# P = alpha * nu + (1 - alpha) R_x * R_y^T 
# such that the marginals are correct (normweights1 and normweights2)
# and there's more chance to sample a "matched" pair of indices
# (than if we were doing sampling from the independent coupling normweights1 * normweights2^T)
#'@export
indexmatching_resampling <- function(xparticles1, xparticles2, normweights1, normweights2, ...){
  nparticles <- nrow(xparticles1)
  # common measure  
  nu <- pmin(normweights1, normweights2)
  alpha <- sum(nu)
  # check if the weight vectors are equal, in which case we don't need to sweat too much
  if (alpha > 1-1e-20){
    ancestors1 <- systematic_resampling_n(normalized_weights = normweights1, N = nparticles, u = runif(1))
    ancestors2 <- ancestors1
    return(cbind(ancestors1, ancestors2))
  }
  mu <- nu / alpha 
  # residuals
  R1 <- normweights1 - nu
  R1 <- R1 / (1 - alpha)
  R2 <- normweights2 -  nu
  R2 <- R2 / (1 - alpha)
  #
  coupled <- (runif(nparticles) < alpha)
  ncoupled <- sum(coupled)
  ancestors1 <- rep(0, nparticles)
  ancestors2 <- rep(0, nparticles)
  if (ncoupled > 0){
    ancestors1[coupled] <- systematic_resampling_n(mu, ncoupled, runif(1))
    ancestors2[coupled] <- ancestors1[coupled]
    if (ncoupled < nparticles){
      u <- runif(1)
      ancestors1[!coupled] <- systematic_resampling_n(R1, nparticles - ncoupled, u)
      ancestors2[!coupled] <- systematic_resampling_n(R2, nparticles - ncoupled, u)
    }
  } else {
    u <- runif(1)
    ancestors1 <- systematic_resampling_n(R1, nparticles, u)
    ancestors2 <- systematic_resampling_n(R2, nparticles, u)
  }
  return(cbind(ancestors1, ancestors2))
}
  

  

  
  
#   nparticles <- nrow(xparticles1)
#   matching_vector <- 1:nparticles
#   nsamples <- nparticles
#   M <- matrix(0, ncol = nparticles, nrow = nparticles)
#   for (i in 1:nparticles){
#     M[i,matching_vector[i]] <- 1
#   }
#   
#   a <- normweights1
#   for (i in 1:nparticles){
#     match_i <- which(M[,i] == 1)
#     if (length(match_i) > 0){
#       b <- sum(a[match_i])
#       c_i <- max(1, b / normweights2[i])
#       for (j in match_i){
#         a[j] <- a[j] / c_i
#       }
#     }
#   }
#   
#   alpha <- sum(a)
#   omega_M <- a / alpha
#   R1 <- normweights1 - a
#   R1[R1 < 0] <- 0
#   R1 <- R1 / sum(R1)
#   R2 <- t(t(normweights2) - t(a) %*% M)
#   R2[R2 < 0] <- 0
#   R2 <- R2 / sum(R2)
#   coupled <- (runif(nsamples) < alpha)
#   ncoupled <- sum(coupled)
#   indx <- rep(0, nsamples)
#   indy <- rep(0, nsamples)
#   if (ncoupled > 0){
#     indx[coupled] <- systematic_resampling_n(omega_M, ncoupled, runif(1))
#     indy[coupled] <- matching_vector[indx[coupled]]
#     if (ncoupled < nsamples){
#       indx[!coupled] <- systematic_resampling_n(R1, nsamples - ncoupled, runif(1))
#       indy[!coupled] <- systematic_resampling_n(R2, nsamples - ncoupled, runif(1))
#       
#     }
#   } else {
#     indx <- systematic_resampling_n(R1, nsamples, runif(1))
#     indy <- systematic_resampling_n(R2, nsamples, runif(1))
#     
#   }
#   return(cbind(indx, indy))
# }
