# algorithm to construct the decomposition
# P = alpha * mu + (1 - alpha) (R * tilde{R}^T)
# such that the marginals are correct (normweights1 and normweights2)
# and there's more chance to sample a "matched" pair of indices
# (than if we were doing sampling from the independent coupling normweights1 * normweights2^T)
#'@export
index_then_transport_resampling <- function(xparticles1, xparticles2, normweights1, normweights2, 
                                            parameters = list(epsilon = 0.1, niterations = 100, MMax = 1e10)){
  nparticles <- nrow(xparticles1)
  xdim <- ncol(xparticles1)
  
  epsilon <- parameters$epsilon
  niterations <- parameters$niterations
  
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
  
  R1_where_zeros <- (R1 < 1e-10)  
  R2_where_zeros <- (R2 < 1e-10)  
  R1_nonzeros <- R1[!R1_where_zeros]
  R2_nonzeros <- R2[!R2_where_zeros]
  M_reduced <- cost_matrix(matrix(xparticles1[!R1_where_zeros,], ncol = xdim),
                           matrix(xparticles2[!R2_where_zeros,], ncol = xdim))
  # super hacky
  M_reduced[M_reduced > parameters$MMax] <- parameters$MMax
  couplingmatrix <- matrix(0, nrow = nparticles, ncol = nparticles)
  
  # flip coins
  coupled <- (runif(nparticles) < alpha)
  ncoupled <- sum(coupled)
  ancestors1 <- rep(0, nparticles)
  ancestors2 <- rep(0, nparticles)
  if (ncoupled > 0){
    ancestors1[coupled] <- systematic_resampling_n(mu, ncoupled, runif(1))
    ancestors2[coupled] <- ancestors1[coupled]
    if (ncoupled < nparticles){
      wd <- wasserstein(R1_nonzeros, R2_nonzeros, M_reduced, epsilon * median(M_reduced), niterations)
      
      couplingmatrix[!R1_where_zeros, !R2_where_zeros] <- wd$transportmatrix
      joint_as_vec <- as.numeric(couplingmatrix) 
      ancestors_as_vec <- systematic_resampling_n(normalized_weights = joint_as_vec, N = nparticles - ncoupled,
                                                  u = runif(1))
      a1 <- 1 + ((ancestors_as_vec-1) %% nparticles)
      a2 <- (ancestors_as_vec - a1) / nparticles + 1
      ancestors1[!coupled] <- a1
      ancestors2[!coupled] <- a2
    }
  } else {
    wd <- wasserstein(p = matrix(R1_nonzeros, ncol = 1), qs = matrix(R2_nonzeros, ncol = 1),
                      M_reduced, epsilon * median(M_reduced), niterations)
    couplingmatrix[!R1_where_zeros, !R2_where_zeros] <- wd$transportmatrix
    joint_as_vec <- as.numeric(couplingmatrix) 
    ancestors_as_vec <- systematic_resampling_n(normalized_weights = joint_as_vec, N = nparticles, u = runif(1))
    ancestors1 <- 1 + ((ancestors_as_vec-1) %% nparticles)
    ancestors2 <- (ancestors_as_vec - ancestors1) / nparticles + 1
  }
  return(cbind(ancestors1, ancestors2))
}
