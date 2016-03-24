#'@export
CR_transport <- function(xparticles1, xparticles2, normweights1, normweights2, 
                                 parameters = list(epsilon = 0.1, desired_alpha = 0.9)){
  nparticles <- nrow(xparticles1)
  M <- cost_matrix(xparticles1, xparticles2)
  epsilon <- parameters$epsilon
  wd <- wasserstein_auto(p = normweights1, q = normweights2, M, epsilon * median(M), parameters$desired_alpha)
  couplingmatrix <- wd$transportmatrix
  # fix the marginals
  u_1 <- rowSums(couplingmatrix)
  u_2 <- colSums(couplingmatrix)
  #  
  nu <- pmin(normweights1 / u_1, normweights2 / u_2)
  alpha <- min(nu)
  R1 <- (normweights1 - alpha * u_1) / (1 - alpha)
  R2 <- (normweights2 - alpha * u_2) / (1 - alpha)
  coupled <- (runif(nparticles) < alpha)
  ncoupled <- sum(coupled)
  
  ancestors1 <- rep(0, nparticles)
  ancestors2 <- rep(0, nparticles)
  if (ncoupled > 0){
    joint_matrix_as_vec <- as.numeric(couplingmatrix)
    ancestors_as_vec <- systematic_resampling_n(normalized_weights = joint_matrix_as_vec, N = ncoupled, u = runif(1))
    ancestors1_ <- 1 + ((ancestors_as_vec-1) %% nparticles)
    ancestors2_ <- (ancestors_as_vec - ancestors1_) / nparticles + 1
    ancestors1[coupled] <- ancestors1_
    ancestors2[coupled] <- ancestors2_
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

#' #'@export
#' transport_resampling <- function(xparticles1, xparticles2, normweights1, normweights2, 
#'                                        parameters = list(epsilon = 0.1, niterations = 100)){
#'   nparticles <- nrow(xparticles1)
#'   M <- cost_matrix(xparticles1, xparticles2)
#'   epsilon <- parameters$epsilon
#'   niterations <- parameters$niterations
#'   wd <- wasserstein(p = normweights1, q = normweights2, M, epsilon * median(M), niterations)
#'   couplingmatrix <- wd$transportmatrix
#'   # fix the marginals
#'   u_1 <- rowSums(couplingmatrix)
#'   u_2 <- colSums(couplingmatrix)
#'   #  
#'   nu <- pmin(normweights1 / u_1, normweights2 / u_2)
#'   alpha <- min(nu)
#'   R1 <- (normweights1 - alpha * u_1) / (1 - alpha)
#'   R2 <- (normweights2 - alpha * u_2) / (1 - alpha)
#'   coupled <- (runif(nparticles) < alpha)
#'   ncoupled <- sum(coupled)
#'   
#'   ancestors1 <- rep(0, nparticles)
#'   ancestors2 <- rep(0, nparticles)
#'   if (ncoupled > 0){
#'     joint_matrix_as_vec <- as.numeric(couplingmatrix)
#'     ancestors_as_vec <- systematic_resampling_n(normalized_weights = joint_matrix_as_vec, N = ncoupled, u = runif(1))
#'     ancestors1_ <- 1 + ((ancestors_as_vec-1) %% nparticles)
#'     ancestors2_ <- (ancestors_as_vec - ancestors1_) / nparticles + 1
#'     ancestors1[coupled] <- ancestors1_
#'     ancestors2[coupled] <- ancestors2_
#'     if (ncoupled < nparticles){
#'       u <- runif(1)
#'       ancestors1[!coupled] <- systematic_resampling_n(R1, nparticles - ncoupled, u)
#'       ancestors2[!coupled] <- systematic_resampling_n(R2, nparticles - ncoupled, u)
#'     }
#'   } else {
#'     u <- runif(1)
#'     ancestors1 <- systematic_resampling_n(R1, nparticles, u)
#'     ancestors2 <- systematic_resampling_n(R2, nparticles, u)
#'   }
#'   return(cbind(ancestors1, ancestors2))
#' }