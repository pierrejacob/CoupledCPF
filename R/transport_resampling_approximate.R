#' #'@export
#' transport_resampling_approximate <- function(xparticles1, xparticles2, normweights1, normweights2, 
#'                                  parameters = list(epsilon = 0.1, niterations = 100, MMax = 1e10)){
#'   nparticles <- nrow(xparticles1)
#'   M <- cost_matrix(xparticles1, xparticles2)
#'   epsilon <- parameters$epsilon
#'   niterations <- parameters$niterations
#'   # super hacky
#'   M[M > parameters$MMax] <- parameters$MMax
#'   # print(M)
#'   wd <- wasserstein(p = matrix(normweights1, ncol = 1), qs = matrix(normweights2, ncol = 1),
#'                     M, epsilon * median(M), niterations)
#'   couplingmatrix <- wd$transportmatrix
#'   joint_matrix_as_vec <- as.numeric(couplingmatrix)
#'   ancestors_as_vec <- systematic_resampling_n(normalized_weights = joint_matrix_as_vec, N = nparticles, u = runif(1))
#'   ancestors1 <- 1 + ((ancestors_as_vec-1) %% nparticles)
#'   ancestors2 <- (ancestors_as_vec - ancestors1) / nparticles + 1
#'   return(cbind(ancestors1, ancestors2))
#' }
