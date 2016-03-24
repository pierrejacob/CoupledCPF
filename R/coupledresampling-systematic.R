#'@export
CR_systematic <- function(xparticles1, xparticles2, normweights1, normweights2, ...){
  u <- runif(1)
  nparticles <- nrow(xparticles1)
  ancestors1 <- systematic_resampling_n(normweights1, nparticles, u)
  ancestors2 <- systematic_resampling_n(normweights2, nparticles, u)
  # ## does systematic resampling on joint probability matrix
  # joint_matrix <- matrix(normweights1, ncol = 1) %*% matrix(normweights2, nrow = 1)
  # joint_matrix_as_vec <- as.numeric(joint_matrix) 
  # ancestors_as_vec <- systematic_resampling_n(normalized_weights = joint_matrix_as_vec, N = nparticles, u = runif(1))
  # ancestors1 <- 1 + ((ancestors_as_vec-1) %% nparticles)
  # ancestors2 <- (ancestors_as_vec - ancestors1) / nparticles + 1
  return(cbind(ancestors1, ancestors2))
}
