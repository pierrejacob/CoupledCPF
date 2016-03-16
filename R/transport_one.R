#'@export
transport_one <- function(normweights1, normweights2,
                          resampling_parameters, distance_matrix){
  if (sum(diag(distance_matrix)) < 1e-200){
    u <- runif(1)
    ancestors1 <- systematic_resampling_n(normweights1, 1, u)
    ancestors2 <- systematic_resampling_n(normweights2, 1, u)
    return(c(ancestors1, ancestors2))
  }
  nparticles <- length(normweights1)
  epsilon <- resampling_parameters$epsilon
  desired_alpha <- resampling_parameters$desired_alpha
  wd <- wasserstein_auto(normweights1, normweights2, distance_matrix, epsilon * median(distance_matrix), desired_alpha)
  couplingmatrix <- wd$transportmatrix
  # fix the marginals
  u_1 <- rowSums(couplingmatrix)
  u_2 <- colSums(couplingmatrix)
  #  
  nu <- pmin(normweights1 / u_1, normweights2 / u_2)
  alpha <- min(nu)
  cat("in transport_one: total distance between two sets of paths = ", sum(diag(distance_matrix)), "\n")
  R1 <- (normweights1 - alpha * u_1) / (1 - alpha)
  R2 <- (normweights2 - alpha * u_2) / (1 - alpha)
  
  coupled <- (runif(1) < alpha)
  if (coupled > 0){
    joint_matrix_as_vec <- as.numeric(couplingmatrix)
    ancestors_as_vec <- systematic_resampling_n(normalized_weights = joint_matrix_as_vec, N = 1, u = runif(1))
    ancestors1_ <- 1 + ((ancestors_as_vec-1) %% nparticles)
    ancestors2_ <- (ancestors_as_vec - ancestors1_) / nparticles + 1
    ancestors1 <- ancestors1_
    ancestors2 <- ancestors2_
  } else {
    u <- runif(1)
    ancestors1 <- systematic_resampling_n(R1, 1, u)
    ancestors2 <- systematic_resampling_n(R2, 1, u)
  }
  return(c(ancestors1, ancestors2))
}
