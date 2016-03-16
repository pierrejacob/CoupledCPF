#'@export
indexmatching_one <- function(normweights1, normweights2, ...){
  ##
  nu <- pmin(normweights1, normweights2)
  alpha <- sum(nu)
  # check if the weight vectors are equal, in which case we don't need to sweat too much
  cat("in indexmatching_one: alpha = ", alpha, "\n")
  if (alpha > 1-1e-20){
    ancestors1 <- systematic_resampling_n(normalized_weights = normweights1, N = 1, u = runif(1))
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
  coupled <- (runif(1) < alpha)
  if (coupled){
    k_path1 <- systematic_resampling_n(mu, 1, runif(1))
    k_path2 <- k_path1
  } else {
    u <- runif(1)
    k_path1 <- systematic_resampling_n(R1, 1, u)
    k_path2 <- systematic_resampling_n(R2, 1, u)
  }
  return(c(k_path1, k_path2))
}
