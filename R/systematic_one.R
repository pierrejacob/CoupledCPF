#'@export
systematic_one <- function(normweights1, normweights2, ...){
  u <- runif(1)
  k_path1 <- systematic_resampling_n(normweights1, 1, u)
  k_path2 <- systematic_resampling_n(normweights2, 1, u)
  return(c(k_path1, k_path2))
}
