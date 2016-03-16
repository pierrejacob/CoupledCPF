#'@export
#'
wasserstein <- function(p, q, cost_matrix, epsilon, niterations){
  return(wasserstein_(p, q, cost_matrix, epsilon, niterations))
}
