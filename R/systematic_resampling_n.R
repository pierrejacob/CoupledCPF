#' @export

systematic_resampling_n <- function(normalized_weights, N, u){
  return(systematic_resampling_n_(normalized_weights, N, u) + 1)
}