#'@rdname kalman_filtering_means
#'@title kalman_filtering_means
#'@description kalman_filtering_means
#'@export

kalman_filtering_means <- function(parameters, observations){
  return(kalman_filtering_means_(parameters, observations))
}
