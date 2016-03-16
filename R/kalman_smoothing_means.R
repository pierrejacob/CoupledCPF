#'@rdname kalman_smoothing_means
#'@title kalman_smoothing_means
#'@description kalman_smoothing_means
#'@export

kalman_smoothing_means <- function(parameters, observations){
  return(kalman_smoothing_means_(parameters, observations))
}
