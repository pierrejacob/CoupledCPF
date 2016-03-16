#'@rdname kalman_loglikelihood
#'@title kalman_loglikelihood
#'@description kalman_loglikelihood
#'@export

kalman_loglikelihood <- function(parameters, observations){
  return(kalman_loglikelihood_(parameters, observations))
}
