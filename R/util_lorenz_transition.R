#'@rdname lorenz_transition
#'@title lorenz_transition
#'@description Solve lorenz ODE for each particle, given each alpha, from time to time + 1,
#' and given the parameters (c, e, ml, mq).
#'@export
#'
lorenz_transition <- function(xparticles, time_start, time_end, dt, parameters){
  return(one_step_lorenz_vector(xparticles, time_start, time_end, dt, parameters))
}
