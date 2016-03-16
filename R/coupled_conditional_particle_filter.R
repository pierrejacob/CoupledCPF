#'@rdname coupled_conditional_particle_filter
#'@title coupled_conditional_particle_filter
#'@description runs two coupled conditional particle filters, with or without ancestor sampling
#'@export

coupled_conditional_particle_filter <- function(nparticles, model, theta, observations, randomness, 
                                                ref_trajectory1, ref_trajectory2, 
                                                coupled_resampling, final_resampling, resampling_parameters,
                                                with_as, TreeClass){
  given_args <- as.list(match.call())
  distance_matrix <- NULL
  compute_distance <- FALSE
  if (identical(as.character(given_args$final_resampling), "transport_one")){
    compute_distance <- TRUE
    distance_matrix <- matrix(0, nparticles, nparticles)
  }
  datalength <- nrow(observations)
  # create tree representation of the trajectories
  Tree1 <- new(TreeClass, nparticles, 10*nparticles, model$dimension)
  Tree2 <- new(TreeClass, nparticles, 10*nparticles, model$dimension)
  # initialization
  xparticles1 <- model$rinit(nparticles, theta, randomness$init)
  xparticles1[nparticles,] <- ref_trajectory1[1,]
  Tree1$init(xparticles1)
  normweights1 <- rep(1/nparticles, nparticles)
  #  
  xparticles2 <- model$rinit(nparticles, theta, randomness$init)
  xparticles2[nparticles,] <- ref_trajectory2[1,]
  Tree2$init(xparticles2)
  normweights2 <- rep(1/nparticles, nparticles)
  #
  if (compute_distance){
    distance_matrix <- distance_matrix + square_cost_matrix(xparticles1, xparticles2)
  }
  #
  # if ancestor sampling, needs to keep the last generation of particles at hand
  if (with_as){
    x_last1 <- xparticles1
    x_last2 <- xparticles2
  }
  # step t > 1
  for (time in 1:datalength){
    ancestors <- coupled_resampling(xparticles1, xparticles2, normweights1, normweights2, resampling_parameters)
    ancestors1 <- ancestors[,1]
    ancestors2 <- ancestors[,2]
    #
    if (compute_distance){
      distance_matrix <- distance_matrix[ancestors1,]
      distance_matrix <- distance_matrix[,ancestors2]
    }
    
    #
    xparticles1 <- xparticles1[ancestors1,]
    xparticles2 <- xparticles2[ancestors2,]
    #
    xparticles1 <- model$rtransition(xparticles1, theta, time, randomness$transition[,time])
    xparticles2 <- model$rtransition(xparticles2, theta, time, randomness$transition[,time])
    #
    if (compute_distance){
      distance_matrix <- distance_matrix + square_cost_matrix(xparticles1, xparticles2) 
    }
    #
    xparticles1[nparticles,] <- ref_trajectory1[time+1,]
    xparticles2[nparticles,] <- ref_trajectory2[time+1,]
    if (with_as){
      # % Ancestor sampling
      logm1 <- model$dtransition(ref_trajectory1[time+1,], x_last1, theta, time)
      logm1 <- log(normweights1) + logm1
      w_as1 <- exp(logm1 - max(logm1))
      w_as1 <- w_as1 / sum(w_as1)
      ancestors1[nparticles] = systematic_resampling_n(w_as1, 1, randomness$resampling[time,2])
      x_last1 <- xparticles1
      #
      logm2 <- model$dtransition(ref_trajectory2[time+1,], x_last2, theta, time)
      logm2 <- log(normweights2) + logm2
      w_as2 <- exp(logm2 - max(logm2))
      w_as2 <- w_as2 / sum(w_as2)
      ancestors2[nparticles] = systematic_resampling_n(w_as2, 1, randomness$resampling[time,2])
      x_last2 <- xparticles2
    } else {
      ancestors1[nparticles] <- nparticles
      ancestors2[nparticles] <- nparticles
    }
    #    
    logw1 <- model$dmeasurement(xparticles1, theta, observations[time,])
    logw2 <- model$dmeasurement(xparticles2, theta, observations[time,])
    #
    maxlw1 <- max(logw1)
    w1 <- exp(logw1 - maxlw1)
    normweights1 <- w1 / sum(w1)
    #
    maxlw2 <- max(logw2)
    w2 <- exp(logw2 - maxlw2)
    normweights2 <- w2 / sum(w2)    
    #
    Tree1$update(xparticles1, ancestors1 - 1)    
    Tree2$update(xparticles2, ancestors2 - 1)    
  }
  if (compute_distance){
    distance_matrix <- sqrt(distance_matrix)
  }
  path_indices <- final_resampling(normweights1, normweights2, resampling_parameters, distance_matrix)
  
  k_path1 <- path_indices[1]
  k_path2 <- path_indices[2]
  ##
  new_trajectory1 <- Tree1$get_path(k_path1 - 1)
  new_trajectory2 <- Tree2$get_path(k_path2 - 1)
  return(list(new_trajectory1 = new_trajectory1, new_trajectory2 = new_trajectory2, distance_matrix = distance_matrix))
}


# coupled_conditional_particle_filter <- function(nparticles, model, theta, observations, randomness, 
#                                                 TreeClass, ref_trajectory1, ref_trajectory2, 
#                                                 coupled_resampling, resampling_parameters = list(),
#                                                 with_as = FALSE){
#   datalength <- nrow(observations)
#   # create tree representation of the trajectories
#   Tree1 <- new(TreeClass, nparticles, 10*nparticles, model$dimension)
#   Tree2 <- new(TreeClass, nparticles, 10*nparticles, model$dimension)
#   # initialization
#   xparticles1 <- model$rinit(nparticles, theta, randomness$init)
#   xparticles1[nparticles,] <- ref_trajectory1[1,]
#   Tree1$init(xparticles1)
#   normweights1 <- rep(1/nparticles, nparticles)
#   #  
#   xparticles2 <- model$rinit(nparticles, theta, randomness$init)
#   xparticles2[nparticles,] <- ref_trajectory2[1,]
#   Tree2$init(xparticles2)
#   normweights2 <- rep(1/nparticles, nparticles)
#   
#   # if ancestor sampling, needs to keep the last generation of particles at hand
#   if (with_as){
#     x_last1 <- xparticles1
#     x_last2 <- xparticles2
#   }
#   # step t > 1
#   for (time in 1:datalength){
#     ancestors <- coupled_resampling(xparticles1, xparticles2, normweights1, normweights2, resampling_parameters)
#     ancestors1 <- ancestors[,1]
#     ancestors2 <- ancestors[,2]
#     #
#     xparticles1 <- xparticles1[ancestors1,]
#     xparticles2 <- xparticles2[ancestors2,]
#     #
#     xparticles1 <- model$rtransition(xparticles1, theta, time, randomness$transition[,time])
#     xparticles2 <- model$rtransition(xparticles2, theta, time, randomness$transition[,time])
#     #
#     xparticles1[nparticles,] <- ref_trajectory1[time+1,]
#     xparticles2[nparticles,] <- ref_trajectory2[time+1,]
#     if (with_as){
#       # % Ancestor sampling
#       logm1 <- model$dtransition(ref_trajectory1[time+1,], x_last1, theta, time)
#       logm1 <- log(normweights1) + logm1
#       w_as1 <- exp(logm1 - max(logm1))
#       w_as1 <- w_as1 / sum(w_as1)
#       ancestors1[nparticles] = systematic_resampling_n(w_as1, 1, randomness$resampling[time,2])
#       x_last1 <- xparticles1
#       #
#       logm2 <- model$dtransition(ref_trajectory2[time+1,], x_last2, theta, time)
#       logm2 <- log(normweights2) + logm2
#       w_as2 <- exp(logm2 - max(logm2))
#       w_as2 <- w_as2 / sum(w_as2)
#       ancestors2[nparticles] = systematic_resampling_n(w_as2, 1, randomness$resampling[time,2])
#       x_last2 <- xparticles2
#     } else {
#       ancestors1[nparticles] <- nparticles
#       ancestors2[nparticles] <- nparticles
#     }
#     #    
#     logw1 <- model$dmeasurement(xparticles1, theta, observations[time,])
#     logw2 <- model$dmeasurement(xparticles2, theta, observations[time,])
#     #
#     maxlw1 <- max(logw1)
#     w1 <- exp(logw1 - maxlw1)
#     normweights1 <- w1 / sum(w1)
#     #
#     maxlw2 <- max(logw2)
#     w2 <- exp(logw2 - maxlw2)
#     normweights2 <- w2 / sum(w2)    
#     #
#     Tree1$update(xparticles1, ancestors1 - 1)    
#     Tree2$update(xparticles2, ancestors2 - 1)    
#   }
#   k_path1 <- systematic_resampling_n(normweights1, 1, randomness$resampling[datalength+1,1])
#   k_path2 <- systematic_resampling_n(normweights2, 1, randomness$resampling[datalength+1,1])
#   new_trajectory1 <- Tree1$get_path(k_path1 - 1)
#   new_trajectory2 <- Tree2$get_path(k_path2 - 1)
#   return(list(new_trajectory1 = new_trajectory1, new_trajectory2 = new_trajectory2))
# }
