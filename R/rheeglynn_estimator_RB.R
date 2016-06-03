#'@export
CPF_RB <- function(nparticles, model, theta, observations, ref_trajectory = NULL, with_as = FALSE,
                   test_function = function(trajectory){ return(trajectory) }){
  datalength <- nrow(observations)
  # create tree representation of the trajectories
  Tree <- new(TreeClass, nparticles, 10*nparticles*model$dimension, model$dimension)
  # initialization
  model_precomputed <- model$precompute(theta)
  xparticles <- model$rinit(nparticles, theta, model$rinit_rand(nparticles, theta), model_precomputed)
  if (!is.null(ref_trajectory)){
    xparticles[,nparticles] <- ref_trajectory[,1]
  }
  Tree$init(xparticles)
  #  
  normweights <- rep(1/nparticles, nparticles)
  # if ancestor sampling, needs to keep the last generation of particles at hand
  if (with_as){
    x_last <- xparticles
  }
  # step t > 1
  for (time in 1:datalength){
    ancestors <- multinomial_resampling_n(normweights, nparticles)
    # if no observation or first time, no resampling
    if (time == 1 || (time > 1 && is.na(observations[time-1,1]))){
      ancestors <- 1:nparticles
    }
    xparticles <- xparticles[,ancestors]
    if (is.null(dim(xparticles))) xparticles <- matrix(xparticles, nrow = dimension)
    xparticles <- model$rtransition(xparticles, theta, time, model$rtransition_rand(nparticles, theta), model_precomputed)
    if (!is.null(ref_trajectory)){
      xparticles[,nparticles] <- ref_trajectory[,time+1]
      if (with_as){
        # Ancestor sampling
        logm <- model$dtransition(ref_trajectory[,time+1], x_last, theta, time, model_precomputed)
        logm <- log(normweights) + logm
        w_as <- exp(logm - max(logm))
        w_as <- w_as / sum(w_as)
        ancestors[nparticles] <- systematic_resampling_n(w_as, 1, runif(1))
        x_last <- xparticles
      } else {
        ancestors[nparticles] <- nparticles
      }
    }
    #    
    logw <- model$dmeasurement(xparticles, theta, observations[time,], model_precomputed)
    maxlw <- max(logw)
    w <- exp(logw - maxlw)
    normweights <- w / sum(w)
    #
    Tree$update(xparticles, ancestors - 1)    
  }
  trajectories <- array(dim = c(model$dimension, datalength + 1, nparticles))
  estimate <- 0
  for (k in 0:(nparticles-1)){
    trajectories[,,k+1] <- Tree$get_path(k)
    estimate <- estimate + normweights[k+1] * test_function(trajectories[,,k+1])
  }
  new_trajectory <- matrix(trajectories[,,systematic_resampling_n(normweights, 1, runif(1))], nrow = model$dimension)
  return(list(new_trajectory = new_trajectory, estimate = estimate))
}
#'@export
CPF_RB_coupled <- function(nparticles, model, theta, observations, ref_trajectory1, ref_trajectory2, 
                           coupled_resampling, with_as = FALSE, 
                           test_function = function(trajectory){ return(trajectory) }){
  #
  datalength <- nrow(observations)
  # create tree representation of the trajectories
  Tree1 <- new(TreeClass, nparticles, 10*nparticles*model$dimension, model$dimension)
  Tree2 <- new(TreeClass, nparticles, 10*nparticles*model$dimension, model$dimension)
  # initialization
  model_precomputed <- model$precompute(theta)
  init_rand <- model$rinit_rand(nparticles, theta)
  xparticles1 <- model$rinit(nparticles, theta, init_rand, model_precomputed)
  xparticles1[,nparticles] <- ref_trajectory1[,1]
  Tree1$init(xparticles1)
  normweights1 <- rep(1/nparticles, nparticles)
  #  
  xparticles2 <- model$rinit(nparticles, theta, init_rand, model_precomputed)
  xparticles2[,nparticles] <- ref_trajectory2[,1]
  Tree2$init(xparticles2)
  normweights2 <- rep(1/nparticles, nparticles)
  #
  # if ancestor sampling, needs to keep the last generation of particles at hand
  if (with_as){
    x_last1 <- xparticles1
    x_last2 <- xparticles2
  }
  # step t > 1
  for (time in 1:datalength){
    ancestors <- coupled_resampling(xparticles1, xparticles2, normweights1, normweights2)
    ancestors1 <- ancestors[,1]
    ancestors2 <- ancestors[,2]
    # if no observation or first time, no resampling
    if (time == 1 || (time > 1 && is.na(observations[time-1,1]))){
      ancestors1 <- 1:nparticles
      ancestors2 <- 1:nparticles
    }
    #
    xparticles1 <- xparticles1[,ancestors1]
    xparticles2 <- xparticles2[,ancestors2]
    
    if (is.null(dim(xparticles1))) xparticles1 <- matrix(xparticles1, nrow = dimension)
    if (is.null(dim(xparticles2))) xparticles2 <- matrix(xparticles2, nrow = dimension)
    #
    transition_rand <- model$rtransition_rand(nparticles, theta)
    xparticles1 <- model$rtransition(xparticles1, theta, time, transition_rand, model_precomputed)
    xparticles2 <- model$rtransition(xparticles2, theta, time, transition_rand, model_precomputed)
    if (is.null(dim(xparticles1))) xparticles1 <- matrix(xparticles1, nrow = dimension)
    if (is.null(dim(xparticles2))) xparticles2 <- matrix(xparticles2, nrow = dimension)
    #
    xparticles1[,nparticles] <- ref_trajectory1[,time+1]
    xparticles2[,nparticles] <- ref_trajectory2[,time+1]
    if (with_as){
      # % Ancestor sampling
      logm1 <- model$dtransition(ref_trajectory1[,time+1], x_last1, theta, time, model_precomputed)
      logm1 <- log(normweights1) + logm1
      w_as1 <- exp(logm1 - max(logm1))
      w_as1 <- w_as1 / sum(w_as1)
      unif_resampling_as <- runif(1)
      ancestors1[nparticles] = systematic_resampling_n(w_as1, 1, unif_resampling_as)
      x_last1 <- xparticles1
      #
      logm2 <- model$dtransition(ref_trajectory2[,time+1], x_last2, theta, time, model_precomputed)
      logm2 <- log(normweights2) + logm2
      w_as2 <- exp(logm2 - max(logm2))
      w_as2 <- w_as2 / sum(w_as2)
      ancestors2[nparticles] = systematic_resampling_n(w_as2, 1, unif_resampling_as)
      x_last2 <- xparticles2
    } else {
      ancestors1[nparticles] <- nparticles
      ancestors2[nparticles] <- nparticles
    }
    #    
    logw1 <- model$dmeasurement(xparticles1, theta, observations[time,], model_precomputed)
    logw2 <- model$dmeasurement(xparticles2, theta, observations[time,], model_precomputed)
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
  u <- runif(1)
  k_path1 <- systematic_resampling_n(normweights1, 1, u)
  k_path2 <- systematic_resampling_n(normweights2, 1, u)
  ##
  trajectories1 <- array(dim = c(model$dimension, datalength + 1, nparticles))
  trajectories2 <- array(dim = c(model$dimension, datalength + 1, nparticles))
  estimate1 <- 0
  estimate2 <- 0
  for (k in 0:(nparticles-1)){
    trajectories1[,,k+1] <- Tree1$get_path(k)
    trajectories2[,,k+1] <- Tree2$get_path(k)
    estimate1 <- estimate1 + normweights1[k+1] * test_function(trajectories1[,,k+1])
    estimate2 <- estimate2 + normweights2[k+1] * test_function(trajectories2[,,k+1])
  }
  
  new_trajectory1 <- matrix(trajectories1[,,k_path1], nrow = model$dimension)
  new_trajectory2 <- matrix(trajectories2[,,k_path2], nrow = model$dimension)
  return(list(new_trajectory1 = new_trajectory1, new_trajectory2 = new_trajectory2,
              estimate1 = estimate1, estimate2 = estimate2))
}

#'@export
rheeglynn_estimator_RB <- function(observations, model, theta, algoparameters,
                                   test_function = function(trajectory){ return(trajectory) }){
  # number of particles
  nparticles <- algoparameters$nparticles
  # ancestor sampling
  with_as <- algoparameters$with_as
  # coupled resampling scheme
  coupled_resampling <- algoparameters$coupled_resampling
  #
  CPF_RB_res <- CPF_RB(nparticles, model, theta, observations, test_function = test_function)
  xref <- CPF_RB_res$new_trajectory
  CPF_RB_res_tilde <- CPF_RB(nparticles, model, theta, observations, test_function = test_function)
  xref_tilde <- CPF_RB_res_tilde$new_trajectory
  iteration <- 0
  estimate <- CPF_RB_res$estimate #test_function(xref)
  iteration <- iteration + 1
  CPF_RB_res <- CPF_RB(nparticles, model, theta, observations, xref, with_as = with_as)
  xref <- CPF_RB_res$new_trajectory
  # estimate <- estimate + (test_function(xref) - test_function(xref_tilde))
  estimate <- estimate + CPF_RB_res$estimate - CPF_RB_res_tilde$estimate
  meeting <- FALSE  
  while (!meeting){
    iteration <- iteration + 1
    res <- CPF_RB_coupled(nparticles, model, theta, observations, xref, xref_tilde, coupled_resampling, with_as = with_as,
                          test_function = test_function)
    xref <- res$new_trajectory1
    xref_tilde <- res$new_trajectory2
    estimate <- estimate + (res$estimate1 - res$estimate2)
    # estimate <- estimate + (test_function(xref) - test_function(xref_tilde))
    if (isTRUE(all.equal(xref, xref_tilde))){
      meeting <- TRUE
      break
    }
  }
  return(list(estimate = estimate, iteration = iteration))
}

#'@export
rheeglynn_estimator_RB2 <- function(observations, model, theta, algoparameters,
                                    test_function = function(trajectory){ return(trajectory) }){
  # number of particles
  nparticles <- algoparameters$nparticles
  # ancestor sampling
  with_as <- algoparameters$with_as
  # coupled resampling scheme
  coupled_resampling <- algoparameters$coupled_resampling
  # step at which to start computing the estimate
  m <- algoparameters$m
  if (is.null(algoparameters$m)){
    m <- 0
  }
  #
  CPF_RB_res <- CPF_RB(nparticles, model, theta, observations, test_function = test_function)
  xref <- CPF_RB_res$new_trajectory
  CPF_RB_res_tilde <- CPF_RB(nparticles, model, theta, observations, test_function = test_function)
  xref_tilde <- CPF_RB_res_tilde$new_trajectory
  iteration <- 0
  estimate <- 0
  if (m == 0){
    estimate <- CPF_RB_res$estimate 
  }
  iteration <- 1
  CPF_RB_res <- CPF_RB(nparticles, model, theta, observations, xref, with_as = with_as, test_function)
  xref <- CPF_RB_res$new_trajectory
  if (m == 1){
    estimate <- CPF_RB_res$estimate
  } else {
    estimate <- estimate + CPF_RB_res$estimate - CPF_RB_res_tilde$estimate
  }
  meeting <- FALSE
  meetingtime <- Inf
  finished <- FALSE
  while (!finished){
    # while (!meeting && iteration < m){
    iteration <- iteration + 1
    res <- CPF_RB_coupled(nparticles, model, theta, observations, xref, xref_tilde, coupled_resampling, with_as, test_function)
    if (!meeting && isTRUE(all.equal(xref, xref_tilde))){
      meeting <- TRUE
      meetingtime <- iteration
    }
    xref <- res$new_trajectory1
    xref_tilde <- res$new_trajectory2
    if (m == iteration){
      estimate <- res$estimate1
    } else {
      estimate <- estimate + (res$estimate1 - res$estimate2)
    }
    if (iteration >= meetingtime && iteration >= m){
      finished <- TRUE
      if (m+1 >= meetingtime){
        estimate <- res$estimate1
      }
    }
  }
  return(list(estimate = estimate, iteration = iteration, meetingtime = meetingtime))
}
