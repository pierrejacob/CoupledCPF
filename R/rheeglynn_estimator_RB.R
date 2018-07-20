
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
