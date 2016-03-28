#'@rdname coupled_pf_given
#'@title coupled_pf_given
#'@description runs a coupled particle filter, given the results of another one
#'@export
coupled_pf_given <- function(nparticles, model, theta2, observations, randomness, 
                                    coupled_resampling_given, resampling_parameters = list(),
                                    system_ref){
  datalength <- nrow(observations)
  # initialization
  xparticles1 <- system_ref$xhistory[[1]]
  normweights1 <- system_ref$whistory[[1]]
  #
  xparticles2 <- model$rinit(nparticles, theta2, randomness$init)
  normweights2 <- rep(1/nparticles, nparticles)
  ll2 <- 0
  #
  xhistory <- rep(list(matrix(nrow = nparticles, ncol = ncol(xparticles2))), datalength + 1)
  whistory <- rep(list(rep(0, nparticles)), datalength + 1)
  ahistory <- rep(list(rep(0, nparticles)), datalength)
  
  xhistory[[1]] <- xparticles2
  whistory[[1]] <- normweights2
  
  # step t > 1
  for (time in 1:datalength){
    ancestors1 <- system_ref$ahistory[[time]]
    ancestors2 <- coupled_resampling_given(xparticles1, xparticles2, normweights1, normweights2, 
                                     resampling_parameters, ancestors1, runif(2*nparticles + 1))
    ahistory[[time]] <- ancestors2
    xparticles2 <- xparticles2[ancestors2,]
    #
    xparticles2 <- model$rtransition(xparticles2, theta2, time, randomness$transition[,time])
    logw2 <- model$dmeasurement(xparticles2, theta2, observations[time,])
    maxlw2 <- max(logw2)
    w2 <- exp(logw2 - maxlw2)
    ll2 <- ll2 + maxlw2 + log(mean(w2))
    normweights2 <- w2 / sum(w2)
    #
    xparticles1 <- system_ref$xhistory[[time + 1]]
    normweights1 <- system_ref$whistory[[time + 1]]
    #
    xhistory[[time+1]] <- xparticles2
    whistory[[time+1]] <- normweights2
  }
  return(list(ll = ll2, xhistory = xhistory, whistory = whistory, ahistory = ahistory))
}

