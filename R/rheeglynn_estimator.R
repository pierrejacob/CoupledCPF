#'@rdname unbiasedestimator
#'@title Rhee--Glynn smoothing estimator
#'@description Estimates the expectation of a function h with respect to the invariant distribution of a Markov chain
#'@param single_kernel function to sample the Markov chain
#'@param coupled_kernel function to sample the coupled Markov chain
#'@param rinit function to sample from the initial distribution
#'@param h test function to estimate (default to h(x) = x)
#'@param k iteration at which to start averaging (default to 0)
#'@param m iteration at which to stop averaging (default to 1)
#'@param max_niterations iteration at which to stop the while loop (default to Infinity)
#'@return a list with the estimator of the smoothing expectation (uestimator) and the meeting time (meetingtime)
#'@export
#'
unbiasedestimator <- function(single_kernel, coupled_kernel, rinit, h = function(x) x, k = 0, m = 1, max_iterations = Inf){
  chain_state1 <- rinit()
  chain_state2 <- rinit()
  # mcmcestimator computes the sum of h(X_t) for t=k,...,m
  mcmcestimator <- h(chain_state1)
  dimh <- length(mcmcestimator)
  if (k > 0){
    mcmcestimator <- rep(0, dimh)
  }
  # correction computes the sum of min(1, (t - k + 1) / (m - k + 1)) * (h(X_{t+1}) - h(X_t)) for t=k,...,max(m, tau - 1)
  correction <- rep(0, dimh)
  chain_state1 <- single_kernel(chain_state1)
  if (k == 0){
    correction <- correction + (min(1, (0 - k + 1)/(m - k + 1))) * (h(chain_state1) - h(chain_state2))
  }
  if (k <= 1 && m >= 1){
    mcmcestimator <- mcmcestimator + h(chain_state1)
  }
  iter <- 1
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){
      chain_state1 <- single_kernel(chain_state1)
      chain_state2 <- chain_state1
      if (k <= iter && iter <= m){
        mcmcestimator <- mcmcestimator + h(chain_state1)
      }
    } else {
      res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2)
      chain_state1 <- res_coupled_kernel$new_trajectory1
      chain_state2 <- res_coupled_kernel$new_trajectory2
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
      if (k <= iter){
        if (iter <= m){
          mcmcestimator <- mcmcestimator + h(chain_state1)
        }
        correction <- correction + (min(1, (iter-1 - k + 1)/(m - k + 1))) * (h(chain_state1) - h(chain_state2))
      }
    }
    # stop after max(m, tau) steps
    if (iter >= max(meetingtime, m)){
      finished <- TRUE
    }
  }
  mcmcestimator <- mcmcestimator / (m - k + 1)
  uestimator <- mcmcestimator + correction
  return(list(mcmcestimator = mcmcestimator, correction = correction, uestimator = uestimator,
              meetingtime = meetingtime, iteration = iter, finished = finished))
}

#'@rdname unbiasedestimator_RB
#'@title Rhee--Glynn smoothing estimator with Rao-Blackwellization
#'@description Estimates the expectation of a function h with respect to the invariant distribution of a Markov chain
#'@param single_kernel_RB function to sample the Markov chain and return an estimate of integral {h(x) pi(x)}
#'@param coupled_kernel function to sample the coupled Markov chain and return two estimates of integral {h(x) pi(x)}
#'@param rinit function to sample from the initial distribution
#'@param h test function to estimate (default to h(x) = x)
#'@param k iteration at which to start averaging (default to 0)
#'@param m iteration at which to stop averaging (default to 1)
#'@param max_niterations iteration at which to stop the while loop (default to Infinity)
#'@return a list with the estimator of the smoothing expectation (uestimator) and the meeting time (meetingtime)
#'@export
unbiasedestimator_RB <- function(single_kernel_RB, coupled_kernel_RB, rinit, h = function(x) x, k = 0, m = 1, max_iterations = Inf){
  chain_state1 <- rinit()
  chain_state2 <- rinit()
  # mcmcestimator computes the sum of h(X_t) for t=k,...,m
  mcmcestimator <- h(chain_state1)
  dimh <- length(mcmcestimator)
  if (k > 0){
    mcmcestimator <- rep(0, dimh)
  }
  # correction computes the sum of min(1, (t - k + 1) / (m - k + 1)) * (h(X_{t+1}) - h(X_t)) for t=k,...,max(m, tau - 1)
  correction <- rep(0, dimh)
  res <- single_kernel_RB(chain_state1, h)
  chain_state1 <- res$new_trajectory
  if (k == 0){
    correction <- correction + (min(1, (0 - k + 1)/(m - k + 1))) * (h(chain_state1) - h(chain_state2))
  }
  if (k <= 1 && m >= 1){
    mcmcestimator <- mcmcestimator + res$estimate # h(chain_state1)
  }
  iter <- 1
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){
      res <- single_kernel_RB(chain_state1, h)
      chain_state1 <- res$new_trajectory
      chain_state2 <- chain_state1
      if (k <= iter && iter <= m){
        mcmcestimator <- mcmcestimator + res$estimate # h(chain_state1)
      }
    } else {
      res_coupled_kernel <- coupled_kernel_RB(chain_state1, chain_state2, h)
      chain_state1 <- res_coupled_kernel$new_trajectory1
      chain_state2 <- res_coupled_kernel$new_trajectory2
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
      if (k <= iter){
        if (iter <= m){
          mcmcestimator <- mcmcestimator + res_coupled_kernel$estimate1
        }
        correction <- correction + (min(1, (iter-1 - k + 1)/(m - k + 1))) * (res_coupled_kernel$estimate1 - res_coupled_kernel$estimate2)
      }
    }
    # stop after max(m, tau) steps
    if (iter >= max(meetingtime, m)){
      finished <- TRUE
    }
  }
  mcmcestimator <- mcmcestimator / (m - k + 1)
  uestimator <- mcmcestimator + correction
  return(list(mcmcestimator = mcmcestimator, correction = correction, uestimator = uestimator,
              meetingtime = meetingtime, iteration = iter, finished = finished))
}

