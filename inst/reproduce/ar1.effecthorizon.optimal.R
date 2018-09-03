# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(doRNG)
library(doParallel)
registerDoParallel(cores = detectCores() - 2)
# load custom theme for the plots
# it requires the packages ggplot2, gridExtra, reshape
setmytheme()
# fix the random seed
set.seed(17)
module_tree <<- Module("module_tree", PACKAGE = "CoupledCPF")
TreeClass <<- module_tree$Tree

dimension <- 1
ar1 <- get_ar(dimension)
# Load synthetic dataset
load("ar1data.RData")
original_observations <- observations
coupled_resampling <- CR_indexmatching
theta <- 0.9

## functions to perform auxiliary particle filter instead of bootstrap particle filter
rinit_model <- function(nparticles){
  return(matrix(rnorm(nparticles), nrow = 1))
}
dtransition <- function(x, x_last, theta, time, model_precomputed){
  return(dnorm(x, mean = theta*x_last, sd = 1, log = TRUE))
}
rtransition_optimal <- function(xparticles, y){
  return(theta * xparticles / 2 + y / 2 + matrix(rnorm(ncol(xparticles), mean = 0, sd = sqrt(1/2)), nrow = 1))
}
dmeasurement_optimal <- function(xparticles, xparticles_next, y){
  return(dnorm(y, mean = theta * xparticles[1,], sd = sqrt(2), log = TRUE))
}
# CPF with APF
CPF_RB_optimal <- function(nparticles, observations, ref_trajectory = NULL, with_as = FALSE,
                           h = function(trajectory){ return(trajectory) }){
  datalength <- nrow(observations)
  # create tree representation of the trajectories
  Tree <- new(TreeClass, nparticles, 10*nparticles*dimension, dimension)
  # initialization
  xparticles <- rinit_model(nparticles)
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
  filtering_means <- rep(0, datalength+1)
  filtering_means[1] <- mean(xparticles[1,])
  # step t > 1
  for (time in 1:datalength){
    ancestors <- multinomial_resampling_n(normweights, nparticles)
    # if no observation or first time, no resampling
    if (time == 1 || (time > 1 && is.na(observations[time-1,1]))){
      ancestors <- 1:nparticles
    }
    #
    xparticles <- xparticles[,ancestors,drop=F]
    xparticles_next <- xparticles
    if (is.null(dim(xparticles_next))) xparticles_next <- matrix(xparticles_next, nrow = dimension)
    xparticles_next <- rtransition_optimal(xparticles_next, observations[time,1])
    if (!is.null(ref_trajectory)){
      xparticles_next[,nparticles] <- ref_trajectory[,time+1]
      xparticles[,nparticles] <- ref_trajectory[,time]
      if (with_as){
        # Ancestor sampling
        logm <- dtransition(ref_trajectory[,time+1], x_last, theta, time, NULL)
        logm <- log(normweights) + logm
        w_as <- exp(logm - max(logm))
        w_as <- w_as / sum(w_as)
        ancestors[nparticles] <- systematic_resampling_n(w_as, 1, runif(1))
        x_last <- xparticles_next
      } else {
        ancestors[nparticles] <- nparticles
      }
    }
    #    
    logw <- dmeasurement_optimal(xparticles, xparticles_next, observations[time,1])
    maxlw <- max(logw)
    w <- exp(logw - maxlw)
    normweights <- w / sum(w)
    #
    xparticles <- xparticles_next
    filtering_means[time+1] <- sum(normweights * xparticles[1,])
    Tree$update(xparticles, ancestors - 1)    
  }
  trajectories <- array(dim = c(dimension, datalength + 1, nparticles))
  estimate <- 0
  for (k in 0:(nparticles-1)){
    trajectories[,,k+1] <- Tree$get_path(k)
    estimate <- estimate + normweights[k+1] * h(trajectories[,,k+1])
  }
  new_trajectory <- matrix(trajectories[,,systematic_resampling_n(normweights, 1, runif(1))], nrow = dimension)
  return(list(new_trajectory = new_trajectory, estimate = estimate, filtering_means = filtering_means))
}

coupled_rinit_model <- function(nparticles){
  xp1 <- matrix(rnorm(nparticles), nrow = 1)
  return(list(xp1 = xp1, xp2 = xp1))
}
coupled_rtransition_optimal <- function(xparticles1, xparticles2, y){
  noises <- matrix(rnorm(ncol(xparticles1), sd = sqrt(1/2)), nrow = 1)
  return(list(xp1 = theta * xparticles1 / 2 + y / 2 + noises, xp2 = theta * xparticles2 / 2 + y / 2 + noises))
}

CPF_coupled_RB_optimal <- function(nparticles, observations, ref_trajectory1, ref_trajectory2, 
                                   with_as = FALSE, h = function(trajectory){ return(trajectory) }){
  #
  datalength <- nrow(observations)
  # create tree representation of the trajectories
  Tree1 <- new(TreeClass, nparticles, 10*nparticles*dimension, dimension)
  Tree2 <- new(TreeClass, nparticles, 10*nparticles*dimension, dimension)
  # initialization
  init_res <- coupled_rinit_model(nparticles)
  xparticles1 <- init_res$xp1
  xparticles1[,nparticles] <- ref_trajectory1[,1]
  Tree1$init(xparticles1)
  normweights1 <- rep(1/nparticles, nparticles)
  #  
  xparticles2 <- init_res$xp2
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
    ancestors <- CR_indexmatching(xparticles1, xparticles2, normweights1, normweights2)
    ancestors1 <- ancestors[,1]
    ancestors2 <- ancestors[,2]
    # if no observation or first time, no resampling
    if (time == 1 || (time > 1 && is.na(observations[time-1,1]))){
      ancestors1 <- 1:nparticles
      ancestors2 <- 1:nparticles
    }
    #
    xparticles1 <- xparticles1[,ancestors1,drop=F]
    xparticles2 <- xparticles2[,ancestors2,drop=F]
    
    #
    xparticles_next1 <- xparticles1
    xparticles_next2 <- xparticles2
    #
    coupled_transition_res <- coupled_rtransition_optimal(xparticles_next1, xparticles_next2, observations[time,1])
    xparticles_next1 <- coupled_transition_res$xp1
    xparticles_next2 <- coupled_transition_res$xp2
    if (is.null(dim(xparticles_next1))) xparticles1 <- matrix(xparticles_next1, nrow = dimension)
    if (is.null(dim(xparticles_next2))) xparticles2 <- matrix(xparticles_next2, nrow = dimension)
    #
    xparticles1[,nparticles] <- ref_trajectory1[,time]
    xparticles2[,nparticles] <- ref_trajectory2[,time]
    xparticles_next1[,nparticles] <- ref_trajectory1[,time+1]
    xparticles_next2[,nparticles] <- ref_trajectory2[,time+1]
    if (with_as){
      # % Ancestor sampling
      logm1 <- dtransition(ref_trajectory1[,time+1], x_last1, theta, time, NULL)
      logm1 <- log(normweights1) + logm1
      w_as1 <- exp(logm1 - max(logm1))
      w_as1 <- w_as1 / sum(w_as1)
      # unif_resampling_as <- runif(1)
      # ancestors1[nparticles] = systematic_resampling_n(w_as1, 1, unif_resampling_as)
      x_last1 <- xparticles_next1
      #
      logm2 <- dtransition(ref_trajectory2[,time+1], x_last2, theta, time, NULL)
      logm2 <- log(normweights2) + logm2
      w_as2 <- exp(logm2 - max(logm2))
      w_as2 <- w_as2 / sum(w_as2)
      # ancestors2[nparticles] = systematic_resampling_n(w_as2, 1, unif_resampling_as)
      x_last2 <- xparticles_next2
      ## sample ancestor indices with index-coupled resampling
      ancestors_ <- CR_indexmatching(NULL, NULL, normweights1 = w_as1, normweights2 = w_as2, ntrials = 1)
      ancestors1[nparticles] <- ancestors_[1,1]
      ancestors2[nparticles] <- ancestors_[1,2]
    } else {
      ancestors1[nparticles] <- nparticles
      ancestors2[nparticles] <- nparticles
    }
    #    
    logw1 <- dmeasurement_optimal(xparticles1, xparticles_next1, observations[time,1])
    logw2 <- dmeasurement_optimal(xparticles2, xparticles_next2, observations[time,1])
    #
    maxlw1 <- max(logw1)
    w1 <- exp(logw1 - maxlw1)
    normweights1 <- w1 / sum(w1)
    #
    maxlw2 <- max(logw2)
    w2 <- exp(logw2 - maxlw2)
    normweights2 <- w2 / sum(w2)    
    #
    xparticles1 <- xparticles_next1
    xparticles2 <- xparticles_next2
    Tree1$update(xparticles1, ancestors1 - 1)    
    Tree2$update(xparticles2, ancestors2 - 1)    
  }
  # u <- runif(1)
  # k_path1 <- systematic_resampling_n(normweights1, 1, u)
  # k_path2 <- systematic_resampling_n(normweights2, 1, u)
  ancestors_ <- CR_indexmatching(NULL, NULL, normweights1 = normweights1, normweights2 = normweights2, ntrials = 1)
  k_path1 <- ancestors_[1,1]
  k_path2 <- ancestors_[1,2]
  ##
  trajectories1 <- array(dim = c(dimension, datalength + 1, nparticles))
  trajectories2 <- array(dim = c(dimension, datalength + 1, nparticles))
  estimate1 <- 0
  estimate2 <- 0
  for (k in 0:(nparticles-1)){
    trajectories1[,,k+1] <- Tree1$get_path(k)
    trajectories2[,,k+1] <- Tree2$get_path(k)
    estimate1 <- estimate1 + normweights1[k+1] * h(trajectories1[,,k+1])
    estimate2 <- estimate2 + normweights2[k+1] * h(trajectories2[,,k+1])
  }
  
  new_trajectory1 <- matrix(trajectories1[,,k_path1], nrow = dimension)
  new_trajectory2 <- matrix(trajectories2[,,k_path2], nrow = dimension)
  return(list(new_trajectory1 = new_trajectory1, new_trajectory2 = new_trajectory2,
              estimate1 = estimate1, estimate2 = estimate2))
}

## number of independent replications
nrep <- 500
# time horizons
# datalengths <- c(100, 200, 400, 800, 1600)
datalengths <- c(50, 100, 200, 400, 800)
# numbers of particles
seq_nparticles <- c(128, 256, 512, 1024, 2048)
# filename to store the results
savefilename <- paste0("ar1.effecttimehorizon.optimal.R", nrep, ".RData")

#### Generate unbiased estimators with index-matching resampling
#### and various numbers of particles
#### and with and without ancestor sampling
estimates_was.df <- data.frame()
with_as <- TRUE
h <- function(x) x

for (i in 1:length(datalengths)){
  datalength <- datalengths[i]
  N <- seq_nparticles[i]
  observations <- matrix(original_observations[1:datalength,], ncol = 1)
  print(N)
  rinit <- function(){
    return(CPF_RB_optimal(N, observations, ref_trajectory = NULL, with_as = with_as, h = h)$new_trajectory)
  }
  single_kernel_RB <- function(chain_state, h){
    return(CPF_RB_optimal(N, observations, ref_trajectory = chain_state, with_as = with_as, h = h))
  }
  coupled_kernel_RB <- function(chain_state1, chain_state2, h){
    return(CPF_coupled_RB_optimal(N, observations, ref_trajectory1 = chain_state1, ref_trajectory2 = chain_state2, 
                                  with_as = with_as, h = h))
  }
  estimates_was.df_ <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
    res <- unbiasedestimator_RB(single_kernel_RB, coupled_kernel_RB, rinit, h = h, k = 0, m = 0, max_iterations = 1e4)
    estimate <- res$uestimator
    if (is.null(dim(estimate))){
      estimate <- matrix(estimate, nrow = 1)
    }
    data.frame(irep = irep, time = 0:datalength, estimate = t(estimate), iteration = res$iteration,
               nparticles = N, datalength = datalength, with_as = with_as)
  }
  estimates_was.df <- rbind(estimates_was.df, estimates_was.df_)
  save(observations, nrep, datalength, dimension, theta, estimates_was.df, file = savefilename)
}

# without ancestor sampling
estimates_woas.df <- data.frame()
with_as <- FALSE
for (i in 1:length(datalengths)){
  datalength <- datalengths[i]
  N <- seq_nparticles[i]
  observations <- matrix(original_observations[1:datalength,], ncol = 1)
  print(N)
  rinit <- function(){
    return(CPF_RB_optimal(N, observations, ref_trajectory = NULL, with_as = with_as, h = h)$new_trajectory)
  }
  single_kernel_RB <- function(chain_state, h){
    return(CPF_RB_optimal(N, observations, ref_trajectory = chain_state, with_as = with_as, h = h))
  }
  coupled_kernel_RB <- function(chain_state1, chain_state2, h){
    return(CPF_coupled_RB_optimal(N, observations, ref_trajectory1 = chain_state1, ref_trajectory2 = chain_state2, 
                                  with_as = with_as, h = h))
  }
  estimates_woas.df_ <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
    res <- unbiasedestimator_RB(single_kernel_RB, coupled_kernel_RB, rinit, h = h, k = 0, m = 0, max_iterations = 1e4)
    estimate <- res$uestimator
    if (is.null(dim(estimate))){
      estimate <- matrix(estimate, nrow = 1)
    }
    data.frame(irep = irep, time = 0:datalength, estimate = t(estimate), iteration = res$iteration,
               nparticles = N, datalength = datalength, with_as = with_as)
  }
  estimates_woas.df <- rbind(estimates_woas.df, estimates_woas.df_)
  save(observations, nrep, datalength, dimension, theta, estimates_was.df, estimates_woas.df, file = savefilename)
}


