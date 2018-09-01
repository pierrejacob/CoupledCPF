# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(doParallel)
library(doRNG)
library(dplyr)
library(ggplot2)
registerDoParallel(cores = detectCores()-2)
# fix the random seed
set.seed(17)
datalength <- 100
# get hidden autoregressive model
get_model <- function(){
  #
  rinit <- function(nparticles, theta, rand, precomputed, ...){
    return(matrix(rand, nrow = 1))
  }
  rinit_rand <- function(nparticles, theta){
    return(rnorm(nparticles))
  }
  #
  rtransition <- function(xparticles, theta, time, rand, precomputed, ...){
    return(theta * xparticles + rand)
  }
  #
  rtransition_rand <- function(nparticles, theta){
    return(rnorm(nparticles))
  }
  #
  dtransition <- function(next_x, xparticles, theta, time, precomputed, ...){
    return(dnorm(x = theta * xparticles[1,], mean = next_x, sd = 1, log = TRUE))
  }
  
  dmeasurement <- function(xparticles, theta, observation, precomputed, ...){
    return(dnorm(x = xparticles[1,], mean = observation, sd = 1, log = TRUE))
  }
  
  precompute <- function(theta){
    return(list())
  }
  ar_model <- list(rinit = rinit, rinit_rand = rinit_rand, rtransition = rtransition,
                   rtransition_rand = rtransition_rand,
                   dtransition = dtransition, 
                   dmeasurement = dmeasurement, precompute = precompute, dimension = 1)
  return(ar_model)
}

# generate observations
theta <- 0.95

observations <- matrix(nrow = datalength, ncol = 1)
x_t <- fast_rmvnorm(1, rep(0, 1), diag(1, nrow = 1, ncol = 1))
for (time in 1:datalength){
  x_t <- t(theta %*% t(x_t)) + fast_rmvnorm(1, rep(0, 1), diag(1, nrow = 1, ncol = 1))
  observations[time,] <- x_t + fast_rmvnorm(1, rep(0, 1), diag(1, nrow = 1, ncol = 1))
}
Y <- observations[,1]

library(dlm)
d <- dlm(list(m0 = 0, C0 = 1, FF = 1, V = 1, GG = theta, W = 1))
dsres <- dlmSmooth(observations[,1], d)


plot(x = 1:datalength, y = observations[,1], type = "l")
lines(x = 0:datalength, y = dsres$s, col = "red")
ks.df <- data.frame(time = 0:datalength, ksm = dsres$s)

# unbiased estimation
ar1 <- get_model()
module_tree <<- Module("module_tree", PACKAGE = "CoupledCPF")
TreeClass <<- module_tree$Tree
# with ancestor sampling
with_as <- TRUE
# coupled resampling scheme
coupled_resampling <- CR_indexmatching
# number of particles
N <- 512
# initial distribution of the chains
rinit <- function(){
  return(CPF(N, ar1, theta, observations, ref_trajectory = NULL, with_as = with_as))
}
# single kernel of CPF, with Rao-Blackwellization
single_kernel_RB <- function(chain_state, h){
  return(CPF_RB(N, ar1, theta, observations, ref_trajectory = chain_state, with_as = with_as, h = h))
}
# coupled kernel of CPF, with Rao-Blackwellization
coupled_kernel_RB <- function(chain_state1, chain_state2, h){
  return(CPF_coupled_RB(N, ar1, theta, observations, ref_trajectory1 = chain_state1, ref_trajectory2 = chain_state2,
                        coupled_resampling = coupled_resampling, with_as = with_as, h = h))
}
# test function
h <- function(x) x
# run one unbiased estimator
res <- unbiasedestimator_RB(single_kernel_RB, coupled_kernel_RB, rinit, h = h, k = 0, m = 1)
res$meetingtime
# now run 100 to get the distribution of meeting times
nrep <- 10
results.df <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  res <- unbiasedestimator_RB(single_kernel_RB, coupled_kernel_RB, rinit, h = h, k = 0, m = 1)
  data.frame(irep = irep, meeting = res$iteration, time = 0:datalength, estimate = as.numeric(res$uestimator))
}
summary(results.df$meeting)
hist(results.df$meeting)

# so we can choose k = 25 and m = 100, say
nrep <- 10
results.kandm.df <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  res <- unbiasedestimator_RB(single_kernel_RB, coupled_kernel_RB, rinit, h = h, k = 25, m = 100)
  data.frame(irep = irep, meeting = res$iteration, time = 0:datalength, estimate = as.numeric(res$uestimator))
}

# smoothing means per time 
means.df <- results.df %>% group_by(time) %>% summarise(m = mean(estimate), s = sd(estimate))
means.kandm.df <- results.kandm.df %>% group_by(time) %>% summarise(m = mean(estimate), s = sd(estimate))

# with k = 0 and m = 1, inaccurate results
g <- ggplot(means.df, aes(x = time)) + geom_errorbar(aes(ymin = m - 2 * s / sqrt(nrep), ymax = m + 2 * s/ sqrt(nrep)))
g <- g + geom_line(data=ks.df, aes(x = time, y = ksm, colour = NULL, group = NULL), colour = "red")
g <- g + ylab("latent process means") + ylim(-10, 5) + xlab("time")
g

# with k = 25 and m = 100, precise results
g <- ggplot(means.kandm.df, aes(x = time)) + geom_errorbar(aes(ymin = m - 2 * s / sqrt(nrep), ymax = m + 2 * s/ sqrt(nrep)))
g <- g + geom_line(data=ks.df, aes(x = time, y = ksm, colour = NULL, group = NULL), colour = "red")
g <- g + ylab("latent process means") + ylim(-10, 5) + xlab("time")
g
