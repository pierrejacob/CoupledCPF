#### This scripts briefly shows how this package works
# First, some initialization
# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(doParallel)
library(doRNG)
library(dplyr)
registerDoParallel(cores = detectCores()-2)
setmytheme()
set.seed(17)
module_tree <<- Module("module_tree", PACKAGE = "CoupledCPF")
TreeClass <<- module_tree$Tree

###
# We are going to work with the nonlinear growth model
# first we import it
gss <- get_gss()
# let us generate some data from the model
dimension <- 1
theta <- c()
datalength <- 50
x <- matrix(0, nrow = datalength + 1, ncol = 1)
observations <- matrix(0, nrow = datalength, ncol = 1)
x[1,] <- gss$rinit(1, theta, gss$rinit_rand(1, theta))
# equivalent to
#x[1,] <- rnorm(1, 0, sd = sqrt(2))
for (i in 1:datalength){
  x[i+1,] <- gss$rtransition(x[i,], theta, i, gss$rtransition_rand(1, theta))
  # equivalent to 
  # x[i+1,] <- 0.5 * x[i,] + 25 * x[i,] / (1 + x[i,]^2) + 8 * cos(1.2*(i - 1)) + sqrt(10) * rnorm(1)
  observations[i,] <- (x[i+1,]^2) / 20 + rnorm(1)
}
# plot the latent process
plot(0:datalength, x, type = "l")
# plot the observations
plot(1:datalength, observations, type = "l")
#

# now let us estimate the smoothing mean with the proposed method
# number of particles
N <- 512
# with ancestor sampling
with_as <- TRUE
# coupled resampling scheme
coupled_resampling <- CR_indexmatching
# initial distribution of the chains: draw a trajectory from a particle filter
rinit <- function(){
  return(CPF(N, gss, theta, observations, ref_trajectory = NULL, with_as = FALSE))
}
# single kernel of CPF, with Rao-Blackwellization and ancestor sampling
single_kernel_RB <- function(chain_state, h){
  return(CPF_RB(N, gss, theta, observations, ref_trajectory = chain_state, with_as = with_as, h = h))
}
# coupled kernel of CPF, with Rao-Blackwellization and ancestor sampling
coupled_kernel_RB <- function(chain_state1, chain_state2, h){
  return(CPF_coupled_RB(N, gss, theta, observations, ref_trajectory1 = chain_state1, ref_trajectory2 = chain_state2,
                        coupled_resampling = coupled_resampling, with_as = with_as, h = h))
}
# test function: identity, so that we estimate smoothing means
h <- function(x) x
# run one unbiased estimator
res <- unbiasedestimator_RB(single_kernel_RB, coupled_kernel_RB, rinit, h = h, k = 0, m = 1)
res$meetingtime

# how many estimators do we want?
R <- 10
# the following should take a few seconds if you have a few processing units 
estimates.df <- foreach(r = 1:R, .combine = rbind) %dorng% {
  res <- unbiasedestimator_RB(single_kernel_RB, coupled_kernel_RB, rinit, h = h, k = 0, m = 1)
  data.frame(irep = r, time = 0:datalength, estimate = as.numeric(res$uestimator), meetingtime = res$meetingtime)
}

# distribution of meeting times
hist(estimates.df$meetingtime)
# from this we decide to choose k such that the meeting time is less than k with large probability
# e.g. k = 20
# and we choose m larger than k, e.g. m = 40

estimates.kandm.df <- foreach(r = 1:R, .combine = rbind) %dorng% {
  res <- unbiasedestimator_RB(single_kernel_RB, coupled_kernel_RB, rinit, h = h, k = 20, m = 40)
  data.frame(irep = r, time = 0:datalength, estimate = as.numeric(res$uestimator), meetingtime = res$meetingtime)
}

means.kandm.df <- estimates.kandm.df %>% group_by(time) %>% summarise(m = mean(estimate), s = sd(estimate))

# with k = 0 and m = 1, inaccurate results
g <- ggplot(means.kandm.df, aes(x = time)) + geom_errorbar(aes(ymin = m - 2 * s / sqrt(R), ymax = m + 2 * s/ sqrt(R))) + geom_line(aes(y = m))
g <- g + ylab("latent process means") + xlab("time")
g

# we can compare to the result of a fixed lag smoother with 2^14 particles and a lag of 20
lag <- 20
nparticles <- 2^14
fixedlagestimates <- fixedlagsmoother(nparticles, lag, gss, theta, observations)
g <- g + geom_line(data = NULL, aes(x = 0:datalength, y = fixedlagestimates$smoothing_means, ymin = NULL, ymax = NULL), colour = "red", linetype = 2, size = 2)
print(g)
# we see an agreement between the two procedures.




