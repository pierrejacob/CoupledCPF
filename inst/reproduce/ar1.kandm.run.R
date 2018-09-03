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
# get model
dimension <- 1
ar1 <- get_ar(dimension)
# Load synthetic dataset
load("ar1data.RData")
datalength <- 100
observations <- matrix(observations[1:datalength,], ncol = 1)
# choice of coupled resampling scheme
coupled_resampling <- CR_indexmatching
# parameter
theta <- 0.9
# number of independent replications (set to e.g. 1000)
nrep <- 1000
# number of particles
N <- 256
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
# with ancestor sampling
with_as <- TRUE
# test function
h <- function(x) x

# try different values of k 
ks <- c(10, 20)
# m multiplier, i.e. m = k * mfactor
mfactors <- c(1,2,3)
# produce estimates
estimates.df <- data.frame()
for (ik in 1:length(ks)){
  for (im in 1:length(mfactors)){
    k <- ks[ik]
    m <- k * mfactors[im]
    cat("k =", k, ", m =", m, "\n")
    # stop chains early if they take max_niterations to couple (shouldn't happen!)
    max_niterations <- 1e5
    # generate estimates in parallel
    estimates <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
      res <- unbiasedestimator_RB(single_kernel_RB, coupled_kernel_RB, rinit, h = h, k = k, m = m, max_iterations = max_niterations)
      estimate <- res$uestimator
      if (is.null(dim(estimate))){
        estimate <- matrix(estimate, nrow = 1)
      }
      data.frame(irep = irep, iteration = res$iteration, meetingtime = res$meetingtime, time = 0:datalength,
                 estimate = t(estimate), k = k, m = m)
    }
    estimates.df <- rbind(estimates.df, estimates)
    save(estimates.df, file = "ar.kandm.RData")
  }
}
load("ar.kandm.RData")

# function to compute computational cost of each run
cost <- function(tau, m) 3 + 2 * (tau-1) + pmax(0, m - tau)
library(dplyr)
df <- estimates.df %>% mutate(c = cost(meetingtime, iteration))
df <- df %>% group_by(time,k,m) %>% summarise(meanestim = mean(estimate), varestim = var(estimate), meancost = mean(c))
df <- df %>% mutate(inef = varestim * meancost)
df <- df %>% mutate(ymin = meanestim - 2 * sqrt(varestim/nrep), ymax = meanestim + 2 * sqrt(varestim/nrep))
# get a peek at the results for time 0
# i.e. the task is the estimation of E[X_0|Y_{1:T}] with T = datalength
df %>% filter(time == 0)
# histogram of meeting times 
hist(estimates.df$meetingtime)
