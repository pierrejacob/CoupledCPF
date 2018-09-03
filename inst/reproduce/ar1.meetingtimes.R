# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
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
coupled_resampling <- CR_indexmatching
theta <- 0.9
# number of independent replicates (set to e.g. 1000)
nrep <- 1000
# different numbers of particles
seq_nparticles <- c(16, 128, 256, 512, 1024)
#### Generate unbiased estimators with index-matching resampling
#### and various numbers of particles
#### and with ancestor sampling
meetings.df <- data.frame()
with_as <- TRUE
h <- function(x) x
# 
for (N in seq_nparticles){
  print(N)
  rinit <- function(){
    return(CPF(N, ar1, theta, observations, ref_trajectory = NULL, with_as = with_as))
  }
  single_kernel_RB <- function(chain_state, h){
    return(CPF_RB(N, ar1, theta, observations, ref_trajectory = chain_state, with_as = with_as, h = h))
  }
  coupled_kernel_RB <- function(chain_state1, chain_state2, h){
    return(CPF_coupled_RB(N, ar1, theta, observations, ref_trajectory1 = chain_state1, ref_trajectory2 = chain_state2,
                          coupled_resampling = coupled_resampling, with_as = with_as, h = h))
  }
  meetings_ <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
    res <- unbiasedestimator_RB(single_kernel_RB, coupled_kernel_RB, rinit, h = h, k = 0, m = 0, max_iterations = 1e4)
    estimate <- res$uestimator
    if (is.null(dim(estimate))){
      estimate <- matrix(estimate, nrow = 1)
    }
    data.frame(irep = irep, iteration = res$iteration, meetingtime = res$meetingtime,
               nparticles = N, with_as = with_as)
  }
  meetings.df <- rbind(meetings.df, meetings_)
  save(observations, nrep, datalength, dimension, theta, meetings.df, file = "ar1.meetingtime.RData")
}

load("ar1.meetingtime.RData")
head(meetings.df)
