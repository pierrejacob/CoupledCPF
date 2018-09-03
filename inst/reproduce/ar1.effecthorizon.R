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
# algoparameters <- list(coupled_resampling = CR_indexmatching, lambda = 0)
theta <- 0.9
nrep <- 500
# time horizons
# datalengths <- c(100, 200, 400, 800, 1600)
datalengths <- c(50, 100, 200, 400, 800)
# numbers of particles
seq_nparticles <- c(128, 256, 512, 1024, 2048)
# filename for results
savefilename <- paste0("ar1.effecttimehorizon.R", nrep, ".RData")

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
    return(CPF(N, ar1, theta, observations, ref_trajectory = NULL, with_as = with_as))
  }
  single_kernel_RB <- function(chain_state, h){
    return(CPF_RB(N, ar1, theta, observations, ref_trajectory = chain_state, with_as = with_as, h = h))
  }
  coupled_kernel_RB <- function(chain_state1, chain_state2, h){
    return(CPF_coupled_RB(N, ar1, theta, observations, ref_trajectory1 = chain_state1, ref_trajectory2 = chain_state2, 
                          coupled_resampling = coupled_resampling, with_as = with_as, h = h))
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
# Without ancestor sampling
estimates_woas.df <- data.frame()
with_as <- FALSE
for (i in 1:length(datalengths)){
  datalength <- datalengths[i]
  N <- seq_nparticles[i]
  observations <- matrix(original_observations[1:datalength,], ncol = 1)
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
