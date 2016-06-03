# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(doRNG)
ncores <- 10
registerDoMC(cores = ncores)
# load custom theme for the plots
# it requires the packages ggplot2, gridExtra, reshape
setmytheme()
# fix the random seed
set.seed(17)
module_tree <<- Module("module_tree", PACKAGE = "CoupledCPF")
TreeClass <<- module_tree$Tree

dimension <- 1
ar1 <- get_ar(dimension)
load("ar1data.RData")
datalength <- 500
observations <- matrix(observations[1:datalength,], ncol = 1)

nparticles <- 2^9
coupled_resampling <- CR_indexmatching
with_as <- FALSE
algoparameters <- list(nparticles = nparticles, coupled_resampling = coupled_resampling, with_as = with_as)

theta <- 0.9

nrep <- 1000

estimates.df <- data.frame()
seq_lambdas <- c(0, 0.025, 0.05)
for (lambda in seq_lambdas){
  print(lambda)
  algoparameters$lambda <- lambda
  estimates.df_ <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
    res <- rheeglynn_estimator(observations, ar1, theta, algoparameters)
    data.frame(irep = irep, time = 0:datalength, estimate = t(res$estimate), iteration = res$iteration, truncation = res$truncation,
               lambda = lambda)
  }
  estimates.df <- rbind(estimates.df, estimates.df_)
  save(estimates.df, theta, dimension, observations, nrep, file = "ar1.effecttruncation.RData")
}
