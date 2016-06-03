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
# Load synthetic dataset
load("ar1data.RData")
datalength <- 500
observations <- matrix(observations[1:datalength,], ncol = 1)

algoparameters <- list(coupled_resampling = CR_indexmatching, lambda = 0)

theta <- 0.9

nrep <- 1000
seq_nparticles <- c(256, 512, 1024, 2048, 4096)
#### Generate unbiased estimators with index-matching resampling
#### and various numbers of particles
#### and with and without ancestor sampling
estimates_was.df <- data.frame()
with_as <- TRUE
for (N in seq_nparticles){
  print(N)
  algoparameters$nparticles <- N
  algoparameters$with_as <- with_as
  estimates_was.df_ <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
    res <- rheeglynn_estimator(observations, ar1, theta, algoparameters)
    data.frame(irep = irep, time = 0:datalength, estimate = t(res$estimate), iteration = res$iteration,
               nparticles = N, with_as = with_as)
  }
  estimates_was.df <- rbind(estimates_was.df, estimates_was.df_)
  save(observations, nrep, datalength, dimension, theta, estimates_was.df, file = "ar1.effectnparticles.RData")
}
estimates_was.df %>% tail

estimates_woas.df <- data.frame()
with_as <- FALSE
for (N in seq_nparticles){
  print(N)
  algoparameters$nparticles <- N
  algoparameters$with_as <- with_as
  estimates_woas.df_ <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
    res <- rheeglynn_estimator(observations, ar1, theta, algoparameters)
    data.frame(irep = irep, time = 0:datalength, estimate = t(res$estimate), iteration = res$iteration,
               nparticles = N, with_as = with_as)
  }
  estimates_woas.df <- rbind(estimates_woas.df, estimates_woas.df_)
  save(observations, nrep, datalength, dimension, theta, estimates_was.df, estimates_woas.df, file = "ar1.effectnparticles.RData")
}
