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


load("gssdata.RData")

dimension <- 1

gss <- get_gss()
theta <- c()
datalength <- 50
observations <- matrix(observations[1:datalength,], ncol = 1)


nparticles <- 2^10
coupled_resampling <- CR_indexmatching
with_as <- TRUE
lambda <- 0
algoparameters <- list(nparticles = nparticles, coupled_resampling = coupled_resampling, lambda = lambda, with_as = with_as)


# # see how many iterations we need to get full coalescence
max_niterations <- 100
nrep <- 10

times_coal <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  res <- rheeglynn_estimator(observations, gss, theta, algoparameters, max_niterations = max_niterations)
  iteration <-  res$iteration
  data.frame(irep = irep, iteration = iteration)
}
#
ggplot(times_coal, aes(x = iteration)) + geom_histogram()

nrep <- 1000
estimates <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  res <- rheeglynn_estimator(observations, gss, theta, algoparameters)
  data.frame(irep = irep, time = 0:datalength, estimate = t(res$estimate), iteration = res$iteration,
             truncation = res$truncation)
}

# names(res)
save(nrep, theta, algoparameters, observations, datalength, dimension, estimates, 
     file = "gss_unbiased.autostop.RData")

