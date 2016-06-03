# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(doRNG)
library(ggthemes)
ncores <- 5
registerDoMC(cores = ncores)
# load custom theme for the plots
# it requires the packages ggplot2, gridExtra, reshape
setmytheme()
# fix the random seed
set.seed(17)
module_tree <<- Module("module_tree", PACKAGE = "CoupledCPF")
TreeClass <<- module_tree$Tree

logit <- function(z) log(z / (1 - z))
expit <- function(z) 1 / (1 + exp(-z))


load("pzdata.RData")

logit <- function(z) log(z / (1 - z))
expit <- function(z) 1 / (1 + exp(-z))

theta_dgp <- c(0.7, 0.5, 0.25, 0.3, 0.1, 0.1)
# generate dataset
datalength <- 365
observations <- matrix(observations[1:datalength,], ncol = 1)

theta <- c(theta_dgp[1], log(theta_dgp[2]), logit(theta_dgp[3]), logit(theta_dgp[4]), logit(theta_dgp[5]), logit(theta_dgp[6]))

dimension <- 2
pz <- get_pz()


nparticles <- 2^12
algoparameters <- list(nparticles = nparticles, coupled_resampling = CR_indexmatching, lambda = 0, with_as = FALSE)

nrep <- 1000

estimates <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  res <- rheeglynn_estimator_RB(observations, pz, theta, algoparameters)
  data.frame(irep = irep, time = 0:datalength, estimate = t(res$estimate), iteration = res$iteration)
}

filename <- paste0("pz.RB.R", nrep, "N", nparticles, "T", datalength, ".RData")
save(nrep, theta, algoparameters, observations, datalength, dimension, estimates, file = filename)

