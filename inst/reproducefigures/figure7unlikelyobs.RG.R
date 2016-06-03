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
datalength <- 10

model <- get_unlikely()
observations <- matrix(NA, nrow = datalength, ncol = 1)
observations[datalength,1] <- 1

#parameter = (rho, tau, sigma, tau_0)
# where x = rho * x + tau epsilon_t
# and   y = x + sigma eta_t at the final step
theta <- c(0.9, 0.1, 0.1, 0.1)

#####
nrep <- 10000
algoparameters <- list(coupled_resampling = CR_indexmatching, lambda = 0, with_as = TRUE)

algoparameters$nparticles <- 2^7
RGestimates <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  res <- rheeglynn_estimator(observations, model, theta, algoparameters)
  data.frame(irep = irep, estimate = res$estimate, iteration = res$iteration)
}
filename <- paste0("unlikely.RG.R", nrep, "N", algoparameters$nparticles, "T", datalength, ".RData")
save(algoparameters, RGestimates, nrep, file = filename)

algoparameters$nparticles <- 2^8
RGestimates <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  res <- rheeglynn_estimator(observations, model, theta, algoparameters)
  data.frame(irep = irep, estimate = res$estimate, iteration = res$iteration)
}
filename <- paste0("unlikely.RG.R", nrep, "N", algoparameters$nparticles, "T", datalength, ".RData")
save(algoparameters, RGestimates, nrep, file = filename)

algoparameters$nparticles <- 2^9
RGestimates <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  res <- rheeglynn_estimator(observations, model, theta, algoparameters)
  data.frame(irep = irep, estimate = res$estimate, iteration = res$iteration)
}
filename <- paste0("unlikely.RG.R", nrep, "N", algoparameters$nparticles, "T", datalength, ".RData")
save(algoparameters, RGestimates, nrep, file = filename)

algoparameters$nparticles <- 2^10
RGestimates <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  res <- rheeglynn_estimator(observations, model, theta, algoparameters)
  data.frame(irep = irep, estimate = res$estimate, iteration = res$iteration)
}
filename <- paste0("unlikely.RG.R", nrep, "N", algoparameters$nparticles, "T", datalength, ".RData")
save(algoparameters, RGestimates, nrep, file = filename)
# 





