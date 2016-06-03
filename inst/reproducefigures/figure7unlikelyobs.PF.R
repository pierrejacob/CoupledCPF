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
nrep <- 10000

#####
nparticles <- 2^10
PFestimates <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  res <- particlefilter(nparticles, model, theta, observations)
  x <- res$trajectories[1,,]
  estimate <- apply(X = x, MARGIN = 1, function(v) sum(v * res$weights))
  data.frame(irep = irep, estimate = t(estimate))
}
filename <- paste0("unlikely.PF.R", nrep, "N", nparticles, "T", datalength, ".RData")
save(nparticles, PFestimates, nrep, file = filename)

#####
nparticles <- 2^11
PFestimates <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  res <- particlefilter(nparticles, model, theta, observations)
  x <- res$trajectories[1,,]
  estimate <- apply(X = x, MARGIN = 1, function(v) sum(v * res$weights))
  data.frame(irep = irep, estimate = t(estimate))
}
filename <- paste0("unlikely.PF.R", nrep, "N", nparticles, "T", datalength, ".RData")
save(nparticles, PFestimates, nrep, file = filename)

#####
nparticles <- 2^12
PFestimates <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  res <- particlefilter(nparticles, model, theta, observations)
  x <- res$trajectories[1,,]
  estimate <- apply(X = x, MARGIN = 1, function(v) sum(v * res$weights))
  data.frame(irep = irep, estimate = t(estimate))
}
filename <- paste0("unlikely.PF.R", nrep, "N", nparticles, "T", datalength, ".RData")
save(nparticles, PFestimates, nrep, file = filename)


# #####
nparticles <- 2^13
PFestimates <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  res <- particlefilter(nparticles, model, theta, observations)
  x <- res$trajectories[1,,]
  estimate <- apply(X = x, MARGIN = 1, function(v) sum(v * res$weights))
  data.frame(irep = irep, estimate = t(estimate))
}
filename <- paste0("unlikely.PF.R", nrep, "N", nparticles, "T", datalength, ".RData")
save(nparticles, PFestimates, nrep, file = filename)
#
#####
nparticles <- 2^14
PFestimates <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  res <- particlefilter(nparticles, model, theta, observations)
  x <- res$trajectories[1,,]
  estimate <- apply(X = x, MARGIN = 1, function(v) sum(v * res$weights))
  data.frame(irep = irep, estimate = t(estimate))
}
filename <- paste0("unlikely.PF.R", nrep, "N", nparticles, "T", datalength, ".RData")
save(nparticles, PFestimates, nrep, file = filename)


