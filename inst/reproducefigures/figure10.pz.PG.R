# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(doRNG)
library(ggthemes)
ncores <- 2
registerDoMC(cores = ncores)
# load custom theme for the plots
# it requires the packages ggplot2, gridExtra, reshape
setmytheme()
# fix the random seed
set.seed(17)
module_tree <<- Module("module_tree", PACKAGE = "CoupledCPF")
TreeClass <<- module_tree$Tree

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

# 
nparticles <- 2^12
xref <- CPF(nparticles, pz, theta, observations)
niterations <- 50000
trajectories <- array(0, dim = c(niterations, dimension, datalength + 1))
for (i in 1:niterations){
  if (i %% 100 == 1) cat("iteration ", i, "/", niterations, "\n")
  xref <- CPF(nparticles, pz, theta, observations, xref)
  trajectories[i,,] <- xref
}

PG.df <- melt(trajectories)
names(PG.df) <- c("iteration", "component", "time", "chain")
PG.df$time <- PG.df$time - 1

PG.mean.df <- PG.df %>% filter(iteration > (niterations/2)) %>% group_by(time, component) %>% summarise(m = mean(chain), v = var(chain))

filename <- paste0("pz.PG.M", niterations, "N", nparticles, "T", datalength, ".RData")
save(observations, nparticles, niterations, dimension, PG.mean.df, file = filename)

