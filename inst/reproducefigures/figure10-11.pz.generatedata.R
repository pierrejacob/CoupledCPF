# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(doRNG)
library(ggthemes)
ncores <- 10
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

# theta = (mu_alpha, sd_alpha, c, e, ml, mq)
theta_dgp <- c(0.7, 0.5, 0.25, 0.3, 0.1, 0.1)
# generate dataset
datalength <- 10000
state <- matrix(exp(log(2) + rnorm(2, mean = 0, sd = 1)), nrow = 2)
states <- matrix(nrow = 2, ncol = datalength+1)
states[,1] <- state
log_obs <- rep(0, datalength)
for (t in 1:datalength){
  alpha <- rnorm(n = 1, mean = theta_dgp[1], sd = theta_dgp[2])
  state <- pz_transition(state, alpha, t-1, theta_dgp[3:6])
  states[,t+1] <- state
  log_obs[t] <- rnorm(1, mean = log(state[1,1]), sd = 0.2)
}

pz <- get_pz()
dimension <- 2
observations <- matrix(log_obs, ncol = 1)
# matplot(t(states), type = "l")

save(observations, datalength, file = "pzdata.RData")
