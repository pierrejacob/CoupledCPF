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

# generate observations
alpha_star <- 0.95
A_star <- create_A(alpha_star, dimension)
datalength <- 10000
observations <- matrix(nrow = datalength, ncol = dimension)
x_t <- fast_rmvnorm(1, rep(0, dimension), diag(1, nrow = dimension, ncol = dimension))
for (time in 1:datalength){
  x_t <- t(A_star %*% t(x_t)) + fast_rmvnorm(1, rep(0, dimension), diag(1, nrow = dimension, ncol = dimension))
  observations[time,] <- x_t + fast_rmvnorm(1, rep(0, dimension), diag(1, nrow = dimension, ncol = dimension))
}

# plot(observations, type = "l")
save(observations, alpha_star, datalength, 
     file = "ar1data.RData")
