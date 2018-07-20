# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
# fix the random seed
set.seed(17)
# Below is a very obfuscated way of generating data from a hidden AR(1) model
# get AR(1) model 
dimension <- 1
ar1 <- get_ar(dimension)

# generate observations
# data-generating parameter
alpha_star <- 0.95
A_star <- create_A(alpha_star, dimension)
# data length T
datalength <- 10000
# matrix storing observations T x d
observations <- matrix(nrow = datalength, ncol = dimension)
# latent process initialization
x_t <- fast_rmvnorm(1, rep(0, dimension), diag(1, nrow = dimension, ncol = dimension))
for (time in 1:datalength){
  # latent process transition
  x_t <- t(A_star %*% t(x_t)) + fast_rmvnorm(1, rep(0, dimension), diag(1, nrow = dimension, ncol = dimension))
  # generate data given current latent state
  observations[time,] <- x_t + fast_rmvnorm(1, rep(0, dimension), diag(1, nrow = dimension, ncol = dimension))
}

# plot(observations, type = "l")
# store data in a file in the current working directory
save(observations, alpha_star, datalength, file = "ar1data.RData")



