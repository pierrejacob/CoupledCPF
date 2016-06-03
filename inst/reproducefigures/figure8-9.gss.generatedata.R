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

gss <- get_gss()
theta <- c()
datalength <- 10000
x <- matrix(0, nrow = datalength + 1, ncol = 1)
observations <- matrix(0, nrow = datalength, ncol = 1)
x[1,] <- 0.1
for (i in 1:datalength){
  x[i+1,] <- 0.5 * x[i,] + 25 * x[i,] / (1 + x[i,]^2) + 8 * cos(1.2*(i - 1)) + sqrt(10) * rnorm(1)
  observations[i,] <- (x[i+1,]^2) / 20 + rnorm(1)
}

# plot(observations[1:100,], type = "l")

save(observations, datalength, file = "gssdata.RData")