# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(doRNG)
library(reshape2)
ncores <- 4
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

#  x_0 \sim N(0, tau_0^2), x_t = rho * x_t-1 + tau epsilon_t
# and y = x + sigma eta_t at the final step
# parameter = (rho, tau, sigma, tau_0)
theta <- c(0.9, 0.1, 0.1, 0.1)

joint_mean <- rep(0, datalength+2)
joint_covariance <- matrix(0, nrow = datalength+2, ncol = datalength+2)
joint_covariance[1,1] <- theta[4]^2
for (i in 2:(datalength+1)){
  joint_covariance[i,i] <- joint_covariance[i-1,i-1] * (theta[1]^2) + theta[2]^2
}
joint_covariance[datalength+2,datalength+2] <- joint_covariance[datalength+1,datalength+1] + theta[3]^2
joint_covariance[datalength+2, datalength+1] <- joint_covariance[datalength+1, datalength+2] <- joint_covariance[datalength+1,datalength+1]
for (i in 1:(datalength+1)){
  for (j in (i+1):(datalength+1)){
    joint_covariance[j, i] <- joint_covariance[i, j] <- theta[1]^(abs(j-i)) * joint_covariance[i,i]
  }
}

joint_covariance[datalength+2,1:(datalength+1)] <- joint_covariance[datalength+1,1:(datalength+1)]
joint_covariance[1:(datalength+1),datalength+2] <- joint_covariance[1:(datalength+1),datalength+1]

mu_1 <- matrix(joint_mean[1:(datalength+1)], ncol = 1)
mu_2 <- matrix(joint_mean[datalength+2], ncol = 1)
Sigma_11 <- matrix(joint_covariance[1:(datalength+1),1:(datalength+1)], ncol = datalength+1)
Sigma_21 <- matrix(joint_covariance[datalength+2,1:(datalength+1)], ncol = datalength+1)
Sigma_12 <- t(Sigma_21)
Sigma_11 <- matrix(joint_covariance[1:(datalength+1),1:(datalength+1)], ncol = datalength+1)
Sigma_22 <- matrix(joint_covariance[datalength+2,datalength+2], ncol = 1)

mu_bar <- mu_1 + Sigma_12 %*% solve(Sigma_22) %*% (matrix(1, 1, 1) - mu_2)
Sigma_bar <- Sigma_11 - Sigma_12 %*% solve(Sigma_22) %*% Sigma_21
mu_bar
sqrt(diag(Sigma_bar))
save(mu_bar, Sigma_bar, file = "unlikely.truth.RData")
