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
nrep <- 1000
lag <- 10
print(lag)
fixedlag.df <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  res <- fixedlagsmoother(nparticles, lag, pz, theta, observations)
  data.frame(irep = irep, time = 0:365, Z = res$smoothing_means[,2])
}
fixedlag.summary.df <- fixedlag.df %>% select(irep, time, Z) %>% group_by(time) %>%
  summarise(m.2 = mean(Z), s.2 = sd(Z) / sqrt(nrep))
fixedlag.variance.dfL10 <- fixedlag.summary.df %>% mutate(var = (s.2/m.2)^2)
fixedlag.variance.dfL10$method = "fixed lag L=10  "
# 
lag <- datalength
print(lag)
fixedlag.df <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  res <- fixedlagsmoother(nparticles, lag, pz, theta, observations)
  data.frame(irep = irep, time = 0:365, Z = res$smoothing_means[,2])
}
fixedlag.summary.df <- fixedlag.df %>% select(irep, time, Z) %>% group_by(time) %>%
  summarise(m.2 = mean(Z), s.2 = sd(Z) / sqrt(nrep))
fixedlag.variance.dfL365 <- fixedlag.summary.df %>% mutate(var = (s.2/m.2)^2)
fixedlag.variance.dfL365$method = "fixed lag L=365  "

save(fixedlag.variance.dfL10, fixedlag.variance.dfL365, 
     file = "pz.fxdlagR1000N4096T365.RData")
