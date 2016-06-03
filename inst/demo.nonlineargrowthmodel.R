#### This scripts briefly shows how this package works
# First, some initialization
# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(doRNG)
ncores <- 10
registerDoMC(cores = ncores)
setmytheme()
set.seed(17)
module_tree <<- Module("module_tree", PACKAGE = "CoupledCPF")
TreeClass <<- module_tree$Tree

###
# We are going to work with the nonlinear growth model
# first we import it
gss <- get_gss()
# let us generate some data from the model
dimension <- 1
theta <- c()
datalength <- 50
x <- matrix(0, nrow = datalength + 1, ncol = 1)
observations <- matrix(0, nrow = datalength, ncol = 1)
x[1,] <- gss$rinit(1, theta, gss$rinit_rand(1, theta))
# equivalent to
#x[1,] <- rnorm(1, 0, sd = sqrt(2))
for (i in 1:datalength){
  x[i+1,] <- gss$rtransition(x[i,], theta, i, gss$rtransition_rand(1, theta))
  # equivalent to 
  # x[i+1,] <- 0.5 * x[i,] + 25 * x[i,] / (1 + x[i,]^2) + 8 * cos(1.2*(i - 1)) + sqrt(10) * rnorm(1)
  observations[i,] <- (x[i+1,]^2) / 20 + rnorm(1)
}
# plot the latent process
plot(0:datalength, x, type = "l")
# plot the observations
plot(1:datalength, observations, type = "l")
#

# now let us estimate the smoothing mean with the proposed method
# choose a number of particles
nparticles <- 2^10
# choose a coupled resampling scheme 
coupled_resampling <- CR_indexmatching
# ancestor sampling or not
with_as <- TRUE
# truncation variable parameter between 0 and 1 (0 means that we don't use any truncation variable)
lambda <- 0
# put all the parameters in a list
algorithmic_parameters <- list(nparticles = nparticles, coupled_resampling = coupled_resampling, lambda = lambda, with_as = with_as)
# how many estimators do we want?
R <- 100
# the following should take a few seconds if you have 10 cores 
estimates <- foreach(r = 1:R, .combine = rbind) %dorng% {
  res <- rheeglynn_estimator(observations, gss, theta, algorithmic_parameters)
  data.frame(irep = r, time = 0:datalength, estimate = t(res$estimate), iteration = res$iteration)
}
#
# how many iterations did the estimators take?
summary(estimates$iteration)
# now plot the confidence intervals on the smoothing mean
estimates.df <- estimates %>% select(irep, time, starts_with("estimate")) %>% group_by(time) %>%
  summarise(m = mean(estimate), s = sd(estimate) / sqrt(R))


g <- ggplot(estimates.df, aes(x = time)) + geom_errorbar(aes(ymin = m - 2 * s, ymax = m + 2 * s)) + geom_point(aes(y = m))
g <- g + geom_line(aes(y = m))
g <- g + ylab("smoothing means")
g <- g + xlab("time")
print(g)

# we can compare to the result of a fixed lag smoother with 2^14 particles and a lag of 20
lag <- 20
nparticles <- 2^14
fixedlagestimates <- fixedlagsmoother(nparticles, lag, gss, theta, observations)
g <- g + geom_line(data = NULL, aes(x = 0:datalength, y = fixedlagestimates$smoothing_means, ymin = NULL, ymax = NULL), colour = "red", linetype = 2)
print(g)
# we see an agreement between the two procedures.




