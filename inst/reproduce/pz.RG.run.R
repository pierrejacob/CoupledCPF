# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(doRNG)
library(doParallel)
library(dplyr)

registerDoParallel(cores = detectCores() - 2)

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

nrep <- 100
seq_nparticles <- c(1024, 2048, 4096)

# nparticles <- 2^12
# algoparameters <- list(nparticles = nparticles, coupled_resampling = CR_indexmatching, lambda = 0, with_as = FALSE)
coupled_resampling <- CR_indexmatching
with_as <- FALSE
# nrep <- 1000
h <- function(x) x
meetings.df <- data.frame()

for (N in seq_nparticles){
  print(N)
  rinit <- function(){
    return(CPF(N, pz, theta, observations, ref_trajectory = NULL, with_as = with_as))
  }
  single_kernel_RB <- function(chain_state, h){
    return(CPF_RB(N, pz, theta, observations, ref_trajectory = chain_state, with_as = with_as, h = h))
  }
  coupled_kernel_RB <- function(chain_state1, chain_state2, h){
    return(CPF_coupled_RB(N, pz, theta, observations, ref_trajectory1 = chain_state1, ref_trajectory2 = chain_state2,
                          coupled_resampling = coupled_resampling, with_as = with_as, h = h))
  }
  meetings_ <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
    # res <- rheeglynn_estimator(observations, pz, theta, algoparameters)
    res <- unbiasedestimator_RB(single_kernel_RB, coupled_kernel_RB, rinit, h = h, k = 0, m = 0, max_iterations = 1e4)
    estimate <- res$uestimator
    if (is.null(dim(estimate))){
      estimate <- matrix(estimate, nrow = 1)
    }
    data.frame(irep = irep, iteration = res$iteration, meetingtime = res$meetingtime,
               nparticles = N, with_as = with_as)
  }
  meetings.df <- rbind(meetings.df, meetings_)
  save(meetings.df, file = "pz.meetingtime.RData")
}

load("pz.meetingtime.RData")
meetings.df %>% tail
g <- ggplot(meetings.df, aes(x = meetingtime, group = factor(nparticles), fill = factor(nparticles)))
g <- g + geom_histogram(aes(y = ..density..), position = position_dodge())
g
meetings.df %>% group_by(nparticles) %>% summarise(m = mean(meetingtime), s = sd(meetingtime), c = mean(nparticles * meetingtime))
ggplot(meetings.df %>% group_by(nparticles) %>% summarise(m = mean(meetingtime)),
       aes(x = nparticles, y = m)) + geom_point()

summ.df <- meetings.df %>% group_by(nparticles) %>% summarise(k = floor(quantile(meetingtime, probs = 0.9)))
summ.df$k

nrep <- 1000
estimates.df <- data.frame()
for (iN in 1:length(seq_nparticles)){
  N <- seq_nparticles[iN]
  k <- summ.df$k[iN]
  m <- 2*k
  print(N)
  rinit <- function(){
    return(CPF(N, pz, theta, observations, ref_trajectory = NULL, with_as = with_as))
  }
  single_kernel_RB <- function(chain_state, h){
    return(CPF_RB(N, pz, theta, observations, ref_trajectory = chain_state, with_as = with_as, h = h))
  }
  coupled_kernel_RB <- function(chain_state1, chain_state2, h){
    return(CPF_coupled_RB(N, pz, theta, observations, ref_trajectory1 = chain_state1, ref_trajectory2 = chain_state2,
                          coupled_resampling = coupled_resampling, with_as = with_as, h = h))
  }
  estimates_ <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
    res <- unbiasedestimator_RB(single_kernel_RB, coupled_kernel_RB, rinit, h = h, k = k, m = m, max_iterations = 1e4)
    estimate <- res$uestimator
    if (is.null(dim(estimate))){
      estimate <- matrix(estimate, nrow = 1)
    }
    data.frame(irep = irep, iteration = res$iteration, meetingtime = res$meetingtime,
               time = 0:datalength, estimate = t(estimate), k = k, m = m,
               nparticles = N, with_as = with_as)
  }
  estimates.df <- rbind(estimates.df, estimates_)
  save(observations, nrep, datalength, dimension, theta, estimates.df, file = "pz.estimate.RData")
}

load("pz.estimate.RData")
cost <- function(tau, m) 2 * tau + pmax(1, m + 1 - tau)
estimates.df <- estimates.df %>% mutate(c = cost(meetingtime, iteration) * nparticles)
estimates.df %>% head
estimates.df$irep %>% summary
summary.df <- estimates.df %>% group_by(nparticles, time) %>% summarise(m = mean(estimate.2), v = var(estimate.2), inef = mean(c) * var(estimate.1),
                                                                        meancost = mean(c))
summary.df %>% filter(time == 0) %>% tail
ggplot(summary.df, aes(x = time, y = v, colour = factor(nparticles))) + geom_point() + geom_line() + scale_y_log10()
ggplot(summary.df, aes(x = time, y = inef, colour = factor(nparticles))) + geom_point() + geom_line() + scale_y_log10()

