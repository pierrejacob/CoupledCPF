# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(doParallel)
library(doRNG)
library(dplyr)
registerDoParallel(cores = detectCores() - 2)
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

#####
nrep <- 100
coupled_resampling <- CR_indexmatching
with_as <- FALSE
# seq_nparticles <- c(2^7, 2^8, 2^9, 2^10)
seq_nparticles <- c(2^7, 2^8, 2^9, 2^10)

meetings.df <- data.frame()

for (iN in 1:length(seq_nparticles)){
  N <- seq_nparticles[iN]
  print(N)
  rinit <- function(){
    return(CPF(N, model, theta, observations, ref_trajectory = NULL, with_as = with_as))
  }
  single_kernel_RB <- function(chain_state, h){
    return(CPF_RB(N, model, theta, observations, ref_trajectory = chain_state, with_as = with_as, h = h))
  }
  coupled_kernel_RB <- function(chain_state1, chain_state2, h){
    return(CPF_coupled_RB(N, model, theta, observations, ref_trajectory1 = chain_state1, ref_trajectory2 = chain_state2,
                          coupled_resampling = coupled_resampling, with_as = with_as, h = h))
  }

  h <- function(x) x
  max_niterations <- 10^4
  # res <- unbiasedestimator_RB(single_kernel_RB, coupled_kernel_RB, rinit, h = h, k = 0, m = 0, max_iterations = max_niterations)
  meetings <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
    res <- unbiasedestimator_RB(single_kernel_RB, coupled_kernel_RB, rinit, h = h, k = 0, m = 0, max_iterations = max_niterations)
    data.frame(irep = irep, iteration = res$iteration, meetingtime = res$meetingtime, N = N)
  }
  meetings.df <- rbind(meetings.df, meetings)
  save(meetings.df, file = paste0("unlikely.RG.AS", with_as, ".meetings.RData"))
}

load(file = paste0("unlikely.RG.AS", with_as, ".meetings.RData"))
meetings.df %>% head
ggplot(meetings.df, aes(x = N, y = iteration)) + geom_point()
summary.meetings.df <- meetings.df %>% group_by(N) %>% summarise(m = mean(iteration), q95 = as.numeric(quantile(iteration, probs = 0.95)))
summary.meetings.df
q95 <- floor(summary.meetings.df$q95)
meantau <- floor(summary.meetings.df$m)

nrep <- 10000
estimates.df <- data.frame()
for (iN in 1:length(seq_nparticles)){
  N <- seq_nparticles[iN]
  k <- meantau[iN]
  m <- 1 * k
  print(N)
  rinit <- function(){
    return(CPF(N, model, theta, observations, ref_trajectory = NULL, with_as = with_as))
  }
  single_kernel_RB <- function(chain_state, h){
    return(CPF_RB(N, model, theta, observations, ref_trajectory = chain_state, with_as = with_as, h = h))
  }
  coupled_kernel_RB <- function(chain_state1, chain_state2, h){
    return(CPF_coupled_RB(N, model, theta, observations, ref_trajectory1 = chain_state1, ref_trajectory2 = chain_state2, 
                          coupled_resampling = coupled_resampling, with_as = with_as, h = h))
  }
  h <- function(x) x
  max_niterations <- 10^4
  # res <- unbiasedestimator_RB(single_kernel_RB, coupled_kernel_RB, rinit, h = h, k = 0, m = 0, max_iterations = max_niterations)
  estimates <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
    res <- unbiasedestimator_RB(single_kernel_RB, coupled_kernel_RB, rinit, h = h, k = k, m = m, max_iterations = max_niterations)
    estimate <- res$uestimator
    if (is.null(dim(estimate))){
      estimate <- matrix(estimate, nrow = 1)
    }
    data.frame(irep = irep, iteration = res$iteration, meetingtime = res$meetingtime, time = 0:datalength, estimate = t(estimate), N = N,
               k = k, m = m)
  }
  estimates.df <- rbind(estimates.df, estimates)
  filename <- paste0("unlikely.RG.AS", with_as, ".estimates.RData")
  save(estimates.df, file = filename)
}

estimates.df %>% head
