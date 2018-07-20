# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(ggthemes)
library(dplyr)
library(reshape2)
# load custom theme for the plots
# it requires the packages ggplot2, gridExtra, reshape
setmytheme()
# fix the random seed
set.seed(17)
load("ar1data.RData")
datalength <- 100
observations <- matrix(observations[1:datalength,], ncol = 1)
## get exact smoothing means
## smoothing to exact means
phi_star <- 0.9
sigma_w_star <- 1
sigma_v_star <- 1
sigma_v_star2 <- sigma_v_star^2
# number of observations
### Kalman filter
T <- datalength
Y <- observations[,1]
kf <- function(phi, sigmaw){
  sigmaw2 <- sigmaw^2
  kf_means <- rep(0, T+1)
  kf_variances <- rep(0, T+1)
  m_current <- 0
  kf_means[1] <- m_current
  V_current <- 1
  kf_variances[1] <- V_current
  loglikelihood <- 0
  # also store predictive means and variances
  kf_predmeans <- rep(0, T)
  kf_predvariances <- rep(0, T)
  for (t in 1:T){
    # prediction step
    m_next <- phi * m_current
    V_next <- phi * V_current * phi + sigmaw2
    loglikelihood <- loglikelihood + dnorm(Y[t], mean = m_next, sd = sqrt(V_next + sigma_v_star2), log = TRUE)
    kf_predmeans[t] <- m_next
    kf_predvariances[t] <- V_next
    # update step
    K <- V_next  /  (V_next + sigma_v_star2)
    m_current <- m_next + K * (Y[t] - m_next)
    V_current <- (1 - K) * V_next
    kf_means[t+1] <- m_current
    kf_variances[t+1] <- V_current
  }
  return(list(kf_means = kf_means, kf_variances = kf_variances, loglikelihood = loglikelihood,
              kf_predmeans = kf_predmeans, kf_predvariances = kf_predvariances))
}

kf_results <- kf(phi_star, sigma_w_star)
### Kalman smoother
ks <- function(kf_results, phi, sigmaw){
  kf_means <- kf_results$kf_means
  kf_variances <- kf_results$kf_variances
  kf_predmeans <- kf_results$kf_predmeans
  kf_predvariances <- kf_results$kf_predvariances
  # smoothing backward pass
  ks_means <- rep(0, T+1)
  ks_means[T+1] <- kf_means[T+1]
  ks_variances <- rep(0, T+1)
  ks_variances[T+1] <- kf_variances[T+1]
  for (t in (T-1):0){
    L <- kf_variances[t+1] * phi * 1/(kf_predvariances[t+1])
    m <- kf_means[t+1] + L * (ks_means[t+2] - kf_predmeans[t+1])
    V <- kf_variances[t+1] + L^2 * (ks_variances[t+2] - kf_predvariances[t+1])
    ks_means[t+1] <- m
    ks_variances[t+1] <- V
  }
  return(list(ks_means = ks_means, ks_variances = ks_variances))
}

ks_results <- ks(kf_results, phi_star, sigma_w_star)
ksmeans.df <- data.frame(time = 0:datalength, truemeans = ks_results$ks_means)

# First, meeting times for N=128 particles
load("ar.kandm.RData")
estimates.df %>% nrow
g <- ggplot(estimates.df, aes(x = meetingtime))
g <- g + geom_histogram(aes(y = ..density..), bins = 15) + xlab(expression(tau))
g
# ggsave(filename = "ar1.meetingtime.pdf", plot = g, width = 5, height = 4)

cost <- function(tau, m) 3 + 2 * (tau-1) + pmax(0, m - tau)
estimates.df <- estimates.df %>% mutate(c = cost(meetingtime, iteration))
df <- estimates.df %>% group_by(time,k,m) %>% summarise(meanestim = mean(estimate), varestim = var(estimate), meancost = mean(c), varcost = var(c), nrep = n())
df <- df %>% mutate(inef = varestim * meancost)
df <- df %>% mutate(ymin = meanestim - 2 * sqrt(varestim/nrep), ymax = meanestim + 2 * sqrt(varestim/nrep))

# For k = 10, m = 20, the average cost is 
df %>% filter(k == 10, m == 20) %>% group_by(k,m) %>% summarise(cost = mean(meancost), sdcost = sqrt(mean(varcost)))
# the cost is 28.151... and the standard deviation 3.7, which is about 13%
df %>% head

df.100 <- estimates.df %>% filter(irep <= 100) %>% group_by(time,k,m) %>% summarise(meanestim = mean(estimate), varestim = var(estimate), meancost = mean(c), varcost = var(c), nrep = n())
df.100 <- df.100 %>% mutate(inef = varestim * meancost)
df.100 <- df.100 %>% mutate(ymin = meanestim - 2 * sqrt(varestim/nrep), ymax = meanestim + 2 * sqrt(varestim/nrep))

g <- ggplot(df.100 %>% filter(time >=0, time <= 100, k == 10, m == 20), aes(x = time))
g <- g + geom_errorbar(aes(ymin = ymin, ymax = ymax), position = position_dodge())
g <- g + geom_line(data=ksmeans.df %>% filter(time >=0, time <= 100), aes(x = time, y = truemeans, group = NULL, colour = NULL), colour = "red", alpha = 0.5)
g <- g + xlab("time") + ylab("smoothing means")
g
# ggsave(filename = "ar1.smoothingmeans.pdf", plot = g, width = 10, height = 4)


## Meeting times for different numbers of particles, fixed time horizon
load("ar1.meetingtime.RData")
meetings.df %>% group_by(nparticles, with_as) %>% summarise(meanmeeting = mean(meetingtime), nrep = n()) %>% mutate(costmeeting = nparticles * meanmeeting)
# from there we extract the mean meeting times, and the mean cost (last column)

