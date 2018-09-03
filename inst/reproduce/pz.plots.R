# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(ggthemes)
library(dplyr)
library(reshape2)
library(latex2exp)

# load custom theme for the plots
# it requires the packages ggplot2, gridExtra, reshape
setmytheme()
# fix the random seed
set.seed(17)

load("pzdata.RData")

theta_dgp <- c(0.7, 0.5, 0.25, 0.3, 0.1, 0.1)
datalength <- 365

### load smoothing means obtained from long run of PG 
load("pz.PG.M1e+05N4096T365.RData")
### obtained from long run of Particle Gibbs and after discarding some burn-in 
truemeanZ.df <- PG.mean.df %>% filter(component == 2) %>% select(-component,-v) %>% rename(truemeanZ = m)

### load results from fixed lag smoother
nparticles <- 2^16
nrep <- 1000
load(file = paste0("pz.fxdlagR", nrep, "N", nparticles, "T365.RData"))
# 
fixedlag.summary.dfL10 <- fixedlag.dfL10 %>% group_by(time) %>% summarize(meanZ = mean(Z), varZ = var(Z), nrep = n())
fixedlag.summary.dfL10$method <- "fixed-lag"
fixedlag.summary.dfL10 <- merge(fixedlag.summary.dfL10, truemeanZ.df, by = "time")
fixedlag.summary.dfL10 <- fixedlag.summary.dfL10 %>% mutate(relvarZ = varZ / (truemeanZ^2), bias2 = (meanZ-truemeanZ)^2)
head(fixedlag.summary.dfL10)

fixedlag.summary.dfL365 <- fixedlag.dfL365 %>% group_by(time) %>% summarize(meanZ = mean(Z), varZ = var(Z), nrep = n())
fixedlag.summary.dfL365$method <- "particle filter"
fixedlag.summary.dfL365 <- merge(fixedlag.summary.dfL365, truemeanZ.df, by = "time")
fixedlag.summary.dfL365 <- fixedlag.summary.dfL365 %>% mutate(relvarZ = varZ / (truemeanZ^2), bias2 = (meanZ-truemeanZ)^2)
head(fixedlag.summary.dfL365)

### Check that we get the smoothing means right 
# g <- ggplot(fixedlag.summary.dfL10 %>% filter(time < 20), aes(x = time, y = meanZ)) + geom_line() +
#   geom_line(aes(y = truemeanZ), colour = "blue", linetype = 2, alpha = 0.5, size = 3)
# g

### Now load RG estimators, corresponding to different values of N
load("pz.estimate.RData")
estimates.df %>% head
# Histogram of the meeting times, when using 4096 particles
g <- ggplot(estimates.df %>% filter(nparticles == 4096, time == 0), aes(x = meetingtime)) +
  geom_histogram(aes(y = ..density..), bins = 10) + 
  scale_x_continuous(breaks = c(1,5,10,15)) + xlab(expression(tau))
g
ggsave(filename = "pz.meetingtime.pdf", plot = g, width = 5, height = 4)

# compute cost 
cost <- function(tau, m) 2 * tau + pmax(1, m + 1 - tau)
estimates.df <- estimates.df %>% mutate(c = cost(meetingtime, iteration) * nparticles)
estimates.df %>% head
RG.summary.df <- estimates.df %>% group_by(nparticles, time,k,m) %>% summarise(meanZ = mean(estimate.2), varZ = var(estimate.2), inef = mean(c) * var(estimate.2),
                                                                        meancost = mean(c), nrep = n())
RG.summary.df$method <- paste0("RG N=", RG.summary.df$nparticles)
RG.summary.df <- merge(RG.summary.df, truemeanZ.df, by = "time")
RG.summary.df <- RG.summary.df %>% mutate(relvarZ = varZ / (truemeanZ^2), bias2 = (meanZ-truemeanZ)^2)### Cost ?
RG.summary.df %>% head
### Cost?
RG.summary.df %>% filter(time == 0) %>% select(nparticles,k,m,meancost)
### Check that we get the right smoothing means
# g <- ggplot(RG.summary.df %>% filter(time < 20), aes(x = time, y = meanZ, group = nparticles, colour = factor(nparticles))) + geom_line() +
#   geom_line(aes(y = truemeanZ, group = NULL, colour = NULL), colour = "blue", linetype = 2, alpha = 0.5)
# g

# ### Now compare relative variance for each smoothing mean
v.df <- rbind(fixedlag.summary.dfL10 %>% select(time, varZ, relvarZ, bias2, method),
              fixedlag.summary.dfL365 %>% select(time, varZ, relvarZ, bias2, method),
              RG.summary.df %>% ungroup() %>% filter(method == "RG N=1024") %>% select(time, varZ, relvarZ, bias2, method),
              RG.summary.df %>% ungroup() %>% filter(method == "RG N=2048") %>% select(time, varZ, relvarZ, bias2, method),
              RG.summary.df %>% ungroup() %>% filter(method == "RG N=4096") %>% select(time, varZ, relvarZ, bias2, method))

# Focusing on N = 4096 for the RG estimators, for clarity
g <- ggplot(v.df %>% filter(method != "RG N=1024", method != "RG N=2048"), aes(x = time, y = relvarZ, group = method)) + geom_line() + scale_y_log10()
g <- g + scale_x_continuous(breaks = c(0, 91, 182, 274, 365), limits = c(-60, 365))
g <- g + ylab("relative variance")
label.df <- data.frame(time = c(-35,-35,-40),
                       y = c(10^{-3.5}, 10^{-2.5},10^{-1.5}),
                       method = c("fixed-lag", "particle filter", "unbiased"))
g <- g + geom_label(data = label.df, aes(x = time, y = y, label = method), size = 6) + theme(legend.position = "none")
g
ggsave(filename = "pz.relativevariance.pdf", plot = g, width = 10, height = 4)

