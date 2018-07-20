# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(dplyr)
library(reshape2)
library(ggthemes)
library(latex2exp)
# load custom theme for the plots
# it requires the packages ggplot2, gridExtra, reshape
setmytheme()
# fix the random seed
set.seed(17)
module_tree <<- Module("module_tree", PACKAGE = "CoupledCPF")
TreeClass <<- module_tree$Tree

dimension <- 1
datalength <- 10

theta <- c(0.9, 0.1, 0.1, 0.1)
load("unlikely.truth.RData")
truemeans <- data.frame(time = 0:datalength, truemean = mu_bar)
#####

load("unlikely.RG.ASFALSE.estimates.RData")
estimates.df %>% head

cost <- function(tau, m) 2 * tau + pmax(1, m + 1 - tau)
estimates.df <- merge(truemeans, estimates.df, by = "time")
summary.df <- estimates.df %>% group_by(k, m, N, time) %>% summarise(truemean = mean(truemean), meanestimate = mean(estimate), v = var(estimate),
                                                               iteration = mean(iteration), 
                                                               meetingtime = mean(meetingtime), 
                                                               c = mean(N * cost(meetingtime, iteration)), nrep = n())
summary.df <- summary.df %>% mutate(ymin = meanestimate - 2*sqrt(v/nrep), ymax = meanestimate + 2*sqrt(v/nrep), inefficiency = c * v)
summary.df <- summary.df %>% mutate(bias2 = (meanestimate - truemean)^2, mse = bias2 + v)
summary.df %>% head

# cost?
summary.df %>% group_by(N, k, m) %>% summarise(meancost = mean(c), meanmeeting = mean(meetingtime), meaniteration = mean(iteration))

timeindex <- 9
summary.df$Nfactor <- factor(summary.df$N, levels = c(128, 256, 512, 1024), 
                             labels = c("N = 128", "N = 256", "N = 512", "N = 1024"))
g <- ggplot(summary.df %>% filter(time == timeindex), aes(x = Nfactor, group = N, ymin = ymin, ymax = ymax)) + geom_errorbar()
g <- g + geom_hline(yintercept = truemeans$truemean[timeindex+1], linetype = 2)
g <- g + coord_flip()
g <- g + ylab(TeX('estimate')) + xlab("")
g

# ggsave(plot = g, filename = "unlikely.rg.variousN.pdf", width = 7, height = 5)


#####
# Particle filter estimates
Ns <- c(2^8, 2^10, 2^12, 2^14)
PFresults.df <- data.frame()
for (N in Ns){
  load(paste0("unlikely.PF.R10000N", N, "T10.RData"))
  nrep <- max(PFestimates$irep)
  meanPF.df <- melt(PFestimates, id = "irep")
  names(meanPF.df) <- c("irep", "time", "estimate")
  meanPF.df$time <- rep(0:datalength, each = nrep)
  meanPF.df <- merge(meanPF.df, truemeans, by = "time")
  meanPF.df$N <- N
  PFresults.df <- rbind(PFresults.df, meanPF.df)
}
PFresults.df$N <- factor(PFresults.df$N, levels = Ns, labels = paste0(Ns, "  "))
PFresults.df %>% head

timeindex <- 9
g <- ggplot(PFresults.df %>% filter(time == timeindex), aes(x = estimate, group = N)) + geom_density(aes(y = ..density..), alpha = 0.5) +
  geom_vline(data=NULL, xintercept = truemeans$truemean[timeindex+1], linetype = 2) + scale_fill_colorblind()
g <- g + xlab("estimate") + theme(legend.position = "none")
g

label.df <- data.frame(x = c(0.3,0.37,0.44,0.51), y = c(3,4,5,6), txt = c("N=256", "N=1024", "N=4096", "N=16384"))
g <- g + geom_label(data=label.df, aes(x = x, y = y, label = txt, colour = NULL, group = NULL, fill = NULL),
                    size = 7)
g
# ggsave(plot = g, filename = "unlikely.pf.variousN.pdf", width = 7, height = 5)



## backup
# PFresults.df %>% filter(time == timeindex) %>% mutate(sqerror = (estimate - truemean)^2) %>%
#   group_by(N) %>% summarize(mse = mean(sqerror), v = var(estimate), bias2 = (mean(estimate) - mean(truemean))^2)
# 
# # For N = 16,384 we would get
# # v / R + bias2
# 0.003993505 / nrep + 0.0003829163
# 
# summary.df %>% filter(time == timeindex) %>% select(N, v)
# # For N = 1024 we would get
# 1.212709 / nrep
