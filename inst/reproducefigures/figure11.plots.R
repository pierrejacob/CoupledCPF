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
load(file = "pz.fxdlagR1000N4096T365.RData")
tail(fixedlag.variance.dfL10)
fixedlag.variance.dfL10$method <- "fixed-lag L=10"
# fixedlag.variance.dfL100$method <- "fixed-lag L=100"
fixedlag.variance.dfL365$method <- "particle filter"

load("pz.RB2.R1000N4096T365.RData")
# summary(pmax(4, estimates$iteration))
summary(estimates$iteration)
summary(estimates$meetingtime)
estimates %>% tail
RB2estimates.df <- estimates %>% filter(irep < nrep) %>% select(irep, time, starts_with("estimate")) %>% group_by(time) %>%
  summarise(m.2 = mean(estimate.2), s.2 = sd(estimate.2) / sqrt(nrep))

RB2estimator.variance.df <- RB2estimates.df %>% mutate(var = (s.2/m.2)^2) %>% select(time, var)
RB2estimator.variance.df$method = "unbiased+RB+m"

load("pz.RB.R1000N4096T365.RData")
# summary(estimates$iteration)
estimates %>% tail
RBestimates.df <- estimates %>% filter(irep < nrep) %>% select(irep, time, starts_with("estimate")) %>% group_by(time) %>%
  summarise(m.2 = mean(estimate.2), s.2 = sd(estimate.2) / sqrt(nrep))

RBestimator.variance.df <- RBestimates.df %>% mutate(var = (s.2/m.2)^2) %>% select(time, var)
RBestimator.variance.df$method = "unbiased+RB"

# estimates.df %>% tail
# RBestimates.df %>% tail
# RB2estimates.df %>% tail

load("pz.unbiased.R1000N4096T365.RData")

estimates.df <- estimates %>% filter(irep < nrep) %>% select(irep, time, starts_with("estimate")) %>% group_by(time) %>%
  summarise(m.2 = mean(estimate.2), s.2 = sd(estimate.2) / sqrt(nrep))
summary(estimates$iteration)
estimator.variance.df <- estimates.df %>% mutate(var = (s.2/m.2)^2) %>% select(time, var)
estimator.variance.df$method = "unbiased"

# v.df <- rbind(fixedlag.variance.dfL10, fixedlag.variance.dfL100, fixedlag.variance.dfL365,
              # estimator.variance.df, RBestimator.variance.df)
v.df <- rbind(fixedlag.variance.dfL10 %>% select(time, var,method), 
              fixedlag.variance.dfL365 %>% select(time, var,method),
              estimator.variance.df, RBestimator.variance.df, RB2estimator.variance.df)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#0072B2", "#F0E442")

g <- ggplot(v.df, aes(x = time, y = var * 1000, group = method, colour = method)) + geom_line() + scale_y_log10()
g <- g + scale_x_continuous(breaks = c(0, 91, 182, 274, 365))
g <- g + scale_colour_manual(name = "", values = cbbPalette) + ylab("relative variance")
label.df <- data.frame(time = c(70,70,70,70, 70),
                       y = 1000*c(10^{-7.8}, 10^{-6},10^{-3.5},10^{-4.6}, 10^{-5.5}),
                       method = c("fixed-lag L=10", "particle filter", "unbiased", "unbiased+RB", "unbiased+RB+m"))
g <- g + geom_label(data = label.df, aes(x = time, y = y, colour = method, label = method), size = 6) + theme(legend.position = "none")
g

ylabels <- (v.df %>% filter(time == 1))$var
ylabels[3] <- 5*1e-3
ylabels[4] <- 1.5*1e-3
g <- ggplot(v.df, aes(x = time, y = var * 1000, group = method, colour = method)) + geom_line() + scale_y_log10()
g <- g + scale_x_continuous(breaks = c(0, 91, 182, 274, 365), limits = c(-60, 365))
g <- g + scale_colour_manual(name = "", values = cbbPalette) + ylab("relative variance")
label.df <- data.frame(time = c(-30,-20,-15,-25,-30),
                       y = 1000*ylabels,
                       method = c("fixed-lag L=10", "particle filter", "unbiased", "unbiased+RB", "unbiased+RB+m"))
g <- g + geom_label(data = label.df, aes(x = time, y = y, colour = method, label = method), size = 6) + theme(legend.position = "none")
g



ggsave(filename = "pzfixedlag.pdf", plot = g, height = 5, width = 15)

