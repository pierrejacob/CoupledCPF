# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(doRNG)
library(ggthemes)

ncores <- 8
registerDoMC(cores = ncores)
# load custom theme for the plots
# it requires the packages ggplot2, gridExtra, reshape
setmytheme()
# fix the random seed
set.seed(17)

load("pz.unbiased.R1000N4096T365.RData")
load("pz.PG.M50000N4096T365.RData")

head(estimates)

cost <- estimates %>% group_by(irep) %>% summarise(iteration = mean(iteration))
summary(cost$iteration)
sum(cost$iteration)

estimates.df <- estimates %>% select(irep, time, starts_with("estimate")) %>% group_by(time) %>%
  summarise(m.1 = mean(estimate.1), s.1 = sd(estimate.1) / sqrt(nrep),
            m.2 = mean(estimate.2), s.2 = sd(estimate.2) / sqrt(nrep))

g <- ggplot(estimates.df %>% filter(time < 30), aes(x = time)) + geom_errorbar(aes(ymin = m.2 - 2 * s.2, ymax = m.2 + 2 * s.2))
g <- g + ylab("smoothing means of Z") + ylim(0,10)
g <- g + geom_line(data = PG.mean.df %>% filter(component == 2, time < 30), aes(x = time, y = m, ymin = NULL, ymax = NULL, group = NULL), colour = "red")
g

# ggsave(filename = "pzsmoothingZbeginning.pdf", plot = g, height = 5, width = 8)

g <- ggplot(estimates.df %>% filter(time > 335), aes(x = time)) + geom_errorbar(aes(ymin = m.2 - 2 * s.2, ymax = m.2 + 2 * s.2))
g <- g + ylab("smoothing means of Z") + ylim(0,10)
g <- g + geom_line(data = PG.mean.df %>% filter(component == 2, time > 335), aes(x = time, y = m, ymin = NULL, ymax = NULL, group = NULL), colour = "red")
g

# ggsave(filename = "pzsmoothingZend.pdf", plot = g, height = 5, width = 8)

