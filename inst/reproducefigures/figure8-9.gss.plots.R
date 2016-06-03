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
load("gss_unbiased.autostop.RData")
load("gss.PGAS.RData")

tail(PGAS.df)
g <- qplot(x = PGAS.df %>% filter(iteration > 25000, time == 36) %>% select(chain), geom = "blank") + geom_histogram(aes(y = ..density..), binwidth = 0.5)
g <- g + xlab(expression(x[36]))
g

ggsave(filename = "gssPGASx36.pdf", plot = g)

meanx36 <- PGAS.df %>% filter(iteration > 25000, time == 36) %>% summarise(m = mean(chain))

g <- ggplot(estimates %>% filter(time == 36), aes(x = estimate)) + geom_histogram(aes(y = ..density..), binwidth = 1)
g <- g + xlab(expression(hat(x)[36]))
g <- g + geom_vline(xintercept = as.numeric(meanx36[1]), colour = "red", linetype = 2, size = 3)
g
ggsave(filename = "gssRGx36.pdf", plot = g)


head(estimates)
cost <- estimates %>% group_by(irep) %>% summarise(iteration = mean(iteration))
summary(cost$iteration)
sum(cost$iteration)

estimates.df <- estimates %>% select(irep, time, starts_with("estimate")) %>% group_by(time) %>%
  summarise(m = mean(estimate), s = sd(estimate) / sqrt(nrep))


PGAS.mean.df <- PGAS.df %>% filter(iteration > 25000) %>% group_by(time) %>% summarise(m = mean(chain))
g <- ggplot(estimates.df, aes(x = time)) + geom_errorbar(aes(ymin = m - 2 * s, ymax = m + 2 * s)) + geom_point(aes(y = m))
g <- g + ylab("smoothing means")
g <- g + geom_line(data = PGAS.mean.df, aes(x = time, y = m, ymin = NULL, ymax = NULL, group = NULL), colour = "red", linetype = 1)
g <- g + xlab("time")
g

ggsave(filename = "gssunbiasedsmoothing.pdf", plot = g, height = 5, width = 12)
