# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(doRNG)
library(ggthemes)
library(xtable)
ncores <- 10
registerDoMC(cores = ncores)
# load custom theme for the plots
# it requires the packages ggplot2, gridExtra, reshape
setmytheme()
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#0072B2", "#F0E442")

# fix the random seed
set.seed(17)
module_tree <<- Module("module_tree", PACKAGE = "CoupledCPF")
TreeClass <<- module_tree$Tree

load(file = "ar1.effecttruncation.RData")

iter.df <- estimates.df %>% group_by(irep, lambda) %>% summarise(iteration = mean(iteration), truncation = mean(truncation))

m <- iter.df %>% group_by(lambda) %>% summarise(m = mean(iteration), s = sd(iteration)) %>% 
  mutate(Meet = paste0(round(m,2), " (", round(s,2), ")")) %>% select(lambda, Meet)
m$lambda <- paste0("p =  ", m$lambda)
names(m) <- c("", "meeting time")
cap <- "Average meeting time, as a function of the Geometric parameter $p$.
Standard deviations are between brackets. Results obtained in the hidden auto-regressive model with $T = 500$.
\\label{table:effecttruncation}"
print.xtable(xtable(m, caption = cap), include.rownames =FALSE)
print.xtable(xtable(m, caption = cap), include.rownames =FALSE,  file = "ar1effecttruncation.tex")


RGSmean.df <- estimates.df %>% group_by(time, lambda) %>% summarise(m = mean(estimate), variance = var(estimate), cost = mean(iteration))

mean.df <- RGSmean.df
mean.df <- mean.df %>% mutate(efficiency = 1/(variance * cost))

seq_p <- c(0, 0.01, 0.025, 0.05)
label.df <- data.frame(time = c(50,50,50,200),
                       lambda = seq_p)
label.df <- merge(label.df, mean.df %>% select(time, lambda, variance, efficiency), by = c('time', 'lambda'))
mean.df$lambda <- factor(mean.df$lambda, levels = seq_p, labels = paste0("p = ", seq_p))
label.df$lambda <- factor(label.df$lambda, levels = seq_p, labels = paste0("p = ", seq_p))

gvar <- ggplot(mean.df, aes(x = time, y = variance, group = lambda, colour = lambda)) 
gvar <- gvar+ geom_line() + scale_color_manual(name = "p:", values = cbbPalette) + scale_y_log10()
gvar <- gvar + xlab("time") + ylab("variance")
gvar <- gvar + geom_label(data = label.df,
                    aes(x = time, y = variance, label = lambda, group = lambda, 
                        colour = lambda),
                    size = 6) + theme(legend.position = "none")
gvar

ggsave(filename = "ar1.effecttruncation.var.pdf", plot = gvar)

label.df[,4] <- c(label.df[1,4], 1e-4, 4.5e-05, 2e-5)
geff <- ggplot(mean.df, aes(x = time, y = efficiency, group = lambda, colour = lambda)) + geom_line() +
  scale_colour_manual(name = "p:", values = cbbPalette) + scale_y_log10()
geff <- geff + xlab("time") + ylab("efficiency")
geff <- geff + geom_label(data = label.df,
                          aes(x = time, y = efficiency, label = lambda, group = lambda, 
                              colour = lambda),
                          size = 6) + theme(legend.position = "none")

geff

ggsave(filename = "ar1.effecttruncation.eff.pdf", plot = geff)
