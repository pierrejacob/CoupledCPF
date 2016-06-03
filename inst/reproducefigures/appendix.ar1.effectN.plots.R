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
# fix the random seed
set.seed(17)
module_tree <<- Module("module_tree", PACKAGE = "CoupledCPF")
TreeClass <<- module_tree$Tree

load(file = "ar1.effectnparticles.RData")

estimates.df <- rbind(estimates_was.df, estimates_woas.df)


times.df <- estimates.df %>% filter(with_as == FALSE) %>% group_by(irep, nparticles) %>% summarise(iteration = mean(iteration))
times.df <- times.df %>% mutate(cost = nparticles * iteration) %>% select(irep, nparticles, iteration, cost)
m <- times.df %>% group_by(nparticles) %>%  summarise(meancost = mean(cost), sdcost = sd(cost), mean_iter = mean(iteration),
                                                      sd_iter = sd(iteration))

m <- m %>% mutate(cost = paste0(round(meancost,0), " (", round(sdcost,0), ")"),
                  iter = paste0(round(mean_iter,2), " (", round(sd_iter,2), ")")) %>% select(nparticles, cost, iter) %>% ungroup()
tm <- t(m)
colnames(tm) <- paste0("N =  ", tm[1,])
tm <- tm[2:3,]
rownames(tm) <- c("cost", "meeting time")

xtable(t(tm))
cap <- "Average cost and meeting time, as a function of the number of particles $N$.
Standard deviations are between brackets. Results obtained in the hidden auto-regressive model with $T = 500$.
\\label{table:effectnparticles}"

print.xtable(xtable(t(tm), caption = cap), file = "ar1effectnparticles.tex")

RGSmean.df <- estimates.df %>% filter(with_as == FALSE) %>% group_by(time, with_as, nparticles) %>%
  summarise(m = mean(estimate), v= var(estimate), iteration = mean(iteration))
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#0072B2", "#F0E442")

#
mean.df <- RGSmean.df
mean.df <- mean.df %>% mutate(efficiency = 1 / (v * iteration * nparticles)) %>% ungroup() %>% select(time, nparticles, v, efficiency)
mean.df$nparticles <- factor(mean.df$nparticles)

label.df <- data.frame(time = c(1,1,1,1,1),
                       nparticles = c(256, 512, 1024, 4096, 2048))
label.df <- merge(label.df, mean.df %>% select(time, nparticles, v), by = c('time', 'nparticles'))
label.df$nparticles <- factor(label.df$nparticles)

g <- ggplot(mean.df, aes(x = time, y = v, group = nparticles, colour = nparticles))
g <- g + geom_line() + ylab("variance") + xlab("time")
g <- g + scale_y_log10(breaks = c(1, 1e2, 1e4)) 
g <- g + scale_colour_manual(name = "N:", values = cbbPalette) #+ scale_linetype(name = "N:")
g <- g + geom_label(data = label.df,
          aes(x = time, y = v, label = nparticles, group =  nparticles, colour =  nparticles),
          size = 6) + theme(legend.position = "none")
g

ggsave(filename = "ar1.effectnparticles.woas.pdf", plot = g)


label.df <- data.frame(time = c(1,1,1, 1, 1),
                       nparticles = c(256, 512, 1024, 4096, 2048))
label.df <- merge(label.df, mean.df %>% select(time, nparticles, efficiency), by = c('time', 'nparticles'))
label.df$nparticles <- factor(label.df$nparticles)
# label.df[3,3] <- 5e-5
label.df[4,3] <- 2.5e-5
g <- ggplot(mean.df, aes(x = time, y = efficiency, group = nparticles, colour = nparticles))
g <- g + geom_line() + ylab("efficiency") + xlab("time")
g <- g + geom_label(data = label.df,
                    aes(x = time, y = efficiency, label = nparticles, group =  nparticles, colour =  nparticles),
                    size = 6) + theme(legend.position = "none")
g <- g + scale_y_log10(breaks = c(1e-8, 1e-6, 1e-4)) + scale_colour_manual(name = "N:", values = cbbPalette) #+ scale_linetype(name = "N:")
g

ggsave(filename = "ar1.effectnparticles.effi.woas.pdf", plot = g)


