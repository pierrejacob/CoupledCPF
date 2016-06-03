# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(doRNG)
library(ggthemes)
ncores <- 4
registerDoMC(cores = ncores)
# load custom theme for the plots
# it requires the packages ggplot2, gridExtra, reshape
setmytheme()
# fix the random seed
set.seed(17)
module_tree <<- Module("module_tree", PACKAGE = "CoupledCPF")
TreeClass <<- module_tree$Tree

dimension <- 1
datalength <- 10
nparticles <- 2^10

theta <- c(0.9, 0.1, 0.1, 0.1)
load("unlikely.truth.RData")
truemeans <- data.frame(time = 0:datalength, truemean = mu_bar)
#####
Ns <- c(128, 256, 512, 1024)
RGresults.df <- data.frame()
iterations.df <- data.frame()

for (N in Ns){
  load(paste0("unlikely.RG.R10000N", N, "T10.RData"))
  nrep <- max(RGestimates$irep)
  print(summary(RGestimates$iteration))
  iterations.df <- rbind(iterations.df, data.frame(N = N, iteration = RGestimates$iteration))
  meanRG.df <- melt(RGestimates %>% select(-iteration), id = "irep")
  names(meanRG.df) <- c("irep", "time", "estimate")
  meanRG.df$time <- rep(0:datalength, each = nrep)
  meanRG.df <- meanRG.df %>% group_by(time) %>% summarise(m = mean(estimate), s = sd(estimate), 
                ymin = mean(estimate) - 2*sd(estimate)/sqrt(nrep), ymax = mean(estimate) + 2*sd(estimate)/sqrt(nrep))
  meanRG.df$N <- N
  RGresults.df <- rbind(RGresults.df, meanRG.df)
}
# RGresults.df %>% tail
RGresults.df$N <- factor(RGresults.df$N, levels = Ns, labels = paste0(Ns, "  "))
g <- ggplot(RGresults.df, aes(x = time, ymin = ymin, ymax = ymax, group = N, colour = N)) + geom_errorbar(position = "dodge")
g <- g + geom_point(data = truemeans, aes(x = time, y = truemean, ymin = NULL, ymax = NULL, colour = NULL), size = 3)
g <- g + scale_colour_colorblind(name = "N=")
g <- g  + scale_x_continuous(breaks = 0:datalength) + ylab("smoothing means") + ylim(0, 1)
g

# ggsave(filename = "unlikelyobsRG.pdf", plot = g, width = 7, height = 7)


iterations.df$N <- factor(iterations.df$N)
# iterations.df %>% tail
iteration.table <- iterations.df %>% group_by(N) %>% summarise(m = mean(iteration), s = sd(iteration), sum = sum(iteration))

library(xtable)
iteration.table <- iteration.table %>% mutate(meeting = paste0(round(m,2), " (", round(s,2), ")")) %>% select(N, meeting)
names(iteration.table) <- c("N", "meeting time")
iteration.table$N <- paste0("N =  ", iteration.table$N)
names(iteration.table) <- c("", "meeting time")

xtable(iteration.table)
cap <- "Average meeting time as a function of the number of particles $N$. 
Standard deviations are between brackets. \\label{table:meetingtimeunlikelyobs}"

# save table
# print.xtable(xtable(iteration.table, caption = cap), include.rownames = FALSE, include.colnames = TRUE, file = "iterationunlikely.tex")

##
Ns <- c(1024, 2048, 4096, 8192, 16384)
PFresults.df <- data.frame()
for (N in Ns){
  load(paste0("unlikely.PF.R10000N", N, "T10.RData"))
  nrep <- max(PFestimates$irep)
  meanPF.df <- melt(PFestimates, id = "irep")
  names(meanPF.df) <- c("irep", "time", "estimate")
  meanPF.df$time <- rep(0:datalength, each = nrep)
  meanPF.df <- meanPF.df %>% group_by(time) %>% summarise(m = mean(estimate), s = sd(estimate), 
                                                          ymin = mean(estimate) - 2*sd(estimate)/sqrt(nrep), ymax = mean(estimate) + 2*sd(estimate)/sqrt(nrep))
  meanPF.df$N <- N
  PFresults.df <- rbind(PFresults.df, meanPF.df)
}
PFresults.df$N <- factor(PFresults.df$N, levels = Ns, labels = paste0(Ns, "  "))
g <- ggplot(PFresults.df, aes(x = time, ymin = ymin, ymax = ymax, group = N, colour = N)) + geom_errorbar(position = "dodge")
g <- g + geom_point(data = truemeans, aes(x = time, y = truemean, ymin = NULL, ymax = NULL, colour = NULL), size = 3)
g <- g + scale_colour_colorblind(name = "N=", guide = guide_legend(title.position = "left", nrow = 1)) 
g <- g  + scale_x_continuous(breaks = 0:datalength) + ylab("smoothing means") + ylim(0, 1)
g

# ggsave(filename = "unlikelyobsPF.pdf", plot = g, width = 7, height = 7)





