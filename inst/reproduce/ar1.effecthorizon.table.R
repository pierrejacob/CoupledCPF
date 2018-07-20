# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(ggthemes)
library(dplyr)
library(reshape2)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#0072B2", "#F0E442")
# load custom theme for the plots
# it requires the packages ggplot2, gridExtra, reshape
setmytheme()
# fix the random seed
set.seed(17)
module_tree <<- Module("module_tree", PACKAGE = "CoupledCPF")
TreeClass <<- module_tree$Tree

# load results with BPF
savefilename <- paste0("ar1.effecttimehorizon.R", 1000, ".RData")
load(savefilename)
estimates.df <- rbind(estimates_was.df, estimates_woas.df)
estimates.df$pf <- "BPF"

# load results with APF
savefilename <- paste0("ar1.effecttimehorizon.optimal.R", 1000, ".RData")
load(savefilename)

estimates.optimal.df <- rbind(estimates_was.df, estimates_woas.df)
estimates.optimal.df$pf <- "APF"
estimates.df <- rbind(estimates.df, estimates.optimal.df)
estimates.df %>% tail
times.df <- estimates.df %>%  group_by(irep, with_as, nparticles, datalength, pf) %>% summarise(iteration = mean(iteration))
times.df <- times.df %>% mutate(cost = nparticles * iteration) %>% select(irep, with_as, nparticles, datalength, iteration, pf)
m <- times.df %>% group_by(with_as, nparticles, datalength, pf) %>%  summarise(mean_iter = mean(iteration),sd_iter = sd(iteration))
m <- m %>% mutate(iter = paste0(round(mean_iter,2), " (", round(sd_iter,2), ")"), nt = paste0("N = ", nparticles, ", T = ", datalength))
m <- m %>% select(nt, with_as, iter, pf) %>% ungroup()
m

tm <- dcast(m, nt ~ with_as + pf, value.var = "iter")
colnames(tm) <- c("", "APF without ancestor sampling", "BPF without ancestor sampling", "APF with ancestor sampling", "BPF with ancestor sampling")
tm <- tm[c(2,4,5,1,3),c(1,3,5,2,4)]

library(xtable)

cap <- "Average meeting time, as a function of the number of particles $N$ and the time
horizon $T$, with bootstrap particle filters and auxiliary particle filters, 
with and without ancestor sampling, computed over $R=1,000$ experiments.
Standard deviations are between brackets. Results obtained in the hidden auto-regressive model
of Section \\ref{sec:numerics:hiddenar}.
\\label{table:effecthorizon}"
print.xtable(xtable(tm, caption = cap), include.rownames = F, file = "tablear1horizon.tex")

 
