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
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#0072B2", "#F0E442")
# fix the random seed
set.seed(17)
module_tree <<- Module("module_tree", PACKAGE = "CoupledCPF")
TreeClass <<- module_tree$Tree

load(file = "ar1.effectnparticles.RData")


estimates.df <- rbind(estimates_was.df, estimates_woas.df)
# estimates.df <- estimates.df %>% filter(nparticles > 256)

times.df <- estimates.df %>%  group_by(irep, with_as, nparticles) %>% summarise(iteration = mean(iteration))
times.df <- times.df %>% mutate(cost = nparticles * iteration) %>% select(irep, with_as, nparticles, iteration, cost)
m <- times.df %>% group_by(with_as, nparticles) %>%  summarise(meancost = mean(cost), sdcost = sd(cost), mean_iter = mean(iteration),
                                                      sd_iter = sd(iteration))

m <- m %>% mutate(cost = paste0(round(meancost,0), " (", round(sdcost,0), ")"),
                  iter = paste0(round(mean_iter,2), " (", round(sd_iter,2), ")")) %>% select(nparticles, cost, iter) %>% ungroup()

tm <- t(dcast(m, nparticles ~ with_as, value.var = "iter"))
colnames(tm) <- paste0("N = ", tm[1,])
tm <- tm[2:3,]
rownames(tm) <- c("without ancestor sampling", "with ancestor sampling")
tm

# dcast(m,  cost  ~  nparticles)

# names(tableau) <- c("Resampling", paste0("N = ", seq_nparticles))
library(xtable)
xtable(t(tm))

cap <- "Average meeting time, as a function of the number of particles $N$, with and without ancestor sampling.
Standard deviations are between brackets. Results obtained in the hidden auto-regressive model with $T = 500$.
\\label{table:effectancestor}"
print.xtable(xtable(t(tm), caption = cap), file = "ar1effectancestor.tex")

RGSmean.df <- estimates.df %>% group_by(time, with_as, nparticles) %>% summarise(m = mean(estimate), v= var(estimate), iteration = mean(iteration))
mean.df <- RGSmean.df # merge(RGSmean.df, smoothing_means_kf, by = "time")
mean.df <- mean.df %>% mutate(efficiency = 1/(v * iteration * nparticles)) %>% ungroup()


mean.df$nparticles <- factor(mean.df$nparticles)
label.df <- data.frame(time = c(1,1,1,1,1),
                       nparticles = c(256, 512, 1024, 4096, 2048))
label.df <- merge(label.df, mean.df %>% filter(with_as == TRUE) %>% select(time, nparticles, v, efficiency), by = c('time', 'nparticles'))
label.df$nparticles <- factor(label.df$nparticles)
label.df[3,3] <- 3.3
label.df[4,3] <- 1.9
label.df[2,3] <- 2.2

g <- ggplot(mean.df %>% filter(with_as == TRUE), aes(x = time, y = v, group = nparticles,
                                                      colour = nparticles))
g <- g + geom_line() + ylab("variance") + xlab("time")
g <- g + scale_y_log10(breaks = c(2, 4, 8), limits = c(1.5,10))
g <- g + scale_colour_manual(name = "N:", values = cbbPalette)
g <- g + geom_label(data = label.df,
                    aes(x = time, y = v, label = nparticles, group =  nparticles, colour =  nparticles),
                    size = 6) + theme(legend.position = "none")
g

ggsave(filename = "ar1.effectnparticles.var.was.pdf", plot = g)

g <- ggplot(mean.df %>% filter(with_as == TRUE), aes(x = time, y = efficiency, group = nparticles,
                                                     colour = nparticles))
g <- g + geom_line() + ylab("efficiency") + xlab("time")
g <- g + scale_y_log10() 
g <- g + scale_colour_manual(name = "N:", values = cbbPalette)
g <- g + geom_label(data = label.df,
                    aes(x = time, y = efficiency, label = nparticles, group =  nparticles, colour =  nparticles),
                    size = 6) + theme(legend.position = "none")
g
ggsave(filename = "ar1.effectnparticles.effi.was.pdf", plot = g)


