# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(doRNG)
library(ggthemes)
library(xtable)
library(dplyr)
library(reshape2)

ncores <- 10
registerDoMC(cores = ncores)
# load custom theme for the plots
# it requires the packages ggplot2, gridExtra, reshape
setmytheme()
# fix the random seed
set.seed(17)

load("ar1.effectresampling.RData")

# distribution of meeting times with systematic resampling
ggplot(meeting_times.systematic.df, aes(x = meeting, group = nparticles, fill = factor(nparticles))) + geom_histogram(position = "dodge") +
  scale_fill_colorblind(name = "N:")

# distribution of meeting times with index-matching resampling
ggplot(meeting_times.indexmatching.df, aes(x = meeting, group = nparticles, fill = factor(nparticles))) + geom_histogram(position = "dodge") +
  scale_fill_colorblind(name = "N:")

# The mean average times can be represented in a plot
meeting_times.df <- rbind(meeting_times.systematic.df, meeting_times.indexmatching.df)
meeting_times.df %>% tail

mean.df <- meeting_times.df %>% group_by(resampling, nparticles) %>% summarise(m = mean(meeting), s = sd(meeting))
ggplot(mean.df, aes(x = nparticles, y = m, group = resampling, colour = resampling)) + geom_point() + geom_line() + 
  facet_wrap(~ resampling)

ggplot(mean.df %>% filter(resampling == "systematic"), aes(x = nparticles, y = m, group = resampling, colour = resampling)) + geom_point() + geom_line() 
ggplot(mean.df %>% filter(resampling == "indexmatching"), aes(x = nparticles, y = m, group = resampling, colour = resampling)) + geom_point() + geom_line() 

# Or in a table, as in the article
m <- mean.df %>% mutate(v = paste0(round(m,2), " (", round(s,2), ")")) %>% select(resampling, nparticles, v) %>% ungroup()
m$resampling <- factor(m$resampling, levels = c("systematic", "indexmatching"), labels = c("systematic", "index-coupled"))
tableau <- dcast(m, resampling ~ nparticles)
names(tableau) <- c("Resampling", paste0("N = ", seq_nparticles))

m <- t(tableau)
colnames(m) <- m[1,]
m <- m[2:nrow(m),]
cap <- "Average meeting time as a function of the number of particles $N$ and of the resampling scheme. 
Standard deviations are between brackets. Results obtained in the hidden auto-regressive model with $T = 20$. \\label{table:effectresampling}"
print.xtable(xtable(m, caption = cap), include.rownames = TRUE, include.colnames = FALSE)

# save table
print.xtable(xtable(m, caption = cap), include.rownames = TRUE, file = "ar1effectresampling.tex")

