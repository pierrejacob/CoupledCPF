# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledCPF)
library(doRNG)
ncores <- 10
registerDoMC(cores = ncores)
# load custom theme for the plots
# it requires the packages ggplot2, gridExtra, reshape
setmytheme()
# fix the random seed
set.seed(17)
#
module_tree <<- Module("module_tree", PACKAGE = "CoupledCPF")
TreeClass <<- module_tree$Tree

dimension <- 1
ar1 <- get_ar(dimension)

# Load synthetic dataset
load("ar1data.RData")
datalength <- 20
observations <- matrix(observations[1:datalength,], ncol = 1)


algoparameters <- list(lambda = 0, with_as = FALSE)

theta <- 0.9
nrep <- 1000
seq_nparticles <- c(50, 100, 150, 200)
#### Generate unbiased estimators with systematic resampling
algoparameters$coupled_resampling <- CR_systematic
meeting_times.systematic.df <- data.frame()
for (N in seq_nparticles){
  print(N)
  algoparameters$nparticles <- N
  meeting_times.df_ <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
    res <- rheeglynn_estimator(observations, ar1, theta, algoparameters)
    data.frame(irep = irep, meeting = res$iteration, nparticles = N, resampling = "systematic")
  }
  meeting_times.systematic.df <- rbind(meeting_times.systematic.df, meeting_times.df_)
  save(nrep, seq_nparticles, datalength, observations, theta, algoparameters, meeting_times.systematic.df,
       file = "ar1.effectresampling.RData")
}
meeting_times.systematic.df %>% tail

#### Generate unbiased estimators with index-matching resampling
algoparameters$coupled_resampling <- CR_indexmatching
meeting_times.indexmatching.df <- data.frame()
for (N in seq_nparticles){
  print(N)
  algoparameters$nparticles <- N
  meeting_times.df_ <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
    res <- rheeglynn_estimator(observations, ar1, theta, algoparameters)
    data.frame(irep = irep, meeting = res$iteration, nparticles = N, resampling = "indexmatching")
  }
  meeting_times.indexmatching.df <- rbind(meeting_times.indexmatching.df, meeting_times.df_)
  save(nrep, seq_nparticles, datalength, observations, theta, algoparameters,
       meeting_times.systematic.df, meeting_times.indexmatching.df,
       file = "ar1.effectresampling.RData")
 }

