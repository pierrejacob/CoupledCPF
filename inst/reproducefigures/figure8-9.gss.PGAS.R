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
module_tree <<- Module("module_tree", PACKAGE = "CoupledCPF")
TreeClass <<- module_tree$Tree

load("gssdata.RData")

dimension <- 1

gss <- get_gss()
theta <- c()
datalength <- 50
observations <- matrix(observations[1:datalength,], ncol = 1)

nparticles <- 2^10
with_as <- TRUE

xref <- CPF(nparticles, gss, theta, observations, with_as = with_as)
niterations <- 50000
trajectories <- matrix(ncol = datalength + 1, nrow = niterations)
for (i in 1:niterations){
  if (i %% 100 == 1) cat("iteration ", i, "/", niterations, "\n")
  xref <- CPF(nparticles, gss, theta, observations, xref, with_as = TRUE)
  trajectories[i,] <- xref
}

PGAS.df <- melt(trajectories)
names(PGAS.df) <- c("iteration", "time", "chain")
PGAS.df$time <- PGAS.df$time - 1

save(observations, nparticles, with_as, dimension, PGAS.df, 
     file = "gss.PGAS.RData")
