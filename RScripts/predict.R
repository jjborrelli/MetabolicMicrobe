library(sybil)
library(sybilSBML)
library(BacArena)
library(deSolve)
library(parallel)

# Functions to find growth rate in individual org simulations

grw <- function(t, N, parms){
  with(as.list(parms), { # extract parameters from parms vector
    dN = N * parms * (1 - N/400)
    return(list(dN)) # return dn/dt as a list
  })
}

grwMin = function(parms, nTrue){
  n0 <- nTrue[1]
  times <- 1:length(nTrue)
  out = as.data.frame(lsoda(n0, times, grw, parms)) # run ode
  mse = mean((out$`1`-nTrue)^2) # calculate mean squared error between simulated and original data
  return(mse) # return mean squared error
}

grwMin_sse = function(parms, nTrue){
  n0 <- nTrue[1]
  times <- 1:length(nTrue)
  out = as.data.frame(lsoda(n0, times, grw, parms)) # run ode
  mse = sum((out$`1`-nTrue)^2) # calculate sum of squared error between simulated and original data
  return(mse) # return sum of squared error
}

# try fitting model with nls???


# Simulation function (parallel)

single_sim <- function(mod, subs, iter = 100, size = 20, init = 10, t = 12, cores = 1){
  
  cl <- makeCluster(cores)
  clusterExport(cl, varlist = c("mod", "size", "init", "t", "subs"), envir = environment())
  
  t0 <- Sys.time()
  simlist <- parLapply(cl, 1:iter, function(x){
    b1 <- BacArena::Bac(mod)
    arena <- BacArena::Arena(n = size, m = size)
    arena <- BacArena::addOrg(arena, b1, amount = init)
    arena <- BacArena::addSubs(arena, smax = subs$concentration.im.mM, 
                               mediac = subs$Exchange, difspeed = subs$diffconstant, addAnyway = TRUE)
    sim <- BacArena::simEnv(arena, time = t)
  })
  stopCluster(cl)
  t1 <- Sys.time()
  print(t1-t0)
  
  return(simlist)
}


# get fitted growth rate parameter

single_R <- function(simlst){
  n1 <- lapply(simlst, function(x) sapply(x@simlist, nrow))
  Rpars <- lapply(n1, function(x){optimize(f = grwMin, interval = c(-3,3), nTrue = x)})
  R <- sapply(Rpars, "[[", 1)
  mse <- sapply(Rpars, "[[", 2)
  
  return(list(r = R, mse = mse))
}


# plot predicted values

pred_plots <- function(simlst){
  n1 <- sapply(sim1, function(x) sapply(x@simlist, nrow))
  matplot(n1, pch = 20)
  
  Rpars <- apply(n1, 2, function(x){optimize(f = grwMin, interval = c(-3,3), nTrue = x)})
  
  pred <- sapply(Rpars, function(x){ode(c(N = 10), times = 1:13, func = grw, parms = x$minimum)[,-1]})
  points(apply(pred, 1, median), typ = "l", lty = 2, col = "darkgreen")
  
  points(apply(pred, 1, quantile, prob = .025), typ = "l", lty = 2, col = "darkgreen")
  points(apply(pred, 1, quantile, prob = .975), typ = "l", lty = 2, col = "darkgreen")
}



# Testing

## fake data
sdat <- ode(10, times = 1:13, func = grw, parms = .8)[,-1]
## add noise
smat <- matrix(0, nrow = 13, ncol = 100)
for(i in 1:100){
  smat[,i] <- abs(rnorm(1:length(sdat), sdat, sdat*.01))
}

## fit model
fit1 <- apply(smat, 2, function(x) optimize(f = grwMin, interval = c(-3,3), nTrue = x))
fit2 <- apply(smat, 2, function(x) optimize(f = grwMin_sse, interval = c(-3,3), nTrue = x))




# Real Data
# filenames of metabolic models
recfiles <- list.files("D:/MetabolicMicrobe/AGORA-1.01-Reconstructions/")
fnames <- paste("D:/MetabolicMicrobe/AGORA-1.01-Reconstructions/", recfiles, sep = "")

# matrix of media concentrations
submat <- read.csv("./Data/media-list.csv", stringsAsFactors = FALSE)

# filepath to dump simulation results
datapath <- "D:/MetabolicMicrobe/SingleMicrobeSimulations/sim_4-16-2018/"

# how many cores to use in parallel sim
ncore <- detectCores()-1

t0 <- Sys.time()
for(i in seq_along(fnames)){
  b1 <- readSBMLmod(fnames[i])
  
  sim.lst <- single_sim(mod = b1, subs = submat, iter = 7, size = 20, init = 10, t = 12, cores = ncore)
  params.lst <- single_R(sim.lst)
  
  saveRDS(sim.lst, paste0(datapath, "simlst_", i, "_", Sys.Date()))
  saveRDS(params.lst, paste0(datapath, "parlst_", i, "_", Sys.Date()))
}
t1 <- Sys.time()
t1=t0
