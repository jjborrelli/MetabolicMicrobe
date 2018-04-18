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



# Simulation function (parallel)

single_sim <- function(mod, subs, iter = 100, size = 20, init = 10, t = 12, cores = 1, ...){
  
  cl <- makeCluster(cores)
  clusterExport(cl, varlist = c("mod", "size", "init", "t", "subs"), envir = environment())
  
  t0 <- Sys.time()
  simlist <- parLapply(cl, 1:iter, function(x){
    b1 <- BacArena::Bac(mod, growtype = "exponential", speed = 10)
    arena <- BacArena::Arena(n = size, m = size, ...)
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


pair_sim <- function(mods, subs, iter = 100, size = 20, init = 10, t = 12, cores = 1, ...){
  
  cl <- makeCluster(cores)
  clusterExport(cl, varlist = c("mods", "size", "init", "t", "subs"), envir = environment())
  
  t0 <- Sys.time()
  simlist <- parLapply(cl, 1:iter, function(x){
    b1 <- BacArena::Bac(mods[[1]], growtype = "exponential", speed = 10)
    b2 <- BacArena::Bac(mods[[2]], growtype = "exponential", speed = 10)
    arena <- BacArena::Arena(n = size, m = size, ...)
    arena <- BacArena::addOrg(arena, b1, amount = init)
    arena <- BacArena::addOrg(arena, b2, amount = init)
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


# pairwise competition coefficients

pair_int <- function(t, state, parms){
  with(as.list(state, parms), { # extract parameters from parms vector
    dN1 = N1 * parms$r1 * (1 - N1/400) + parms$a21 * N2
    dN2 = N2 * parms$r2 * (1 - N2/400) + parms$a12 * N1
    return(list(c(dN1, dN2))) # return dn/dt as a list
  })
}

lvcomp = function(t, state, parms) {
  with(as.list(state, parms), { # extract parameters from parms vector
    
    dN1 = parms$r1*N1*(parms$K1-N1-parms$a1*N2)/parms$K1 # species 1
    dN2 = parms$r2*N2*(parms$K2-N2-parms$a2*N1)/parms$K2 # species 2
    
    return(list(c(dN1, dN2))) # return dn/dt as a list
  })
}

prMin = function(parms, growth, nTrue){
  n0 <- c(nTrue[1,1], nTrue[2,1])
  names(n0) <- c("N1", "N2")
  times <- 1:ncol(nTrue)
  parms <- append(growth, parms)
  out <- as.data.frame(lsoda(n0, times = times, func = lvcomp, parms = parms)) # run ode
  mse <- sum((out$N1-nTrue[1,])^2)+sum((out$N2 - nTrue[2,])^2) # calculate mean squared error between simulated and original data
  return(mse) # return mean squared error
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
sdat <- ode(10, times = 1:13, func = grw, parms = 1.9)[,-1]
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
t1-t0


simres <- list.files("D:/MetabolicMicrobe/SingleMicrobeSimulations/sim_4-16-2018/")
rlist <- list()
for(i in 1:773){
  res1 <- grep(paste0("parlst_", i, "_"), simres)
  plist <- readRDS(paste0("D:/MetabolicMicrobe/SingleMicrobeSimulations/sim_4-16-2018/",simres[res1]))
  rlist[[i]] <- plist$r
}

hist(sapply(rlist, mean))

b1 <- readSBMLmod(fnames[1])
b2 <- readSBMLmod(fnames[3])

simb1 <- single_sim(mod = b1, subs = submat, iter = 7, size = 20, init = 10, t = 12, cores = ncore)
simb2 <- single_sim(mod = b2, subs = submat, iter = 7, size = 20, init = 10, t = 12, cores = ncore)

p <- plotGrowthCurve(simb1)
p[[2]]
p <- plotGrowthCurve(simb2)
p[[2]]

psim <- pair_sim(list(b1, b2), subs = submat, iter = 7, size = 20, init = 10, t = 12, cores = ncore)

p <- plotGrowthCurve(psim)
p[[2]]

n1 <- sapply(psim[[1]]@simlist, function(x) c(sum(x$type == 1), sum(x$type == 2)))
rb1 <- mean(single_R(simb1)$r)
optimize(grwMin, interval = c(-3,3), nTrue = n1[1,])$minimum
rb2 <- mean(single_R(simb2)$r)
optimize(grwMin, interval = c(-3,3), nTrue = n1[2,])$minimum


optpar <- optim(par = list(a1 = .5, a2 = .5), fn = prMin, growth = list(r1 = rb1, r2 = rb2, K1 = 400, K2 = 400), nTrue = n1)
optpar
optpar$par



pl <- list(r1 = rb1, r2 = rb2, a21 = optpar$par[1], a12 = optpar$par[2], K1 = 400, K2 = 400)
n0 <- c(nTrue[1,1], nTrue[2,1])
names(n0) <- c("N1", "N2")
times <- 1:ncol(nTrue)
out <- as.data.frame(ode(n0, times = times, func = lvcomp, parms = pl))
out

matplot(out[,-1], typ = "l", ylim = c(0,400))
points(n1[1,])
points(n1[2,])
