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

# Simulation function (parallel)

single_sim <- function(mod, iter = 100, size = 20, init = 10, t = 12, cores = 1){
  
  cl <- makeCluster(cores)
  clusterExport(cl, varlist = c("mod", "size", "init", "t"), envir = environment())
  
  t0 <- Sys.time()
  simlist <- parLapply(cl, 1:iter, function(x){
    b1 <- BacArena::Bac(mod)
    arena <- BacArena::Arena(n = size, m = size)
    arena <- BacArena::addOrg(arena, b1, amount = init)
    arena <- BacArena::addDefaultMed(arena, b1)
    sim <- BacArena::simEnv(arena, time = t)
  })
  stopCluster(cl)
  t1 <- Sys.time()
  print(t1-t0)
  
  return(simlist)
}

recfiles <- list.files("D:/MetabolicMicrobe/AGORA-1.01-Reconstructions/")
fnames <- paste("D:/MetabolicMicrobe/AGORA-1.01-Reconstructions/", recfiles, sep = "")

b1 <- readSBMLmod(fnames[1])

ncore <- detectCores()-1
sim1 <- single_sim(mod = b1, iter = 7, size = 20, init = 10, t = 12, cores = ncore)

cores <- detectCores()-1
cl <- makeCluster(cores)
clusterExport(cl, "b1")
t0 <- Sys.time()
simlist <- parLapply(cl, 1:100, function(x){
  b1 <- BacArena::Bac(b1)
  arena <- BacArena::Arena(n = 20, m = 20)
  arena <- BacArena::addOrg(arena, b1, amount = 10)
  arena <- BacArena::addDefaultMed(arena, b1)
  sim <- BacArena::simEnv(arena, time = 12)
})
stopCluster(cl)
t1 <- Sys.time()
t1-t0

p <- plotGrowthCurve(simlist)
p[[2]]



n1 <- sapply(sim1, function(x) sapply(x@simlist, nrow))
matplot(n1, pch = 20)

Rpars <- apply(n1, 2, function(x){optimize(f = grwMin, interval = c(-3,3), nTrue = x)})

pred <- sapply(Rpars, function(x){ode(c(N = 10), times = 1:13, func = grw, parms = x$minimum)[,-1]})
points(apply(pred, 1, mean), typ = "l", lty = 2)

points(apply(pred, 1, quantile, prob = .025), typ = "l", lty = 2, col = "blue")
points(apply(pred, 1, quantile, prob = .975), typ = "l", lty = 2, col = "blue")


rs <- quantile(sapply(Rpars, function(x) x$minimum), prob = c(.025, .5, .975))
pred <- sapply(rs, function(x){ode(c(N = 10), times = 1:25, func = grw, parms = x)[,-1]})
points(pred[,1], typ = "l", col = "darkgreen")
points(pred[,3], typ = "l", col = "darkgreen")
