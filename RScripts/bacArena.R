# Load BacArena package
library(BacArena)
####
# Single organism
####

# Set up E. coli metabolic model
data("Ec_core")
str(Ec_core)

# Construct model with proper class
bac <- Bac(Ec_core)

# Set up environment
arena <- Arena(n = 20, m = 20)
## 20x20 grid
arena

# put the organism (E. coli) into the environment
arena <- addOrg(arena, bac, amount = 20)
## adds 20 individuals of E. coli
arena

# add substances to the environment
arena <- addDefaultMed(arena, bac)
arena
arena <- addSubs(arena, smax = 0.5, mediac <- "EX_glc(e)", unit = "mM")
arena

# Start in silico experiment
## 10 hours of growth
eval <- simEnv(arena, time = 12)

#### 
## Check out the results
####

# look at substances with high variation
getVarSubs(eval)
## Time series for one substance
getSubHist(eval, "EX_glc(e)")

# Growth curve with substance variations
par(mfrow = c(2,1))
plotCurves2(eval, legendpos = "topleft")

# Spatio-temporal changes
par(mfrow= c(5,2))
evalArena(eval, show_legend = F, time = 1:10)

eval@simlist[[2]]

grw <- function(t, N, parms){
  with(as.list(parms), { # extract parameters from parms vector
    dN = N * R * (1 - N/K)
    return(list(dN)) # return dn/dt as a list
  })
}

grw2 <- function(t, N, parms){
  with(as.list(parms), { # extract parameters from parms vector
    dN = R * N * (K - N - A * N[c(2,1)]) / K
    return(list(dN)) # return dn/dt as a list
  })
}


grw3 <- function(t, N, parms){
  with(as.list(parms), { # extract parameters from parms vector
    dN = N * R[ceiling(t)] * (1 - N/K)
    return(list(dN)) # return dn/dt as a list
  })
}


plot(ode(c(N = 20), times = 1:10, func = grw, parms = list(R = .3, K = 400))[,-1])



grwMin = function(parms, nTrue){
  n0 <- nTrue[1]
  times <- 1:length(nTrue)
  out = as.data.frame(lsoda(n0, times, grw, parms)) # run ode
  mse = mean((out$`1`-nTrue)^2) # calculate mean squared error between simulated and original data
  return(mse) # return mean squared error
}

grwMin2 <- function(parms, nTrue){
  n0 <- nTrue[1]
  times <- length(nTrue)
  parms <- append(parms, c(K = 900))
  out <- grw2(times, n0, parms)
  mse <- mean((out-nTrue)^2)
  return(mse)
}

n <- t(sapply(eval@simlist, function(x) as.vector(table(x$type))))
optout2 <- optim(par = list(R = 1.4, K = 400), fn = grwMin, nTrue = n)

plot(ode(c(N = 20), times = 1:13, func = grw, parms = optout2$par)[,-1], typ= "o")
points(n, typ = "o", pch = 20)
plot(grw2(200, 20, list(R = 4, K = 500)), typ = "l")

pts <- matrix(nrow = 13, ncol = 2000)
for(i in 1:1000){
  pts[,i] <- (ode(c(N = 20), times = 1:13, func = grw3, parms = list(R = rnorm(15, .7768, .2), K = 400))[,-1])
  print(i)
}
matplot(pts, typ = "l", col = "gray", lty = 1)
points(n, typ = "o", col = "blue")

Rs <- rnorm(15, .7768, .1)
n1 <- (ode(c(N = 20), times = 1:13, func = grw3, parms = list(R = Rs, K = 400))[,-1])
plot(n1)
optout2 <- optim(par = list(R = 1.4, K = 400), fn = grwMin, nTrue = n1)

####
# Multiple Organisms
####
# Create two E. coli types
bac1 <- Bac(Ec_core, type = "ecoli_wt") # "wild type" E. coli

ecore_aux <- changeBounds(Ec_core, "EX_o2(e)", lb=0)
bac2 <- Bac(ecore_aux, type = "ecoli_aux", setExInf = FALSE) # auxotrophic mutant that can't use aerobic resp

# Set up environment with orgs and subs
arena <- Arena(n = 30, m = 30)
arena <- addOrg(arena, bac1, amount = 20)
arena <- addOrg(arena, bac2, amount = 20)
arena <- addDefaultMed(arena, bac1)
arena <- addSubs(arena, smax = 0.5, mediac = "EX_glc(e)", unit = "mM")
arena

# Simulate
eval <- simEnv(arena, time = 15)

# Plotting results
par(mfrow = c(2,1))
plotCurves2(eval, legendpos = "topleft")[2]


grwMin1 = function(parms, nTrue, K){
  n0 <- nTrue[1]
  parms <- list(R = parms, K = K)
  times <- 1:length(nTrue)
  out = as.data.frame(lsoda(n0, times, grw, parms)) # run ode
  mse = mean((out$`1`-nTrue)^2) # calculate mean squared error between simulated and original data
  return(mse) # return mean squared error
}

wboth <- t(sapply(eval@simlist, function(x) as.vector(table(x$type))))
n.aux <- sapply(eval@simlist, nrow)
n.aux2 <- wboth[,2]
opt_aux <- optim(par = list(R = 1.4), fn = grwMin1, nTrue = n.aux, K = 900)
opt_aux2 <- optim(par = list(R = 1.4), fn = grwMin1, nTrue = n.aux2, K = 900)
n1 <- (ode(c(N = 20), times = 1:20, func = grw, parms = list(R = opt_aux$par[1], K = 900))[,-1])
n1.2 <- (ode(c(N = 20), times = 1:20, func = grw, parms = list(R = opt_aux2$par[1], K = 900))[,-1])

n.wt <- sapply(eval@simlist, nrow)
n.wt2 <- wboth[,1]
opt_wt <- optim(par = list(R = 1.4), fn = grwMin, nTrue = n.wt, K = 900)
opt_wt2 <- optim(par = list(R = 1.4), fn = grwMin, nTrue = n.wt2, K = 900)
n2 <- (ode(c(N = 20), times = 1:20, func = grw, parms = list(R = opt_wt2$par[1], K = 900))[,-1])
n2.2 <- (ode(c(N = 20), times = 1:20, func = grw, parms = list(R = opt_wt$par[1], K = 900))[,-1])


dgr.aux <- opt_aux2[[1]]-opt_aux[[1]]/opt_aux[[1]]
dgr.wt <- opt_wt2[[1]]-opt_wt[[1]]/opt_wt[[1]]

intMin <- function(parms, nTrue, K, R){
  n0 <- nTrue[1,]
  parms <- list(R = R, K = K, A = parms)
  times <- 1:length(nTrue)
  out = as.data.frame(lsoda(n0, times, grw2, parms)) # run ode
  mse = mean((out-nTrue)^2) # calculate mean squared error between simulated and original data
}

lvm <- function(t, N, parms){
  with(as.list(parms), { # extract parameters from parms vector
    dN = N * R * (1 - N/K) + m %*% N * N
    return(list(dN)) # return dn/dt as a list
  })
}

lvm2 <- function(t, N, parms){
  with(as.list(parms), { # extract parameters from parms vector
    dN <- c(0,0)
    dN[1] = N[1] * R[1] * (K[1] - N[1] - A[1]*N[2]) / K[1]
    dN[2] = N[2] * R[2] * (K[2] - N[2] - A[2]*N[1]) / K[2]
    return(list(dN)) # return dn/dt as a list
  })
}

opt.both <- optim(par = c(A = c(.8,.3)), fn = intMin, nTrue = wboth, K = c(900, 900), R = c(.775, .372))
out1 <- lsoda(c(20,20), 1:20, grw2, parms = list(R = c(.775, .372), K = c(900,900), A = opt.both$par))[,-1]
matplot(out1, typ = "l")
points(n[,1])
points(n[,2])
lines(1:16, out1[1:16,1],col = "red", lty = 2)
lines(1:16, out1[1:16,2],col = "blue", lty = 2)
####
# Communities
####

# sample run of simEnv on SIHUMI data
data("sihumi_test")
eval <- sihumi_test
plotAbundance(eval)

plotCurves2(eval, legendpos = "topleft")
g <- findFeeding3(eval, time = 10, mets = c("EX_acald(e)", "EX_lac_D(e)", "EX_glc(e)","EX_co2(e)", "EX_h(e)",   "EX_h2o(e)", "EX_nh4(e)", "EX_o2(e)",  "EX_pi(e)"))

g <- findFeeding3(eval, time = 10, mets = names(eval@subchange[eval@subchange > 0])[-1])
g1 <- frame2webs(data.frame(lower = g[[1]][,3], higher = g[[1]][,2], webID = 1, freq = 1))[[1]]
g2 <- frame2webs(data.frame(higher = g[[1]][,3], lower = g[[1]][,1], webID = 1, freq = 1))[[1]]

plotweb(g1)
plotweb(g2)
plotweb2(g1, g2)

plot(graph.adjacency(get.adjacency(g[[2]], sparse = F), weighted = T), edge.width = E(graph.adjacency(get.adjacency(g[[2]], sparse = F), weighted = T))$weight, edge.arrow.size = 1)
