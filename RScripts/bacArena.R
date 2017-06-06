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


grw <- function(t, N, parms){
  with(as.list(parms), { # extract parameters from parms vector
    dN = N * R * (1 - N/K)
    return(list(dN)) # return dn/dt as a list
  })
}

grw2 <- function(t, N, parms){
  parms <- as.list(parms)
  for(i in 2:t){
    N[i] <- round(N[i-1] * parms$R * (1 - N[i-1]/parms$K))
  }
  return(N)
}

ode(c(N = 10), times = 1:10, func = grw, parms = list(R = .8, K = 400))

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

n <- sapply(eval@simlist, nrow)
optout2 <- optim(par = list(R = 1.4), fn = grwMin2, nTrue = n)

plot(ode(c(N = 20), times = 1:20, func = grw, parms = optout2$par)[,-1], typ= "o")
plot(grw2(200, 20, list(R = 4, K = 500)), typ = "l")


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
eval <- simEnv(arena, time = 20)

# Plotting results
par(mfrow = c(2,1))
plotCurves2(eval, legendpos = "topleft")




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
