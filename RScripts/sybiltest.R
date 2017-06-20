library(sybil)
library(sybilSBML)
library(BacArena)

recfiles <- list.files("D:/MetabolicMicrobe/AGORA-1.01-Reconstructions/")
f1 <- sample(1:length(recfiles), 2)
fnames <- paste("D:/MetabolicMicrobe/AGORA-1.01-Reconstructions/", recfiles[f1], sep = "")
grep("Clostridium", recfiles)

b1 <- readSBMLmod(fnames[1])
b1 <- Bac(b1)
b2 <- readSBMLmod(fnames[2])
b2 <- Bac(b2)


test <- readSBMLmod("C:/Users/jjborrelli/Downloads/Clostridium_difficile_NAP07.xml")
data("Ec_core")
cdif <- Bac(test)
ecoli <- Bac(Ec_core)

arena <- Arena(n=20, m=20)
arena <- addOrg(arena, b2,amount=20)
arena <- addOrg(arena, b1, amount = 20)
arena <- addDefaultMed(arena, b2)
arena <- addDefaultMed(arena, b1)
arena <- addSubs(arena, smax=0.5, mediac="EX_glc(e)", unit="mM")

eval <- simEnv(arena,time=12)
par(mfrow=c(1,2))
plotCurves2(eval, legendpos = "bottomright")
par(mfrow=c(1,1))
sum(eval@simlist[[12]]$biomass)
sapply(eval@simlist, function(x) sum(x$biomass)) %>% plot(typ = "l")

pmat <- getPhenoMat(eval)
pmat[,which(colSums(pmat)>0)]

strt <- Sys.time()
cdf <- grep("Clostridium", recfiles)
dbio <- list()
for(i in 1:length(grep("Clostridium", recfiles))){
  fnames <- paste("D:/MetabolicMicrobe/AGORA-1.01-Reconstructions/", recfiles[cdf[i]], sep = "")
  b1 <- Bac(readSBMLmod(fnames))
  arena <- Arena(n=20, m=20)
  arena <- addOrg(arena, b1,amount=20)
  arena <- addDefaultMed(arena, b1)
  eval <- simEnv(arena,time=10)
  dbio[[i]] <- sapply(eval@simlist, function(x) sum(x$biomass))

}

fin <- Sys.time()
fin-strt



recfiles[cdf]
grw <- function(t, N, parms){
  with(as.list(parms), { # extract parameters from parms vector
     dN = N * R * (1 - N/K)
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

optout2 <- lapply(dbio, function(x) optim(par = list(R = 1.4, K = 400), fn = grwMin, nTrue = x))
