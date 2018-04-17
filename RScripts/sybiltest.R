library(sybil)
library(sybilSBML)
library(BacArena)

recfiles <- list.files("D:/MetabolicMicrobe/AGORA-1.01-Reconstructions/")
#f1 <- sample(1:length(recfiles), 2)
#fnames <- paste("D:/MetabolicMicrobe/AGORA-1.01-Reconstructions/", recfiles[f1], sep = "")
fnames <- paste("D:/MetabolicMicrobe/AGORA-1.01-Reconstructions/", recfiles, sep = "")
grep("Clostridium", recfiles)
grep("Lachnospiracea", recfiles)
grep("Barnesiella", recfiles)

grep("Lactobacillus_rhamnos", recfiles)[1]
grep("Bacteroides_cacca", recfiles)
b1 <- readSBMLmod(fnames[77]) #bacteroides
b1 <- Bac(b1)
b2 <- readSBMLmod(fnames[490]) #lactobacillus
b2 <- Bac(b2)


test <- readSBMLmod("C:/Users/jjborrelli/Downloads/Clostridium_difficile_NAP07.xml")
data("Ec_core")
cdif <- Bac(test)
ecoli <- Bac(Ec_core)

arena <- Arena(n=20, m=20)
arena <- addOrg(arena, b2, amount = 20)
arena <- addOrg(arena, b1, amount = 20)
arena <- addDefaultMed(arena, b2)
arena <- addDefaultMed(arena, b1)
#arena <- addSubs(arena, smax=0.5, mediac="EX_glc(e)", unit="mM")

eval <- simEnv(arena,time=12)
par(mfrow=c(1,2))
plotCurves2(eval, legendpos = "bottomright")
par(mfrow=c(1,1))
sum(eval@simlist[[12]]$biomass)
sapply(eval@simlist, function(x) sum(x$biomass)) %>% plot(typ = "l")


arena <- Arena(n=20, m=20)
arena2 <- addOrg(arena, b2,amount=20)
arena1 <- addOrg(arena, b1, amount = 20)
arena2 <- addDefaultMed(arena2, b2)
arena1 <- addDefaultMed(arena1, b1)

eval1 <- simEnv(arena1,time=12)
eval2 <- simEnv(arena2,time=12)



b.1 <- sapply(eval1@simlist, function(x) sum(x$biomass))# %>% plot(typ = "l")
b.2 <- sapply(eval2@simlist, function(x) sum(x$biomass))# %>% plot(typ = "l")
b.1_2 <- t(sapply(eval@simlist, function(x) aggregate(x$biomass, list(x$type), sum)$x))

sapply(lapply(list(b.1, b.2), function(x) optim(par = list(R = 1.4, K = 800), fn = grwMin, nTrue = x)), function(x) x$par)
sapply(apply(b.1_2, 2, function(x) optim(par = c(R = 1), fn = grwMin, nTrue = x)), function(x) x$par)


ggplot(eval@simlist[[13]], aes(x, y)) + geom_point(aes(alpha = biomass, col = factor(type))) + scale_color_manual(values = c("blue", "darkgreen")) + theme_bw(base_size = 20) + guides(alpha = FALSE, col = FALSE)
ggplot(data.frame(bio = b.1_2, time = 1:13), aes(x = time, y = bio.1)) + geom_line(col = "blue") + geom_point(col = "blue") + geom_line(aes(x = time, y = bio.2), col = "darkgreen") + geom_point(aes(x = time, y = bio.2), col = "darkgreen") + theme_bw(base_size = 20) + labs(y = "Biomass")

library(multipanelfigure)
library(animation)


saveGIF(
  {
    for(i in c(2:13,rep(13, 5))){
      g1 <- ggplot(eval@simlist[[i]], aes(x, y)) + 
        geom_point(aes(size = biomass, col = factor(type))) + scale_color_manual(values = c("blue", "green3")) +
        theme_bw(base_size = 20) + guides(size = FALSE, col = FALSE) 
      g2 <- ggplot(data.frame(bio = b.1_2[1:i,], time = 1:i), aes(x = time, y = bio.1)) + geom_line(col = "blue") + 
        geom_point(col = "blue") + geom_line(aes(x = time, y = bio.2), col = "green3") + 
        geom_point(aes(x = time, y = bio.2), col = "green3") + theme_bw(base_size = 20) + labs(y = "Biomass")
      
      cols <- 2
      rows <- 1
      figure <- multi_panel_figure(
        width = 10,
        columns = cols,
        height = 6,
        rows = rows, unit = "in", figure_name = paste("Time = ", i))
      
      figure %<>% fill_panel(g1, row = 1, column = 1)
      figure %<>% fill_panel(g2, row = 1, column = 2)
      print(figure)
    }
  }, 
  movie.name = "metabolic.gif", interval = 1, ani.width = 1000, ani.height = 600, outdir = "~/Desktop/"
)




lvcomp = function(t, n, parms) {
  with(as.list(parms), { # extract parameters from parms vector
    dn = rep(0, 2) # initialize dn/dt vector
    dn[1] = r1*n[1]*(K1-n[1]-a1*n[2])/K1 # species 1
    dn[2] = r2*n[2]*(K2-n[2]-a2*n[1])/K2 # species 2
    return(list(dn)) # return dn/dt as a list
  })
}

# function to minimize for fitting data to competition model
LVcompMin = function(par, nTrue, times = 1:13) {
  parms = list(r1 = par[1], r2 = par[2], K1 = par[3], K2 = par[4], a1 = par[5], a2 = par[6]) 
  out = as.data.frame(lsoda(nTrue[1,], times, lvcomp, parms)) # run ode
  mse = mean((rbind(out$`1`, out$`2`)-t(nTrue))^2) # calculate mean squared error between simulated and original data
  return(mse) # return mean squared error
}

optim(par = c(1, 1, 400, 400, .2, -.1), fn = LVcompMin, nTrue = b.1_2)

par = c(0.9234324, 0.9390510, 438.8911796, 399.6226466, 0.1202402, 0.2700562)
parms = list(r1 = par[1], r2 = par[2], K1 = par[3], K2 = par[4], a1 = par[5], a2 = par[6]) 
out = ode(y = nTrue[1,], times=1:13, func = lvcomp,parms =  parms)
matplot(out[,-1], typ = "l")
library(rootSolve)
jf <- jacobian.full(y = out[13,-1], func = lvcomp, parms = parms)
max(Re(eigen(jf)))


ggplot(data.frame(out), aes(x = time, y = X1)) + geom_line(col = "blue") + 
  geom_point(col = "blue") + geom_line(aes(x = time, y = X2), col = "green3") + 
  geom_point(aes(x = time, y = X2), col = "green3") + 
  geom_line(data = data.frame(bio = b.1_2[1:13,], time = 1:13), aes(x = time, y = bio.1), alpha = .25) +
  geom_point(data = data.frame(bio = b.1_2[1:13,], time = 1:13), aes(x = time, y = bio.1), alpha = .25, col = "blue") + 
  geom_line(data = data.frame(bio = b.1_2[1:13,], time = 1:13), aes(x = time, y = bio.2), alpha = .25, col = "green3") + 
  geom_point(data = data.frame(bio = b.1_2[1:13,], time = 1:13), aes(x = time, y = bio.2), alpha = .25, col = "green3") + 
  theme_bw(base_size = 20) + labs(y = "Biomass")


pmat <- getPhenoMat(eval)
pmat[,which(colSums(pmat)>0)]

strt <- Sys.time()
#cdf <- grep("Clostridium", recfiles)
dbio <- list()
#for(i in 1:length(grep("Clostridium", recfiles))){
for(i in 1:length(recfiles)){
  #fnames <- paste("D:/MetabolicMicrobe/AGORA-1.01-Reconstructions/", recfiles[cdf[i]], sep = "")
  fnames <- paste("D:/MetabolicMicrobe/AGORA-1.01-Reconstructions/", recfiles[i], sep = "")
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

grwMin = function(parms, nTrue){
  n0 <- nTrue[1]
  times <- 1:length(nTrue)
  out = as.data.frame(lsoda(n0, times, grw, list(R = parms, K = 700))) # run ode
  mse = mean((out$`1`-nTrue)^2) # calculate mean squared error between simulated and original data
  return(mse) # return mean squared error
}

optout2 <- lapply(dbio, function(x) optim(par = list(R = 1.4, K = 400), fn = grwMin, nTrue = x))
