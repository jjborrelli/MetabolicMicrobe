library(deSolve)


lvcomp = function(t, state, parms) {
  with(as.list(state, parms), { # extract parameters from parms vector
    
    dN1 = parms$r1*N1*(parms$K1-N1-parms$a1*N2)/parms$K1 # species 1
    dN2 = parms$r2*N2*(parms$K2-N2-parms$a2*N1)/parms$K2 # species 2
    
    return(list(c(dN1, dN2))) # return dn/dt as a list
  })
}

pl <- list(r1 = .8, r2 = .81, a21 = -.5, a12 = .6, K1 = 400, K2 = 400)
n0 <- c(nTrue[1,1], nTrue[2,1])
names(n0) <- c("N1", "N2")
times <- 1:ncol(nTrue)
out <- as.data.frame(ode(n0, times = times, func = lvcomp, parms = pl))
out

n2<- matrix(c(out$N1, out$N2), nrow = 2, byrow = TRUE)

optpar <- optim(par = list(a21 = -.2, a12 = -.2), fn = prMin, growth = list(r1 = rb1, r2 = rb2, K1 = 400, K2 = 400), nTrue = n2)
optpar


load("D:/MetabolicMicrobe/SIHUMI_models.RData")

arena <- BacArena::Arena(100, 100, stir=F, Lx=0.025, Ly=0.025, tstep=1)
for(i in 1:length(sihumi)){
  model = sihumi[[i]]
  bac = BacArena::Bac(model=model, growtype="exponential", speed=10)
  arena = BacArena::addOrg(arena, bac, amount=200)
}

arena = BacArena::addSubs(arena, smax=submat$concentration.im.mM,difspeed=submat$diffconstant,unit="mM",
                          mediac=as.character(submat$Exchange))

sim = BacArena::simEnv(arena, time=16)

plotGrowthCurve(sim)
