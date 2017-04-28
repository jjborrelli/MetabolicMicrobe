single.gr <- read.csv("~/Downloads/nbt.3703-S6.csv")
pairw.gr <- read.csv("~/Downloads/nbt.3703-S7.csv")

strA <- as.character(single.gr$AGORA.Strain)
strA[grep("=", single.gr$AGORA.Strain)] <- sapply(strA[grep("=", single.gr$AGORA.Strain)], function(x) substr(x, 1, gregexpr("=", x)[[1]][1][1]-2))
strA <- sapply(strA, gsub, pattern = "[[:punct:]]", replacement = " ")
strA <- sapply(strA, function(x) paste(strsplit(x, "  ")[[1]], collapse = " "))
strA <- sapply(strA, trimws, "right")


st1 <- sapply(pairw.gr$Strain1, function(x){paste(strsplit(as.character(x), "_")[[1]], collapse = " ")})
#st1 <- sapply(st1, function(x) paste(strsplit(x, " ")[[1]], collapse = ""))
st2 <- sapply(pairw.gr$Strain2, function(x){paste(strsplit(as.character(x), "_")[[1]], collapse = " ")})

wdi <- matrix(nrow = 773, ncol = 773)
colnames(wdi) <- single.gr$AGORA.Strain
rownames(wdi) <- single.gr$AGORA.Strain

fdi <- matrix(nrow = 773, ncol = 773)
colnames(fdi) <- single.gr$AGORA.Strain
rownames(fdi) <- single.gr$AGORA.Strain


res <- list()
res2 <- list()
for(i in 1:nrow(single.gr)){
  as1 <- (st1 %in% strA[i])
  gch <- ((pairw.gr$WesternDietAnaerobic_Strain1[as1] - single.gr$Western.diet..anaerobic..h.1.[i])/single.gr$Western.diet..anaerobic..h.1.[i])*100
  wdi[i,strA %in% st2[as1]] <- gch
  
  gchF <- ((pairw.gr$HighFiberDietAnaerobic_Strain1[as1] - single.gr$High.fiber.diet..anaerobic..h.1.[i])/single.gr$High.fiber.diet..anaerobic..h.1.[i])*100
  fdi[i,strA %in% st2[as1]] <- gchF
  
  as2 <- st2 %in% strA[i]
  gch2 <- ((pairw.gr$WesternDietAnaerobic_Strain2[as2] - single.gr$Western.diet..anaerobic..h.1.[i])/single.gr$Western.diet..anaerobic..h.1.[i])*100
  wdi[i, strA %in% st1[as2]] <- gch2
  
  gchF2 <- ((pairw.gr$HighFiberDietAnaerobic_Strain2[as2] - single.gr$High.fiber.diet..anaerobic..h.1.[i])/single.gr$High.fiber.diet..anaerobic..h.1.[i])*100
  fdi[i, strA %in% st1[as2]] <- gchF2
  
  print(i)
  
}

wdi[abs(wdi) < 10] <- 0
fdi[abs(wdi) < 10] <- 0
diag(wdi) <- 0
diag(fdi) <- 0
itystr(wdi)
itystr(fdi)

wdg <- graph.adjacency(abs(sign(wdi)))
genera <- sapply(strsplit(colnames(wdi), " "), "[[", 1)
table(genera)
ngen <- (aggregate(genera, list(genera), length)$x)
vsize <- ifelse(ngen > 5, ifelse(ngen > 20, 20, 15), 5)

adj <- get.adjacency(wdg, sparse = F)

library(ggraph)
ggraph(wdg) + geom_edge_link(alpha = .1) + geom_node_point()

lvmodK <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha * (1 - state/(parms$K/sum(state > 10^-10))) + state * parms$m %*% state 
    #dB <- state * parms$alpha * (1 - state/(parms$K)) + state * parms$m %*% state 
    
    list(dB)
  })
}


samp <- sample(1:773, 50)
wdi2 <- apply(wdi, c(1,2), function(x) runif(1)*sign(x))
wdi2 <- apply(wdi, 1, function(x) x/max(x))
diag(wdi2) <- -1
N <- 20
samp <- sample(1:773, N)
par.ag <- list(alpha = single.gr$Western.diet..anaerobic..h.1.[samp], m = wdi2[samp,samp], K = 20)
test <- ode(runif(N, 0, 20/N), times = 1:20, func = lvmodK, parms = par.ag, events = list(func = ext1, time =  1:20))
matplot(test[,-1], typ = "l", lwd = 3)
test[20,-1]
