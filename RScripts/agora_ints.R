#single.gr <- read.csv("~/Downloads/nbt.3703-S6.csv")
#pairw.gr <- read.csv("~/Downloads/nbt.3703-S7.csv")
single.gr <- read.csv("C:/Users/jjborrelli/Desktop/GitHub/cometsUtil/Data/nbt.3703-S6.csv")
pairw.gr <- read.csv("C:/Users/jjborrelli/Desktop/GitHub/cometsUtil/Data/nbt.3703-S7.csv")


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

induced_subgraph(wdg, V(wdg)[genera == "Clostridium"])


adj <- get.adjacency(wdg, sparse = F)
df1 <- split(as.data.frame(adj), genera)
df1 <- lapply(df1, function(x){aggregate(abs(colSums(x)-nrow(x))/nrow(x), list(genera), max)$x})
# (aggregate(colSums(x), list(genera), sum)$x-nrow(x))/nrow(x))
df1 <- do.call(rbind, df1)

diag(df1) <- 0
#df1[df1 <= .65] <- 0 
wdg2 <- graph.adjacency(sign(df1))

l <-layout_with_lgl(wdg2)
plot(wdg2, vertex.label = NA, vertex.size = vsize, edge.arrow.mode = 0, layout = l)

wdg3 <- induced_subgraph(wdg2, V(wdg2)[vsize > 5])
l <-layout_with_kk(wdg3)

png(filename = "C:/Users/jjborrelli/Desktop/common-bac-graph.png", width = 700, height = 700, res = 100)
plot(wdg3, vertex.label.cex = 1, vertex.label = unique(genera)[vsize > 5], vertex.size = ceiling(ngen[vsize > 5]/2), edge.arrow.size = .75, layout = l, margin = c(0,0,0,0), vertex.color = "green3", vertex.frame.color = "black", vertex.label.color = "black")
dev.off()

library(ggraph)
ggraph(wdg) + geom_edge_link(alpha = .1) + geom_node_point()

library(deSolve)
lvmodK <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha * (1 - state/(parms$K/sum(state > 10^-10))) + state * parms$m %*% state 
    #dB <- state * parms$alpha * (1 - state/(parms$K)) + state * parms$m %*% state 
    
    list(dB)
  })
}
ext1 <- function (times, states, parms){
  with(as.list(states), {
    states[states < 10^-10] <- 0 
    
    return(c(states))
  })
}


strt <- Sys.time()
#wdi2 <- apply(wdi/100, 1, function(x) (x)/max(x))
wdi2 <- apply(sqrt(abs(wdi/100))*sign(wdi), 1, function(x) (x)/max(x))
diag(wdi2) <- -1

N <- 5
samp <- sample(1:773, N)
par.ag <- list(alpha = single.gr$Western.diet..anaerobic..h.1., m = wdi2, K = 20)


ist <- rep(0, 773)
ist[samp] <- runif(N, 10^-9, 20/N)
test <- ode(ist, times = 1:25, func = lvmodK, parms = par.ag, events = list(func = ext1, time =  1:25))
matplot(test[,-1], typ = "l", main = paste("i =",1, "\t", ":", "\t", "N =", sum(test[25,-1] > 0)))


istmat <- matrix(0, nrow = 100, ncol = 773)
istmat[1,] <- test[25,-1]
x <- c()
for(i in 1:99){
  cond <- FALSE
  c2 <- 0
  while(!cond){
    x[i] <- sample(1:773, 1)
    istmat[i,x[i]] <- istmat[i,x[i]]+runif(1, 10^-9, 20/sum(istmat[i,] > 0))
    test <- ode(istmat[i,], times = 1:25, func = lvmodK, parms = par.ag, events = list(func = ext1, time =  1:25))
    c1 <- nrow(test) == 25
    c2 <- c2 + 1
    cond <- c1 | c2 == 20
    if(!c1){istmat[i,x[i]] <- 0}
    print(c2)
  }
  if(nrow(test) != 25){print("Unable to find solution");break}
  matplot(test[,-1], typ = "l", main = paste("i =",i, "\t", ":", "\t", "N =", sum(test[25,-1] > 10^-10)))
  test[25, -1][test[25, -1] < 10^-10] <- 0
  istmat[i+1,] <- test[25,-1]
  #cat("\t", i, "\t", sum(test[25,-1] > 0), "\n")
}
fin <- Sys.time()
fin - strt

nsp <- apply(istmat[1:i,], 1, function(x) sum(x > 0))
plot(nsp, typ = "o")
tail(istmat[1:i,])
istmat[which.max(nsp),istmat[which.max(nsp),]>0]

barplot(t(istmat[1:100,]))
dim(istmat)
barplot(apply(istmat[1:100,], 1, function(x) x/sum(x)))


res3 <-t(apply(istmat[1:100,], 1, function(x) x/sum(x)))
res3 <- melt(res3)
res3 <- res3[res3$value > 0,]
ggplot(res3, aes(x = Var1, y = value, fill = factor(Var2))) + geom_area(position = "stack")
??area

l <- layout_with_lgl(wdg)
i = 30
ig <- induced_subgraph(wdg, V(wdg)[istmat[i, ] > 0])
plot(ig,
     layout = l[istmat[i, ] > 0,],
     vertex.size = ceiling(istmat[i,istmat[i, ] > 0]), 
     vertex.label = NA, 
     edge.arrow.size = .7
     )
text(0,-1.25, label = paste("Time =", i), cex = 2)

library(animation)
saveGIF(
  {
    for(i in c(rep(1,5), 2:99, rep(100, 5))){
      ig <- induced_subgraph(wdg, V(wdg)[istmat[i, ] > 0])
      plot(ig,
           layout = l[istmat[i, ] > 0,],
           vertex.size = ceiling(istmat[i,istmat[i, ] > 0]), 
           vertex.label = NA, 
           edge.arrow.size = .7, 
           main = paste(sum(istmat[i,] > 0), "SPECIES")
      )
      text(0,-1.25, label = paste("Time =", i), cex = 2)
    }
  }, 
  movie.name = "net-assem2.gif", interval = 0.2, ani.width = 600, ani.height = 600, outdir = "~/Desktop/"
)

l2 <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
allv <- apply(istmat, 2, function(x) sum(x) > 0)
ig <- induced_subgraph(wdg, V(wdg)[allv])
saveGIF(
  {
    for(i in c(rep(1,5), 2:99, rep(100, 5))){
      ig2 <- induced_subgraph(wdg, V(wdg)[istmat[i, ] > 0])
      active <- allv[allv]
      active[!istmat[i,allv] > 0] <- FALSE 
      V(ig)$color <- ifelse(active, "green4", "gray")
      V(ig)$size <- ceiling(istmat[i,allv])+1
      E(ig)$color = "gray"
      plot(ig,
           layout = l2[allv,]*1,
           vertex.label = NA, 
           edge.arrow.size = .7, 
           vertex.frame.color = "black",
           rescale = F,
           margin = c(0,0,0,0),
           main = paste(sum(istmat[i,] > 0), "SPECIES")
      )
      plot(ig2,
           layout = l2[allv,][active,]*1,
           vertex.size = V(ig)$size[active],
           vertex.color = "green4",
           edge.color = "slategrey",
           edge.arrow.size = .7,
           vertex.label = NA,
           rescale = F,
           margin = c(0,0,0,0),
           add = T)
      text(0,-1.25, label = paste("Time =", i), cex = 2)
    }
  }, 
  movie.name = "net-assem2a.gif", interval = 0.2, ani.width = 600, ani.height = 600, outdir = "~/Desktop/"
)

allgen <- unique(genera[apply(istmat, 2, function(x) sum(x) > 0)])
bar1 <- data.frame(allgen, num = rep(0, length(allgen)))
saveGIF(
  {
    for(i in c(rep(1,5), 2:99, rep(100, 5))){
      bar1$num <- 0
      t1 <- data.frame(table(genera[istmat[i,] > 0]))
      bar1$num[allgen %in% t1$Var1] <- t1$Freq  
      
      p <- ggplot(bar1, aes(x = allgen, y = num)) + 
        geom_bar(stat = "identity") + 
        xlab("Genus") + ylab("Number of Species") + 
        ggtitle(paste("Time = ", i)) +
        theme_bw() +
        theme(text = element_text(size=20),axis.text.x = element_text(size = 10,angle = 90, hjust = 1))
      print(p)
    }
  }, 
  movie.name = "gen-assem.gif", interval = 0.2, ani.width = 600, ani.height = 600, outdir = "~/Desktop/"
)
