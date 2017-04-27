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

pimat <- matrix(nrow = 773, ncol = 773)
colnames(pimat) <- single.gr$AGORA.Strain
rownames(pimat) <- single.gr$AGORA.Strain

res <- list()
res2 <- list()
for(i in 1:nrow(single.gr)){
  as1 <- (st1 %in% strA[i])
  gch <- ((pairw.gr$WesternDietAnaerobic_Strain1[as1] - single.gr$Western.diet..anaerobic..h.1.[i])/single.gr$Western.diet..anaerobic..h.1.[i])*100
  #pimat[i,-i] <- gch
  res[[i]] <- gch
  
  as2 <- st2 %in% strA[i]
  gch2 <- ((pairw.gr$WesternDietAnaerobic_Strain2[as2] - single.gr$Western.diet..anaerobic..h.1.[i])/single.gr$Western.diet..anaerobic..h.1.[i])*100
  res2[[i]] <- gch2
  print(i)
}
plot(sapply(res, length), typ = "l")
plot(sapply(res2, length), typ = "l")
which(sapply(res, length) == 0)

i = 94
strA[i]
unique(st1)[i]
