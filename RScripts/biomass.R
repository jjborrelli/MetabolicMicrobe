outp <- readLines("~/Desktop/GitHub/cometsUtil/Data/biomass.txt")
sapply(outp[1:2], split, "_")
grep(x = outp[1:2], pattern =  "_*_")

sapply(outp[1:2], function(x) grep("_*_", x))

outp[(1+grep("sparse", outp)[1]):(grep("sparse", outp)[2]-1)]
outp[(1+grep("sparse", outp)[2]):(grep("sparse", outp)[3]-1)]
head(outp[(1+grep("sparse", outp)[3]):(grep("sparse", outp)[4]-1)])
plot(c(1e-7, 9.6130185655e-8, sum(9.2645345082e-8, 4.6626846156e-9, 4.6626846156e-9, 4.6626846156e-9, 4.6626846156e-9)))

test <- (outp[(1+grep("sparse", outp)[2]):(grep("sparse", outp)[3]-1)])
strsplit(test, "=") %>% sapply(function(x) as.numeric(strsplit(x[2], ";")[1])) %>% sum()
pos1 <- gregexpr('_', test)
ts <- sapply(1:length(pos1), function(x) substr(test[x], pos1[[x]][1]+1, pos1[[x]][2]-1)) %>% as.numeric # substr(test, pos1[1]+1, pos1[2]-1)
pos2 <- gregexpr("\\(", test)
pos3 <- gregexpr("\\)", test)
posco <- strsplit(substr(test, pos2[[1]]+1, pos3[[1]]-1), ",")[[1]] %>% as.numeric


spatdat <- list()
for(i in 1:(length(grep("sparse", outp))-1)){
  test <- (outp[(1+grep("sparse", outp)[i]):(grep("sparse", outp)[i+1]-1)])
  bio1 <- strsplit(test, "=") %>% sapply(function(x) as.numeric(strsplit(x[2], ";")[1])) 
  pos1 <- gregexpr('_', test)
  ts <- sapply(1:length(pos1), function(x) substr(test[x], pos1[[x]][1]+1, pos1[[x]][2]-1)) %>% as.numeric
  pos2 <- gregexpr("\\(", test)
  pos3 <- gregexpr("\\)", test)
  posco <- sapply(1:length(pos1), function(x) strsplit(substr(test[x], pos2[[x]]+1, pos3[[x]]-1), ",")) %>% sapply(as.numeric) %>% t()
  spatdat[[i]] <- cbind(posco, bio1, ts)
  colnames(spatdat[[i]]) <- c("x", "y", "biomass", "time")
  print(i)
}

bmat <- lapply(1:200, function(x) matrix(0, nrow = 50, ncol = 50))
for(i in 1:length(spatdat)){
  nz <- which(spatdat[[i]][,3] > 0)
  bmat[[i]][spatdat[[i]][nz,1],spatdat[[i]][nz,2]] <- spatdat[[i]][nz,3]
}

plot(sapply(bmat, function(x) max(apply(x, 2, function(y) sum(y > 0))))~sapply(spatdat, function(x) sum(x[,3])))

ggplot(as.data.frame(spatdat[[200]]), aes(x = x, y = y, fill = biomass)) + geom_tile()

outp2 <- read.table("~/Desktop/GitHub/cometsUtil/Data/total_biomass3spp.txt", sep = "\t", col.names = c("Time", "Sp1", "Sp2", "Sp3"))
matplot(outp2[,-1], typ = "l", lwd = 5)
