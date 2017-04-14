outp <- readLines("~/Desktop/GitHub/cometsUtil/biomass.txt")
sapply(outp[1:2], split, "_")
grep(x = outp[1:2], pattern =  "_*_")

sapply(outp[1:2], function(x) grep("_*_", x))

outp[(1+grep("sparse", outp)[1]):(grep("sparse", outp)[2]-1)]
outp[(1+grep("sparse", outp)[2]):(grep("sparse", outp)[3]-1)]
head(outp[(1+grep("sparse", outp)[3]):(grep("sparse", outp)[4]-1)])
plot(c(1e-7, 9.6130185655e-8, sum(9.2645345082e-8, 4.6626846156e-9, 4.6626846156e-9, 4.6626846156e-9, 4.6626846156e-9)))
