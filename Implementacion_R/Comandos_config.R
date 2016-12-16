setwd("C:/Users/Ana/Documents/Usach/Bioinformatica/Proyecto_filogenia/Proyecto 1")

install.packages("RInside")
install.packages("Rcpp")
install.packages("RJSONIO")

install.packages("Rserve")
install.packages("RSclient")

# subir servicio
require("Rserve")
Rserve()
run.Rserve()

#bajar servicio
require("RSclient")
c <- RSconnect()
RSshutdown(c)

source("sources/system/packages.R") 
source("sources/system/import2.R")

setwd('C:/Users/Ana/Documents/Usach/Bioinformatica/Proyecto_filogenia/Proyecto 2')









library("phangorn");library("ape");library("RJSONIO");
data=as.phyDat(read.dna("http://localhost:8081/sequence.seq"));
toJSON(nj(dist.hamming(data)))




