setwd("C:/Users/Ana/Documents/Usach/Bioinformatica/Proyecto_filogenia/Proyecto 1") #Cambiar por tu directorio
source("sources/system/packages.R")         #Carga un conjunto de paquetes
cat("\014")                                 #Limpiar consola
library("RJSONIO")

#===================================
#1 Lectura de secuencia
#===================================
data=as.phyDat(read.dna("data/sequence.seq")) #Aplicar read.aa si es aminoácido


data=as.phyDat(read.dna("http://localhost:8081/sequence.seq"))



#===================================
#2 Estimación de árboles - Métodos basados en distancia
#===================================
# 2.1- NJ
#====
  #Hamming
  tree1=nj(dist.hamming(data))
  response <- toString(tree1)
  plot(tree1)
  edgelabels(round(tree1$edge.length,1),cex=0.6,bg = "white")
  
  #Evolutionary model
  #model puede ser "raw", "N", "TS", "TV", "JC69", "K80" (the default), "F81", "K81", "F84",
  #"BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", or "indelblock"
  tree1=nj(dist.dna(as.DNAbin(data),model="F81")) #Aplicar dist.aa si es aminoácido
  plot(tree1)
  edgelabels(round(tree1$edge.length,1),cex=0.6,bg = "white")
  
# 2.2- BioNJ
#====
  #Hamming
  tree1=bionj(dist.hamming(data))
  plot(tree1)
  edgelabels(round(tree1$edge.length,1),cex=0.6,bg = "white")
  
  #Evolutionary model
  #model puede ser "raw", "N", "TS", "TV", "JC69", "K80" (the default), "F81", "K81", "F84",
  #"BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", or "indelblock"
  tree1=bionj(dist.dna(as.DNAbin(data),model="F81"))
  plot(tree1)
  edgelabels(round(tree1$edge.length,1),cex=0.6,bg = "white")

# 2.3- BioNJ
#====
  #Hamming
  tree1=bionj(dist.hamming(data))
  plot(tree1)
  edgelabels(round(tree1$edge.length,1),cex=0.6,bg = "white")
  
  #Evolutionary model
  #model puede ser "raw", "N", "TS", "TV", "JC69", "K80" (the default), "F81", "K81", "F84",
  #"BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", or "indelblock"
  tree1=bionj(dist.dna(as.DNAbin(data),model="F81"))
  plot(tree1)
  edgelabels(round(tree1$edge.length,1),cex=0.6,bg = "white")

# 2.4- UPGMA
#====
  #Hamming
  tree1=upgma(dist.hamming(data))
  plot(tree1)
  edgelabels(round(tree1$edge.length,1),cex=0.6,bg = "white")
  
  #Evolutionary model
  #model puede ser "raw", "N", "TS", "TV", "JC69", "K80" (the default), "F81", "K81", "F84",
  #"BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", or "indelblock"
  tree1=upgma(dist.dna(as.DNAbin(data),model="F81"))
  plot(tree1)
  edgelabels(round(tree1$edge.length,1),cex=0.6,bg = "white")

# 2.5- WPGMA
#====
  #Hamming
  tree1=wpgma(dist.hamming(data))
  plot(tree1)
  edgelabels(round(tree1$edge.length,1),cex=0.6,bg = "white")
  
  #Evolutionary model
  #model puede ser "raw", "N", "TS", "TV", "JC69", "K80" (the default), "F81", "K81", "F84",
  #"BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", or "indelblock"
  tree1=wpgma(dist.dna(as.DNAbin(data),model="F81"))
  plot(tree1)
  edgelabels(round(tree1$edge.length,1),cex=0.6,bg = "white")
  
# 2.6- Fast Minimum evolution 
#====
  #Balanced
  tree1=fastme.bal(dist.hamming(data), nni = TRUE, spr = TRUE, tbr = TRUE)
  #OLS
  tree1=fastme.ols(dist.hamming(data), nni = TRUE)

# 2.7- Least-Squares 
#====
  #Fitch-Margoliash
  tree1=Rfitch(dist.hamming(data), path="sources/exe",quiet=TRUE,method="fm")
  #Least-squares
  tree1=Rfitch(dist.hamming(data), path="sources/exe",quiet=TRUE,method="ls")
  
  
  
  
  
  
  
  

# 2.8 - Evaluación de calidad métodos de distancia - BOOTSTRAP
#===================
  tree1=wpgma(dist.hamming(data))
  NJtrees = bootstrap.phyDat(data, FUN=function(x)wpgma(dist.hamming(x)), bs=100)
  treeNJ <- plotBS(tree1, NJtrees, "phylogram")
  
  
#===================================
#3 Estimación de árboles - Métodos basados en caracteres
#===================================
  # 3.1 - Parsimony
  #===================
  tree1=optim.parsimony(tree1, data)
  NJtrees = bootstrap.phyDat(data, FUN=function(x)optim.parsimony(tree1,x), bs=100)
  tree1=acctran(tree1,data)         
  treeNJ <- plotBS(tree1, NJtrees, "phylogram")
  add.scale.bar()
  
  # 3.2 - Likelihood
  #===================
  #Evolutionary model estimation to Likelihood treee
  cat("- Including initial likelihood optimisation...\n")
  mT = NULL
  mT = modelTest(data,tree1)
  env <- attr(mT, "env")
  ev_tree <- eval(get(mT$Model[which.min(mT$AIC)], env), env)
  
  #Evolutionary model selection
  if ((mT$Model[which.min(mT$AIC)]=="GTR") | (mT$Model[which.min(mT$AIC)]=="GTR+I") | (mT$Model[which.min(mT$AIC)]=="GTR+G")| (mT$Model[which.min(mT$AIC)]=="GTR+G+I")) {ev_tree$model="GTR"}
  if ((mT$Model[which.min(mT$AIC)]=="HKY") | (mT$Model[which.min(mT$AIC)]=="HKY+I") | (mT$Model[which.min(mT$AIC)]=="HKY+G")| (mT$Model[which.min(mT$AIC)]=="HKY+G+I")) {ev_tree$model="HKY"}
  if ((mT$Model[which.min(mT$AIC)]=="JC") | (mT$Model[which.min(mT$AIC)]=="JC+I") | (mT$Model[which.min(mT$AIC)]=="JC+G")| (mT$Model[which.min(mT$AIC)]=="JC+G+I")) {ev_tree$model="JC"}            #Acá poner los otros hacia abajo

  tree1=pml(tree = ev_tree$tree, data = data, bf = ev_tree$bf, Q = ev_tree$Q, inv = ev_tree$inv,k = ev_tree$k, shape = ev_tree$shape)
  
  #Análisis de calidada Bootstrap para likelihood
  bs <- bootstrap.pml(tree1, bs=100, optNni=TRUE)
  treeBS <- plotBS(tree1$tree,bs,"phylogram")
  
  #===================================
  #6 Guardar árboles
  #===================================
  write.tree(tree1,"tree1_saved.tree")
  