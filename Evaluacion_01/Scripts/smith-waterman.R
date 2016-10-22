install.packages("stringr")
library(stringr)


# Algoritmo Smith-Waterman


smith <- function (s1,s2)
{
  matriz = matrix(nrow=nchar(s2)+1,ncol=nchar(s1) + 1)
  s1 <- paste("",s1)
  s2 <- paste("",s2)
  colnames(matriz) <- c(unlist(strsplit(s1, "")))
  rownames(matriz) <- c(unlist(strsplit(s2, "")))
  matriz <- inicializar_matriz(matriz)
  matriz <- fill_table(matriz)
  vectorScore <- score(matriz)
  scoreResult <- paste("El score obtenido es de : ",vectorScore[1])
  print(scoreResult)
  #print(matriz)
  vectorString <- reconstructionString(matriz,vectorScore[2],vectorScore[3])
  print(vectorString[1])
  print(vectorString[2])
}

score <- function(matriz)
{
  maxValue <- 0
  iMAxPos <- 0
  jMAxPos <- 0
  for(i in 2:nrow(matriz))
  {
    for(j in 2:ncol(matriz))
    {
      maxValue<- pmax(maxValue,matriz[i,j])
      if(maxValue == pmax(maxValue,matriz[i,j]))
      {
        iMAxPos <- i
        jMAxPos <- j
      }
    }
  }
  
  return(c(maxValue,iMAxPos,jMAxPos))
}


fill_table <- function(matriz)
{
  match <- 1
  gap <- -2
  mmatch <- -1
  puntaje <- 0
  for(i in 2:nrow(matriz))
  {
    for(j in 2:ncol(matriz))
    {
      if(rownames(matriz)[i]==colnames(matriz)[j])
      {
        puntaje <- highScore(matriz[i-1,j]+gap,matriz[i,j-1]+gap,matriz[i-1,j-1]+match)
      }
      else
      {
        puntaje <- highScore(matriz[i-1,j]+gap,matriz[i,j-1]+gap,matriz[i-1,j-1]+mmatch)
      }
      matriz[i,j] <- puntaje
    }
  }
  return(matriz)
}


inicializar_matriz <- function (matriz)
{
  for(i in 1:nrow(matriz))
  {
    matriz[i,1] <- 0
  }
  for(j in 1:ncol(matriz))
  {
    matriz[1,j] <- 0
  }
  return(matriz)
}

highScore <- function (v1,v2,v3)
{
  return(pmax(0,v1,v2,v3))
}

reconstructionString <- function(matriz, imax, jmax)
{
  s1 <- ""
  s2 <- ""
  while((imax > 1 && jmax > 1) && matriz[imax,jmax] > 0 )
  {
    if(rownames(matriz)[imax]==colnames(matriz)[jmax])
    {
      s1 <- paste(colnames(matriz)[jmax],s1,sep="")
      s2 <- paste(rownames(matriz)[imax],s2,sep="")
      imax <- imax-1
      jmax <- jmax-1
    }
    else
    {
      if(matriz[imax-1,jmax] > matriz[imax,jmax-1])
      {
        s2 <- paste(rownames(matriz)[imax],s2,sep="")
        s1 <- paste("-",s1,sep="")
        imax <- imax-1
      }
      else
      {
        s2 <- paste("-",s2,sep="")
        s1 <- paste(colnames(matriz)[jmax],s1,sep="")
        jmax <- jmax-1
      }
    }
  }
  return(c(s1,s2))
}




smith('AAAA','AAAA')

smith('GAATTCCTACTACGAAGAATTCCTACTACGAAACTACGAAAATTCCTACTACGA',
      'GAATTCCTACTACGGAATTCCCCTCCCATAATTCCTACTACGA')

matriz = matrix(nrow = nchar("hola") + 1 , ncol = nchar("hola")+1)
s1 <- paste("","hola")
s2 <- paste("","hola")
rownames(matriz) <- c(unlist(strsplit(s1, "")))
colnames(matriz) <- c(unlist(strsplit(s2, "")))
matriz








