install.packages("stringr")
library(stringr)


# Algoritmo Smith-Waterman

prueba <- function (s1,s2)
{
  print(0)
}

smith <- function (s1,s2)
{
  matriz = matrix(nrow = nchar(s1) + 1 , ncol = nchar(s2)+1)
  s1 <- paste("",s1)
  s2 <- paste("",s2)
  rownames(matriz) <- c(unlist(strsplit(s1, "")))
  colnames(matriz) <- c(unlist(strsplit(s2, "")))
  inicializar_matriz(matriz)
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
    print(matriz)
}



highScore <- function (v1,v2,v3)
{
  pmax(0,v1,v2,v3)
}





smith('hola','hola')






