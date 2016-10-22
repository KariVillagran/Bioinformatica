library(Biostrings)

# Actividad 3.2 Aliniamiento con R

setwd("C:\\Users\\Kari\\Documents\\GitHub\\Bioinformatica\\Evaluacion_01\\Scripts")

seq1<-DNAString("GAATTCCTACTACGAAGAATTCCTACTACGAAACTACGAAAATTCCTACTACGA") 
seq2<-DNAString("GAATTCCTACTACGGAATTCCCCTCCCATAATTCCTACTACGA")

# Alineamiento Global

mat = nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)
globalAlign = pairwiseAlignment(seq1, seq2, type="global", substitutionMatrix = mat,gapOpening=0, gapExtension = 0)
print(globalAlign)

# Alineamiento Local


mat = nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)
localAlign = pairwiseAlignment(seq1, seq2, type="local", substitutionMatrix = mat,gapOpening=0, gapExtension = 2)
print(localAlign)