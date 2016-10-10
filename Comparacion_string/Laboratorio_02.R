#biocLite(Biostrings)
library(Biostrings)

# Manipulacion de secuencias
dnaSeq = DNAString("TTCAGATCTAGTTCGTGTGTGACTGATGATCTGTCACACGTTTTTCTGATCTTCTGACTAGTCGAT")

# Ambos entregan el mismo resultado
dnaSeq[1:10]
substr(dnaSeq,1,10)

length(dnaSeq) #Establece el largo de la secuencia
reverseComplement(dnaSeq) #Genera una cadena complementaria reversa
alphabetFrequency(dnaSeq) #Indica la cardinalidad de cada caracter


# Este es solo un string matching (comparativa sencilla)
compareStrings(dnaSeq,dnaSeq)
dnaSeq2 = DNAString("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
compareStrings(dnaSeq,dnaSeq2)


# Carga archivo fasta
dnaSeq = readDNAStringSet("D:\\Proyectos\\Bioinformatica\\primates_14.fasta")

names(dnaSeq) #Conocer el nombre de las especies
names(dnaSeq[1]) #Reconocer una especie
length(dnaSeq) #Saber cuantas especies son
reverseComplement(dnaSeq[1]) #Complemento

alphabetFrequency(dnaSeq) #Contador de caracteres


