library(Biostrings)

setwd("C:\\Users\\Kari\\Documents\\GitHub\\Bioinformatica\\Evaluacion_01\\Scripts")

dnaSeqUndifenid = readDNAStringSet("secuencia.fasta")
dnaSeqRatt = readDNAStringSet("rattusnorvegicuschromosomeY.fasta")


names(dnaSeqUndifenid)
names(dnaSeqRatt)

length(dnaSeqUndifenid$`gi|661902763|gb|CM002824.1| Rattus norvegicus strain partialchromosome Y, partial genome shotgun sequence`)
length(dnaSeqRatt)

reverseComplement(dnaSeqUndifenid)
reverseComplement(dnaSeqRatt)

alphabetFrequency(dnaSeqUndifenid)
alphabetFrequency(dnaSeqRatt)

dinucleotideFrequency(dnaSeqUndifenid)
dinucleotideFrequency(dnaSeqRatt)


# Alineamiento Local

mat = nucleotideSubstitutionMatrix(match=1, mismatch= -1,baseOnly = TRUE)
localAlign = pairwiseAlignment(dnaSeqUndifenid[1],dnaSeqRatt[1],type="local",substitutionMatrix=mat, gapOpening=0, gapExtension=-1)

print(localAlign)
