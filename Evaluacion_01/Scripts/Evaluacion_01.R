library(Biostrings)

setwd("C:\\Users\\Kari\\Documents\\GitHub\\Bioinformatica\\Evaluacion_01\\Scripts")

dnaSeqUndifenid = readDNAStringSet("secuencia.fasta")
dnaSeqRatt = readDNAStringSet("rattusnorvegicuschromosomeY.fasta")

names(dnaSeqUndifenid)
names(dnaSeqRatt)

length(dnaSeqUndifenid)
length(dnaSeqRatt)

reverseComplement(dnaSeqUndifenid)
reverseComplement(dnaSeqRatt)

alphabetFrequency(dnaSeqUndifenid)
alphabetFrequency(dnaSeqRatt)

dinucleotideFrequency(dnaSeqUndifenid)
dinucleotideFrequency(dnaSeqRatt)
