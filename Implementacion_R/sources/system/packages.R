
#====================================================
#A)Configuration
#====================================================
#Packages
cat("- Loading packages...\n")
if (!require("phangorn")) install.packages("phangorn")
library("phangorn")
if (!require("ape")) install.packages("ape")
library("ape")
if (!require("Rphylip")) install.packages("Rphylip")
library("Rphylip")
if (!require("RJSONIO")) install.packages("RJSONIO")
library("RJSONIO")