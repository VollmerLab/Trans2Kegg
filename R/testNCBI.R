#!/usr/bin/env Rscript
# testNCBI.R
library(tidyr)
library(RCurl)
library(XML)
library(Biostrings)
library(utils)
library(KEGGREST)
library(annotate)

source("getKegg.R")
source("annotateNCBI.R")

filepath <- system.file("extdata", "aiptasia.fa", package="Trans2Kegg")
annotateNCBI(c("KXJ29331.1"), filepath,"annot.csv")

