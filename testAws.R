#!/usr/bin/env Rscript
# test.R
library(tidyr)
library(RCurl)
library(XML)
library(Biostrings)
library(utils)
library(KEGGREST)
library(Trans2Kegg)

instance <- 'i-07da948c2d85b7388'
dns <- 'ec2-54-175-9-203.compute-1.amazonaws.com'
filepath <- system.file("extdata", "aiptasia.fa", package="Trans2Kegg")
annotateAws(c("KXJ29331.1"), filepath,"annot.csv", instance=instance, 
            dns=dns, threads=2)

