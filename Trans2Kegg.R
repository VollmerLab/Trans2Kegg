## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(Trans2Kegg)
library(Biostrings)

## ------------------------------------------------------------------------
transFile <- system.file("extdata", "transcriptSample.fa", package="Trans2Kegg")
transIds <- c("rna3 gene=rfx6")
library(Trans2Kegg)
annotateTranscripts(transIds, transFile, "annot.csv")

