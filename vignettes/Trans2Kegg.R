#!/usr/bin/env Rscript

## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----message=FALSE, warning=FALSE----------------------------------------
# Load the Trans2Kegg library
library(dplyr)
library(Trans2Kegg)
library(knitr)
library(tidyr)

## ----AnnotateDE, eval=FALSE----------------------------------------------
## annotateDE(ddsAll)

## ----Annotate, eval=FALSE------------------------------------------------
## deAll <- read.csv("dfAll.csv", row.names=1)
## # Get the IDs (rownames) from deAll.
## ids <- rownames(deAll)
## # Select a subset of rows for testing purposes
## ids2 <- head(ids, n=10)
## annotateNCBI(ids2, "aiptasia.fa")

## ----AnnotateAWS, eval=TRUE, message=FALSE, warning=FALSE----------------
deAll <- read.csv("dfAll.csv", row.names=1)
# Get the IDs (rownames) from deAll.
ids <- rownames(deAll)
# Select a subset of rows for testing purposes
ids2 <- head(ids, n=200)
# NOTE: instance and DNS show below are examples only. 
# Users must launch their own NCBI BLAST AMI
instance <- 'i-07da948c2d85b7388'
dns <- 'ec2-54-175-9-203.compute-1.amazonaws.com'
annotateAWS(ids2, "aiptasia.fa", instance=instance, dns=dns, threads=2)

## ----mergeAnnotations----------------------------------------------------
mergeAnnotations()

## ----getPaths------------------------------------------------------------
getPathways()

## ----mergePaths----------------------------------------------------------
mergePaths() 

## ----Table1--------------------------------------------------------------
deAll <- read.csv("dfAll.csv", row.names=1)
kable(head(deAll), caption = "Sample of deAll.csv, the output of function 
    AnnotateDE, which combines all comparisons from a DESeqDataSet
    into a single file to facilitate annotation.", booktabs=TRUE)

## ----Table2--------------------------------------------------------------
outFile <- "annot.csv"
annot <- read.csv(outFile, stringsAsFactors=FALSE)
kable(head(annot, n= 20L), 
    caption="Sample annot.csv output of annotateTranscript and annotateAws 
    functions.", booktab=TRUE)

## ----Table3--------------------------------------------------------------
dfCovCount <- read.csv("cvCnt.csv", stringsAsFactors=FALSE)
dfCovCount$Factor <- gsub("_", " ", dfCovCount$Factor)
kable(head(dfCovCount, n=20L), 
    caption="Sample cvCnt.csv output of mergePathways function.
    This table summarizes the annotation results. qCov is the mean query
    coverage for the species BLAST hits matching the KEGG ortholog.
    n is the number of species BLAST hits matching the KEGG ortholog.", 
    booktab=TRUE, row.names=FALSE)

## ----Table4--------------------------------------------------------------
dfPathsKos <- read.csv("pthKo.csv", stringsAsFactors=FALSE)
kable(head(dfPathsKos, n=20L), 
    caption="Sample pthKo.csv output of getPathways function.", 
    booktab=TRUE, row.names=FALSE)

## ----Table5--------------------------------------------------------------
dfPaths <- read.csv("path.csv", stringsAsFactors=FALSE)
kable(head(dfPaths, n=20L), 
    caption="Sample path.csv output of getPathways function.", 
    booktab=TRUE, row.names=FALSE)

## ----Table6--------------------------------------------------------------
countByPath <- read.csv("countByPath.csv", stringsAsFactors=FALSE)
# Replace _ with spaces in Factor to allow word-wrap within the column.
countByPath$Factor <- gsub("_", " ", countByPath$Factor)
# Subset to include only "Organismal Systems" category
countByPathOrg <- subset(countByPath, category %in% 
    c("Organismal Systems"))
countByPathOrg <- arrange(countByPathOrg, class, path, direction)
kable(head(countByPathOrg, n=15L), 
    caption="Sample countByPath.csv output of mergePathways function.", 
    booktab=TRUE, row.names=FALSE)


## ----Table7--------------------------------------------------------------
# Read in countByClass.csv
countByClass <- read.csv("countByClass.csv", stringsAsFactors=FALSE)
# Replace _ with spaces in Factor to allow word-wrap within the column.
countByClass$Factor <- gsub("_", " ", countByClass$Factor)
# Subset to include only "Organismal Systems" category
countByClassOrg <- subset(countByClass, category %in% 
    c("Organismal Systems"))
countByClassOrg <- arrange(countByClassOrg, class, direction)
kable(head(countByClassOrg, n=20L), 
    caption="Sample countByClass.csv output of mergePathways function.", 
    booktab=TRUE, row.names=FALSE)

## ----Table8--------------------------------------------------------------
# Read in dePathsDetails.csv
deDetails <- read.csv("dePathsDetails.csv", stringsAsFactors=FALSE)
# Replace _ with spaces in Factor to allow word-wrap within the column.
deDetails$Factor <- gsub("_", " ", deDetails$Factor)
# Subset genes with Immune system pathways
deDetailsImmune <- subset(deDetails, class %in% c(" Immune system"))
kable(head(deDetailsImmune, n=10L), 
    caption="Sample dePathsDetails.csv output of mergePathways function, 
    subset to show the first ten genes within Immune system.
    The combinations of DE and KEGG information in dePathsDetails.csv 
    facilitates subsetting and summarizing DE results", 
    booktab=TRUE, row.names=FALSE)

