#!/bin/bash
Rscript rox.R
cd ..
R CMD build Trans2Kegg 1>Trans2Kegg/build.log 2>Trans2Kegg/build.err
R CMD check --no-build-vignettes Trans2Kegg_0.99.0.tar.gz 1>>Trans2Kegg/build.log 2>>Trans2Kegg/build.err
R CMD BiocCheck Trans2Kegg_0.99.0.tar.gz 1>>Trans2Kegg/build.log 2>>Trans2Kegg/build.err
#Rscript install.R
cd Trans2Kegg
