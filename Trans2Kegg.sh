#!/bin/bash
rm Icon?
Rscript rox.R
cd ..
R CMD build Trans2Kegg
R CMD check --no-build-vignettes Trans2Kegg_0.99.0.tar.gz
R CMD BiocCheck Trans2Kegg_0.99.0.tar.gz
Rscript install.R
cd Trans2Kegg
