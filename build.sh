#!/bin/bash
Rscript rox.R
cd ..
R CMD build Trans2Kegg
cd Trans2Kegg
