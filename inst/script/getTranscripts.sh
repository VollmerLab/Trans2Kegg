#!/bin/bash
# The test transcriptome GCA_001417965.1_Aiptasia_transcripts.fa was generated
# using gffread (http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread)
# The genomic and GFF files were retrieved from
# http://aiptasia.reefgenomics.org/download/
# Symbiotic sea anemone genome
# Sebastian Baumgarten, Oleg Simakov, Lisl Y. Esherick, Yi Jin Liew,
# Erik M. Lehnert, Craig T. Michell, Yong Li, Elizabeth A. Hambleton,
# Annika Guse, Matt E. Oates, Julian Gough, Virginia M. Weis, Manuel Aranda,
# John R. Pringle, Christian R. Voolstra
# Proceedings of the National Academy of Sciences Aug 2015,
# 201513318; DOI: 10.1073/pnas.1513318112

# Get the Aiptasia genome GFF file
wget http://aiptasia.reefgenomics.org/\
download/GCA_001417965.1_Aiptasia_genome_1.1_genomic.gff.gz

# Get the Aiptasia genome scaffolds
wget http://aiptasia.reefgenomics.org/\
download/GCA_001417965.1_Aiptasia_genome_1.1_genomic.fna.gz

# Unzip the files
gunzip GCA_001417965.1_Aiptasia_genome_1.1_genomic.gff.gz
gunzip GCA_001417965.1_Aiptasia_genome_1.1_genomic.fna.gz

# Extract CDS as FASTA file
gffread GCA_001417965.1_Aiptasia_genome_1.1_genomic.gff \
-g GCA_001417965.1_Aiptasia_genome_1.1_genomic.fna \
-x GCA_001417965.1_Aiptasia_transcripts.fa -C

# Take a subset to keep file size within package limits
head -n30000 GCA_001417965.1_Aiptasia_transcripts.fa \
> ../extdata/transcriptSample.fa

