#!/bin/bash
nice -n 19 blastx -db ~/blastDB/uniprot_sprot.fasta \
-outfmt "6 qseqid sseqid pident evalue qcovs" \
-query OGS2_20140407_transNoSpaceHeader2.fa \
-num_threads 32 -evalue .0001 \
-out blastSP.tsv 1>blastx.log 2>blastx.err &
