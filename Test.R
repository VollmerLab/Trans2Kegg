library(Biostrings)
library(Trans2Kegg)
out <- "annot.csv"
fasta <- "inst/extdata/aiptasia.fa"
transDone <- data.frame()
# Check to see if BLAST already started.
if(file.exists(out)) {
  transDone <-
    read.csv(out, stringsAsFactors = FALSE)$Iteration_query.def
}
ids <- c("KXJ29331.1")
# If BLAST already started, what's left?
transLeft <- setdiff(ids, transDone)
# Read transcriptome FASTA file
refTrans <- readDNAStringSet(fasta)
for(accession in transLeft){
  # Get sequences for transcripts remaining to be BLASTed
  transSeqs <- refTrans[accession]
  print(transSeqs)
}

instance <- 'i-07da948c2d85b7388'
dns <- 'ec2-54-175-9-203.compute-1.amazonaws.com'
filepath <- system.file("extdata", "aiptasia.fa", package="Trans2Kegg")

annotateAws(c("KXJ29331.1"), filepath,"annot.csv", instance=instance, dns=dns, threads=2)

