blastTrans <- function(transFasta){
  spFile1 <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/"
  spFile2 <- "knowledgebase/complete/uniprot_sprot.fasta.gz"
  spFile <- paste0(spFile1, spFile2)
  method<-"auto"
  filename <- "sp.fasta.gz"
  download.file(spFile, filename, method=method, quiet = FALSE, mode = "w",
    cacheOK = TRUE)
  gunzip(filename, destname = gsub("[.]gz$", "", filename), overwrite = TRUE,
    remove = TRUE, BFR.SIZE = 1e+07)
  makeDbCmd <- "makeblastdb -dbtype 'prot' -in sp.fasta"
  system(makeDbCmd)
  #transFasta <- 'OGS2_20140407_transNoSpaceHeader2.fa'
  blastCmd1 <- 'blastx -db sp.fasta -outfmt '
  blastCmd2 <- '"6 qseqid sseqid pident evalue qcovs" -query '
  blastCmd3 <- ' -num_threads 4 -evalue .00001 -out blastSP.tsv '
  blastCmd4 <- '1>blastx.log 2>blastx.err'
  blastCmd <- paste0(blastCmd1,blastCmd2, transFasta, blastCmd3, blastCmd4)
  system(blastCmd)
}
