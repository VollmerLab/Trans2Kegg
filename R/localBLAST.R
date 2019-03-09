#' Perform a local BLAST of DE transcripts
#' @name localBlast
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom utils read.csv write.csv write.table
#' @param accessions A character vector of accessions from FASTA file to
#' be annotated
#' @param refTransFile FASTA file of transcripts to be annotated
#' @param fastaOut Output file for DE sequences
#' @return DE Sequences are written to specified FASTA file
#' @examples
#' \dontrun{
#' filepath <- system.file("extdata", "aiptasia.fa", package="Trans2Kegg")
#' localBlast(c("KXJ29317.1"), filepath)
#' }
#' @export
localBlast <- function(accessions, refTransFile, fastaOut="deTrans.fa"){ 
    refTrans <- readDNAStringSet(refTransFile)
    # Create a FASTA file of selected IDs to BLAST locally
    for (accession in accessions){
        transSeq <- refTrans[accession]
        writeXStringSet(transSeq, fastaOut, append=TRUE, compress=FALSE, 
            compression_level=NA, format="fasta")
    }
}
