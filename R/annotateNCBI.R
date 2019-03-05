#' Annotate selected transcripts from FASTA file
#' @name annotateNCBI
#' @import tidyr
#' @import KEGGREST
#' @importFrom Biostrings readDNAStringSet
#' @import annotate
#' @importFrom utils read.csv write.csv write.table
#' @param ids A character vector of accessions from FASTA file to
#' be annotated
#' @param fasta FASTA file of transcripts to be annotated
#' @param out csv output file for annotation results
#' @param evalue e-value cutoff for BLAST results
#' @return Annotation results are written to the csv file specified by outFile
#' @examples
#' filepath <- system.file("extdata", "aiptasia.fa", package="Trans2Kegg")
#' annotateNCBI(c("KXJ29331.1"), filepath, "annot.csv")
#' @export
annotateNCBI <- function(ids, fasta, out="annot.csv", evalue="1e-5"){
    refTrans <- readDNAStringSet(fasta)
    rowNum <- 0
    accDone <- c()
    dfUniKegg <- data.frame()
    # Check if BLAST already started. If already started, pick up
    # where you left off.
    if (file.exists(out)) {
        annot <- read.csv(out, stringsAsFactors = FALSE)
        rowNum <- nrow(annot)
        accDone <- unique(annot$Iteration_query.def)
    }
    # Determine which sequences have not yet been BLASTed.
    accRemaining <- setdiff (ids, accDone)
    for (accession in accRemaining){
        transSeq <- refTrans[accession]
        seq <- as.character(transSeq)
        # BLAST using NCBI servers
        blastResult <- tryCatch(blastSequences(paste0(">",accession,"\n", seq),
            database="swissprot", program="blastx",as="data.frame",
            expect=evalue), error=function(e) print(accession))
        if(is.data.frame(blastResult) & length(blastResult) > 0){
            blastResult <- subset(blastResult, select=c("Iteration_query-def",
            "Iteration_query-len", "Hit_accession", "Hit_len",
            "Hsp_evalue", "Hsp_identity", "Hsp_gaps", "Hsp_align-len"))
            write.csv(blastResult, file="blastResult.csv")
            uniprots <- unique(blastResult["Hit_accession"])
            # Loop through SwissProt (uniprot) IDs to get KEGG /IDs
            for (uniprot in uniprots){
                getKegg(blastResult, out)
            }
        }else{
            # If no BLAST hits append a default row for the id to prevent
            # re-BLASTing on restart.
            emptyRow <- data.frame(uniprot=NA, kegg = NA, ko=NA, desc=NA,
            Iteration_query.def=accession, Iteration_query.len=NA, Hit_len=NA,
            Hsp_evalue=NA, Hsp_identity=NA, Hsp_gaps = NA, Hsp_align.len=NA)
            write.table(emptyRow, file=out,
            col.names=(rowNum ==1),
            append=(rowNum !=1), sep=',', row.names=FALSE)
        }
    }
    return(0)
}
