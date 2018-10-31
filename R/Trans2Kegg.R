#' Annotate selected transcripts from FASTA file
#' @name annotateTranscripts
#' @import tidyr
#' @import KEGGREST
#' @import Biostrings
#' @import annotate
#' @importFrom utils read.csv write.csv write.table
#' @param accessions A character vector of accessions from FASTA file to
#' be annotated
#' @param refTransFile FASTA file of transcripts to be annotated
#' @param outFile csv output file for annotation results
#' @examples
#' filepath <- system.file("extdata", "aiptasia.fa", package="Trans2Kegg")
#' annotateTranscripts(c("KXJ29317.1", "KXJ29331.1"), filepath,
#'  "annot.csv")
#' @export
annotateTranscripts <- function(accessions, refTransFile, outFile){
  refTrans <- readDNAStringSet(refTransFile)
  rowNum <- 0
  accDone <- c()
  dfUniKegg <- data.frame()
  if (file.exists(outFile)) {
    annot <- read.csv(outFile, stringsAsFactors = FALSE)
    rowNum <- nrow(annot)
    accDone <- unique(annot$Iteration_query.def)
  }
  accRemaining <- setdiff (accessions, accDone)
  for (accession in accRemaining){
    transSeq <- refTrans[accession]
    seq <- as.character(transSeq)
    blastResult <- tryCatch(blastSequences(paste0(">",accession,"\n", seq),
      database="swissprot",
      program="blastx",as="data.frame", expect=1e-5),
      error=function(e) print(accession))
    if(is.data.frame(blastResult) & length(blastResult) > 0){
      blastResult <- subset(blastResult, select=c("Iteration_query-def",
        "Iteration_query-len", "Hit_accession", "Hit_len", "Hsp_evalue",
        "Hsp_identity", "Hsp_gaps", "Hsp_align-len"))
      write.csv(blastResult, file="blastResult.csv")
      uniprots <- unique(blastResult$Hit_accession)
      for (uniprot in uniprots){
        uniKegg <- keggConv("genes", paste0("up:", uniprot))
        if(length(uniKegg) > 0){
          newRow <- data.frame("uniprot"=uniprot, "kegg"=uniKegg)
          dfUniKegg <- rbind(dfUniKegg, newRow)
          for(kegg in newRow$kegg){
            koDetail <- tryCatch(keggGet(kegg),
              error=function(e) print(kegg))
            ortho <- koDetail[[1]]$ORTHOLOGY
            if(length(names(ortho)) > 0){
              rowNum <- rowNum + 1
              newOrtho <- data.frame(row.names=kegg,
                "kegg"=kegg, "ko"=names(ortho),
                "desc"=ortho[names(ortho)])
              uniToKO <- merge(newRow, newOrtho)
              write.csv(uniToKO, file="uniToKO.csv")
              #print(uniToKO)
              blastToKO <- merge(uniToKO,blastResult, by.x="uniprot",
                by.y="Hit_accession")
              write.table(blastToKO, file=outFile,
                col.names=(rowNum ==1),
                append=(rowNum !=1), sep=',', row.names=FALSE)
            }
          }
        }
      }
    }else{
      emptyRow <- data.frame(uniprot=NA, kegg = NA, ko=NA, desc=NA,
        Iteration_query.def=accession, Iteration_query.len=NA, Hit_len=NA,
        Hsp_evalue=NA, Hsp_identity=NA, Hsp_gaps = NA, Hsp_align.len=NA)
      write.table(emptyRow, file=outFile,
        col.names=(rowNum ==1),
        append=(rowNum !=1), sep=',', row.names=FALSE)
    }
  }
}

