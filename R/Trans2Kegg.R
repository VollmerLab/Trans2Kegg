rm(list = ls())

#' @import tidyr
#' @import KEGGREST
#' @import Biostrings
#' @import annotate
#' @export
annotateTranscripts <- function(accessions, refTransFile){
  refTrans <- readDNAStringSet(refTransFile)
  dfUniKegg <- data.frame()
  rowNum <- 0
  for (accession in accessions){
    transSeq <- refTrans[accession]
    seq <- as.character(transSeq)
    blastResult <- blastSequences(paste0(">",accession,"\n", seq),
      database="swissprot",
      program="blastx",as="data.frame", expect=1e-5)
    if(length(blastResult) > 0){
      blastResult <- subset(blastResult, select=c("Iteration_query-def",
        "Iteration_query-len", "Hit_accession", "Hit_len", "Hsp_evalue",
        "Hsp_identity", "Hsp_gaps", "Hsp_align-len"))
      print(blastResult)
      uniprots <- unique(blastResult$Hit_accession)
      for (uniprot in uniprots){
        uniKegg <- keggConv("genes", paste0("up:", uniprot))
        if(length(uniKegg) > 0){
          newRow <- data.frame("uniprot"=uniprot, "kegg"=uniKegg)
          dfUniKegg <- rbind(dfUniKegg, newRow)
          for(kegg in newRow$kegg){
            koDetail <- tryCatch(keggGet(kegg), error=function(e) print(kegg))
            ortho <- koDetail[[1]]$ORTHOLOGY
            if(length(names(ortho)) > 0){
              rowNum <- rowNum + 1
              newOrtho <- data.frame(row.names=kegg, "kegg"=kegg, "ko"=names(ortho),
                "desc"=ortho[names(ortho)])
              uniToKO <- merge(newRow, newOrtho)
              #print(uniToKO)
              blastToKO <- merge(uniToKO,blastResult, by.x="uniprot",
                by.y="Hit_accession")
              write.table(blastToKO, file="blastToKO.csv", col.names=(rowNum ==1),
                append=(rowNum != 1), row.names=FALSE)
            }
          }
        }
      }
    }
  }
}

#accessions <- c("KXJ29317.1", "KXJ29331.1", "KXJ29349.1", "KXJ29355.1", "KXJ29295.1", "KXJ29270.1")

#annotateTranscriptss(accessions, "kTranscripts.fa")

