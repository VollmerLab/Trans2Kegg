#' Get KEGG info for uniprot ID
#' @name getKegg
#' @param uniprot Uniprot ID to use in KEGG query
#' @param blastResult BLAST result to merge with KEGG
#' @param out Output file
#' @return KEGG results are written to the csv file specified by outFile
getKegg <- function(uniprot, blastResult, out){
  # Get KEGG ID for SwissProt ID
  print(uniprot)
  uniKegg <- keggConv("genes", paste0("up:", uniprot))
  if(length(uniKegg) > 0 & (length(uniprot) == length(uniKegg))){
    newRow <- data.frame("uniprot"=uniprot, "kegg"=uniKegg)
    for(kegg in newRow["kegg"]){
      # Get KEGG orthologs for species-specific KEGG ID
      koDetail <- tryCatch(keggGet(kegg),
                           error=function(e) print(kegg))
      ortho <- koDetail[[1]]["ORTHOLOGY"]
      message(ortho)
      if(length(names(ortho)) > 0){
        if(!is.na(kegg)){
          rowNum <- rowNum + 1
          newOrtho <- tryCatch(data.frame(row.names=kegg,
                                          "kegg"=kegg, "ko"=names(ortho),
                                          "desc"=ortho[names(ortho)]),
                               error=function(e) print(ortho))
          uniToKO <- merge(newRow, newOrtho)
          write.csv(uniToKO, file="uniToKO.csv")
          # Merge othologs and BLAST rows
          blastToKO <- merge(uniToKO,blastResult,
                             by.x="uniprot", by.y="Hit_accession")
          # Append BLAST and KEGG info to out file.
          write.table(blastToKO, file=out,
                      col.names=(rowNum ==1), append=(rowNum !=1),
                      sep=',', row.names=FALSE)
        }else{
          message("kegg was NA")
        }
      }
    }
  }
}
#' Annotate selected transcripts from FASTA file
#' @name annotateTranscripts
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
#' annotateTranscripts(c("KXJ29331.1"), filepath,
#'     "annot.csv")
#' @export
annotateTranscripts <- function(ids, fasta, out="annot.csv", evalue="1e-5"){
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
    # Write BLAST results for debug purposes
    if(is.data.frame(blastResult) & length(blastResult) > 0){
      blastResult <- subset(blastResult, select=c("Iteration_query-def",
                                                  "Iteration_query-len", "Hit_accession", "Hit_len",
                                                  "Hsp_evalue", "Hsp_identity", "Hsp_gaps", "Hsp_align-len"))
      write.csv(blastResult, file="blastResult.csv")
      uniprots <- unique(blastResult["Hit_accession"])
      # Loop through SwissProt (uniprot) IDs to get KEGG /IDs
      for (uniprot in uniprots){
        getKegg(uniprot, blastResult, out)
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
