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
#' \dontrun{
#' filepath <- system.file("extdata", "aiptasia.fa", package="Trans2Kegg")
#' annotateTranscripts(c("KXJ29317.1"), filepath,
#'     "annot.csv")
#' }
#' @export
annotateTranscripts <- function(ids, fasta, out="annot.csv", evalue="1e-5"){
    refTrans <- readDNAStringSet(fasta)
    rowNum <- 0
    accDone <- c()
    dfUniKegg <- data.frame()
    if (file.exists(out)) {
        annot <- read.csv(out, stringsAsFactors = FALSE)
        rowNum <- nrow(annot)
        accDone <- unique(annot$Iteration_query.def)
    }
    accRemaining <- setdiff (ids, accDone)
    for (accession in accRemaining){
        transSeq <- refTrans[accession]
        seq <- as.character(transSeq)
        #blastResult <- blastSequences(paste0(">",accession,"\n", seq),
        #        database="swissprot", program="blastx",as="data.frame", 
        #        expect=evalue)
        blastResult <- tryCatch(blastSequences(paste0(">",accession,"\n", seq),
                database="swissprot", program="blastx",as="data.frame", 
                expect=evalue), error=function(e) print(accession))
        if(is.data.frame(blastResult) & length(blastResult) > 0){
            blastResult <- subset(blastResult, select=c("Iteration_query-def",
                "Iteration_query-len", "Hit_accession", "Hit_len", 
                "Hsp_evalue", "Hsp_identity", "Hsp_gaps", "Hsp_align-len"))
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
                            if(!is.na(kegg)){
                            rowNum <- rowNum + 1
                            newOrtho <- tryCatch(data.frame(row.names=kegg,
                                "kegg"=kegg, "ko"=names(ortho),
                                "desc"=ortho[names(ortho)]), 
                                error=function(e) print(ortho))
                            uniToKO <- merge(newRow, newOrtho)
                            write.csv(uniToKO, file="uniToKO.csv")
                            blastToKO <- merge(uniToKO,blastResult, 
                                by.x="uniprot", by.y="Hit_accession")
                            write.table(blastToKO, file=out,
                                col.names=(rowNum ==1), append=(rowNum !=1), 
                                sep=',', row.names=FALSE)
                            }else{
                                print(ortho)
                                return(-1)
                            }
                        }
                    }
                }
            }
        }else{
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

