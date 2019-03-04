#' Get KEGG Ortholog info for BLAST Results
#' @name getKegg
#' @import tidyr
#' @import KEGGREST
#' @importFrom utils read.csv write.csv write.table URLencode
#' @param blastResult BLAST results from blastSequencesAws
#' @param outFile csv output file for annotation results
#' @return Annotation results are written to the csv file specified by outFile
getKegg <- function(blastResult,
    outFile="annot.csv"){
    rowNum <- 0
    accDone <- c()
    #dfUniKegg <- data.frame()
    uniprots <- unique(blastResult["Hit_accession"])
    for (uniprot in uniprots){
        uniKegg <- keggConv("genes", paste0("up:", uniprot))
        if(length(uniKegg) > 0 & (length(uniKegg) == length(uniprot))){
            newRow <- data.frame("uniprot"=uniprot, "kegg"=uniKegg)
            #dfUniKegg <- rbind(dfUniKegg, newRow)
            # Get KEGG orthologs using KEGGREST
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
                        # Start new file with headers if first,
                        # otherwise append.
                        write.table(blastToKO, file=outFile,
                            col.names=!file.exists(outFile),
                            append=file.exists(outFile),
                            sep=',', row.names=FALSE)
                    }else{
                        message("kegg was NA")
                    }
                }
            }
        }
    }
}

