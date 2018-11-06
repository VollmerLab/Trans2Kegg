#' Get KEGG pathway information for differentially expressed (DE) transcripts
#' @name annotateDE
#' @import DESeq2
#' @importFrom utils write.csv
#' @param ddsAll DESeqDataSet to annotate
#' @param transcriptFasta FASTA file of transcripts use for abundance 
#' estimation
#' @return Annotation results are written to dfAll.csv, annot.csv, dfPaths.csv, 
#' dfPathsKos.csv, dfRefs.csv
#' @examples
#' \dontrun{
#' annotateDE(dds, "transcripts.fa")
#' }
#' @export
dfAll <- data.frame()
annotateDE <- function(ddsAll, transcriptFasta){
    for (result in resultsNames(ddsAll)){
        if(result != 'Intercept'){
            res <- results(ddsAll, alpha=.05, name=result)
            dfRes <- as.data.frame(res)
            dfRes <- subset(subset(dfRes, padj < .01 & abs(log2FoldChange) > 1,
                           select=c(log2FoldChange, padj)))
            dfRes$Factor <- result
            dfAll <- rbind(dfAll, dfRes)
        }
    }
    write.csv(dfAll, file="dfAll.csv")
    ids <- rownames(dfAll)
    annotateTranscripts(ids, transcriptFasta, "annot.csv")
}
