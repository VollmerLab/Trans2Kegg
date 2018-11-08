#' Get KEGG information for DE transcripts
#' @name annotateDE
#' @import DESeq2
#' @importFrom utils write.csv
#' @param ddsAll DESeqDataSet to annotate
#' @param transcriptFasta FASTA file of transcripts use for abundance 
#' estimation
#' @return Annotation results are written to dfAll.csv and annot.csv
#' @examples
#' library(DESeq2)
#' library(repmis)
#' source_data(paste0("https://github.com/croesel/",
#'     "SampleData/blob/master/salmon.RData?raw=true"))
#' samples <- system.file("extdata", "Samples.csv", package="Trans2Kegg")
#' transcripts <- system.file("extdata", "aiptasia.fa", package="Trans2Kegg")
#' colData <- read.csv(samples, row.names=1)
#' colData$Menthol <- gsub("Control","NoMenthol",colData$Menthol)
#' colData$Vibrio <- gsub("Control","NoVibrio",colData$Vibrio)
#' colData <- colData[with(colData, order(row.names(colData))), ]
#' ddsAll <- DESeqDataSetFromTximport(counts.salmon,
#'     colData = colData,
#'     design = ~ Vibrio + Menthol + Menthol:Vibrio)
#' ddsAll$Menthol <- relevel(ddsAll$Menthol, ref="NoMenthol")
#' ddsAll$Vibrio <- relevel(ddsAll$Vibrio, ref="NoVibrio")
#' ddsAll <- ddsAll[ rowSums(counts(ddsAll)) > 10, ]
#' ddsAll <- DESeq(ddsAll)
#' annotateDE(ddsAll,transcripts) 
#' @export

annotateDE <- function(ddsAll, transcriptFasta){
    dfAll <- data.frame()
    for (result in resultsNames(ddsAll)){
        print(result)
        if(result != 'Intercept'){
            res <- results(ddsAll, alpha=.05, name=result)
            print("got res")
            dfRes <- as.data.frame(res)
            print("got dfRes")
            dfRes <- subset(subset(dfRes, padj < .01 & abs(log2FoldChange) > 1,
                select=c(log2FoldChange, padj)))
            print("got subset")
            dfRes$Factor <- result
            print("set factor")
            dfAll <- rbind(dfAll, dfRes)
            print("rbind")
        }
    }
    write.csv(dfAll, file="dfAll.csv")
    #print("wrote csv")
    #ids <- rownames(dfAll)
    #print("got ids")
    #annotateTranscripts(ids, transcriptFasta, "annot.csv")
    #print("called annotateTranscripts")
}
