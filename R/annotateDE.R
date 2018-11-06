#' Get KEGG pathway information for differentially expressed (DE) transcripts
#' @name annotateDE
#' @import DESeq2
#' @importFrom utils write.csv
#' @param ddsAll DESeqDataSet to annotate
#' @param fasta FASTA file of transcripts use for abundance 
#' estimation
#' @return Annotation results are written to dfAll.csv and annot.csv
#' @examples
#' library(DESeq2)
#' library(tximport)
#' trans2gene <- system.file("extdata", "tx2gene.csv", package="Trans2Kegg")
#' samples <- system.file("extdata", "Samples.csv", package="Trans2Kegg")
#' quants <- system.file("extdata", "quants/", package="Trans2Kegg")
#' tx2gene <- read.csv(trans2gene)
#' samplesList <- list.files(quants)
#' files <- file.path(quants, samplesList, "quant.sf")
#' colData <- read.csv(samples, row.names=1)
#' colData$Menthol <- gsub("Control","NoMenthol",colData$Menthol)
#' colData$Vibrio <- gsub("Control","NoVibrio",colData$Vibrio)
#' colData <- colData[with(colData, order(row.names(colData))), ]
#' counts.salmon <- tximport(files, type = "salmon", tx2gene=tx2gene)
#' txi.salmon <- tximport(files, type = "salmon", tx2gene=tx2gene)
#' ddsAll <- DESeqDataSetFromTximport(txi.salmon,
#'     colData = colData,
#'     design = ~ Vibrio + Menthol + Menthol:Vibrio)
#' ddsAll$Menthol <- relevel(ddsAll$Menthol, ref="NoMenthol")
#' ddsAll$Vibrio <- relevel(ddsAll$Vibrio, ref="NoVibrio")
#' ddsAll <- ddsAll[ rowSums(counts(ddsAll)) > 10, ]
#' ddsAll <- DESeq(ddsAll)
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
