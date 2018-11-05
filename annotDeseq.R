library(Trans2Kegg)
library(DESeq2)
load("ddsAll.RData")
dfAll <- data.frame()
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
annotateTranscripts(ids, "kTranscripts.fa", "annot.csv")

