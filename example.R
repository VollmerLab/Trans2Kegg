library(DESeq2)
library(repmis)
source_data("https://github.com/croesel/SampleData/blob/master/salmon.RData?raw=true")
samples <- system.file("extdata", "Samples.csv", package="Trans2Kegg")
colData <- read.csv(samples, row.names=1)
colData$Menthol <- gsub("Control","NoMenthol",colData$Menthol)
colData$Vibrio <- gsub("Control","NoVibrio",colData$Vibrio)
colData <- colData[with(colData, order(row.names(colData))), ]
ddsAll <- DESeqDataSetFromTximport(counts.salmon,
     colData = colData,
     design = ~ Vibrio + Menthol + Menthol:Vibrio)
ddsAll$Menthol <- relevel(ddsAll$Menthol, ref="NoMenthol")
ddsAll$Vibrio <- relevel(ddsAll$Vibrio, ref="NoVibrio")
ddsAll <- ddsAll[ rowSums(counts(ddsAll)) > 10, ]
ddsAll <- DESeq(ddsAll)
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
