library(DESeq2)
library(tximport)
library(Biostrings)

dna <- readDNAStringSet("kTranscripts.fa")
tx2gene <- data.frame(Transcript=names(dna), Gene=names(dna)) 

dir <- 'quants/'
samples <- list.files(dir)
files <- file.path(dir, samples, "quant.sf")

colData <- read.csv("Samples.csv", row.names=1)
colData$Menthol <- gsub("Control","NoMenthol",colData$Menthol)
colData$Vibrio <- gsub("Control","NoVibrio",colData$Vibrio)
colData <- colData[with(colData, order(row.names(colData))), ]

counts.salmon <- tximport(files, type = "salmon", tx2gene=tx2gene)

txi.salmon <- tximport(files, type = "salmon", tx2gene=tx2gene)
ddsAll <- DESeqDataSetFromTximport(txi.salmon,
                                     colData = colData,
                                     design = ~ Vibrio + Menthol + Menthol:Vibrio)
ddsAll$Menthol <- relevel(ddsAll$Menthol, ref="NoMenthol")
ddsAll$Vibrio <- relevel(ddsAll$Vibrio, ref="NoVibrio")
##Remove low counts
ddsAll <- ddsAll[ rowSums(counts(ddsAll)) > 10, ]
##Run DESeq
ddsAll <- DESeq(ddsAll)
save(ddsAll, file="ddsAll.RData")
dfAll <- data.frame()
for (result in resultsNames(ddsAll)){
    if(result != 'Intercept'){
    res <- results(ddsAll, alpha=.05, name=result)
    dfRes <- as.data.frame(res)
    dfRes <- subset(subset(dfRes, padj < .05, select=c(log2FoldChange, padj)))
    dfRes$Factor <- result
    dfAll <- rbind(dfAll, dfRes)
    }
}
write.csv(dfAll, file="dfAll.csv")
