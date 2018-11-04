library(Trans2Kegg)
dfVibrio <- read.csv("dfVibrio.csv")
head(dfVibrio)

annotateTranscripts(unique(dfVibrio$X), "kTranscripts.fa", "annot.csv")
