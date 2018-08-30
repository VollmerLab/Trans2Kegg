rm(list = ls())
library(rjson)
library(stringr)
library(KEGGREST)
library(tidyr)
library(splitstackshape)
library(Biostrings)
library(annotate)
library(DESeq2)
library(tximport)

refTransFile <- "kTranscripts.fa"
salmonOutDir <- 'quants/'
refTrans <- readDNAStringSet(refTransFile)

dfUniKegg <- data.frame()
blastResults <- data.frame()
accession <- "KXJ08417.1"
transSeq <- refTrans[accession]
seq <- as.character(transSeq)
blastResult <- blastSequences(paste0(">",accession,"\n", seq),
  database="swissprot",
  program="blastx",as="data.frame", expect=1e-5)
blastResult <- subset(blastResult, select=c("Iteration_query-def",
  "Iteration_query-len", "Hit_accession", "Hit_len", "Hsp_evalue",
  "Hsp_identity", "Hsp_gaps", "Hsp_align-len"))
blastResults <- rbind(blastResults, blastResult)
uniprots <- unique(blastResult$Hit_accession)
for (uniprot in uniprots){
  uniKegg <- keggConv("genes", paste0("up:", uniprot))
  if(length(uniKegg) > 0){
    newRow <- data.frame("uniprot"=uniprot, "kegg"=uniKegg)
    dfUniKegg <- rbind(dfUniKegg, newRow)
    for(kegg in newRow$kegg){
      koDetail <- tryCatch(keggGet(kegg), error=function(e) print(kegg))
      ortho <- koDetail[[1]]$ORTHOLOGY
      names(ortho)
      newOrtho <- data.frame(row.names=kegg, "kegg"=kegg, "ko"=names(ortho),
        "desc"=ortho[names(ortho)])
      keggToKO <- rbind(keggToKO, newOrtho)
      uniToKO <- merge(newRow, newOrtho)
      print(uniToKO)
    }
  }
}

keggToKO <- data.frame()
kegg <- "mmu:70790"

tx2gene <- data.frame(Transcript=names(refTrans), Gene=names(refTrans))
samples <- list.files(salmonOutDir)
files <- file.path(dir, samples, "quant.sf")
colData <- read.csv("Samples.csv", row.names=1)
colData$Menthol <- gsub("Control","NoMenthol",colData$Menthol)
colData$Vibrio <- gsub("Control","NoVibrio",colData$Vibrio)
colData <- colData[with(colData, order(row.names(colData))), ]
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
resultsNames(ddsAll)
res1 <- as.data.frame(results(ddsAll, alpha=.05, name="Vibrio_Vibrio_vs_NoVibrio"))
highlyDE <- subset(res1, abs(log2FoldChange) > 1 & pvalue <= .05)
highlyDEgenes <- unique(rownames(highlyDE))
head(highlyDE)


if (!file.exists("cache.RData")) {
  data1 <- fromJSON(file="http://www.kegg.jp/kegg-bin/download_htext?htext=br08901.keg&format=json&filedir=")
  geneToPath <- keggLink("pathway", "ko")
  koDetailList <- list()
  save.image(file="cache.RData")
}else{
  load("cache.RData")
}
pathCategories <- data.frame(path=character(),PathDesc=character(), Category1=character(), Category=character())
for (node1 in data1$children){
  for (node2 in node1$children){
    for (node3 in node2$children){
      pathName <- paste0("ko", node3$name)
      IdName <- str_split(pathName, "\\s\\s")
      newRow <- data.frame(path=unlist(IdName)[1],
        PathDesc=unlist(IdName)[2],Category1=node1$name, Category=node2$name)
      pathCategories <- rbind(pathCategories,newRow)
    }
  }
}
formatEntry <- function(i){
  entry <- i$ENTRY
  genes <-  i$GENES
  dfEntry <- data.frame("ko"=NA, "kegg"=genes, "name"=NA, "desc"=NA)
  dfEntry$ko <- entry
  dfEntry$name <- i$NAME
  dfEntry$desc <- i$DEFINITION
  dfEntry <- separate(dfEntry, kegg, c("species", "kegg"), sep=":")
  dfEntry <- cSplit(dfEntry, "kegg", sep=" ", direction="long", type.convert=F)
  dfEntry <- cSplit(dfEntry, "kegg", sep="(",fixed=T ,direction="wide",
    type.convert=F)
  dfEntry$kegg <- paste(tolower(dfEntry$species), dfEntry$kegg_1, sep=":")
  dfEntry <- subset(dfEntry, select=c("kegg", "ko", "name", "desc", "species"))
  return(dfEntry)
}
selCat <- subset(pathCategories, PathDesc=="Arachidonic acid metabolism")
dfGeneToPath <- data.frame(path=gsub("path:","",geneToPath),
  gene=gsub("ko:","",names(geneToPath)))
selCatGenes <- merge (selCat, dfGeneToPath, by.x="path", by.y="path")
selGenes <- unique(selCatGenes$gene)
for (ko in selGenes){
  if(!(ko %in% dfKeggToKO$ko)){
    print(ko)
    koDetail <- tryCatch(keggGet(ko), error=function(e) print(ko))
    koDetailList <- c(koDetailList, koDetail)
  }
}
keggToKoList <- lapply(koDetailList, formatEntry)
dfKeggToKO <- do.call("rbind", keggToKoList)
dfKeggToKO <- separate(dfKeggToKO, "kegg", c("org", "gene"),":", remove = F)
save.image(file="cache.RData")
write.csv(dfKeggToKO, file="dfKeggToKO.csv")
unique(dfKeggToKO$kegg)

dfKeggUni<-data.frame()
names(dfKeggUni)<-c("kegg","uniprot")
uniKegg <- keggConv("uniprot", kegg)
dfKeggUni <- rbind(dfKeggUni, uniKegg)


for (kegg in unique(dfKeggToKO$kegg)){
  print(kegg)
  uniKegg <- keggConv("uniprot", kegg)
  dfKeggUni <- rbind(dfKeggUni, uniKegg)
}

unique(dfKeggToKO$kegg)
library(annotate)
blastResult <- blastSequences("agataaaataccggcggccatgttggaggaatgaacaatagatactaatgaaaaatcttttctagaagtt
cctccaagatggccgctgtgacgtaacttgaaaacgaagaattgatcgtttttttcaaccattcccgtga
  gattcaaaagataaaataccggcggccatgttggaggaatgaataATAGTTACTATCttgaaatcctttg
  ttgaagttcctccaaagatggccgctgtgaccTATTaggaaaacgaagaatatgaacaatgaaaaattat
  aaagaattcATTTGCATGATTTTCTTTAGAGAAGGAACCTTCGCATCACGGAAGTACTCTCTAAGCGCCA
  AGACAGGAACACTGCTACCCGACTTTCCCAATGCGTTATTACTGGAATTACCACCAGAAGTATCGCTAGA
  TAAGGTCCATACTTTTATAATCATGTACAAGACGCACTGTCAATGTATATTAGATACTGCCATCAATGCT
  AACTTTGATGAAATCAAaaattttcttcttcatttctGGCAAGGAATGCCAGACCATCTCCTTCCGCTCA
  TGGACACTGAGGTTGTATCTGACGCCATCTGCGTGTGTGATTCTATACTTTATAAGGTATTGATGGATGT
  CCTGCTACCATCATCAATGCAAGACGTCCCGGAGGGtgttggTGACGGCACATGGAATCCAAGCTAATTA
  TTCACAGAAATCGATCGAGATTTTTTATACACTTCTTTTTCGCATTTGTCAAATGGCTTGATCATAACCA
  AATGTCTTGGTAGATAATTCAACTGGATTTGTTTGAATTTAGGTTGCTTGCGGATATAAGAACTTTTGCG
  CGTCACCTAGAAACCTGGGTAACCACGGCAACAGAGCAACTACCGGAGTACATAAAAGTCGGCAAAGTGC
  AAAGTAAGAGGGCCTTTTGTTAATTAATAGACCAATYCCgcattcgcttggatgaacacgaaacaatcac
  ttcatttcgaggctggagttttcgacaaaagataaaactccagcctcgaaatgaagccattgttccATGT
  TAATCCAAGCGAAtgcgggattggtctatttagacatcaattaaaaatatgcaaagacgctattttgtaa
  aacaaaacaacatggGATTTGTTTTCCTCGCACCTTAAAAAATCTSCCTAAAATCTTTTTTGCGGTCGAC
  GCATCATTTTTAACCACTCACTAATCGAGAAATAACGCATTCCCCATYGTCTCTTGAATTtcgctggatg
  gcaaaagtgaGAACTCTGGGTATGAGATTGGGACAGGGCCTTACACTTCTTTTCGTTATTCAGTTGTCAG
  AAGATTCGCTCAATCCATCAGAAGACAGACATCGTTCCTTCATCTTTCGCAGGTATGAAATGCAACTGCT
  ATCATCGAGATCCgggttgttgttttatttgttgatgCTGCTGCTGCagcaataatccagtttcgaattc
  tcccatgctctgtgttagcctcatttctGCCCAAAATATTcgtaattgtgagaataaaagtgtgagattt
  caaatcaatttaattgtttttcttcgcaaattccatgcgttttcttgggctatcaataaaacaaaaaagt
  attttacttggtttttgaggctaacacagagcaaatcagaactcgaaattggagtattattattgttgct
  gctgttgttgttgttgatgctgctgctgtttttgtggttgttgttgttgatgttgttgttgctgctgctg
  ttgttgttgttgttgttgttgattttgctttattgttgttgttgctgctgctactgttgttgttgttgtt
  gttaattttgctatattgttgttgttgttgtttttgctatattgttgctgttgctgttgttgtttttgct
  attgtcgttgttgttgttgttgttgttgttgttgtcgttgttgttgttgctgttaatTTTGctatattgt
  tgttgttgttgtttttgctatattgttgctgttgctgttgttgtttttgctattgtcgttgttgttgttg
  ttgttgttgttgttgtcgttgttgttgttgctgttaatTTTGctatattgttgttgttgtNAATGTTTAG
  TGTTGAGACTTGAAACTATGTcttcttacttattttgtTCCCCTCTTTCCATTTTCRTTAAAATTTATAT
  CcataattcaagaaaaaggACAATTTTGTACGGAAAATCTATGGTTATTCACATTTGAATACGTTTTCAT
  CGGCGAGAGATAGATTTCCtgctttttttgctttttttgttgttttgttttggtttcctgTTCAATTACA
  TTCGTCCTTCATTCATTCCCTGATACAATGCGCAAAAGCAAACGAATGCATCGGTTCGGATATACTAGGG
  AAATTATTTCACCTGATTTAGTAGTTGATCGAGTGAAGCATTTAGGACAGTTTCAAACTCctaaagaata
  ataaataatgaaaaaaatcaaaaatatcataaagaataataaataccAGGACAATATATacagaataata
  aataatagcaaaaatatataatacaatgaATAATGGGGGAAAAATAAGACACGAATAAAGWATAATTGGC",
  database="swissprot",
  program="blastx",as="data.frame", expect=1e-5)
