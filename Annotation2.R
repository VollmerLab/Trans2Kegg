#library(dplyr)
#library(tidyr)
#library(knitr)
#library(kableExtra)
#library(pander)
#library(UpSetR)
#library(pathview)
#library(igraph)
#library(data.table)

blastFile <- "blastSP.tsv"

blastHits <- read.table(blastFile, sep="\t", header=F, quote="")
head(blastHits)


blastHits$cov <- blastHits$V4/blastHits$V3
blastHits <- subset(blastHits, select=c("V1", "V2", "V6", "cov"))
head(blastHits)
colnames(blastHits) <- c("Gene", "SP", "e", "cov")
spBlast <- separate(blastHits, SP, c("spPrefix", "sp", "symbol"), sep = "\\|")
spBlast <- separate(spBlast, Gene, c("Gene", "Tran"), sep = "-")
head(spBlast)
spBlast <- separate(spBlast, symbol, c("symbol", "species"), sep = "_")
spBlast <- subset(spBlast, select=c("Gene", "sp", "symbol", "e", "cov"))
head(spBlast)
uniToKegg <- read.table("uniprot.tsv", header=F, sep="\t", stringsAsFactors = F, quote="")
head(uniToKegg)
uniToKegg <- separate(uniToKegg, V1, c("source", "sp"), sep = ":")
head(uniToKegg)
uniToKegg <- subset(uniToKegg, select=c("sp", "V2"))
colnames(uniToKegg) <- c("sp", "kegg")
head(uniToKegg)
keggToKo <- read.table("keggToKO.tsv", sep="\t", header=F, quote="", stringsAsFactors = F)
head(keggToKo)
keggToKo <- separate(keggToKo, V2, c("prefix", "KO"), sep = ":")
head(keggToKo)
colnames(keggToKo) <- c("kegg", "prefix", "KO", "Desc")
keggToKo <- subset(keggToKo, select=c("kegg", "KO", "Desc"))
head(keggToKo)
spToKo <- merge(uniToKegg, keggToKo)
head(spToKo)
write.csv(spToKo, file="spToKo.csv")
kxToKo <- merge(spBlast, spToKo, by="sp")
head(kxToKo)
write.csv(kxToKo, file="kxToKo.csv")
keggPathways <- read.table("koToPathway.tsv", sep='\t', header=TRUE, quote='')
head(keggPathways)
#every pathway and catagory in all of KEGG
pathCategories <- data.frame(PathDesc=character(), Category1=character(), Category=character())
data1 <- fromJSON(file="http://www.kegg.jp/kegg-bin/download_htext?htext=br08901.keg&format=json&filedir=")
for (node1 in data1$children){
  for (node2 in node1$children){
    for (node3 in node2$children){
      pathName <- paste0("ko", node3$name)
      pathName <- str_replace(pathName, "  ", " ")
      newRow <- data.frame(PathDesc=pathName,Category1=node1$name, Category=node2$name)
      pathCategories <- rbind(pathCategories,newRow)
    }
  }
}
pathCategories2 <- separate(pathCategories, PathDesc, c("Path", "PathDesc"), sep=8 )
pathCategories2$Path <- gsub(' ', '',pathCategories2$Path)
head(pathCategories2)
write.csv(pathCategories2, file="pathCategories2.csv")
pathAndCategory <- merge(keggPathways, pathCategories2)
colnames(pathAndCategory) <- c("Path", "PathDesc", "KO", "Category1", "Category")
head(pathAndCategory)
head(kxToKo)
kxToKo$lenDiff <- abs(kxToKo$cov - 1)
head(kxToKo)
kxToKo <- kxToKo[order(kxToKo$lenDiff, kxToKo$e),]
head(kxToKo)
DT <- as.data.table(kxToKo)
testKO <- DT[order(lenDiff, e)]
test2KO <- testKO[, head(.SD,1), by = .(KO)]
test3KO <- test2KO[, head(.SD,1), by = .(Gene)]
test3KO
genePathCat <- merge(test3KO, pathAndCategory)
write.csv(genePathCat, file="msexAnnot.csv")







