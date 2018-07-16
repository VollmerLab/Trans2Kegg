load("dfUniKegg.RData")
load("koDetailList.RData")
library(tidyr)
library(dplyr)
library(splitstackshape)
formatEntry <- function(i){
  entry <- i$ENTRY
  genes <-  i$GENES
  dfEntry <- data.frame("ko"=NA, "kegg"=genes, "name"=NA, "desc"=NA)
  dfEntry$ko <- entry
  dfEntry$name <- i$NAME
  dfEntry$desc <- i$DEFINITION
  dfEntry <- separate(dfEntry, kegg, c("species", "kegg"), sep=":")
  dfEntry <- cSplit(dfEntry, "kegg", sep=" ", direction="long", type.convert=F)
  dfEntry <- cSplit(dfEntry, "kegg", sep="(",fixed=T ,direction="wide", type.convert=F)
  dfEntry$kegg <- paste(tolower(dfEntry$species), dfEntry$kegg_1, sep=":")
  dfEntry <- subset(dfEntry, select=c("kegg", "ko", "name", "desc"))
  return(dfEntry)
}
keggToKoList <- lapply(koDetailList, formatEntry)
dfKeggToKO <- do.call("rbind", keggToKoList)
dfKeggToKO <- separate(dfKeggToKO, "kegg", c("org", "gene"),":", remove = F)
head(dfKeggToKO)
org1 <- unique(dfKeggToKO$org)
library(KEGGREST)
org <- keggList("organism")
org2 <- unique(org[,2])
length(org1)
length(org2)
readBLAST <- function(blastFile){
  load("dfUniKegg.RData")
  load("koDetailList.RData")
  #test <- mclapply(koDetailList, formatEntry, mc.cores=4L)
  dfUniKegg <- separate(dfUniKegg, "uni", c("prefix", "uni"), sep=":")
  getPathways <- function(i){
    entry <- i$ENTRY
    pathway <-  i$PATHWAY
    if(!is.null(pathway)){
      #print(pathway)
      dfEntry <- data.frame("ko" = NA,"pathdesc"=pathway)
      dfEntry$ko <- entry
      dfEntry$pathId <- rownames(dfEntry)
      return(dfEntry)
    }
  }
  pathList <- lapply(koDetailList, getPathways)
  dfKoPathway <- do.call("rbind", pathList)

  blastHits <- read.table(blastFile, sep="\t", header=F, quote="")
  colnames(blastHits) <- c("gene", "sp", "ident", "evalue", "cov")
  spBlast <- separate(blastHits, sp, c("spPrefix", "sp", "symbol"), sep = "\\|")
  spBlast <- separate(spBlast, gene, c("gene", "Tran"), sep = "-", remove=F)
  spBlast2 <- subset(spBlast, evalue < 1e-10 & cov > 50 & ident > 50,
    select=-c(spPrefix, Tran))
  spBlast3 <- separate(spBlast2, symbol, c("symbol", "species"), sep = "_")

  head(dfUniKegg)
  head(spBlast3)
  blastToKegg <- merge(spBlast3, dfUniKegg, by.x="sp", by.y="uni")
  head(blastToKegg)
  blastToKO <- merge(blastToKegg, dfKeggToKO, all.x=T)
  write.csv(blastToKO, file="blastToKO.csv")
  koSummCount <- aggregate(sp~gene+ko, blastToKO,FUN=NROW)
  koSummIdent <- aggregate(ident~gene+ko, blastToKO,FUN=mean)
  koSummCov <- aggregate(cov~gene+ko, blastToKO,FUN=mean)
  koSumm <- merge(koSummCount, koSummIdent)
  koSumm <- merge(koSumm, koSummCov)
  head(koSumm)
  write.csv(koSumm, file="koSumm.csv")
  koToDesc <- unique(subset(blastToKO, select=c("ko", "name", "desc")))
  koSummDesc<- merge(koSumm, koToDesc)
  koSummDesc<- merge(koSummDesc, dfKoPathway)
  write.csv(koSummDesc, file="koSummDesc.csv")
}
