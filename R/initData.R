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

getOrg <- function(){
  sink(file="dfUniKegg.log")
  org <- keggList("organism")
  orgs <- unique(org[,2])
  dfUniKegg <- data.frame("uni"=NA, "kegg"=NA)
  for (org in orgs){
    print(org)
    conv <- tryCatch(keggConv(org, "uniprot"), error=function(e) print(org))
    if(length(conv) > 10){
      dfConv <- data.frame("uni"=names(conv), "kegg"=conv)
      dfUniKegg <- rbind(dfUniKegg, dfConv)
    }
  }
  save(dfUniKegg, file="dfUniKegg.RData")
}

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

#' Initialize KEGG data
#'
#' @param x A number
#' @return The sum of \code{x} and \code{y}
#' @examples
#' init_data()

getKO <- function(){
  sink(file="koDetailList.log")
  koList <- keggList("ko")
  koDetailList <- list()
  for (ko in names(koList)){
    print(ko)
    koDetail <- tryCatch(keggGet(ko), error=function(e) print(ko))
    #koDetail <- keggGet(ko)
    koDetailList <- c(koDetailList, koDetail)
  }
  keggToKoList <- lapply(koDetailList, formatEntry)
  save(keggToKoList, file="keggToKoList.RData")
  dfKeggToKO <- do.call("rbind", keggToKoList)
  dfKeggToKO <- separate(dfKeggToKO, "kegg", c("org", "gene"),":", remove = F)
  save(dfKeggToKO, file="dfKeggToKO.RData")
  pathList <- lapply(koDetailList, getPathways)
  dfKoPathway <- do.call("rbind", pathList)
  save(dfKoPathway, file="dfKoPathway.RData")
}

#' Initialize KEGG data
#'
#' @param x A number
#' @return The sum of \code{x} and \code{y}
#' @examples
#' init_data()
init_data <- function(){
  plan(multiprocess)
  #future(blastTrans('OGS2_20140407_transNoSpaceHeader2.fa'))
  future(getOrg())
  future(getKO())
  #getKO()
}

