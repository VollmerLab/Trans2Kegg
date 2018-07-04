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
  save.image("dfUniKegg.RData")
}

getKO <- function(){
  sink(file="koDetailList.log")
  koList <- keggList("ko")
  koDetailList <- list()
  for (ko in names(koList)){
    print(ko)
    koDetail <- keggGet(ko)
    koDetailList <- c(koDetailList, koDetail)
  }
  save.image("koDetailList.RData")
}

init_data1 <- function(){
  plan(multiprocess)
  future(getKO())
}
init_data2 <- function(){
  plan(multiprocess)
  future(getOrg)
}

