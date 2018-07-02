init_data <- function(){
  org <- keggList("organism")
  orgs <- unique(org[,2])
  dfUniKegg <- data.frame("uni"=NA, "kegg"=NA)
  for (org in orgs){
    conv <- tryCatch(keggConv(org, "uniprot"), error=function(e) print(org))
    if(length(conv) > 10){
      dfConv <- data.frame("uni"=names(conv), "kegg"=conv)
      dfUniKegg <- rbind(dfUniKegg, dfConv)
    }
  }
  save.image("dfUniKegg.RData")

  koList <- keggList("ko")
  koDetailList <- list()
  for (ko in names(koList)){
    koDetail <- keggGet(ko)
    koDetailList <- c(koDetailList, koDetail)
  }
  save.image("keggOrthologs.RData")
}
