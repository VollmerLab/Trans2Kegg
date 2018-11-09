#' Get KEGG pathway information for annotated transcripts
#' @name getPathways
#' @import KEGGREST
#' @importFrom utils read.csv write.csv
#' @importFrom stats na.omit
#' @param annotFile csv output file from annotateTranscripts
#' @return Annotation results are written to dfPaths.csv, dfPathsKos.csv
#' @examples
#' filepath <- system.file("extdata", "deCovAndCountDesc.csv", 
#' package="Trans2Kegg")
#' getPathways(filepath)
#' @export
getPathways <- function(annotFile){
    annot <- read.csv(annotFile, stringsAsFactors = FALSE)
    annotNoNA <- na.omit(annot)
    koList <- unique(annotNoNA$ko)
    dfPaths <- data.frame( "id" = c(), "class" = c(), "path" = c())
    dfDone <- data.frame("id" =c(), "ko" = c())
    dfPathsKos <- data.frame("id" = c(), "ko" = c())
    if (file.exists("dfPaths.csv")) {
        dfPaths <- read.csv("dfPaths.csv", stringsAsFactors = FALSE)
    }else{
        write.csv(dfPaths, file="dfPaths.csv", row.names=TRUE)
    }
    if (file.exists("dfDone.csv")) {
        dfDone <- read.csv("dfDone.csv", stringsAsFactors = FALSE)
    }else{
        write.csv(dfDone, file="dfDone.csv", row.names=TRUE)
    } 
    if (file.exists("dfPathsKos.csv")){
        dfPathsKos <- read.csv("dfPathsKos.csv", stringsAsFactors = FALSE)
    }else{
        write.csv(dfPathsKos, file="dfPathsKos.csv", row.names=TRUE)
    }
    kosLeft <- setdiff(koList, unique(dfDone$ko))
    for (deKo in kosLeft){
        pathways <- keggLink("pathway", c(deKo))
        uniquePathways <- unique(pathways)
        if(length(uniquePathways) > 0){
            pathsLeft <- setdiff(uniquePathways, dfPaths$id)
            for(uniquePathway in pathsLeft){
                pathwayDetails <- keggGet(uniquePathway)
                for (pathwayDetail in pathwayDetails){
                    #print(pathwayDetail)
                    class <- pathwayDetail$CLASS
                    name <- pathwayDetail$NAME
                    id <- names(pathwayDetail$PATHWAY_MAP)
                    kos <- pathwayDetail$ORTHOLOGY
                    koDescriptions <- unname(kos)
                    koIds <- names(kos)
                    for (ko in koIds){
                        dfPathKo <- data.frame("id"=id, "ko"=koIds)
                        write.table(dfPathKo, file="dfPathsKos.csv", 
                                    col.names=FALSE, append=TRUE, 
                                    sep=',', row.names=FALSE)
                        #dfPathsKos <- rbind(dfPathKos, dfPathsKo)
                    }
                }
                if(length(class) > 0){  
                    dfPath <- data.frame("id" = id, "class" = class,
                        "path"=name)
                    write.table(dfPath, file="dfPaths.csv", 
                                col.names=FALSE, append=TRUE, 
                                sep=',', row.names=FALSE)
                    #dfPaths <- rbind(dfPath, dfPaths)
                }
            }
        }
        dfKoDone<- data.frame("id" =deKo, "ko"=deKo)
        write.table(dfKoDone, file="dfDone.csv",
            col.names=FALSE, append=TRUE,
            sep=',', row.names=FALSE)
    }
}
#library(KEGGREST)
#getPathways("deCovAndCountDesc.csv")
