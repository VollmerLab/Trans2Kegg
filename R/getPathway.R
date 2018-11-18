#' Get KEGG pathway information for annotated transcripts
#' @name getPathways
#' @import KEGGREST
#' @importFrom utils read.csv write.csv
#' @importFrom stats na.omit
#' @param annotFile csv output file from annotateTranscripts
#' @return Annotation results are written to path.csv, pthKo.csv
#' @examples
#' filepath <- system.file("extdata", "cvCnt.csv", 
#' package="Trans2Kegg")
#' getPathways(filepath)
#' @export
getPathways <- function(annotFile = "cvCnt.csv"){
    pthKo <- "pthKo.csv"
    path <- "path.csv"
    done <- "dfDone.csv"
    annot <- read.csv(annotFile, stringsAsFactors = FALSE)
    annotNoNA <- na.omit(annot)
    koList <- unique(annotNoNA$ko)
    dfPaths <- data.frame( "id" = c(), "class" = c(), "path" = c())
    dfDone <- data.frame("id" =c(), "ko" = c())
    dfPathsKos <- data.frame("id" = c(), "ko" = c())
    if (file.exists(path)) {
        dfPaths <- read.csv(path, stringsAsFactors = FALSE)
    }    
    if (file.exists(done)) {
        dfDone <- read.csv(done, stringsAsFactors = FALSE)
    }
    if (file.exists(pthKo)){
        dfPathsKos <- read.csv(pthKo, stringsAsFactors = FALSE)
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
                        write.table(dfPathKo, file=pthKo, 
                                    col.names=!file.exists(pthKo), 
                                    append=file.exists(pthKo), 
                                    sep=',', row.names=FALSE)
                    }
                }
                if(length(class) > 0){  
                    dfPath <- data.frame("id" = id, "class" = class,
                        "path"=name)
                    write.table(dfPath, file=path, 
                                col.names=!file.exists(path), 
                                append=file.exists(path),
                                sep=',', row.names=FALSE)
                }
            }
        }
        dfKoDone<- data.frame("id" =deKo, "ko"=deKo)
        write.table(dfKoDone, file=done,
            col.names=!file.exists(done), 
            append=file.exists(done),
            sep=',', row.names=FALSE)
    }
}
