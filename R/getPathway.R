#' Get KEGG pathway information for annotated transcripts
#' @name getPathways
#' @import KEGGREST
#' @importFrom utils read.csv write.csv
#' @importFrom stats na.omit
#' @param annotFile csv output file from annotateTranscripts
#' @return Annotation results are written to dfPaths.csv, dfPathsKos.csv,
#' dfRefs.csv
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
    dfPathsKos <- data.frame("id" = c(), "ko" = c())
    dfRefs <- data.frame("id" = c(), "title" = c(), "authors" = c())
    for (deKo in koList){
        pathways <- keggLink("pathway", c(deKo))
        uniquePathways <- unique(pathways)
        if(length(uniquePathways) > 0){
            pathsLeft <- setdiff(uniquePathways, dfPaths$id)
            for(uniquePathway in pathsLeft){
                pathwayDetails <- keggGet(uniquePathway)
                for (pathwayDetail in pathwayDetails){
                    #print(pathwayDetail)
                    references <- pathwayDetail$REFERENCE
                    class <- pathwayDetail$CLASS
                    name <- pathwayDetail$NAME
                    id <- names(pathwayDetail$PATHWAY_MAP)
                    kos <- pathwayDetail$ORTHOLOGY
                    koDescriptions <- unname(kos)
                    koIds <- names(kos)
                    for (reference in references){
                        dfRef <- data.frame("id" = id, 
                            "title" = reference$TITLE,
                            "authors" = reference$AUTHORS)
                        dfRefs <- rbind(dfRef, dfRefs)
                    }
                    for (ko in koIds){
                        dfPathKos <- data.frame("id"=id, "ko"=koIds)
                        dfPathsKos <- rbind(dfPathKos, dfPathsKos)
                    }
                    if(length(class) > 0){  
                        dfPath <- data.frame("id" = id, "class" = class,
                                             "path"=name)
                        dfPaths <- rbind(dfPath, dfPaths)
                    }
                }
            }
        }
    }
    write.csv(dfPaths, file="dfPaths.csv")
    write.csv(dfPathsKos, file="dfPathsKos.csv")
    write.csv(dfRefs, file="dfRefs.csv")
}
