#' Get KEGG pathway information for annotated transcripts
#' @name getPathways
#' @import KEGGREST
#' @importFrom utils read.csv write.csv
#' @param annotFile csv output file from annotateTranscripts
#' @examples
#' filepath <- system.file("extdata", "annot.csv", package="Trans2Kegg")
#' getPathways(filepath)
#' @export
getPathways <- function(annotFile){
    annot <- read.csv(annotFile, stringsAsFactors = FALSE)
    annotNoNA <- na.omit(annot)
    koList <- unique(annotNoNA$ko)
    pathways <- keggLink("pathway", koList)
    uniquePathways <- unique(pathways)
    pathwayDetails <- keggGet(uniquePathways)
    dfPaths <- data.frame( "id" = c(), "class" = c(), "path" = c())
    dfPathsKos <- data.frame("id" = c(), "ko" = c())
    dfRefs <- data.frame("id" = c(), "title" = c(), "authors" = c())
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
        dfPath <- data.frame("id" = id, "class" = class,"path"=name)
        dfPaths <- rbind(dfPath, dfPaths)
    }
    write.csv(dfPaths, file="dfPaths.csv")
    write.csv(dfPathsKos, file="dfPathsKos.csv")
    write.csv(dfRefs, file="dfRefs.csv")
}
