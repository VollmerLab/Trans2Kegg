#' Get KEGG pathway information for annotated transcripts
#' @name getPathways
#' @import KEGGREST
#' @importFrom utils read.csv write.csv
#' @importFrom stats na.omit
#' @param annotFile csv output file from annotateTranscripts
#' @return Annotation results are written to path.csv, pthKo.csv
#' @examples
#' \dontrun{
#' filepath <- system.file("extdata", "cvCnt.csv", 
#' package="Trans2Kegg")
#' getPathways(filepath)
#' }
#' @export
getPathways <- function(annotFile = "cvCnt.csv"){
    # These were going to be params, but Roxygen kept
    # putting non-standard indents and causing NOTEs
    pthKo <- "pthKo.csv"
    path <- "path.csv"
    # in BiocCheck because param line longer than 80
    done <- "dfDone.csv"
    # Read file of match counts and mean query coverage for each
    # query/ko pair.
    annot <- read.csv(annotFile, stringsAsFactors = FALSE)
    # Remove the "no BLAST hit" rows.
    annotNoNA <- na.omit(annot)
    # Make sure there are no duplicates
    koList <- unique(annotNoNA$ko)
    # Create empty dataframes to rbind to
    dfPaths <- data.frame( "id" = c(), "class" = c(), "path" = c())
    dfDone <- data.frame("id" =c(), "ko" = c())
    dfPathsKos <- data.frame("id" = c(), "ko" = c())
    # Check and see if any paths already retrieved during
    # previous run
    if (file.exists(path)) {
        dfPaths <- read.csv(path, stringsAsFactors = FALSE)
    }    
    if (file.exists(done)) {
        dfDone <- read.csv(done, stringsAsFactors = FALSE)
    }
    if (file.exists(pthKo)){
        dfPathsKos <- read.csv(pthKo, stringsAsFactors = FALSE)
    }
    # Which paths remain to be queried
    kosLeft <- setdiff(koList, unique(dfDone$ko))
    for (deKo in kosLeft){
        # Get pathways for each ko
        pathways <- keggLink("pathway", c(deKo))
        uniquePathways <- unique(pathways)
        if(length(uniquePathways) > 0){
            # Eliminate paths already retrieved due to association 
            # with other kos
            pathsLeft <- setdiff(uniquePathways, dfPaths$id)
            for(uniquePathway in pathsLeft){
                # Get pathway details
                getPathDetails(uniquePathway, pthKo, path)
            }
        }
        # Write to "done" list in case of restart
        dfKoDone<- data.frame("id" =deKo, "ko"=deKo)
        write.table(dfKoDone, file=done, col.names=!file.exists(done), 
            append=file.exists(done), sep=',', row.names=FALSE)
    }
}                

#' Get KEGG pathway details for one path
#' @name getPathDetails
#' @import KEGGREST
#' @importFrom utils write.csv
#' @importFrom stats na.omit
#' @param uniquePathway KEGG pathway ID
#' @param pthKo Output file for path to KO mapping
#' @param path Output file for path details
#' @return Annotation results are written to path.csv, pthKo.csv
#' @examples
#' getPathDetails("ko04621", "pthKo.csv", "path.csv")
#' @export
getPathDetails <- function(uniquePathway, pthKo, path) {
    pathwayDetails <- keggGet(uniquePathway)
    for (pathwayDetail in pathwayDetails){
        # Get fields we want
        class <- pathwayDetail$CLASS
        name <- pathwayDetail$NAME
        id <- names(pathwayDetail$PATHWAY_MAP)
        kos <- pathwayDetail$ORTHOLOGY
        koDescriptions <- unname(kos)
        koIds <- names(kos)
        for (ko in koIds){
            # Append all ko/pathway associations
            dfPathKo <- data.frame("id"=id, "ko"=koIds)
            write.table(dfPathKo, file=pthKo, 
                col.names=!file.exists(pthKo), append=file.exists(pthKo), 
                sep=',', row.names=FALSE)
        }
    }
    if(length(class) > 0){ 
        # Append pathway details 
        dfPath <- data.frame("id" = id, "class" = class, "path"=name)
        write.table(dfPath, file=path, col.names=!file.exists(path), 
            append=file.exists(path), sep=',', row.names=FALSE)
    }
} 
