#' Merge pathway information and summarize DE counts
#' @name mergePaths
#' @import tidyr
#' @importFrom utils read.csv write.csv
#' @importFrom stats aggregate
#' @param covCount csv output file from
#' @param pathKo csv output file from
#' @param paths csv output file from
#' @return Annotation results are written to dePathsDetails.csv,
#' countByClass.csv, countByPath.csv
#' @examples
#' covCount <- system.file("extdata", "deCovAndCountDesc.csv", 
#'     package="Trans2Kegg")
#' pathKo <- system.file("extdata", "dfPathsKos.csv", package="Trans2Kegg")
#' paths <- system.file("extdata", "dfPaths.csv", package="Trans2Kegg")
#' mergePaths(covCount, pathKo, paths)
#' @export

mergePaths <- function(covCount, pathKo, paths){
    deCovAndCountDesc <- read.csv(covCount)
    deCovAndCountDesc <- subset(deCovAndCountDesc, select=-c(X))
    dfPathsKos <- read.csv(pathKo)
    dfPathsKos <- subset(dfPathsKos, select=-c(X))
    deWithPaths <- merge(deCovAndCountDesc, dfPathsKos, by="ko")
    deWithPaths <- unique(deWithPaths)
    pathDetails <- read.csv(paths)
    pathDetails <- subset(pathDetails, select=-c(X))
    dePathsDetails <- merge(deWithPaths, pathDetails, by="id")
    dePathsDetails$direction[dePathsDetails$log2FoldChange > 0] = 'Up'
    dePathsDetails$direction[dePathsDetails$log2FoldChange < 0] = 'Down'
    dePathsDetails <- separate(data = dePathsDetails, 
        col = class, into = c("category",
        "class"), sep = ";")
    write.csv(dePathsDetails, file="dePathsDetails.csv", 
        row.names=FALSE)
    countByClass <- aggregate(ko ~ Factor + category + class + direction, 
        data=dePathsDetails, NROW)
    write.csv(countByClass, file="countByClass.csv", row.names=FALSE)
    countByPath <- aggregate(ko ~ Factor + category + class + path + direction, 
        data=dePathsDetails, NROW)
    write.csv(countByPath, file="countByPath.csv", row.names=FALSE)
}
