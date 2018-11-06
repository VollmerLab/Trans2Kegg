#' Merge pathway information and summarize DE counts
#' @name mergePaths
#' @import tidyr
#' @importFrom utils read.csv write.csv
#' @importFrom stats aggregate
#' @return Annotation results are written to deWithPathsAndDetails.csv,
#' countByClass.csv, countByPath.csv
#' @examples
#' mergePaths()
#' @export

mergePaths <- function(){
    deCovAndCountDesc <- read.csv("deCovAndCountDesc.csv")
    deCovAndCountDesc <- subset(deCovAndCountDesc, select=-c(X))
    dfPathsKos <- read.csv("dfPathsKos.csv")
    dfPathsKos <- subset(dfPathsKos, select=-c(X))
    deWithPaths <- merge(deCovAndCountDesc, dfPathsKos, by="ko")
    deWithPaths <- unique(deWithPaths)
    pathDetails <- read.csv("dfPaths.csv")
    pathDetails <- subset(pathDetails, select=-c(X))
    deWithPathsAndDetails <- merge(deWithPaths, pathDetails, by="id")
    deWithPathsAndDetails$direction[deWithPathsAndDetails$log2FoldChange > 0] = 'Up'
    deWithPathsAndDetails$direction[deWithPathsAndDetails$log2FoldChange < 0] = 'Down'
    deWithPathsAndDetails <- separate(data = deWithPathsAndDetails, col = class, into = c("category", "class"), sep = ";")
    write.csv(deWithPathsAndDetails, file="deWithPathsAndDetails.csv", row.names=FALSE)
    countByClass <- aggregate(ko ~ Factor + category + class + direction, data=deWithPathsAndDetails, NROW)
    write.csv(countByClass, file="countByClass.csv", row.names=FALSE)
    countByPath <- aggregate(ko ~ Factor + category + class + path + direction, data=deWithPathsAndDetails, NROW)
    write.csv(countByPath, file="countByPath.csv", row.names=FALSE)
}
