#' Merge pathway information and summarize DE counts
#' @name mergePaths
#' @import tidyr
#' @importFrom utils read.csv write.csv
#' @importFrom stats aggregate
#' @return Annotation results are written to dePathsDetails.csv,
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
