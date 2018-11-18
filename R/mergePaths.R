#' Merge pathway information and summarize DE counts
#' @name mergePaths
#' @import tidyr
#' @importFrom utils read.csv write.csv
#' @importFrom stats aggregate
#' @param prefix File prefix for input files.
#' @return Annotation results are written to dePathsDetails.csv,
#' countByClass.csv, countByPath.csv
#' @examples
#' covCount <- system.file("extdata", "cvCnt.csv", 
#'     package="Trans2Kegg")
#' prefix <- system.file("extdata/", package="Trans2Kegg")
#' mergePaths(prefix=prefix)
#' @export

mergePaths <- function(prefix = "") {
    cvCnt <- paste0(prefix, "cvCnt.csv")
    pthKo <- paste0(prefix, "pthKo.csv")
    pth <- paste0(prefix, "path.csv")
    deCovAndCountDesc <- unique(read.csv(cvCnt))
    dfPathsKos <- unique(read.csv(pthKo))
    deWithPaths <- merge(deCovAndCountDesc, dfPathsKos, by="ko")
    deWithPaths <- unique(deWithPaths)
    pathDetails <- unique(read.csv(pth))
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
    write.csv(dfPathsKos, file=pthKo, row.names=FALSE)
    write.csv(deCovAndCountDesc, file=cvCnt, row.names=FALSE)
    write.csv(pathDetails, file=pth, row.names=FALSE)
}
