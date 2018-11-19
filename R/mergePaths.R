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
    # load DE, ortholog, and pathway info
    cvCnt <- paste0(prefix, "cvCnt.csv")
    pthKo <- paste0(prefix, "pthKo.csv")
    pth <- paste0(prefix, "path.csv")
    # Make sure no duplicates
    deCovAndCountDesc <- unique(read.csv(cvCnt))
    dfPathsKos <- unique(read.csv(pthKo))
    # Merge DE and pathways
    deWithPaths <- merge(deCovAndCountDesc, dfPathsKos, by="ko")
    deWithPaths <- unique(deWithPaths)
    pathDetails <- unique(read.csv(pth))
    # Merge pathway details
    dePathsDetails <- merge(deWithPaths, pathDetails, by="id")
    # Calculate direction column for up-regulated or down-regulated
    dePathsDetails$direction[dePathsDetails$log2FoldChange > 0] = 'Up'
    dePathsDetails$direction[dePathsDetails$log2FoldChange < 0] = 'Down'
    # Separate the KEGG category and class into separate columns
    dePathsDetails <- separate(data = dePathsDetails, 
        col = class, into = c("category",
        "class"), sep = ";")
    # Write detail file
    write.csv(dePathsDetails, file="dePathsDetails.csv", 
        row.names=FALSE)
    # Count DE by category, class, direction
    countByClass <- aggregate(ko ~ Factor + category + class + direction, 
        data=dePathsDetails, NROW)
    write.csv(countByClass, file="countByClass.csv", row.names=FALSE)
    # Count DE by class, pathway, and direction
    countByPath <- aggregate(ko ~ Factor + category + class + path + direction, 
        data=dePathsDetails, NROW)
    # Write summary files
    write.csv(countByPath, file="countByPath.csv", row.names=FALSE)
    write.csv(dfPathsKos, file=pthKo, row.names=FALSE)
    write.csv(deCovAndCountDesc, file=cvCnt, row.names=FALSE)
    write.csv(pathDetails, file=pth, row.names=FALSE)
}
