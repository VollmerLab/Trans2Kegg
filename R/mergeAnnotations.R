#' Merge DE, BLAST, KEGG Info and filter
#' @name mergeAnnotations
#' @import dplyr
#' @importFrom utils read.csv write.csv
#' @param annotFile csv output file from annotateTranscripts
#' @return Annotation results are written to deCovAndCountDesc.csv.
#' @examples
#' annot <- system.file("extdata", "annot.csv", 
#'     package="Trans2Kegg")
#' mergeAnnotations(annot)
#' @export

mergeAnnotations <- function(annotFile){
    annot <- read.csv(annotFile)
    de <- read.csv("dfAll.csv")
    colnames(de) <- c("Iteration_query.def", "log2FoldChange", "padj", "Factor")
    annot$qCov <- (annot$Hsp_align.len - annot$Hsp_gaps)/annot$Hit_len
    annot <- subset(annot, qCov > .5 & Hsp_evalue < 1e-10, 
        select=c("ko","Iteration_query.def", "Hsp_evalue",
        "qCov", "desc"))  
    annotAggregated <- count(annot, ko, Iteration_query.def)
    annotCov <- subset(annot, select=c("ko", "Iteration_query.def", "qCov"))
    meanCov <- aggregate(qCov ~ ko + Iteration_query.def,data=annotCov, mean)
    koDesc <- unique(subset(annot, select=c("ko", "desc")))
    covAndCount1 <- merge(annotAggregated, meanCov)
    covAndCount2 <- aggregate(n ~ Iteration_query.def,data=covAndCount1, max)
    covAndCount3 <- merge(covAndCount1, covAndCount2)
    covAndCount4 <- aggregate(n ~ ko,data=covAndCount3, max)
    covAndCount5 <- merge(covAndCount1, covAndCount4)
    covAndCount6 <- aggregate(qCov ~ ko,data=covAndCount5, max)
    covAndCount7 <- merge(covAndCount1, covAndCount6)
    covAndCount8 <- aggregate(qCov ~ Iteration_query.def,data=covAndCount7, max)
    covAndCount9 <- merge(covAndCount1, covAndCount8)
    covAndCountDesc <- merge(covAndCount9, koDesc)
    deCovAndCountDesc <- merge(de, covAndCountDesc)
    write.csv(deCovAndCountDesc, file="deCovAndCountDesc.csv", row.names=FALSE)
}
