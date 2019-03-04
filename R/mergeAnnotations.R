#' Merge DE, BLAST, KEGG Info and filter
#' @name mergeAnnotations
#' @importFrom utils read.csv write.csv
#' @importFrom dplyr count
#' @param annotFile csv output file from annotateTranscripts
#' @param cov Query coverage cutoff
#' @param e_val e-value cutoff
#' @return Annotation results are written to cvCnt.csv.
#' @examples
#' annot <- system.file("extdata", "annot.csv", 
#'     package="Trans2Kegg")
#' mergeAnnotations(annot)
#' @export

mergeAnnotations <- function(annotFile="annot.csv", cov=.5, e_val=1e-10){
    #These variables are defined in the header of the csv file read below.
    #Defining here to avoid note.
    ko <- ''
    Iteration_query.def <- ''
    Hsp_evalue <- 1
    qCov <- 0
    desc <- ''
    annot <- read.csv(annotFile)
    # Load DE information
    de <- read.csv("dfAll.csv")
    colnames(de) <- c("Iteration_query.def", "log2FoldChange", "padj", "Factor")
    # Calculate query coverage
    annot$qCov <- (annot$Hsp_align.len - annot$Hsp_gaps)/annot$Hit_len
    # Filter on coverage and e-value
    annot <- subset(annot, qCov > cov & Hsp_evalue < e_val, 
        select=c("ko","Iteration_query.def", "Hsp_evalue",
        "qCov", "desc"))
    # Count number of matches for each query/ortholog pair  
    annotAggregated <- count(annot, ko, Iteration_query.def)
    annotCov <- subset(annot, select=c("ko", "Iteration_query.def", "qCov"))
    # Calculate mean coverage for each query/ortholog pair
    meanCov <- aggregate(qCov ~ ko + Iteration_query.def,data=annotCov, mean)
    koDesc <- unique(subset(annot, select=c("ko", "desc")))
    covAndCount1 <- merge(annotAggregated, meanCov)
    # Get the ko with most matches for each query
    covAndCount2 <- aggregate(n ~ Iteration_query.def,data=covAndCount1, max)
    covAndCount3 <- merge(covAndCount1, covAndCount2)
    # Get the query with most matches for each ko
    covAndCount4 <- aggregate(n ~ ko,data=covAndCount3, max)
    covAndCount5 <- merge(covAndCount1, covAndCount4)
    # Get the max coverage for each ko
    covAndCount6 <- aggregate(qCov ~ ko,data=covAndCount5, max)
    covAndCount7 <- merge(covAndCount1, covAndCount6)
    # Get the max coverage for each query
    covAndCount8 <- aggregate(qCov ~ Iteration_query.def,data=covAndCount7, max)
    covAndCount9 <- merge(covAndCount1, covAndCount8)
    covAndCountDesc <- merge(covAndCount9, koDesc)
    deCovAndCountDesc <- merge(de, covAndCountDesc)
    # Write to output file
    write.csv(deCovAndCountDesc, file="cvCnt.csv", row.names=FALSE)
}
