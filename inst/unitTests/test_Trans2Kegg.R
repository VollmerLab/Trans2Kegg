test_Trans2Kegg <- function(){
    filepath <- system.file("extdata", "aiptasia.fa", package="Trans2Kegg")
    returnVal <- annotateTranscripts(c("KXJ29317.1"), filepath, "annot.csv")
    checkEquals(returnVal, 0)
}

