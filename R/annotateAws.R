#' Annotate DE transcripts using NCBI BLAST AMI
#' @name annotateAws
#' @import tidyr
#' @import KEGGREST
#' @import Rcurl
#' @import XML
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom utils read.csv write.csv write.table
#' @param accessions A character vector of accessions from FASTA file to
#' be annotated
#' @param refTransFile FASTA file of transcripts to be annotated
#' @param outFile csv output file for annotation results
#' @param instance Instance from Amazon EC2 running NCBI BLAST AMI
#' @param dns DNS from Amazon EC2 running NCBI BLAST AMI
#' @param threads Number of threads to use for BLAST
#' @return Annotation results are written to the csv file specified by outFile
#' @examples
#' \dontrun{
#' instance <- 'i-007ad8d6b0e9981ec'
#' dns <- 'ec2-52-67-171-31.sa-east-1.compute.amazonaws.com'
#' filepath <- system.file("extdata", "aiptasia.fa", package="Trans2Kegg")
#' annotateAws(c("KXJ29317.1"), filePath,"annot.csv", instance=instance, dns=dns, threads=4)
#' annotateTranscripts(c("KXJ29317.1"), filepath,
#'     "annot.csv")
#' }
#' @export
annotateAws <- function(accessions, refTransFile,
                     outFile="annot.csv", instance, dns, threads){
    transDone <- data.frame() 
    if(file.exists(outFile)){
        transDone <- read.csv(outFile, stringsAsFactors=FALSE)$Iteration_query.def
    }
    transLeft <- setdiff(accessions, transDone)
    refTrans <- readDNAStringSet(refTransFile)
    for(accession in transLeft){
        print(accession)
        transSeqs <- refTrans[accession]
        writeXStringSet(transSeqs, "deTrans.fa", append=FALSE, compress=FALSE, 
                   compression_level=NA, format="fasta")
        query <- readChar(fastaOut, file.info(fastaOut)$size)
        blastResult <- blastSequencesAws(query,
        database="swissprot",program="blastx", 
        as="data.frame", expect=1e-5,instance=instance,
        dns=dns, threads=threads)
        if(is.data.frame(blastResult) & length(blastResult) > 0){
            blastResult <- subset(blastResult, select=c("Iteration_query-def",
                "Iteration_query-len", "Hit_accession", "Hit_len", 
                "Hsp_evalue", "Hsp_identity", "Hsp_gaps", "Hsp_align-len"))
            write.csv(blastResult, file="blastResult.csv")
            annotateTranscriptsAws(blastResult, outFile)
        }
        emptyRow <- data.frame(uniprot=NA, kegg = NA, ko=NA, desc=NA,
            Iteration_query.def=accession, Iteration_query.len=NA, Hit_len=NA,
            Hsp_evalue=NA, Hsp_identity=NA, Hsp_gaps = NA, Hsp_align.len=NA)
        write.table(emptyRow, file=outFile,
            col.names=!file.exists(outFile),
            append=file.exists(outFile), sep=',', 
            row.names=FALSE)
    }       
}

annotateTranscriptsAws <- function(blastResult,
    outFile="annot.csv"){
    rowNum <- 0
    accDone <- c()
    dfUniKegg <- data.frame()
    #blastResult <- read.csv("blastResult.csv")
    uniprots <- unique(blastResult$Hit_accession)
   for (uniprot in uniprots){
        uniKegg <- keggConv("genes", paste0("up:", uniprot))
        if(length(uniKegg) > 0){
            newRow <- data.frame("uniprot"=uniprot, "kegg"=uniKegg)
            dfUniKegg <- rbind(dfUniKegg, newRow)
            for(kegg in newRow$kegg){
                koDetail <- tryCatch(keggGet(kegg),
                    error=function(e) print(kegg))
                ortho <- koDetail[[1]]$ORTHOLOGY
                if(length(names(ortho)) > 0){
                    if(!is.na(kegg)){
                        rowNum <- rowNum + 1
                        newOrtho <- tryCatch(data.frame(row.names=kegg,
                            "kegg"=kegg, "ko"=names(ortho),
                            "desc"=ortho[names(ortho)]), 
                            error=function(e) print(ortho))
                        uniToKO <- merge(newRow, newOrtho)
                        write.csv(uniToKO, file="uniToKO.csv")
                        blastToKO <- merge(uniToKO,blastResult, 
                            by.x="uniprot", by.y="Hit_accession")
                        write.table(blastToKO, file=outFile,
                            col.names=!file.exists(outFile), 
                            append=file.exists(outFile), 
                            sep=',', row.names=FALSE)
                    }else{
                        print(ortho)
                        return(-1)
                    }
                }
            }
        }
    }
}

# Note: All functions below this line are based on source code
# from annotate.
# This subroutine was taken un-changed from annotate source.
.blastSequencesToDataFrame <- function(xml) {
    if (xpathSApply(xml, "count(//Hit)") == 0L) {
        message("'blastSequences' returned 0 matches")
        return(data.frame())
    }

    iter <- xml["//Iteration"]
    iterlen <- sapply(iter, xpathSApply, "count(.//Hsp)")
    iterdf <- xmlToDataFrame(iter, stringsAsFactors=FALSE)

    hit <- xml["//Hit"]
    hitlen <- sapply(hit, xpathSApply, "count(.//Hsp)")
    hitdf <- xmlToDataFrame(hit, stringsAsFactors=FALSE)
    hitdf <- hitdf[, names(hitdf) != "Hit_hsps", drop=FALSE]

    hsp <- xmlToDataFrame(xml["//Hsp"] , stringsAsFactors=FALSE)

    df <- cbind(
        iterdf[rep(seq_len(nrow(iterdf)), iterlen),, drop=FALSE],
        hitdf[rep(seq_len(nrow(hitdf)), hitlen),, drop=FALSE],
        hsp)
    rownames(df) <- NULL
    df
}

#This subroutine was taken un-changed from annotate source
.tryParseResult <- function(baseUrl, rid, rtoe, timeout) {
    message("estimated response time ", rtoe, " seconds")
    start <- Sys.time()
    end <- Sys.time() + timeout
    url <- sprintf("%s?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=%s",
                   baseUrl, rid)
    Sys.sleep(min(rtoe, timeout))
    repeat {
        elapsed <- as.double(Sys.time() - start, units="secs")
        result <- as(htmlParse(getURL(url, followlocation=TRUE),
                               error = xmlErrorCumulator(immediate=FALSE)),
                     "character")

        if (grepl("Status=FAILED", result))
            stop("BLAST search failed")
        else if  (grepl("Status=UNKNOWN", result))
            stop("BLAST search expired")
        else if (grepl("Status=READY", result)) {
            url <- sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl, rid)
            result <- xmlParse(getURL(url, followlocation=TRUE),
                               error = xmlErrorCumulator(immediate=FALSE))
            return(result)
        } else if (grepl("Status=WAITING", result)) {
            message(sprintf("elapsed time %.0f seconds", elapsed))
            if (Sys.time() > end && interactive()) {
                msg <- sprintf("wait another %d seconds? [y/n] ", timeout)
                repeat {
                    ans <- substr(trimws(tolower(readline(msg))), 1, 1)
                    if (ans %in% c("y", "n"))
                        break
                }
                if (ans == "n")
                    break
                end <- Sys.time() + timeout
            }
            Sys.sleep(10)
        } else
            stop("BLAST search unknown response") 
    }
    msg <- sprintf("'blastSequences' timeout after %.0f seconds",
                   elapsed)
    stop(msg, call.=FALSE)
}

# The majority of this subroutine was taken from annotate source, 
# and modifications were added to support NCBI BLAST AMI.
## Using the REST-ish API described at
## http://www.ncbi.nlm.nih.gov/blast/Doc/node2.html
blastSequencesAws <- function(x, database="nr",
                           hitListSize="10",
                           filter="L",
                           expect="10",
                           program="blastn",
                           timeout=40,
                           as=c("DNAMultipleAlignment", "data.frame", "XML"),
                           instance,
                           dns,
                           threads=8)
{
    PARSE <- switch(match.arg(as),
                    DNAMultipleAlignment=.blastSequencesToDNAMultipleAlignment,
                    data.frame=.blastSequencesToDataFrame,
                    XML=identity)
    ## TODO: lots of argument checking and testing.  Also,
    ## depending on which program string is used we need to make the correct
    ## kind of object at the end (so blastn means DNAMultipleAlignment, and
    ## blastp means AAMultipleAlignment etc.

    ## So:
    ## 1) get online values these parameters can be
    ## 2) document those
    ## 3) restrict their vals in the code here.
    ## 4) for program, use this to determine what object is returned.
    
    ## assemble the query
    #baseUrl <- "https://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
    #instance <- 'i-007ad8d6b0e9981ec'
    #dns <- 'ec2-52-67-171-31.sa-east-1.compute.amazonaws.com'
    baseUrl <- paste0("http://blast:",instance, "@", dns, "/cgi-bin/blast.cgi")
    query <- paste("QUERY=", URLencode(as.character(x)), "&DATABASE=",database,
                   "&HITLIST_SIZE=",hitListSize,"&FILTER=",filter, "&NUM_THREADS=",threads,
                   "&EXPECT=",expect,"&PROGRAM=",program, sep="")
    url0 <- sprintf("%s?%s&CMD=Put", baseUrl, query)
    post <- htmlParse(getURL(url0, followlocation=TRUE))
    
    x <- post[['string(//comment()[contains(., "QBlastInfoBegin")])']]
    print(x)
    rid <- sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
    rtoe <- as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1", x))
    #if(is.na(rtoe)){
    #    rtoe <- 60
    #}
    result <- .tryParseResult(baseUrl, rid, rtoe, timeout)
    PARSE(result)
}

