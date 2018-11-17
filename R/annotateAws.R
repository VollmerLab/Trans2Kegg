#' Annotate DE transcripts using NCBI BLAST AMI
#' @name annotateAws
#' @import tidyr
#' @importFrom RCurl getURL
#' @import XML
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom utils read.csv write.csv write.table URLencode
#' @param accessions A character vector of accessions from FASTA file to
#' be annotated
#' @param refTransFile FASTA file of transcripts to be annotated
#' @param outFile csv output file for annotation results
#' @param instance Instance from Amazon EC2 running NCBI BLAST AMI
#' @param dns DNS from Amazon EC2 running NCBI BLAST AMI
#' @param threads Number of threads to use for BLAST
#' @return Annotation results are written to the csv file specified by outFile
#' @examples
#' instance <- 'i-0eef3748ecac0a146'
#' dns <- 'ec2-35-153-135-164.compute-1.amazonaws.com'
#' filepath <- system.file("extdata", "aiptasia.fa", package="Trans2Kegg")
#' annotateAws(c("KXJ29317.1"), filepath,"annot.csv", instance=instance, 
#' dns=dns, threads=4)
#' annotateTranscripts(c("KXJ29317.1"), filepath,
#'     "annot.csv")
#' @export
annotateAws <- function(accessions, 
                        refTransFile,
                        outFile="annot.csv", 
                        instance, 
                        dns, 
                        threads){
    transDone <- data.frame() 
    fastaOut <- "deTrans.fa"
    if(file.exists(outFile)){
        transDone <- read.csv(outFile, 
            stringsAsFactors=FALSE)$Iteration_query.def
    }
    transLeft <- setdiff(accessions, transDone)
    refTrans <- readDNAStringSet(refTransFile)
    for(accession in transLeft){
        #print(accession)
        transSeqs <- refTrans[accession]
        writeXStringSet(transSeqs, fastaOut, append=FALSE, compress=FALSE, 
            compression_level=NA, format="fasta")
        query <- readChar(fastaOut, file.info(fastaOut)$size)
        blastResult <- blastSequencesAws(query,
        database="swissprot",program="blastx", 
        fmt="data.frame", expect=1e-5,instance=instance,
        dns=dns, threads=threads)
        if(is.data.frame(blastResult) & length(blastResult) > 0){
            blastResult <- subset(blastResult, select=c("Iteration_query-def",
                "Iteration_query-len", "Hit_accession", "Hit_len", 
                "Hsp_evalue", "Hsp_identity", "Hsp_gaps", "Hsp_align-len"))
            write.csv(blastResult, file="blastResult.csv")
            getKegg(blastResult, outFile)
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

#' Get KEGG Ortholog info for BLAST Results
#' @name getKegg
#' @import tidyr
#' @import KEGGREST
#' @importFrom utils read.csv write.csv write.table URLencode
#' @param blastResult BLAST results from blastSequencesAws
#' @param outFile csv output file for annotation results
#' @return Annotation results are written to the csv file specified by outFile
#' @examples
#' filepath <- system.file("extdata", "blastResult.csv", package="Trans2Kegg")
#' blastResult <- read.csv(filepath, stringsAsFactors=FALSE)
#' getKegg(blastResult, "annot.csv")
#' @export
getKegg <- function(blastResult,
    outFile="annot.csv"){
    rowNum <- 0
    accDone <- c()
    dfUniKegg <- data.frame()
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
# R. Gentleman (2018). annotate: Annotation for microarrays. R package
# version 1.61.0.
# This subroutine was taken un-changed from annotate source because
# it is not exported in annotate. This can be eliminated if AWS
# option is incorporated into annotate.
# Aws suffix added to avoid conflict with standard version

#' BLAST Sequences to dataframe
#' @name .blastSequencesToDataFrameAws
#' @import XML
#' @param xml XML BLAST results to be parsed
#' @return Returns parsed BLAST results as dataframe
.blastSequencesToDataFrameAws <- function(xml) {
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

# This subroutine was taken un-changed from annotate source
# R. Gentleman (2018). annotate: Annotation for microarrays. R package
# version 1.61.0.
# This subroutine was taken un-changed from annotate source because
# it is not exported in annotate. This can/will be eliminated
# if AWS option is incorporated into annotate.
# Aws suffix added to avoid conflict with standard version.

#' Parse results from AWS BLAST
#' @name .tryParseResultAws 
#' @importFrom methods as
#' @param baseUrl URL from which to retrieve response
#' @param rid ID associated with BLAST response
#' @param rtoe Estimated response time
#' @param timeout Timeout value
#' @return Returns BLAST results as XML
.tryParseResultAws <- function(baseUrl, rid, rtoe, timeout) {
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
# R. Gentleman (2018). annotate: Annotation for microarrays. R package
# version 1.61.0.
# Modifications were added to support NCBI BLAST AMI. Code is freely
# available for incorporation into annotate if the author chooses to
# do so. 
# Using the REST-ish API described at
# http://www.ncbi.nlm.nih.gov/blast/Doc/node2.html
# Added instance and dns parameters to support AWS
# Added Aws suffix to avoid conflict with standard version
# Changed as to fmt to avoid conflict with function as
blastSequencesAws <- function(x, database="nr",
                                hitListSize="10",
                                filter="L",
                                expect="10",
                                program="blastn",
                                timeout=40,
                                fmt=c("DNAMultipleAlignment", "data.frame", 
                                    "XML"),
                                instance,
                                dns,
                                threads=4)
{

    PARSE <- switch(match.arg(fmt),
                    data.frame=.blastSequencesToDataFrameAws,
                    XML=identity)
    # assemble the query from AWS instance and DNS
    baseUrl <- paste0("http://blast:",instance, "@", dns, "/cgi-bin/blast.cgi")
    query <- paste("QUERY=", URLencode(as.character(x)), "&DATABASE=",database,
            "&HITLIST_SIZE=",hitListSize,"&FILTER=",filter,
            "&NUM_THREADS=",threads, "&EXPECT=",expect,"&PROGRAM=",program, 
            sep="")
    url0 <- sprintf("%s?%s&CMD=Put", baseUrl, query)
    post <- htmlParse(getURL(url0, followlocation=TRUE))
    
    x <- post[['string(//comment()[contains(., "QBlastInfoBegin")])']]
    #print(x)
    rid <- sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
    rtoe <- as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1", x))
    if(is.na(rtoe)){
        return(-1)
    }
    result <- .tryParseResultAws(baseUrl, rid, rtoe, timeout)
    PARSE(result)
}

