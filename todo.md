Subdirectory 'R' contains invalid file names:
  ‘annot.csv’ ‘blastResult.csv’ ‘deTrans.fa’ ‘en.utf-8.add’
  ‘en.utf-8.add.spl’ ‘uniToKO.csv’
   * NOTE: Avoid sapply(); use vapply()
      Found in files:
        annotateAws.R (line 89, column 16)
        annotateAws.R (line 93, column 15)
* Checking function lengths................
    * NOTE: Recommended function length <= 50 lines.
      There are 1 functions > 50 lines.
      The longest 5 functions are:
        getPathways() (R/getPathway.R, line 13): 78 lines
        annotateAws() (R/annotateAws.R, line 23): 46 lines
        annotateNCBI() (R/annotateNCBI.R, line 19): 45 lines
        .tryParseResultAws() (R/annotateAws.R, line 122): 43 lines
        mergeAnnotations() (R/mergeAnnotations.R, line 15): 42 lines
* Checking man page documentation...
    * NOTE: Consider adding runnable examples to the following man
      pages which document exported objects:
      getKegg.Rd
   * NOTE: Consider multiples of 4 spaces for line indents, 23
      lines(2%) are not.
    First 6 lines:
      R/annotateNCBI.R:20   refTrans <- readDNAStringSet(fasta)
      R/annotateNCBI.R:21   rowNum <- 0
      R/annotateNCBI.R:22   accDone <- c()
      R/annotateNCBI.R:23   dfUniKegg <- data.frame()
      R/annotateNCBI.R:24   # Check if BLAST already started. If already star...
      R/annotateNCBI.R:25   # where you left off.
    See http://bioconductor.org/developers/how-to/coding-style/
    See FormatR package:
      https://cran.r-project.org/web/packages/formatR/index.html


