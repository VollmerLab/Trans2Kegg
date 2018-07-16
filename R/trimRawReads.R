trimReads <- function(){
  trimmomaticDir <- '/home/chuck/Trimmomatic-0.36'
  rawReadDir <- "/media/chuck/NGS/Manduca"
  ProjName <- rawReadDir
  fastqExt <- '.fastq'
  rawFastaFiles <- list.files(rawReadDir, paste0('*', fastqExt))
  dfFiles <- data.frame("Path"=rawFastaFiles)
  dfFiles <- separate(dfFiles, "Path", c("Sample", "ext"),
    sep=fastqExt)
  dfFiles <- separate(dfFiles, "Sample", c("Sample", "Direction"),
    sep="\\.")
  dfFiles
  suffix <- gsub("\\.","",unique(dfFiles$ext))
  suffix <- paste(gsub("\\.","",fastqExt), suffix, sep=".")
  samples <- unique(dfFiles$Sample)
  direction <- unique(dfFiles$Direction)

  dfSummary <- data.frame(matrix(ncol=4, nrow=0))
  colnames(dfSummary) <- c("Sample", "In", "Out", "Percent")
  pairedOut <- paste0(rawReadDir, "/Paired")
  unpairedOut <- paste0(rawReadDir, "/Unpaired")
  system(paste('mkdir -p', pairedOut))
  system(paste('mkdir -p', unpairedOut))
  i <- 1
  for(i in i:length(samples)){
    sample <- samples[i]
    print(sample)
    leftRead <- paste(ProjName, paste(sample, direction[1], suffix[1],
      sep="."), sep="/")
    rightRead <- paste(ProjName, paste(sample, direction[2], suffix[1],
      sep="."), sep="/")
    leftPaired <- paste(pairedOut, paste(sample, direction[1], suffix[1],
      sep="."), sep="/")
    leftUnpaired <- paste(unpairedOut, paste(sample, direction[1], suffix[1],
      sep="."), sep="/")
    rightPaired <- paste(pairedOut, paste(sample, direction[2], suffix[1],
      sep="."), sep="/")
    rightUnpaired <- paste(unpairedOut, paste(sample, direction[2], suffix[1],
      sep="."), sep="/")
    trimCmd <- paste('nice -n 19 java -jar',
      paste0(trimmomaticDir,'/trimmomatic-0.36.jar PE'),
      '-threads 4 -phred33',
      leftRead,
      rightRead,
      leftPaired,
      leftUnpaired,
      rightPaired,
      rightUnpaired,
      paste0('ILLUMINACLIP:',trimmomaticDir,
        '/adapters/TruSeq3-PE.fa:2:30:10'),
      'LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:36 2>&1')
    print(trimCmd)
    tryCatch({
      trimLog <- system(trimCmd, intern=T)
      print(trimLog)
      input <- as.numeric(str_match(trimLog[5], "Pairs: (.*?) ")[, 2])
      #Get output pairs
      passed <- as.numeric(str_match(trimLog[5], "Surviving: (.*?) ")[, 2])
      pctPassed <- round(100*(passed/input))
      dfSummary[nrow(dfSummary)+1, ] <- c(sample, input, passed, pctPassed)
      }, error=function(e) print("trim error"))
  }
  save(dfSummary, file="dfSummary.RData")
}


