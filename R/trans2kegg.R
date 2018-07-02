library(dplyr)
library(tidyr)
library(KEGGREST)
blastFile <- "blastSP.tsv"
blastHits <- read.table(blastFile, sep="\t", header=F, quote="")
colnames(blastHits) <- c("gene", "sp", "ident", "evalue", "cov")
spBlast <- separate(blastHits, sp, c("spPrefix", "sp", "symbol"), sep = "\\|")
spBlast <- separate(spBlast, gene, c("gene", "Tran"), sep = "-", remove=F)
spBlast2 <- subset(spBlast, evalue < 1e-10 & cov > 50 & ident > 50,
  select=-c(spPrefix, Tran))
spBlast3 <- separate(spBlast2, symbol, c("symbol", "species"), sep = "_")
head(spBlast3)
speciesCount <- aggregate(sp ~ species, data=spBlast3, FUN=NROW)
